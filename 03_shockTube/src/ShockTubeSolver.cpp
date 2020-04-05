#include "ShockTubeSolver.hpp"

#include <cmath>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <sstream>

ShockTubeSolver::ShockTubeSolver(Bounds calc_bounds, double endtime, double dt,
                                 int n_cells)
    : calc_bounds_(calc_bounds), endtime_(endtime), dt_(dt), n_cells_(n_cells)
{
    dx_         = (calc_bounds_[1] - calc_bounds_[0]) / n_cells_;
    n_timestep_ = static_cast<int>(endtime / dt_);

    std::cout << "dx : " << dx_ << std::endl;
    std::cout << "dt : " << dt_ << std::endl;
}

ShockTubeSolver::~ShockTubeSolver() {}

void ShockTubeSolver::SetGrid()
{
    cells.resize(n_cells_, 0);
    for (int i = 0; i < n_cells_; ++i)
    {
        cells[i] = i * dx_;
    }
}

void ShockTubeSolver::SetFlowField(double wall_position, double rhoL, double pL,
                                   double rhoR, double pR)
{
    flowVars.resize(n_cells_);

    for (int i = 0; i < n_cells_; ++i)
    {
        if (cells[i] <= wall_position)
        {
            Eigen::Vector3d prim(3);
            prim(0)     = rhoL;
            prim(1)     = 0.0;
            prim(2)     = pL;
            flowVars[i] = prim;
        }
        else
        {
            Eigen::Vector3d prim(3);
            prim(0)     = rhoR;
            prim(1)     = 0.0;
            prim(2)     = pR;
            flowVars[i] = prim;
        }
    }
}

Eigen::Vector3d ShockTubeSolver::PrimToCons(Eigen::Vector3d prim)
{
    Eigen::Vector3d cons(3);
    cons(0) = prim(0);            // rho
    cons(1) = prim(0) * prim(1);  // rho * u
    double H =
        (GAMMA / (GAMMA - 1)) * prim(2) / prim(0) + 0.5 * prim(1) * prim(1);
    cons(2) = prim(0) * H - prim(2);  // rho * E = rho * H - p

    return cons;
}
Eigen::Vector3d ShockTubeSolver::ConsToPrim(Eigen::Vector3d cons)
{
    Eigen::Vector3d prim(3);
    prim(0)      = cons(0);  // rho
    double r_rho = 1.0 / cons(0);
    prim(1)      = r_rho * cons(1);  // u
    prim(2)      = (GAMMA - 1) *
              (cons(2) - 0.5 * prim(0) * prim(1) *
                             prim(1));  // p = (gamma - 1) * (E - 1/2 rho u^2)

    return prim;
}

Eigen::Vector3d ShockTubeSolver::Roe(Eigen::Vector3d primL,
                                     Eigen::Vector3d primR)
{
    /* --- Set Variables --- */
    double rhoL = primL(0);
    double uL   = primL(1);
    double pL   = primL(2);
    double aL   = std::sqrt(GAMMA * pL / rhoL);  // sound speed
    double HL =
        aL * aL / (GAMMA - 1) + 0.5 * uL * uL;  // Total enthalpy per unit mass

    double rhoR = primR(0);
    double uR   = primR(1);
    double pR   = primR(2);
    double aR   = std::sqrt(GAMMA * pR / rhoR);
    double HR   = aR * aR / (GAMMA - 1) + 0.5 * uR * uR;

    /* --- Roe average --- */
    double r       = std::sqrt(rhoR / rhoL);
    double rho_ave = r * rhoL;
    double u_ave   = (uL + r * uR) / (1.0 + r);
    double H_ave   = (HL + r * HR) / (1.0 + r);
    double c_ave   = std::sqrt((GAMMA - 1) * (H_ave - 0.5 * u_ave * u_ave));

    /* --- Calculate Eigenvectors --- */
    double lambda1 = std::abs(u_ave - c_ave);
    double lambda2 = std::abs(u_ave);
    double lambda3 = std::abs(u_ave + c_ave);
    Eigen::Matrix3d A_diag;
    A_diag       = Eigen::Matrix3d::Zero();
    A_diag(0, 0) = lambda1;
    A_diag(1, 1) = lambda2;
    A_diag(2, 2) = lambda3;

    // std::cout << "lambda1: " << lambda1 << std::endl;

    /* --- Right Eigenvectors --- */
    Eigen::MatrixXd R(3, 3);
    R(0, 0) = 1;
    R(0, 1) = 1;
    R(0, 0) = 1;
    R(1, 0) = u_ave - c_ave;
    R(1, 1) = u_ave;
    R(1, 2) = u_ave + c_ave;
    R(2, 0) = H_ave - u_ave * c_ave;
    R(2, 1) = 0.5 * u_ave * u_ave;
    R(2, 2) = H_ave + u_ave * c_ave;

    /* --- Light Eigenvectors --- */
    double b1 = 0.5 * u_ave * u_ave * (GAMMA - 1) / (c_ave * c_ave);
    double b2 = (GAMMA - 1) / (c_ave * c_ave);
    Eigen::MatrixXd R_inv(3, 3);
    R_inv(0, 0) = 0.5 * (b1 + u_ave / c_ave);
    R_inv(0, 1) = -0.5 * (1.0 / c_ave + b2 * u_ave);
    R_inv(0, 2) = 0.5 * b2;
    R_inv(1, 0) = 1.0 - b1;
    R_inv(1, 1) = b2 * u_ave;
    R_inv(1, 2) = -b2;
    R_inv(2, 0) = 0.5 * (b1 - u_ave / c_ave);
    R_inv(2, 1) = 0.5 * (1.0 / c_ave - b2 * u_ave);
    R_inv(2, 2) = 0.5 * b2;

    /* --- Flux Jacobian Matrix --- */
    Eigen::MatrixXd A_ave(3, 3);
    A_ave = R * A_diag * R_inv;

    /* --- Compute Flux --- */
    Eigen::Vector3d QL          = PrimToCons(primL);
    Eigen::Vector3d QR          = PrimToCons(primR);
    Eigen::Vector3d dissipation = A_ave * (QR - QL);

    Eigen::Vector3d fluxL(3);
    fluxL(0) = rhoL * uL;
    fluxL(1) = pL + rhoL * uL * uL;
    fluxL(2) = (rhoL * HL) * uL;  // rho * H : Total enthalpy per unit volume

    Eigen::Vector3d fluxR(3);
    fluxR(0) = rhoR * uR;
    fluxR(1) = pR + rhoR * uR * uR;
    fluxR(2) = (rhoR * HR) * uR;

    Eigen::Vector3d roe_flux(3);
    roe_flux = 0.5 * (fluxR + fluxL - dissipation);
    return roe_flux;
}

Eigen::Vector3d ShockTubeSolver::vanLeer(Eigen::Vector3d primL,
                                         Eigen::Vector3d primR)
{
}

void ShockTubeSolver::WriteFlowFile(std::string filename)
{
    std::ofstream ofile(filename);
    if (!ofile)
    {
        std::cerr << "Cannot file open..." << std::endl;
        std::exit(1);
    }

    ofile << "x\t"
          << "rho\t"
          << "u\t"
          << "p" << std::endl;

    for (std::size_t i = 0; i < flowVars.size(); ++i)
    {
        ofile << cells[i] << "\t" << flowVars[i](0) << "\t" << flowVars[i](1)
              << "\t" << flowVars[i](2) << std::endl;
    }
}

void ShockTubeSolver::Solve()
{
    std::cout << "\nStart solving flow..." << std::endl;

    // Initialize Flow Convervatives
    flowConservatives.resize(n_cells_);
    for (int i = 0; i < n_cells_; ++i)
    {
        flowConservatives[i] = PrimToCons(flowVars[i]);
    }

    // Start timestep calculations
    for (int tstep = 0; tstep < n_timestep_; ++tstep)
    {
        // // write file
        // bool isOutputDir = std::filesystem::exists("./flow_output");
        // if (!isOutputDir)
        //     std::filesystem::create_directory("./flow_output");
        // else
        // {
        //     std::string ofile_flow;
        //     std::stringstream ss;
        //     ss << "flow_" << std::setfill('0') << std::setw(4) << std::right
        //        << std::to_string(tstep) << ".dat";
        //     ss >> ofile_flow;
        //     WriteFlowFile("./flow_output/" + ofile_flow);
        // }

        std::vector<Eigen::Vector3d> flowConservatives_new(flowConservatives);
        for (int i = 0; i < n_cells_; ++i)
        {
            // BCs.
            if (i == 0 || i == n_cells_ - 1)
            {
                flowConservatives_new[i] = flowConservatives[i];
            }
            else
            {
                Eigen::Vector3d E_right(3);  // E_j+1/2
                Eigen::Vector3d E_left(3);   // E_j-1/2

                E_right = Roe(ConsToPrim(flowConservatives[i]),
                              ConsToPrim(flowConservatives[i + 1]));
                E_left  = Roe(ConsToPrim(flowConservatives[i - 1]),
                             ConsToPrim(flowConservatives[i]));

                flowConservatives_new[i] =
                    flowConservatives[i] - dt_ / dx_ * (E_right - E_left);
                flowVars[i] = ConsToPrim(flowConservatives_new[i]);
                // std::cout << flowConservatives_new[i] << std::endl;
            }
        }

        // upadate flowfield
        flowConservatives = flowConservatives_new;
    }

    std::cout << "Finish solving flow!" << std::endl;
}

int main()
{
    /*--- Setting up Configurations ---*/
    Bounds calc_bounds_ = {0.0, 1.0};
    double endtime      = 0.25;
    double dt           = 1e-4;
    int n_cells         = 100;

    /*--- Configure ---*/
    ShockTubeSolver flow(calc_bounds_, endtime, dt, n_cells);
    flow.SetGrid();
    /*--- Initial Conditions ---*/
    double wall_position = 0.5;
    double rhoL          = 1.0;
    double pL            = 1.0;
    double rhoR          = 0.125;
    double pR            = 0.1;
    flow.SetFlowField(wall_position, rhoL, pL, rhoR, pR);

    /*--- Solve flow ---*/
    flow.Solve();
    flow.WriteFlowFile("flowData.dat");
    return 0;
}
