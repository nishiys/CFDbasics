#include "ShockTubeSolver.hpp"

#include <cmath>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <sstream>





ShockTubeSolver::ShockTubeSolver(Bounds calc_bounds, double endtime, double dt,
                                 int n_cells)
    : calc_bounds_(calc_bounds),
      endtime_(endtime),
      dt_(dt),
      n_cells_(n_cells)
{
    dx_         = (calc_bounds_[1] - calc_bounds_[0]) / n_cells_;
    n_timestep_ = static_cast<int>(endtime / dt_);

    std::cout << "[Configurations]: " << std::endl;
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

void ShockTubeSolver::SetSchemes(std::string convectiveflux_scheme, bool MUSCL, double k)
{
    convectiveflux_scheme_ = convectiveflux_scheme;
    MUSCL_ = MUSCL;
    k_ = k;

    if (MUSCL_ == true)
    {
        std::cout << "MUSCL scheme enabled." << std::endl;
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


double ShockTubeSolver::minmod(double a, double b)
{
    if (a*b <= 0)
        return 0.0;
    else
    {
        if (std::abs(a) < std::abs(b))
            return a;
        else
            return b;
    }
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
    /* --- Set Variables --- */
    double rhoL = primL(0);
    double uL   = primL(1);
    double pL   = primL(2);
    double aL   = std::sqrt(GAMMA * pL / rhoL);  // sound speed
    double HL =
        aL * aL / (GAMMA - 1) + 0.5 * uL * uL;  // Total enthalpy per unit mass
    double ML = uL / aL;                        // Mach No.

    double rhoR = primR(0);
    double uR   = primR(1);
    double pR   = primR(2);
    double aR   = std::sqrt(GAMMA * pR / rhoR);
    double HR   = aR * aR / (GAMMA - 1) + 0.5 * uR * uR;
    double MR   = uR / aR;

    /* --- Calculate Positive & Negative Fluxes --- */
    Eigen::Vector3d PositiveFlux(3);
    Eigen::Vector3d NegativeFlux(3);

    if (ML <= -1.0)
    {
        Eigen::Vector3d fluxR(3);
        fluxR(0) = rhoR * uR;
        fluxR(1) = pR + rhoR * uR * uR;
        fluxR(2) = (rhoR * HR) * uR;

        PositiveFlux = Eigen::VectorXd::Zero(3);
        NegativeFlux = fluxR;
    }

    else if (ML >= 1.0)
    {
        Eigen::Vector3d fluxL(3);
        fluxL(0) = rhoL * uL;
        fluxL(1) = pL + rhoL * uL * uL;
        fluxL(2) =
            (rhoL * HL) * uL;  // rho * H : Total enthalpy per unit volume

        PositiveFlux = fluxL;
        NegativeFlux = Eigen::VectorXd::Zero(3);
    }

    else  // -1 < ML < 1
    {
        // Positive Flux
        double Mp       = 0.25 * (ML + 1) * (ML + 1);
        double MLG1     = 1 + 0.5 * (GAMMA - 1) * ML;
        PositiveFlux(0) = rhoL * aL * Mp;
        PositiveFlux(1) = Mp * (2.0 * rhoL * aL * aL / GAMMA) * MLG1;
        PositiveFlux(2) = Mp *
                          (2.0 * rhoL * aL * aL * aL / (GAMMA * GAMMA - 1)) *
                          (MLG1 * MLG1);
        // Negative Flux
        double Mn       = -0.25 * (MR - 1) * (MR - 1);
        double MRG1     = -1 + 0.5 * (GAMMA - 1) * MR;
        NegativeFlux(0) = rhoR * aR * Mn;
        NegativeFlux(1) = Mn * (2.0 * rhoR * aR * aR / GAMMA) * MRG1;
        NegativeFlux(2) = Mn *
                          (2.0 * rhoR * aR * aR * aR / (GAMMA * GAMMA - 1)) *
                          (MRG1 * MRG1);
    }

    return (PositiveFlux + NegativeFlux);  // (F_j+1/2)^+ + (F_j+1/2)^-
}

Eigen::Vector3d ShockTubeSolver::AUSM(Eigen::Vector3d primL,
                                      Eigen::Vector3d primR)
{
    /* --- Set Variables --- */
    double rhoL = primL(0);
    double uL   = primL(1);
    double pL   = primL(2);
    double aL   = std::sqrt(GAMMA * pL / rhoL);  // sound speed
    double HL =
        aL * aL / (GAMMA - 1) + 0.5 * uL * uL;  // Total enthalpy per unit mass
    double ML = uL / aL;                        // Mach No.

    double rhoR = primR(0);
    double uR   = primR(1);
    double pR   = primR(2);
    double aR   = std::sqrt(GAMMA * pR / rhoR);
    double HR   = aR * aR / (GAMMA - 1) + 0.5 * uR * uR;
    double MR   = uR / aR;

    /* --- Calculate M_1/2 --- */
    double MLp;  // (ML)^+
    if (std::abs(ML) > 1)
    {
        MLp = 0.5 * (ML + std::abs(ML));
    }
    else
    {
        MLp = 0.25 * (ML + 1) * (ML + 1);
    }
    double MRm;  // (MR)^-
    if (std::abs(MR) > 1)
    {
        MRm = 0.5 * (MR - std::abs(MR));
    }
    else
    {
        MRm = -0.25 * (MR - 1) * (MR - 1);
    }
    double M_half = MLp + MRm;

    /* --- Calculate Flux without pressure (F_1/2)^c --- */
    Eigen::Vector3d Fc;
    if (M_half >= 0)
    {
        Fc(0) = M_half * (rhoL * aL);
        Fc(1) = M_half * (rhoL * aL * uL);
        Fc(2) = M_half * (rhoL * aL * HL);
    }
    else
    {
        Fc(0) = M_half * (rhoR * aR);
        Fc(1) = M_half * (rhoR * aR * uR);
        Fc(2) = M_half * (rhoR * aR * HR);
    }
    /* --- Calculate p_half --- */
    double pLp;
    if (std::abs(ML) > 1)
    {
        pLp = 0.5 * pL * (ML + std::abs(ML)) / ML;
    }
    else
    {
        pLp = 0.5 * pL * (1 + ML);
    }

    double pRm;
    if (std::abs(MR) > 1)
    {
        pRm = 0.5 * pR * (MR - std::abs(MR)) / MR;
    }
    else
    {
        pRm = 0.5 * pR * (1 - MR);
    }

    double p_half = pLp + pRm;

    /* --- Return the resulting flux --- */
    Fc(1) = Fc(1) + p_half;
    return Fc;
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
    std::cout << "Convective flux scheme: " << convectiveflux_scheme_
              << std::endl;
    std::cout << "Time accuracy: 2nd order" << std::endl;


    // Initialize Flow Convervatives
    flowConservatives.resize(n_cells_);
    for (int i = 0; i < n_cells_; ++i)
    {
        flowConservatives[i] = PrimToCons(flowVars[i]);
    }

    // Start timestep calculations
    for (int tstep = 0; tstep < n_timestep_; ++tstep)
    {
        std::vector<Eigen::Vector3d> flowConservatives_new(flowConservatives);
        std::vector<Eigen::Vector3d> flowConservatives_halfstep(flowConservatives);
        for (int i = 0; i < n_cells_; ++i)
        {
            // BCs.
            if (i == 0 || i == 1 || i == n_cells_ - 2 || i == n_cells_ - 1)
            {
                flowConservatives_new[i] = flowConservatives[i];
            }
            // Other than Boundaries
            else
            {
                /* --- Scheme of 2nd order accuracy in time is used. --- */
                // <1st step> Calculate Q^n+1/2 with 1st order accurate fluxes.
                Eigen::Vector3d QL_right = flowConservatives[i];
                Eigen::Vector3d QR_right =flowConservatives[i+1];
                Eigen::Vector3d QL_left = flowConservatives[i-1];
                Eigen::Vector3d QR_left = flowConservatives[i];
                Eigen::Vector3d F_right(3);  // F_j+1/2
                Eigen::Vector3d F_left(3);   // F_j-1/2

                /* --- Convective flux evaluation --- */
                if (convectiveflux_scheme_ == "Roe")
                {
                    F_right = Roe(ConsToPrim(QL_right),
                                  ConsToPrim(QR_right));
                    F_left  = Roe(ConsToPrim(QL_left),
                                 ConsToPrim(QR_left));
                }
                else if (convectiveflux_scheme_ == "vanLeer")
                {
                    F_right = vanLeer(ConsToPrim(QL_right),
                                      ConsToPrim(QR_right));
                    F_left  = vanLeer(ConsToPrim(QL_left),
                                     ConsToPrim(QR_left));
                }
                else if (convectiveflux_scheme_ == "AUSM")
                {
                    F_right = AUSM(ConsToPrim(QL_right),
                                   ConsToPrim(QR_right));
                    F_left  = AUSM(ConsToPrim(QL_left),
                                  ConsToPrim(QR_left));
                }
                else
                {
                    std::cerr << convectiveflux_scheme_
                              << " scheme is not supported..." << std::endl;
                    std::exit(1);
                }
                // Update conservative quantities
                flowConservatives_halfstep[i] =
                    flowConservatives[i] - 0.5 * dt_ / dx_ * (F_right - F_left);

                // <2nd step> Calculate Fluxes with Q^n+1/2 we get at the 1st step
                /* --- MUSCL scheme --- */
                if (MUSCL_ == true)
                {
                    // with minmod
                    double b = (3 - k_)/(1 - k_);
                    // fd: fowrard difference, bd: backward difference
                    Eigen::Vector3d fdL_right = flowConservatives[i+1] - flowConservatives[i];
                    Eigen::Vector3d bdL_right = flowConservatives[i] - flowConservatives[i-1];
                    Eigen::Vector3d fdR_right = flowConservatives[i+2] - flowConservatives[i+1];
                    Eigen::Vector3d bdR_right = flowConservatives[i+1] - flowConservatives[i];
                    
                    Eigen::Vector3d fdL_left = flowConservatives[i] - flowConservatives[i-1];
                    Eigen::Vector3d bdL_left = flowConservatives[i-1] - flowConservatives[i-2];
                    Eigen::Vector3d fdR_left = flowConservatives[i+1] - flowConservatives[i];
                    Eigen::Vector3d bdR_left = flowConservatives[i] - flowConservatives[i-1];

                    for (int i = 0; i < 3; ++i)
                    {
                        fdL_right(i) = minmod(fdL_right(i), b * bdL_right(i));
                        bdL_right(i) = minmod(bdL_right(i), b * fdL_right(i));
                        fdR_right(i) = minmod(fdR_right(i), b * bdR_right(i));
                        bdR_right(i) = minmod(bdR_right(i), b * fdR_right(i));

                        fdL_left(i) = minmod(fdL_left(i), b * bdL_left(i));
                        bdL_left(i) = minmod(bdL_left(i), b * fdL_left(i));
                        fdR_left(i) = minmod(fdR_left(i), b * bdR_left(i));
                        bdR_left(i) = minmod(bdR_left(i), b * fdR_left(i));
                    }

                    QL_right = flowConservatives_halfstep[i] + 0.25 * ((1 - k_) * bdL_right + (1 + k_) * fdL_right);
                    QR_right = flowConservatives_halfstep[i+1] - 0.25 * ((1 - k_) * fdR_right + (1 + k_) * bdR_right);
                    QL_left = flowConservatives_halfstep[i-1] + 0.25 * ((1 - k_) * bdL_left + (1 + k_) * fdL_left);
                    QR_left = flowConservatives_halfstep[i] - 0.25 * ((1 - k_) * fdR_left + (1 + k_) * bdR_left);
                }
                else
                {
                    // 1st order upwind
                    QL_right = flowConservatives[i];
                    QR_right =flowConservatives[i+1];
                    QL_left = flowConservatives[i-1];
                    QR_left = flowConservatives[i];
                }

                /* --- Convective flux evaluation --- */
                if (convectiveflux_scheme_ == "Roe")
                {
                    F_right = Roe(ConsToPrim(QL_right),
                                  ConsToPrim(QR_right));
                    F_left  = Roe(ConsToPrim(QL_left),
                                 ConsToPrim(QR_left));
                }
                else if (convectiveflux_scheme_ == "vanLeer")
                {
                    F_right = vanLeer(ConsToPrim(QL_right),
                                      ConsToPrim(QR_right));
                    F_left  = vanLeer(ConsToPrim(QL_left),
                                     ConsToPrim(QR_left));
                }
                else if (convectiveflux_scheme_ == "AUSM")
                {
                    F_right = AUSM(ConsToPrim(QL_right),
                                   ConsToPrim(QR_right));
                    F_left  = AUSM(ConsToPrim(QL_left),
                                  ConsToPrim(QR_left));
                }
                else
                {
                    std::cerr << convectiveflux_scheme_
                              << " scheme is not supported..." << std::endl;
                    std::exit(1);
                }

                flowConservatives_new[i] =
                    flowConservatives[i] - dt_ / dx_ * (F_right - F_left);
                flowVars[i] = ConsToPrim(flowConservatives_new[i]);
            }
        }

        // upadate flowfield
        flowConservatives = flowConservatives_new;
    }

    std::cout << "Finish solving flow!\n" << std::endl;
}


// void ShockTubeSolver::Solve()
// {
//     std::cout << "\nStart solving flow..." << std::endl;
//     std::cout << "Convective flux scheme: " << convectiveflux_scheme_
//               << std::endl;

//     // Initialize Flow Convervatives
//     flowConservatives.resize(n_cells_);
//     for (int i = 0; i < n_cells_; ++i)
//     {
//         flowConservatives[i] = PrimToCons(flowVars[i]);
//     }

//     // Start timestep calculations
//     for (int tstep = 0; tstep < n_timestep_; ++tstep)
//     {
//         std::vector<Eigen::Vector3d> flowConservatives_new(flowConservatives);
//         for (int i = 0; i < n_cells_; ++i)
//         {
//             // BCs.
//             if (i == 0 || i == 1 || i == n_cells_ - 2 || i == n_cells_ - 1)
//             {
//                 flowConservatives_new[i] = flowConservatives[i];
//             }
//             // Other than Boundaries
//             else
//             {
//                 Eigen::Vector3d F_right(3);  // F_j+1/2
//                 Eigen::Vector3d F_left(3);   // F_j-1/2

//                 /* --- MUSCL scheme --- */
//                 Eigen::Vector3d QL_right(3);
//                 Eigen::Vector3d QR_right(3);
//                 Eigen::Vector3d QL_left(3);
//                 Eigen::Vector3d QR_left(3);

//                 if (MUSCL_ == true)
//                 {
//                     // with minmod
//                     double b = (3 - k_)/(1 - k_);
//                     // fd: fowrard difference, bd: backward difference
//                     Eigen::Vector3d fdL_right = flowConservatives[i+1] - flowConservatives[i];
//                     Eigen::Vector3d bdL_right = flowConservatives[i] - flowConservatives[i-1];
//                     Eigen::Vector3d fdR_right = flowConservatives[i+2] - flowConservatives[i+1];
//                     Eigen::Vector3d bdR_right = flowConservatives[i+1] - flowConservatives[i];
                    
//                     Eigen::Vector3d fdL_left = flowConservatives[i] - flowConservatives[i-1];
//                     Eigen::Vector3d bdL_left = flowConservatives[i-1] - flowConservatives[i-2];
//                     Eigen::Vector3d fdR_left = flowConservatives[i+1] - flowConservatives[i];
//                     Eigen::Vector3d bdR_left = flowConservatives[i] - flowConservatives[i-1];

//                     for (int i = 0; i < 3; ++i)
//                     {
//                         fdL_right(i) = minmod(fdL_right(i), b * bdL_right(i));
//                         bdL_right(i) = minmod(bdL_right(i), b * fdL_right(i));
//                         fdR_right(i) = minmod(fdR_right(i), b * bdR_right(i));
//                         bdR_right(i) = minmod(bdR_right(i), b * fdR_right(i));

//                         fdL_left(i) = minmod(fdL_left(i), b * bdL_left(i));
//                         bdL_left(i) = minmod(bdL_left(i), b * fdL_left(i));
//                         fdR_left(i) = minmod(fdR_left(i), b * bdR_left(i));
//                         bdR_left(i) = minmod(bdR_left(i), b * fdR_left(i));
//                     }

//                     QL_right = flowConservatives[i] + 0.25 * ((1 - k_) * bdL_right + (1 + k_) * fdL_right);
//                     QR_right = flowConservatives[i+1] - 0.25 * ((1 - k_) * fdR_right + (1 + k_) * bdR_right);
//                     QL_left = flowConservatives[i-1] + 0.25 * ((1 - k_) * bdL_left + (1 + k_) * fdL_left);
//                     QR_left = flowConservatives[i] - 0.25 * ((1 - k_) * fdR_left + (1 + k_) * bdR_left);
//                 }
//                 else
//                 {
//                     // 1st order upwind
//                     QL_right = flowConservatives[i];
//                     QR_right =flowConservatives[i+1];
//                     QL_left = flowConservatives[i-1];
//                     QR_left = flowConservatives[i];
//                 }

//                 /* --- Convective flux evaluation --- */
//                 if (convectiveflux_scheme_ == "Roe")
//                 {
//                     F_right = Roe(ConsToPrim(QL_right),
//                                   ConsToPrim(QR_right));
//                     F_left  = Roe(ConsToPrim(QL_left),
//                                  ConsToPrim(QR_left));
//                 }
//                 else if (convectiveflux_scheme_ == "vanLeer")
//                 {
//                     F_right = vanLeer(ConsToPrim(QL_right),
//                                       ConsToPrim(QR_right));
//                     F_left  = vanLeer(ConsToPrim(QL_left),
//                                      ConsToPrim(QR_left));
//                 }
//                 else if (convectiveflux_scheme_ == "AUSM")
//                 {
//                     F_right = AUSM(ConsToPrim(QL_right),
//                                    ConsToPrim(QR_right));
//                     F_left  = AUSM(ConsToPrim(QL_left),
//                                   ConsToPrim(QR_left));
//                 }
//                 else
//                 {
//                     std::cerr << convectiveflux_scheme_
//                               << " scheme is not supported..." << std::endl;
//                     std::exit(1);
//                 }

//                 flowConservatives_new[i] =
//                     flowConservatives[i] - dt_ / dx_ * (F_right - F_left);
//                 flowVars[i] = ConsToPrim(flowConservatives_new[i]);
//             }
//         }

//         // upadate flowfield
//         flowConservatives = flowConservatives_new;
//     }

//     std::cout << "Finish solving flow!\n" << std::endl;
// }





int main()
{
    /*--- Setting up Configurations ---*/
    Bounds calc_bounds_ = {0.0, 1.0};
    double endtime      = 0.2;
    double dt           = 1e-4;
    int n_cells         = 201;
    std::string convectiveflux_scheme;
    bool MUSCL;
    double kappa;
    /*--- Initial Conditions ---*/
    double wall_position = 0.5;
    double rhoL          = 1.0;
    double pL            = 1.0;
    double rhoR          = 0.125;
    double pR            = 0.1;

    // Roe
    convectiveflux_scheme = "Roe";
    MUSCL = false;
    kappa = 0;
    ShockTubeSolver flow_roe(calc_bounds_, endtime, dt, n_cells);
    flow_roe.SetGrid();
    flow_roe.SetFlowField(wall_position, rhoL, pL, rhoR, pR);
    flow_roe.SetSchemes(convectiveflux_scheme, MUSCL, kappa);
    flow_roe.Solve();
    flow_roe.WriteFlowFile("flowData_Roe.dat");

    // van Leer
    convectiveflux_scheme = "vanLeer";
    MUSCL = false;
    kappa = 0;
    ShockTubeSolver flow_vanLeer(calc_bounds_, endtime, dt, n_cells);
    flow_vanLeer.SetGrid();
    flow_vanLeer.SetFlowField(wall_position, rhoL, pL, rhoR, pR);
    flow_vanLeer.SetSchemes(convectiveflux_scheme, MUSCL, kappa);
    flow_vanLeer.Solve();
    flow_vanLeer.WriteFlowFile("flowData_vanLeer.dat");

    // AUSM
    convectiveflux_scheme = "AUSM";
    MUSCL = false;
    kappa = 0;
    ShockTubeSolver flow_AUSM(calc_bounds_, endtime, dt, n_cells);
    flow_AUSM.SetGrid();
    flow_AUSM.SetFlowField(wall_position, rhoL, pL, rhoR, pR);
    flow_AUSM.SetSchemes(convectiveflux_scheme, MUSCL, kappa);
    flow_AUSM.Solve();
    flow_AUSM.WriteFlowFile("flowData_AUSM.dat");

    // Roe MUSCL
    convectiveflux_scheme = "Roe";
    MUSCL = true;
    kappa = -1;
    ShockTubeSolver flow_roe_MUSCL(calc_bounds_, endtime, dt, n_cells);
    flow_roe_MUSCL.SetGrid();
    flow_roe_MUSCL.SetFlowField(wall_position, rhoL, pL, rhoR, pR);
    flow_roe_MUSCL.SetSchemes(convectiveflux_scheme, MUSCL, kappa);
    flow_roe_MUSCL.Solve();
    flow_roe_MUSCL.WriteFlowFile("flowData_Roe_MUSCL.dat");

    // van Leer MUSCL
    convectiveflux_scheme = "vanLeer";
    MUSCL = true;
    kappa = -1;
    ShockTubeSolver flow_vanLeer_MUSCL(calc_bounds_, endtime, dt, n_cells);
    flow_vanLeer_MUSCL.SetGrid();
    flow_vanLeer_MUSCL.SetFlowField(wall_position, rhoL, pL, rhoR, pR);
    flow_vanLeer_MUSCL.SetSchemes(convectiveflux_scheme, MUSCL, kappa);
    flow_vanLeer_MUSCL.Solve();
    flow_vanLeer_MUSCL.WriteFlowFile("flowData_vanLeer_MUSCL.dat");

    // AUSM MUSCL
    convectiveflux_scheme = "AUSM";
    MUSCL = true;
    kappa = -1;
    ShockTubeSolver flow_AUSM_MUSCL(calc_bounds_, endtime, dt, n_cells);
    flow_AUSM_MUSCL.SetGrid();
    flow_AUSM_MUSCL.SetFlowField(wall_position, rhoL, pL, rhoR, pR);
    flow_AUSM_MUSCL.SetSchemes(convectiveflux_scheme, MUSCL, kappa);
    flow_AUSM_MUSCL.Solve();
    flow_AUSM_MUSCL.WriteFlowFile("flowData_AUSM_MUSCL.dat");

    return 0;
}
