#include "inviscidBurgersSolver.hpp"

#include <filesystem>
#include <fstream>
#include <iomanip>
#include <sstream>

inviscidBurgersSolver::inviscidBurgersSolver(Bounds calc_bounds, double endtime,
                                             int n_timestep, int n_cells,
                                             double wave_velocity)
    : calc_bounds_(calc_bounds),
      endtime_(endtime),
      n_timestep_(n_timestep),
      n_cells_(n_cells),
      c_(wave_velocity)
{
    dx_  = (calc_bounds_[1] - calc_bounds_[0]) / n_cells_;
    dt_  = endtime_ / static_cast<double>(n_timestep_);
    CFL_ = c_ * (dt_ / dx_);

    std::cout << "dx : " << dx_ << std::endl;
    std::cout << "dt : " << dt_ << std::endl;
    std::cout << "Courant Number : " << CFL_ << std::endl;

    U.resize(n_cells_, 0);
    x.resize(n_cells_, 0);
}

inviscidBurgersSolver::~inviscidBurgersSolver() {}

void inviscidBurgersSolver::writeFlowFile(std::string filename)
{
    std::ofstream ofile(filename);
    if (!ofile)
    {
        std::cerr << "Cannot file open..." << std::endl;
        std::exit(1);
    }

    ofile << "x\t"
          << "U\t" << std::endl;
    for (std::size_t i = 0; i < U.size(); ++i)
    {
        ofile << x[i] << "\t" << U[i] << std::endl;
    }
}

void inviscidBurgersSolver::solve()
{
    std::cout << "\nStart solving flow..." << std::endl;

    for (int tstep = 0; tstep < n_timestep_; ++tstep)
    {
        Vector U_new(U);
        for (int i = 0; i < n_cells_; ++i)
        {
            if (i == 0 || i == n_cells_)
                U_new[i] = U[i];
            else
                U_new[i] = 0.5 * (U[i + 1] + U[i - 1]) -
                           0.5 * CFL_ * (U[i + 1] - U[i - 1]);
        }

        // write file
        bool isOutputDir = std::filesystem::exists("./flow_output");
        if (!isOutputDir)
            std::filesystem::create_directory("./flow_output");
        else
        {
            std::string ofile_flow;
            std::stringstream ss;
            ss << "flow_" << std::setfill('0') << std::setw(4) << std::right
               << std::to_string(tstep) << ".dat";
            ss >> ofile_flow;
            writeFlowFile("./flow_output/" + ofile_flow);
        }
        // upadate flowfield
        U = U_new;
    }

    std::cout << "Finish solving flow!" << std::endl;
}

int main()
{
    /*--- Setting up Configurations ---*/
    Bounds calc_bounds_  = {0.0, 10.0};
    double endtime       = 0.05;
    int n_timestep       = 100;
    int n_cells          = 20;
    double wave_velocity = 343;

    /*--- Configure ---*/
    inviscidBurgersSolver flow(calc_bounds_, endtime, n_timestep, n_cells,
                               wave_velocity);

    int NCells = flow.get_n_cells();
    double dx  = flow.get_dx();

    for (int i = 0; i < NCells; ++i)
    {
        if (i < NCells / 2)
            flow.U[i] = 1;
        else
            flow.U[i] = 0;
    }
    for (int i = 0; i < NCells; ++i)
    {
        if (i == 0)
            flow.x[i] = 0;
        else
            flow.x[i] = flow.x[i - 1] + dx;
    }

    /*--- Solve flow ---*/
    flow.solve();

    return 0;
}
