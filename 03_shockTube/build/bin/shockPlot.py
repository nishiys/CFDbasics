import matplotlib.pyplot as plt
import numpy as np

import Sod_ShockTube


class Flowfield:
    def __init__(self, schemename):
        self.x_array = np.empty(0)
        self.rho_array = np.empty(0)
        self.u_array = np.empty(0)
        self.p_array = np.empty(0)
        self.scheme = schemename

    def getFlowField(self, filename):
        with open(filename, 'r') as file:
            file.readline()
            data_length = 0
            for line in file:
                data_length += 1
            self.x_array = np.zeros(data_length)
            self.rho_array = np.zeros(data_length)
            self.u_array = np.zeros(data_length)
            self.p_array = np.zeros(data_length)

        with open(filename, 'r') as file:
            file.readline()
            i = 0
            for line in file:
                line = line.strip().split("\t")
                x = float(line[0])
                rho = float(line[1])
                u = float(line[2])
                p = float(line[3])
                self.x_array[i] = x
                self.rho_array[i] = rho
                self.u_array[i] = u
                self.p_array[i] = p
                i += 1

    def plotSolution(self, filename):

        fig = plt.figure(figsize=(8, 12))
        # plt.subplots_adjust(hspace=0.3)

        rho_fig = fig.add_subplot(311)
        rho_fig.plot(self.x_array, self.rho_array)
        rho_fig.set_ylabel(r'$\rho$')

        u_fig = fig.add_subplot(312)
        u_fig.plot(self.x_array, self.u_array)
        u_fig.set_ylabel(r'$u$')

        p_fig = fig.add_subplot(313)
        p_fig.plot(self.x_array, self.p_array)
        p_fig.set_xlabel(r'$x [m]$')
        p_fig.set_ylabel(r'$p$')

        plt.savefig(filename)
        plt.show()

    def plotSolutionWithExactSol(self, exactSol: Sod_ShockTube.SodExact):
        fig = plt.figure(figsize=(8, 12))
        # plt.subplots_adjust(hspace=0.3)

        rho_fig = fig.add_subplot(311)
        rho_fig.plot(self.x_array, self.rho_array)
        rho_fig.plot(exactSol.cells,
                     exactSol.solution[0, :])
        rho_fig.set_ylabel(r'$\rho$')

        u_fig = fig.add_subplot(312)
        u_fig.plot(self.x_array, self.u_array)
        u_fig.plot(exactSol.cells,
                   exactSol.solution[1, :])
        u_fig.set_ylabel(r'$u$')

        p_fig = fig.add_subplot(313)
        p_fig.plot(self.x_array, self.p_array)
        p_fig.plot(exactSol.cells,
                   exactSol.solution[2, :])
        p_fig.set_xlabel(r'$x [m]$')
        p_fig.set_ylabel(r'$p$')

        plt.savefig("flowData_withSol.png")
        plt.show()


def plotFlowData(flowfield_list: list,
                 exactSol: Sod_ShockTube.SodExact):
    rho_fig = plt.figure(1)
    ax_rho = rho_fig.add_subplot(111)
    u_fig = plt.figure(2)
    ax_u = u_fig.add_subplot(111)
    p_fig = plt.figure(3)
    ax_p = p_fig.add_subplot(111)

    n = len(flowfield_list)
    for i in range(0, n):
        ax_rho.plot(flowfield_list[i].x_array, flowfield_list[i].rho_array,
                    label=flowfield_list[i].scheme)
        ax_u.plot(flowfield_list[i].x_array, flowfield_list[i].u_array,
                  label=flowfield_list[i].scheme)
        ax_p.plot(flowfield_list[i].x_array, flowfield_list[i].p_array,
                  label=flowfield_list[i].scheme)

    ax_rho.plot(exactSol.cells,
                exactSol.solution[0, :], label="Exact Solution")
    ax_u.plot(exactSol.cells, exactSol.solution[1, :], label="Exact Solution")
    ax_p.plot(exactSol.cells, exactSol.solution[2, :], label="Exact Solution")

    ax_rho.set_ylabel(r'$\rho$')
    ax_rho.set_xlabel(r'$x [m]$')
    ax_u.set_ylabel(r'$u$')
    ax_u.set_xlabel(r'$x [m]$')
    ax_p.set_ylabel(r'$p$')
    ax_p.set_xlabel(r'$x [m]$')

    rho_fig.legend()
    u_fig.legend()
    p_fig.legend()

    plt.show()
    rho_fig.savefig("rho.png")
    u_fig.savefig("u.png")
    p_fig.savefig("p.png")


def main():

    # flowfield_roe = Flowfield("Roe")
    # flowfield_roe.getFlowField("flowData_Roe.dat")
    # # flowfield_roe.plotSolution("flowData_Roe.png")

    # flowfield_vanLeer = Flowfield("van Leer")
    # flowfield_vanLeer.getFlowField("flowData_vanLeer.dat")
    # # flowfield_vanLeer.plotSolution("flowData_vanLeer.png")

    # flowfield_AUSM = Flowfield("AUSM")
    # flowfield_AUSM.getFlowField("flowData_AUSM.dat")
    # # flowfield_vanLeer.plotSolution("flowData_vanLeer.png")

    flowfield_roe_MUSCL = Flowfield("Roe w/ MUSCL")
    flowfield_roe_MUSCL.getFlowField("flowData_Roe_MUSCL.dat")

    flowfield_vanLeer_MUSCL = Flowfield("van Leer w/ MUSCL")
    flowfield_vanLeer_MUSCL.getFlowField("flowData_vanLeer_MUSCL.dat")

    flowfield_AUSM_MUSCL = Flowfield("AUSM w/ MUSCL")
    flowfield_AUSM_MUSCL.getFlowField("flowData_AUSM_MUSCL.dat")

    flowfield_list = []
    # flowfield_list.append(flowfield_roe)
    # flowfield_list.append(flowfield_vanLeer)
    # flowfield_list.append(flowfield_AUSM)
    flowfield_list.append(flowfield_roe_MUSCL)
    flowfield_list.append(flowfield_vanLeer_MUSCL)
    flowfield_list.append(flowfield_AUSM_MUSCL)

    # --- Exact Solution --- #
    # Initial conditions
    rhoL = 1.0
    uL = 0.0
    pL = 1.0
    primL = np.array([rhoL, uL, pL])

    rhoR = 0.125
    uR = 0.0
    pR = 0.1
    primR = np.array([rhoR, uR, pR])

    # calculation region
    x_min = 0.0
    x_max = 1.0
    N = 1000
    cells = np.linspace(x_min, x_max, num=N)
    # wall position
    x_wall = 0.5
    # time settings
    endtime = 0.2

    # Solve flow
    sod = Sod_ShockTube.SodExact(endtime, cells, x_wall, primL, primR)
    sod.Solve()

    # Plot solution
    plotFlowData(flowfield_list, sod)


if __name__ == "__main__":
    main()
