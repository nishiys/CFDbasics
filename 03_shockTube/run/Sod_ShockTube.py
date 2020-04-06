import numpy as np
from scipy.optimize import fsolve
import matplotlib.pyplot as plt


class SodExact:
    # class variables
    GAMMA = 1.4

    # instance variables

    def __init__(self, time, cells: np.array,
                 x_wall, primL: np.array, primR: np.array):
        self.time = time
        self.x_wall = x_wall
        self.cells = cells
        self.N = cells.shape[0]
        self.primL = primL
        self.primR = primR
        self.solution = np.zeros((3, self.N))
        self.aL = self.SoundSpeed(primL[2], primL[0])
        self.aR = self.SoundSpeed(primR[2], primR[0])
        self.a3 = 0

    def SoundSpeed(self, pressure, density):
        a = np.sqrt(self.GAMMA * pressure / density)
        return a

    def ShockTubeEq(self, Ms):
        pL = self.primL[2]
        pR = self.primR[2]

        GR1 = (self.GAMMA + 1) / (self.GAMMA - 1)
        GR2 = (self.GAMMA - 1)/(2 * self.GAMMA)
        GR3 = (2 * self.GAMMA) / (self.GAMMA + 1)
        RHS = self.aL * GR1 * (1 - (pR/pL * (GR3 * Ms**2 - 1/GR1)) ** GR2)
        LHS = Ms - 1 / Ms

        return (LHS - RHS)

    def SolveShockTubeEq(self):
        guess = 1.0
        sol = fsolve(self.ShockTubeEq, guess)
        return sol[0]

    def CalcRegion4(self, Ms: float):
        # Rankine-Hugoniot relations
        GR1 = (self.GAMMA + 1) / (self.GAMMA - 1)
        GR3 = (2 * self.GAMMA) / (self.GAMMA + 1)
        GR4 = 2 / (self.GAMMA + 1)
        pR = self.primR[2]
        p4 = pR * (GR3 * Ms**2 - 1 / GR1)
        rhoR = self.primR[0]
        r_rho4 = 1 / rhoR * (GR4/(Ms**2) + 1/GR1)
        rho4 = 1 / r_rho4
        u4 = GR4 * (Ms - 1 / Ms)

        return rho4, u4, p4

    def CalcRegion3(self, u4, p4):
        u3 = u4
        p3 = p4
        # isotropic flow
        rhoL = self.primL[0]
        pL = self.primL[2]
        rho3 = rhoL * (p3 / pL)**(1 / self.GAMMA)

        return rho3, u3, p3

    def ReturnRegion(self, x: float, Us, u3):
        x12 = self.x_wall - self.aL * self.time
        x23 = self.x_wall + (u3 - self.a3) * self.time
        x34 = self.x_wall + u3 * self.time
        x45 = self.x_wall + Us * self.time

        if (x <= x12):
            region = 1
        elif (x12 < x < x23):
            region = 2
        elif (x23 <= x < x34):
            region = 3
        elif (x34 <= x < x45):
            region = 4
        elif (x45 <= x):
            region = 5

        return region

    # Solve Sod ShockTube problem

    def Solve(self):
        Ms = self.SolveShockTubeEq()
        Us = self.aR * Ms
        print("Ms: ", Ms)
        print("Us: ", Us)

        # Region 4
        rho4, u4, p4 = self.CalcRegion4(Ms)
        print("rho4: ", rho4)
        print("u4: ", u4)
        print("p4: ", p4)
        # Region 3
        rho3, u3, p3 = self.CalcRegion3(u4, p4)
        self.a3 = self.SoundSpeed(p3, rho3)
        print("rho3: ", rho3)
        print("u3: ", u3)
        print("p3: ", p3)

        # Construct Solutions
        for i in range(0, self.N):
            REGION = self.ReturnRegion(self.cells[i], Us, u3)
            if (REGION == 1):
                self.solution[0][i] = self.primL[0]
                self.solution[1][i] = self.primL[1]
                self.solution[2][i] = self.primL[2]
            elif (REGION == 2):
                x = self.cells[i]
                u2 = 2 / (self.GAMMA + 1) * \
                    (self.aL + (x - self.x_wall) / self.time)
                a2 = self.aL - (self.GAMMA - 1) * u2 / 2
                p2 = self.primL[2] * \
                    (a2 / self.aL) ** (2 * self.GAMMA / (self.GAMMA - 1))
                rho2 = self.GAMMA * p2 / (a2**2)
                self.solution[0][i] = rho2
                self.solution[1][i] = u2
                self.solution[2][i] = p2
            elif (REGION == 3):
                self.solution[0][i] = rho3
                self.solution[1][i] = u3
                self.solution[2][i] = p3
            elif (REGION == 4):
                self.solution[0][i] = rho4
                self.solution[1][i] = u4
                self.solution[2][i] = p4
            elif (REGION == 5):
                self.solution[0][i] = self.primR[0]
                self.solution[1][i] = self.primR[1]
                self.solution[2][i] = self.primR[2]
        return 0


def main():
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
    sod = SodExact(endtime, cells, x_wall, primL, primR)
    sod.Solve()

    # Plot solution
    fig = plt.figure(1)
    plt.plot(sod.cells, sod.solution[0, :])
    plt.ylabel("density")
    fig.savefig("density_Sod_Solution.png")

    fig = plt.figure(2)
    plt.plot(sod.cells, sod.solution[1, :])
    plt.ylabel("velocity")
    fig.savefig("velocity_Sod_Solution.png")

    fig = plt.figure(3)
    plt.plot(sod.cells, sod.solution[2, :])
    plt.ylabel("pressure")
    fig.savefig("pressure_Sod_Solution.png")


if __name__ == "__main__":
    main()
