import matplotlib.pyplot as plt
import numpy as np


class Flowfield:
    def __init__(self):
        self.x_array = np.empty(0)
        self.rho_array = np.empty(0)
        self.u_array = np.empty(0)
        self.p_array = np.empty(0)

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
            print(self.x_array.shape[0])

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

    def plot1Ddata(self):
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

        plt.savefig("flowData.png")
        plt.show()


def main():

    flowfield = Flowfield()
    flowfield.getFlowField("flowData.dat")

    flowfield.plot1Ddata()


if __name__ == "__main__":
    main()
