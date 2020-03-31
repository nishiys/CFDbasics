import matplotlib.pyplot as plt
from pathlib import Path
import os

import numpy as np
import matplotlib.animation as animation


def getFileList():
    filepath = Path(__file__)
    print(filepath)
    cdir = Path(__file__).parent.resolve()  # get parent dir path as absolute
    print(cdir)

    os.chdir(cdir)
    print(os.getcwd())

    print(type(cdir.glob('*.dat')))
    file_list = sorted(list(cdir.glob('flow_*.dat')))

    return file_list


class Flowfield:
    def __init__(self):
        self.x_array = np.empty(0)
        self.U_array = np.empty(0)

    def getFlowField(self, file_list):
        with open(file_list[0], 'r') as file1st:
            file1st.readline()
            count = 0
            for line in file1st:
                count += 1
            self.x_array = np.zeros((1, count))
            self.U_array = np.zeros((1, count))
            print(self.x_array)

        i = 0
        for filename in file_list:
            with open(filename, 'r') as file:
                file.readline()
                j = 0
                for line in file:
                    line = line.strip().split("\t")
                    x = float(line[0])
                    U = float(line[1])
                    self.x_array[i][j] = x
                    self.U_array[i][j] = U
                    j += 1

                self.x_array = np.append(self.x_array, np.zeros(
                    (1, self.x_array.shape[1])), axis=0)
                self.U_array = np.append(self.U_array, np.zeros(
                    (1, self.U_array.shape[1])), axis=0)
            i += 1
        self.x_array = self.x_array[0:-1, :]
        self.U_array = self.U_array[0:-1, :]

    def plot1Ddata(self):
        fig = plt.figure()
        ims = []

        for i in range(0, self.x_array.shape[1]):
            im = plt.plot(self.x_array[i], self.U_array[i])
            ims.append(im)
        ani = animation.ArtistAnimation(fig, ims, interval=200)
        plt.show()


def main():

    os.chdir("./flow_output")

    flowfield = Flowfield()
    flowfield.x_array = np.empty([0, 0])
    flowfield.U_array = np.empty([0, 0])

    file_list = getFileList()
    flowfield.getFlowField(file_list)

    print(flowfield.x_array)
    print(flowfield.U_array)

    flowfield.plot1Ddata()


if __name__ == "__main__":
    main()
