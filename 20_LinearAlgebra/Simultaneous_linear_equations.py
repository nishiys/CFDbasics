# -*- coding: utf-8 -*-
"""
Created on Sun Oct 30 17:04:38 2016

@author: Owner
"""


import numpy as np
import time
import matplotlib.pyplot as plt

# ===============  Gauss Elimination =============== #


def Gaussian_elimination(n, A, b):
    x = np.zeros(n)

    # Forward elimination
    for k in range(1, n):  # e.g. for column 1
        for i in range(k + 1, n + 1):  # from line 2 to n
            w = A[i-1, k-1]/A[k-1, k-1]

            for j in range(k, n + 1):  # from line 1 to n in column
                A[i-1, j-1] = A[i-1, j-1] - w*A[k-1, j-1]

            b[i-1] = b[i-1] - w*b[k-1]

    x[n-1] = b[n-1]/A[n-1, n-1]
    # print np.rint(A)#でAを整数(四捨六入で５は偶数になるほうに丸める)に変換して表示（上三角)

    # Back substitution
    for i in range(n - 1, 0, -1):
        for j in range(i + 1, n + 1):
            b[i-1] = b[i-1] - A[i-1, j-1]*x[j-1]

        x[i-1] = b[i-1]/A[i-1, i-1]

    return x
    print(x)
# ==================================================== #

# ============== Gauss Elimination w/ Partial Pivoting ============= #


def Gaussian_elimination_pivot(n, A, b):
    x = np.zeros(n)

    piv = np.zeros(n).astype(int)  # numpyオブジェクト.astype(型)でまとめて型変換
    for i in range(0, n):
        piv[i] = i  # piv = [0,1,2,...,n-1]

    for k in range(0, n):
        amax = A[piv[k], k]  # k列の(k,k)以下の行の要素のなかで最大値を探す
        imax = k
        for i in range(k+1, n):
            if amax < abs(A[piv[i], k]):
                amax = abs(A[piv[i], k])
                imax = i
            else:
                pass
            if imax != k:  # Aの行を直接は入れ替えず、pivで入れ替えを記憶
                a = piv[k]
                piv[k] = piv[imax]
                piv[imax] = a
            else:
                pass
        for i in range(k+1, n):
            w = A[piv[i], k]/A[piv[k], k]
            for j in range(k+1, n):
                A[piv[i], j] = A[piv[i], j] - w*A[piv[k], j]

            b[piv[i]] = b[piv[i]] - w*b[piv[k]]

    x[n-1] = b[piv[n-1]]/A[piv[n-1], n-1]

    for i in range(n-2, -1, -1):
        for j in range(i+1, n):
            b[piv[i]] = b[piv[i]] - A[piv[i], j]*x[j]
        x[i] = b[piv[i]]/A[piv[i], i]

    return x
    print(x)
# ===================================================== #

# ================== Jacobi method ====================== #


def Jacobi_method(n, A, b, eps):
    # A = L + D + Uと分解し、Dx = (-L-U)x + bを収束させる。
    e = 1.0  # 誤差ベクトル（epsより大きければなんでもいい）
    list_e = []

    x = np.zeros(n)
    for i in range(0, n):
        x[i] = b[i]/A[i, i]
        # Aの対角成分が支配的となることがおおいので、解の初期値は何の情報もない場合はこうしておく。
    D = np.diagflat(np.diag(A))
    T = -(A - D)

    while e > eps:
        list_e.append(e)
        x_old = x
        x = np.dot(np.diagflat(1/np.diag(A)), np.dot(T, x_old) + b)
        e = np.linalg.norm(x - x_old)

    return x
    print(x)
    print("Jacobi法ーー誤差:{0}".format(e))

    plt.plot(list_e)
    plt.show()

# ================================================================= #

# ============================= Gauss-Seidel method ======================= #


def Gauss_Seidel(n, A, b, eps):
    # A = L + D + Uと分解し、(D+L)x = -Ux + bを収束させる。
    e = 1.0  # 誤差ベクトル（epsより大きければなんでもいい）
    list_e = []

    x = np.zeros(n)
    for i in range(0, n):
        x[i] = b[i]/A[i, i]
        # Aの対角成分が支配的となることがおおいので、解の初期値は何の情報もない場合はこうしておく。

    S = np.tril(A)  # S = L + D
    T = -(A - S)
    S_inv = np.linalg.inv(S)

    while e > eps:
        list_e.append(e)
        x_old = x
        x = np.dot(S_inv, np.dot(T, x_old) + b)
        e = np.linalg.norm(x-x_old)

    return x
    print(x)
    print("Gauss-Seidel法ーー誤差:{0}".format(e))

    plt.plot(list_e)
    plt.show()

# ================================================================ #

# ====================== SOR method ==================== #


def SOR(n, A, b, eps, omega):
    # A = L + D + Uと分解し、(D + omega*L)x = ((1-omega)D - U)x + omega*bを収束させる。
    e = 1.0  # 誤差ベクトル（epsより大きければなんでもいい）
    list_e = []

    x = np.zeros(n)
    for i in range(0, n):
        x[i] = b[i]/A[i, i]
        # Aの対角成分が支配的となることがおおいので、解の初期値は何の情報もない場合はこうしておく。

    L = np.tril(A, -1)  # 下三角行列を返す関数。parameter<0で対角成分含まず下のみ抽出。
    U = np.triu(A, 1)  # 上三角行列を返す関数。parameter>0で対角成分含まず上のみ抽出。
    D = np.diagflat(np.diag(A))

    S = D + omega*L
    T = (1-omega)*D - omega*U
    S_inv = np.linalg.inv(S)

    while e > eps:
        list_e.append(e)
        x_old = x
        x = np.dot(S_inv, np.dot(T, x_old) + omega*b)
        e = np.linalg.norm(x-x_old)

    return x
    print(x)
    print("SOR法ーー誤差:{0}".format(e))

    plt.plot(list_e)
    plt.show

# =========================================================== #


if __name__ == '__main__':

    # 1 Benchmark solution by numpy.linalg.solve()
    start = time.time()

# initialization
    A = np.array([[5., -3., -1., 0., 2., 1.],
                  [0., -5., 4., 1., 0., 2.],
                  [-1., 3., -5., -2., 3., 3.],
                  [-1., 0., 3., 4., -2., -1.],
                  [0., 3., 3., -1., 3., -4.],
                  [2., 3., 2., 3., 2., -5.]])

    b = np.array([12., 6., 3., 9., 12., 21.])

# ----------------------main routine----------------------- #
    x = np.linalg.solve(A, b)
    print(x)

    elapsed_time = time.time() - start
    print("numpy.linalg.solve()")
    print("  Execution time: {0}".format(elapsed_time))
# --------------------------------------------------------- #


# 2 Gauss Elimination (normal)
    start = time.time()

# initialization
    n = 6

    A = np.array([[5., -3., -1., 0., 2., 1.],
                  [0., -5., 4., 1., 0., 2.],
                  [-1., 3., -5., -2., 3., 3.],
                  [-1., 0., 3., 4., -2., -1.],
                  [0., 3., 3., -1., 3., -4.],
                  [2., 3., 2., 3., 2., -5.]])

    b = np.array([12., 6., 3., 9., 12., 21.])

# -------------------------main routine------------------ #
    Gaussian_elimination(n, A, b)

    elapsed_time = time.time() - start
    print("Gauss Elimination（normal):")
    print("  Execution time: {0}".format(elapsed_time))
# ------------------------------------------------------------- #


# 3 Gauss Elimination（w/ partial pivot)
    start = time.time()

# initialization
    n = 6

    A = np.array([[5., -3., -1., 0., 2., 1.],
                  [0., -5., 4., 1., 0., 2.],
                  [-1., 3., -5., -2., 3., 3.],
                  [-1., 0., 3., 4., -2., -1.],
                  [0., 3., 3., -1., 3., -4.],
                  [2., 3., 2., 3., 2., -5.]])

    b = np.array([12., 6., 3., 9., 12., 21.])

# ------------------main routine----------------------- #
    Gaussian_elimination_pivot(n, A, b)

    elapsed_time = time.time() - start
    print("Gauss Elimination（w/ partial pivot):")
    print("  Execution time: {0}".format(elapsed_time))
# -------------------------------------------------------- #


# 4 Jacobi
    start = time.time()

# initialization
    n = 6

    A = np.array([[5., -3., -1., 0., 2., 1.],
                  [0., -5., 4., 1., 0., 2.],
                  [-1., 3., -5., -2., 3., 3.],
                  [-1., 0., 3., 4., -2., -1.],
                  [0., 3., 3., -1., 3., -4.],
                  [2., 3., 2., 3., 2., -5.]])

    b = np.array([12., 6., 3., 9., 12., 21.])

# ----------------------main routine-------------------------- #
    Jacobi_method(n, A, b, 1e-10)

    elapsed_time = time.time() - start
    print("Jacobi method:")
    print("  Execution time: {0}".format(elapsed_time))

# ----------------------------------------------------------- #


# 5 Gauss-Seidel
    start = time.time()

# initialization
    n = 6

    A = np.array([[5., -3., -1., 0., 2., 1.],
                  [0., -5., 4., 1., 0., 2.],
                  [-1., 3., -5., -2., 3., 3.],
                  [-1., 0., 3., 4., -2., -1.],
                  [0., 3., 3., -1., 3., -4.],
                  [2., 3., 2., 3., 2., -5.]])

    b = np.array([12., 6., 3., 9., 12., 21.])

# ---------------- main routine ----------------- #
    Gauss_Seidel(n, A, b, 1e-10)

    elapsed_time = time.time() - start
    print("Gauss-Seidel")
    print("  Execution time: {0}".format(elapsed_time))

# -------------------------------------------------- #

# 6 SOR
    start = time.time()

# initialization
    n = 6

    A = np.array([[5., -3., -1., 0., 2., 1.],
                  [0., -5., 4., 1., 0., 2.],
                  [-1., 3., -5., -2., 3., 3.],
                  [-1., 0., 3., 4., -2., -1.],
                  [0., 3., 3., -1., 3., -4.],
                  [2., 3., 2., 3., 2., -5.]])

    b = np.array([12., 6., 3., 9., 12., 21.])

# ---------------------- main routine --------------------- #
    SOR(n, A, b, 1e-10, 0.8)

    elapsed_time = time.time() - start
    print("SOR method")
    print("  Execution time: {0}".format(elapsed_time))

# --------------------------------------------------------- #
