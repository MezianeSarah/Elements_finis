import numpy as np
import matplotlib.pyplot as plt
from scipy import integrate

N = 50
i = np.arange(N)
Points = np.cos((N + 1 - i) * np.pi / (2 * (N + 1)))
print(i)
print(Points)

def createMesh(Points, deg):
    N = len(Points)
    if deg == 1:
        Nint = N
        Next = 2
        Nodes = np.zeros((N + 2, 2))
        for i in range(0, N):
            Nodes[i, 0] = Points[i]
        Nodes[N, 0] = 0
        Nodes[N, 1] = 1
        Nodes[N + 1, 0] = 1
        Nodes[N + 1, 1] = 1
        Cells = np.zeros((N + 1, 2))
        for i in range(1, N + 1):
            Cells[i, 0] = i - 1
            Cells[i, 1] = i
        Cells[0, 0] = Nint
        Cells[N, 1] = Nint + 1
        Cells = Cells.astype(int)
    return Nint, Next, Nodes, Cells

Nint, Next, Nodes, Cells = createMesh(Points, 1)
print('Verification', Nint, Next)
print('Tableau de noeuds')
print(Nodes)
print('Tableau de cellules')
print(Cells)

def RefBasisFct(xhat, deg):
    res = np.zeros(deg + 1)
    if deg == 1:
        res[0] = 1 - xhat
        res[1] = xhat
    if deg == 2:
        res[0] = 2 * (xhat ** 2) - 3 * xhat + 1
        res[1] = -4 * (xhat ** 2) + 4 * xhat
        res[2] = 2 * (xhat ** 2) - xhat
    return res

def DerRefBasisFct(xhat, deg):
    res = np.zeros(deg + 1)
    if deg == 1:
        res[0] = -1
        res[1] = 1
    if deg == 2:
        res[0] = 4 * xhat - 3
        res[1] = -8 * xhat + 4
        res[2] = 4 * xhat - 1
    return res

print('Verification')
print('Valeur au point x_hat=0 :', RefBasisFct(0.4, 1), DerRefBasisFct(0.4, 1))

def FK2(xhat, noeud1, noeud2):
    h = noeud2 - noeud1
    x = xhat * h + noeud1
    return x

def jacFK2(xhat, noeud1, noeud2):
    return abs(noeud2 - noeud1)

def BasisFct(xhat, j, noeud1, noeud2, deg):
    vals = RefBasisFct(xhat, deg)
    return vals

def DerBasisFct(xhat, j, noeud1, noeud2, deg):
    h = abs(noeud2 - noeud1)
    vals = DerRefBasisFct(xhat, deg) / h
    return vals

def main(f, deg, Points, ex):
    print("Creation du maillage")
    Nint, Next, Nodes, Cells = createMesh(Points, deg)
    print("Assemblage de A")
    A = np.zeros((Nint, Nint))
    for j in range(0, np.size(Cells, 0)):
        for mu in range(0, np.size(Cells, 1)):
            for lamb in range(0, np.size(Cells, 1)):
                i = Cells[j, lamb]
                k = Cells[j, mu]
                if (i < Nint) and (k < Nint):
                    x_j = Nodes[Cells[j, 0], 0]
                    x_jp1 = Nodes[Cells[j, 1], 0]
                    g = lambda xhat: DerRefBasisFct(xhat, deg)[lamb] * DerRefBasisFct(xhat, deg)[mu] / jacFK2(xhat, x_j, x_jp1)
                    A[i, k] = A[i, k] + integrate.quad(g, 0, 1)[0]
    
    print("Assemblage du second membre")
    F = np.zeros((Nint, 1))
    for j in range(0, np.size(Cells, 0)):
        for lamb in range(0, np.size(Cells, 1)):
            i = Cells[j, lamb]
            if (i < Nint):
                x_j = Nodes[Cells[j, 0], 0]
                x_jp1 = Nodes[Cells[j, 1], 0]
                g = lambda xhat: -RefBasisFct(xhat, deg)[lamb] * f(FK2(xhat, x_j, x_jp1)) * jacFK2(xhat, x_j, x_jp1)
                F[i] = F[i] + integrate.quad(g, 0, 1)[0]
    
    h = 1 - Points[len(Points) - 1]
    F[Nint - 1] = F[Nint - 1] + 1 / h
    
    print(A)
    print(F)
    
    print("Resolution du systeme lineaire")
    U = np.linalg.solve(A, F)
    
    print("Affichage...")
    X = np.zeros(Nint + Next)
    X[0] = 0
    X[Nint + Next - 1] = 1
    X[1:Nint + 1] = Nodes[0:Nint, 0]
    UU = np.zeros(Nint + Next)
    UU[1:Nint + 1] = U[0:Nint, 0]
    UU[0] = 0
    UU[Nint + 1] = 1
    
    plt.figure()
    plt.plot(X, UU, 'b', label='P1-Lag')
    plt.plot(X, ex(X), 'k', label='exacte')
    plt.legend()
    plt.title("solution approchée d'un maillage non uniforme pour N=50")
    plt.show()
    print("That's all, folks!!")

def f(x):
    return -((np.pi ** 2) / 4) * np.sin((np.pi * x) / 2)

def g(x):
    return np.sin((np.pi * x) / 2)

main(f, 1, Points, g)
