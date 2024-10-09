import numpy as np
import matplotlib.pyplot as plt
from numpy.polynomial.legendre import leggauss

N = 50
a_val = -1
b_val = 1
h_val = 1 / (N + 1)

def laplacian_matrix(n):
    h = 1 / (n + 1)
    return (2 * np.eye(n) - np.diag(np.ones(n - 1), -1) - np.diag(np.ones(n - 1), 1)) / h

x_vals = np.linspace(0, 1, N + 2)

def hat_function(x, left, mid, right):
    if left <= x <= mid:
        return (x - left) / (mid - left)
    elif mid <= x <= right:
        return (right - x) / (right - mid)
    else:
        return 0

def f(x):
    if 0.4 < x < 0.6:
        return -1
    else:
        return 0

def calculate_integral(f, left, mid, right, method):
    if method == 'trapeze':
        h = right - left
        return (h / 2) * (f(left) + f(right))
    elif method == 'quadrature':
        nodes, weights = leggauss(N)
        integral = 0
        for i in range(len(nodes)):
            xi = 0.5 * (left + right) + 0.5 * (right - left) * nodes[i]
            integral += weights[i] * f(xi) * hat_function(xi, left, mid, right)
        return integral * 0.5 * (right - left)

F_interior_trapeze = np.zeros(N - 2)
F_left_trapeze = a_val / h_val - calculate_integral(f, 0, x_vals[0], x_vals[1], method='trapeze')
F_right_trapeze = b_val / h_val - calculate_integral(f, x_vals[-2], x_vals[-1], 1, method='trapeze')

for i in range(1, N - 1):
    F_interior_trapeze[i - 1] = -calculate_integral(f, x_vals[i - 1], x_vals[i], x_vals[i + 1], method='trapeze')

F_trapeze = np.concatenate(([F_left_trapeze], F_interior_trapeze, [F_right_trapeze]))

F_interior_quadrature = np.zeros(N - 2)
F_left_quadrature = a_val / h_val - calculate_integral(f, 0, x_vals[0], x_vals[1], method='quadrature')
F_right_quadrature = b_val / h_val - calculate_integral(f, x_vals[-2], x_vals[-1], 1, method='quadrature')

for i in range(1, N - 1):
    F_interior_quadrature[i - 1] = -calculate_integral(f, x_vals[i - 1], x_vals[i], x_vals[i + 1], method='quadrature')

F_quadrature = np.concatenate(([F_left_quadrature], F_interior_quadrature, [F_right_quadrature]))

U_trapeze = np.linalg.solve(laplacian_matrix(N), F_trapeze)
U_quadrature = np.linalg.solve(laplacian_matrix(N), F_quadrature)

U_trapeze = np.concatenate(([a_val], U_trapeze, [b_val]))
U_quadrature = np.concatenate(([a_val], U_quadrature, [b_val]))

plt.plot(x_vals, U_trapeze[:N + 2], label='Approximation (Trapèze)')
plt.plot(x_vals, U_quadrature[:N + 2], label='Approximation (Quadrature)')

plt.xlabel('x')
plt.ylabel('u(x)')
plt.legend()
plt.title("solution approchée pour N=50")
plt.show()

err_trapeze_quadrature = np.sqrt(h_val) * np.linalg.norm(U_trapeze - U_quadrature[:N + 2])

print("Erreur entre U_trapeze et U_quadrature =", err_trapeze_quadrature)