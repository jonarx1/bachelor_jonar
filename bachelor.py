# -*- coding: utf-8 -*-
"""
Spyder Editor

This is a temporary script file.
"""

import numpy as np
import quadpy
import matplotlib.pyplot as plt

def sphere_size(points, radius):
    #Changes the radius of the sphere
    new_radius = points * radius
    return new_radius

def sphere_position(points, new_position):
    # Changes the center position of the sphere from [0,0,0] to new_position
    for i,x in enumerate(points):
        x[0] = x[0]+new_position[0]
        x[1] = x[1]+new_position[1]
        x[2] = x[2]+new_position[2]
    return points

def total_areal(grid_areals):
    result = sum(grid_areals)
    return result

def S(points, weights, k=1.0694):
    # Create a square matrix
    s = np.zeros((points.shape[0], points.shape[0]))
    # Fill the square matrix
    for i, p in enumerate(points):
        for j, q in enumerate(points):
            if i == j:
                s[i, i] = k * np.sqrt( 4 * np.pi / (R ** 2 * weights[i])) 
            else:
                s[i, j] = ( 1 / np.linalg.norm(p-q))
    return s

def D_matrix(points, weights, radius, normals, k=1.0694):
    D = np.zeros((points.shape[0], points.shape[0]))
    for i, p in enumerate(points):
        for j, q in enumerate(points):
            if i == j:
                D[i, i] = k * np.sqrt( 4 * np.pi * (weights[i])) / (2 * radius)
            else:
                D[i, j] = ( np.dot(p-q, normals[j]) / np.linalg.norm(p-q)**3 )
    return D

def A_matrix(weights):
    matrix = np.diagflat(weights)
    return matrix

def potential(points, position_charge, charge):
    # Create array for the potential
    V = np.zeros(points.shape[0]) #matches number of lebedev points
    # Fill array with potential
    for i, x in enumerate(points):
        for A, R in enumerate(position_charge):
            V[i] += charge[A] / np.linalg.norm(x-R)
    return V

def COSMO(S, V):
    # Invert S
    S_matrix_inverted = np.linalg.inv(S)
    #print(S_matrix_inverted)
    # Matrix-vector multiply S^-1 V
    sigma = -np.dot(S_matrix_inverted, V)
    return sigma

def Solvation_energy(V_i, q_i):
    E = np.dot(V_i, q_i)
    E = E * 0.5 
    return E


def IEF(S, D, A, V, epsilon):
    A_matrix_inverted = np.linalg.inv(A)
    #print(A_matrix_inverted)
    R = 2 * np.pi * A_matrix_inverted - D
    G = 2 * np.pi * ((epsilon + 1)/(epsilon-1)) * A_matrix_inverted - D
    T = np.dot(G, S)
    K = np.dot(np.linalg.inv(T), R)
    sigma = -np.dot(K, V)
    return sigma

def testing(Z):
    a = np.zeros((Z.shape[0],3))
    for i,z in enumerate(Z):
        a[i, 2] = z
    return a

def dipole_testing(Z):
    x = [0.1, -0.1]
    a = np.zeros((Z.shape[0], 3))
    b = np.zeros((Z.shape[0], 3))
    #c = np.zeros((Z.shape[0], 1))
    #c = np.array([])
    c = np.zeros((2*Z.shape[0], 3))
    for i,z in enumerate(Z):
        a[i, 0] = x[0]
        b[i, 0] = x[1]
        a[i, 2] = z
        b[i, 2] = z
    return a,b

# TODO (coding)
# - Create a function to generate the D matrix.
# - Rename `solve` to `cosmo`.
# - Create a functioni similar to `solve` but for IEF and call it `ief`.
# - Make routine to compute the solvation energy  E = 0.5 * \sum_i = V_i asc_i

#print("Weights = ",scheme.weights) # Shows the weights for the sphere
#print('Lebedev Degree:', scheme.degree)
##########################
#Paramaters:
Z = np.array([-1.8, -1.6, -1.4, -1.2, -1, -0.8, -0.6, -0.4, -0.2, 0, 0.2, 0.4, 0.6, 0.8, 1, 1.2, 1.4, 1.6, 1.8])

k = 1.0694
R = 1 # Radius of the sphere
q = np.array([+1]) # Charges
xyz_sphere = [0, 0, 0] # Center position of the sphere
xyz_charge = np.array([[0, 0, 0]]) # Position of the charge
epsilon = 80
#Lebedev Quadrature
###
# q1 = 2 q2 = -1
# (0, 0, -1)  (0, 0, +1)
#
# TODO (testing)
# Test 1 (done for different quadratures)
# q=1 position=(0,0,z) z=0 ..... R-0.2
# Collect the total ASC
# Compute the solvation energy
# Plot ASC(z) E(z)
#
# Test 2 
# Similar to test 1 but now with a dipole
# 
# q=1 (a,0,z) q=-1 (-a,0,z)  a=0.1 z=0 .... R-0.2
# GAUSS THEOREM q_ASC = q * (1 - eps) / eps



scheme = quadpy.sphere.lebedev_065 ()    # Which precision of Lebedev
grid_areals = scheme.integrate(lambda x: 1, xyz_sphere, R)
normals = scheme.points

points = sphere_size(scheme.points, R)
points = sphere_position(points, xyz_sphere)
#print("Points on the sphere:")
#print(points)
print("Number of quadrature points:", points.shape[0])
total_areal_sphere = total_areal(grid_areals)
#print ('Total areal:', total_areal_sphere)
#print("Weights = ",grid_areals)
w_i = scheme.weights * total_areal_sphere
#print(w_i.shape[0])
S_matrix = S(points,w_i)





"""

xyz_1, xyz_2 = dipole_testing(Z)
test = testing(Z)
"""
"""
for r, xyz in enumerate(xyz_1):
    new_xyz = np.vstack([xyz, xyz_2[r]])
    #print(new_xyz)
    r_i = potential(points, new_xyz, q)
    Sigma = COSMO(S_matrix, r_i)
    Energy = Solvation_energy(Sigma, r_i)
    #print(f"XYZ (charge) = {new_xyz} Total charge = {np.sum(Sigma)} Energy = {Energy} ")
"""
"""
for xyz in new_xyz:
    #print(xyz)
    r_i = potential(points, new_xyz, q)
    Sigma = COSMO(S_matrix, r_i)
    Energy = Solvation_energy(Sigma, r_i)
    print(f"XYZ (charge) = {xyz} Total charge = {np.sum(Sigma)} Energy = {Energy} ")
"""    


r_i = potential(points, xyz_charge, q)
#print("r_i", r_i)
Sigma = COSMO(S_matrix, r_i)
#print(f"Sigma = {Sigma}")
print(f"Total charge = {np.sum(Sigma)}")


Energy = Solvation_energy(Sigma, r_i)
print(Energy)

A_matrix = A_matrix(w_i)

D_matrix = D_matrix(points, w_i, R, normals)


IEF = IEF(S_matrix, D_matrix, A_matrix, r_i, epsilon)

print(f"Total charge IEF = {(np.sum(IEF))}")


def gauss_theorem(q, eps):
    gauss = np.sum(q) * (1- eps) / eps
    return gauss

Gauss = gauss_theorem(q,epsilon)
print(Gauss)
#print(IEF)


#print(S_matrix)
#print(D_matrix)
#np.testing.assert_allclose(np.sum(Sigma), - q, rtol=1e-2)
