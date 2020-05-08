# -*- coding: utf-8 -*-
"""
Spyder Editor

This is a temporary script file.
"""

import numpy as np
import quadpy

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

def S(points, weights):
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

def potential(points, position_charge, charge):
    # Create array for the potential
    V = np.zeros(points.shape[0]) #matches number of lebedev points
    # Fill array with potential
    for i, x in enumerate(points):
        for A, R in enumerate(position_charge):
            V[i] += charge[A] / np.linalg.norm(x-R)
    return V

def solve(S, V):
    # Invert S
    S_matrix_inverted = np.linalg.inv(S)
    # Matrix-vector multiply S^-1 V
    sigma = -np.dot(S_matrix_inverted, V)
    return sigma

def IEF(S, D, A, V):
    pass

# TODO (coding)
# - Create a function to generate the D matrix.
# - Rename `solve` to `cosmo`.
# - Create a functioni similar to `solve` but for IEF and call it `ief`.
# - Make routine to compute the solvation energy  E = 0.5 * \sum_i = V_i asc_i

#print("Weights = ",scheme.weights) # Shows the weights for the sphere
#print('Lebedev Degree:', scheme.degree)
##########################
#Paramaters:

k = 1.0694
R = 2 # Radius of the sphere
q = np.array([+2, -1]) # Charges
xyz_sphere = [0, 0, 0] # Center position of the sphere
xyz_charge = np.array([[0, 0, -1], [0,0, +1]]) # Position of the charge
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

scheme = quadpy.sphere.lebedev_131 ()    # Which precision of Lebedev
grid_areals = scheme.integrate(lambda x: 1, xyz_sphere, R)
points = sphere_size(scheme.points, R)
points = sphere_position(points, xyz_sphere)
print("Points on the sphere:")
print(points)
print(points.shape[0])
total_areal_sphere = total_areal(grid_areals)
print ('Total areal:', total_areal_sphere)
print("Weights = ",grid_areals)
w_i = scheme.weights * total_areal_sphere
print(w_i.shape[0])
S_matrix = S(points,w_i)
r_i = potential(points, xyz_charge, q)
print("r_i", r_i)
Sigma = solve(S_matrix, r_i)
print(f"Sigma = {Sigma}")
print(f"Total charge = {np.sum(Sigma)}")

#np.testing.assert_allclose(np.sum(Sigma), - q, rtol=1e-2)
