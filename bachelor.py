# -*- coding: utf-8 -*-
"""
Spyder Editor

This is a temporary script file.
"""

import numpy as np
import quadpy


# The lebedev sphere

#scheme.show()

#print("Weights = ",scheme.weights) # Shows the weights for the sphere
#print('Lebedev Degree:', scheme.degree)
##########################
#Paramaters:




k = 1.0694
R = 2 # Radius of the sphere
q = 1 # Charge
xyz_sphere = [1, 1, 1] # Center position of the sphere
xyz_charge = [0, 0, 0] # Position of the charge

#Lebedev Quadrature
###
scheme = quadpy.sphere.lebedev_003a ()    # Which order of Lebedev


val = scheme.integrate(lambda x: 1, xyz_sphere, R)

def sphere_size(points, radius):
    #Changes the radius of the sphere
    new_radius = points * radius
    return new_radius

points = sphere_size(scheme.points, R)

def sphere_position(points, new_position):
    # Changes the center position of the sphere from [0,0,0] to xyz_sphere
    for i,x in enumerate(points):
        x[0] = x[0]+new_position[0]
        x[1] = x[1]+new_position[1]
        x[2] = x[2]+new_position[2]
    return points
    

points = sphere_position(points, xyz_sphere)



print("Points on the sphere:")
print(points)
print(points.shape[0])



def total_areal(val):
    result = sum(val)
    return result

if __name__ == "__main__":
    print ('Total areal: {}'.format(total_areal(val)))

 
    print("Weights = ",val)

###

w_i = scheme.weights * total_areal(val)
print(w_i.shape[0])

##
#Cosmo
################################################


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

#s_ij -> S_matrix
S_matrix = S(points,w_i)


################################################
def potential(points, position_charge, charge):
    # Create array for the potential
    V = np.zeros(points.shape[0]) #matches number of lebedev points
    # Fill array with potential
    for i, x in enumerate(points):
        for A, R in enumerate(position_charge):
            V[i] += charge[A] / np.linalg.norm(x-R)
    return V

r_i = potential(points, [xyz_charge], [q])
print("r_i", potential(points, [xyz_charge], [q]))



# s_ij_inverted -> S_matrix_inverted
def solve(S, V):
    # Invert S
    S_matrix_inverted = np.linalg.inv(S)
    # Matrix-vector multiply S^-1 V
    sigma = np.dot(S_matrix_inverted, V)
    return sigma

Sigma = solve(S_matrix, r_i)

print(Sigma)
print(np.sum(Sigma))
#np.testing.assert_allclose(np.sum(Sigma), - q, rtol=1e-2)
