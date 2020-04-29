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
R = 1 # Radius of the sphere
q = 1 # Charge
xyz_sphere = [1, 1, 1] # Center position of the sphere
xyz_charge = [1, 1, 1] # Position of the charge

#Lebedev Quadrature
###
scheme = quadpy.sphere.lebedev_003a ()    # Which order of Lebedev


val = scheme.integrate(lambda x: 1, xyz_sphere, R)


points = scheme.points * R

for i, x in enumerate(points):   # Changes the center position of the sphere from [0,0,0] to xyz_sphere
    x[0] = x[0]+xyz_sphere[0]
    x[1] = x[1]+xyz_sphere[1]
    x[2] = x[2]+xyz_sphere[2]


print("Points on the sphere:")
print(points)



def total_areal(val):
    result = sum(val)
    return result

if __name__ == "__main__":
    print ('Total areal: {}'.format(total_areal(val)))

 
    print("Weights = ",val)

###
x_i = np.array([x[0] for x in points]) # Array of all x coordinates
y_i = np.array([y[1] for y in points]) # Array of all y coordinates
z_i = np.array([z[2] for z in points]) # Array of all z coordinates



w_i = scheme.weights * total_areal(val)


###
#Cosmo
################################################
s_ii = k * np.sqrt( 4 * np.pi / (R ** 2 * w_i))  
print("s_ii", s_ii)

list1 = []

for i in points:
    for j in points:
        list1.append( 1 / np.linalg.norm(i-j))
        

s_ij0 = np.asarray(list1).reshape(int(len(list1)**0.5),int(len(list1)**0.5))  # Reshapes into a i x j matrix 


for i in s_ii:
    s_ij = np.where(s_ij0==np.inf,i,s_ij0)  # Replaces the 0 in S_ij matrix with S_ii 
#################################################


r_i = np.zeros(len(points)) #matches number of lebedev points


#print(scheme.points[0, :])


for i, x in enumerate(points):
    r_i[i] = q / np.linalg.norm(x-xyz_charge)
    
    print("r_i:", r_i[i])
################################################


print("S_ij:")
print(s_ij)

s_ij_inverted = np.linalg.inv(s_ij)

#print(s_ij_inverted)

sigma = np.dot(-s_ij_inverted, r_i)
print("Sigma:", sigma)

print(np.sum(sigma))


#np.testing.assert_allclose(np.sum(sigma), - q)




