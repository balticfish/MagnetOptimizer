import numpy as np
import scipy.optimize as optimize
from elements import *
from sigma4D import *
from magnets import *

#Distribution upstream of the round-to-flat transformer
fname = 'fort.1'

print "--------------------------------"
gamma = readgamma(fname)
print "Beam gamma factor: ",gamma
print "--------------------------------"

print "--------------------------------"
print "Initial 4D second moments matrix"
sigma4 = readsigma4(fname)
print sigma4
print "--------------------------------"

print "--------------------------------"
print "Correlation matrix"
C = correlation_matrix(sigma4)
print C
print "--------------------------------"

print "--------------------------------"
#Obtain a first guess via analytical solution
q1 = calc_q1(C[0,0],C[0,1],C[1,0],C[1,1])
q2 = calc_q2(q1, C[0,0],C[0,1],C[1,0],C[1,1])
q3 = calc_q3(q1,q2,C[0,0],C[0,1],C[1,0],C[1,1])

print "Magnets q = 1 / f analytical approximation (q1,q2,q3):"
print q1,q2,q3
print "--------------------------------"

print "--------------------------------"
#Calculate beam sigma matrix after the transformer
sigmaF = RTFB(q1,q2,q3,d2,d3).dot(sigma4).dot(RTFB(q1,q2,q3,d2,d3).T)
print "Non-optimized 4D second moments matrix"
print sigmaF
print "--------------------------------"


def obj_function(qs):  #building the function for minimization

	q1 = qs[0]
	q2 = qs[1]
	q3 = qs[2]
	sigmaF = RTFB(q1,q2,q3,d2,d3).dot(sigma4).dot(RTFB(q1,q2,q3,d2,d3).T)
	res = 1.0e6*np.sqrt(sigmaF[0,2]**2 + sigma4[1,2]**2 + sigmaF[0,3]**2 + sigmaF[1,3]**2)

	return res

print "Objective value before optimization: ", obj_function([q1,q2,q3])

result = optimize.minimize(obj_function, [q1,q2,q3],method='Nelder-Mead',tol=1e-5)

print "Objective value after optimization: ", obj_function(result.x)

#updating the values of quadrupoles
q1 = result.x[0]
q2 = result.x[1]
q3 = result.x[2]

print "--------------------------------"
#Calculate beam sigma matrix after the transformer with optimized quads
sigmaF = RTFB(q1,q2,q3,d2,d3).dot(sigma4).dot(RTFB(q1,q2,q3,d2,d3).T)
print "Optimized 4D second moments matrix"
print sigmaF
print "--------------------------------"

print "--------------------------------"
print "Gradients in T/m (Impact-T format): "
g1 = calc_g(q1, gamma)
g2 = calc_g(q2, gamma)
g3 = calc_g(q3, gamma)
print g1, g2, g3
print "--------------------------------"

#Compute quadrupoles K-values
#print "Magnets K :"
K1 = calc_K(q1)
K2 = calc_K(q2)
K3 = calc_K(q3)
#print K1, K2, K3
print "--------------------------------"
print "FAST quadrupole currents (in Amps):"
Iq1 = FASTQIq(K1,gamma)
Iq2 = FASTQIq(K2,gamma)
Iq3 = FASTQIq(K3,gamma)
print "Q106      Q107    Q111"
print ("%.3f   %.3f   %.3f"% (Iq1, Iq2, Iq3))
print "--------------------------------"


