import numpy as np
import sys
import math
import os

def readsigma4(fname):

	disti=np.loadtxt(fname)
	dist=disti.copy()

	#convert from impact format
	dist[:,1]=disti[:,1]/disti[:,5]
	dist[:,3]=disti[:,3]/disti[:,5]

	#calculate beam second moments
	sigma4=np.zeros((4,4))
	for i in range(0,4):
	   for j in range(0,4):
	      sigma4[i,j]=np.mean(dist[:,i]*dist[:,j])

	return sigma4

def readgamma (fname):

	disti=np.loadtxt(fname)
	betagamma=np.sqrt(np.mean(disti[:,1]**2+disti[:,3]**2+disti[:,5]**2))
	gamma=np.sqrt(1.+betagamma**2)
	beta=np.sqrt(1.-1./gamma**2)

	return gamma


def correlation_matrix(sigma4):

	#Correlation matrix  C = <YX>.<XX>^-1
	XX=np.array([[sigma4[0,0],sigma4[0,1]],[sigma4[1,0],sigma4[1,1]]])
	XY=np.array([[sigma4[0,2],sigma4[0,3]],[sigma4[1,2],sigma4[1,3]]])
	YX=np.array([[sigma4[2,0],sigma4[2,1]],[sigma4[3,0],sigma4[3,1]]])
	YY=np.array([[sigma4[2,2],sigma4[2,3]],[sigma4[3,2],sigma4[3,3]]])
	C=np.dot(YX,np.linalg.inv(XX))

	return C
