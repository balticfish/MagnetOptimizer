import numpy as np
from parameters import *

def Rotate(phi):

	R = [[np.cos(phi), 0, np.sin(phi), 0], [0, np.cos(phi), 0,  np.sin(phi)], [-np.sin(phi), 0, np.cos(phi), 0], [0, -np.sin(phi), 0, np.cos(phi)]]

	return np.asarray(R)

def QskewFD(q, L): #First magnet focusing in xx', defocusing in yy'

	Kappa = q / L

	Qs11 = np.cos(L*np.sqrt(Kappa))
	Qs12 = np.sin(L*np.sqrt(Kappa))/np.sqrt(Kappa)
	Qs13 = 0
	Qs14 = 0

	Qs21 = -np.sqrt(Kappa)* np.sin(L*np.sqrt(Kappa))
	Qs22 = np.cos(L*np.sqrt(Kappa))
	Qs23 = 0
	Qs24 = 0

	Qs31 = 0
	Qs32 = 0
	Qs33 = np.cosh(L*np.sqrt(Kappa))
	Qs34 = np.sinh(L*np.sqrt(Kappa))/np.sqrt(Kappa)

	Qs41 = 0
	Qs42 = 0
	Qs43 = np.sqrt(Kappa)*np.sinh(L*np.sqrt(Kappa))
	Qs44 = np.cosh(L*np.sqrt(Kappa))

	Qs = [[Qs11,Qs12,Qs13,Qs14],[Qs21,Qs22,Qs23,Qs24],[Qs31,Qs32,Qs33,Qs34],[Qs41,Qs42,Qs43,Qs44]]

	Qskew = Rotate(-np.pi/4.0).dot(np.asarray(Qs)).dot(Rotate(np.pi/4.0))
	
	return np.asarray(Qskew)


def QskewDF(q, L): #First magnet focusing in yy', defocusing in xx'

	Kappa = q / L

	Qs11 = np.cosh(L*np.sqrt(Kappa))
	Qs12 = np.sinh(L*np.sqrt(Kappa))/np.sqrt(Kappa)
	Qs13 = 0
	Qs14 = 0

	Qs21 = np.sqrt(Kappa)*np.sinh(L*np.sqrt(Kappa))
	Qs22 = np.cosh(L*np.sqrt(Kappa))
	Qs23 = 0
	Qs24 = 0

	Qs31 = 0
	Qs32 = 0
	Qs33 = np.cos(L*np.sqrt(Kappa))
	Qs34 = np.sin(L*np.sqrt(Kappa))/np.sqrt(Kappa)

	Qs41 = 0
	Qs42 = 0
	Qs43 = -np.sqrt(Kappa)* np.sin(L*np.sqrt(Kappa))
	Qs44 = np.cos(L*np.sqrt(Kappa))

	Qs = [[Qs11,Qs12,Qs13,Qs14],[Qs21,Qs22,Qs23,Qs24],[Qs31,Qs32,Qs33,Qs34],[Qs41,Qs42,Qs43,Qs44]]

	Qskew = Rotate(-np.pi/4.0).dot(np.asarray(Qs)).dot(Rotate(np.pi/4.0))
	
	return np.asarray(Qskew)


def Drift(L):

	D = [[1, L, 0, 0], [0, 1, 0, 0], [0, 0, 1, L], [0, 0, 0, 1]]
	
	return np.asarray(D)


def RTFB(q1,q2,q3,d2,d3):

	if (q2<0):
		RTFB = QskewFD(np.abs(q3),L).dot(Drift(d3)).dot(QskewDF(np.abs(q2),L)).dot(Drift(d2)).dot(QskewFD(np.abs(q1),L))
	else:
		RTFB = QskewDF(np.abs(q3),L).dot(Drift(d3)).dot(QskewFD(np.abs(q2),L)).dot(Drift(d2)).dot(QskewDF(np.abs(q1),L))

	return RTFB
