import numpy as np
from parameters import *

def FASTQK ( Iq, gamma ): #FAST Quadrupole calibration

	return (10.135*40*Iq)/(1.8205*gamma*mc)


def FASTQIq ( K, gamma ): #FAST Quadrupole calibration

	return (K * 1.18205 * gamma * mc) / (10.135 * 40.0) 


#Analytical solutions for three-quadrupole round-to-flat beam transformer
#assuming thin lens approximation


def calc_q1(s11,s12,s21,s22):

	return np.sqrt((s12 - d02*(s11 + d0T*s21) + d0T*s22)/(d02*d0T*s12))

def calc_q2 (q1,s11,s12,s21,s22):

	return (-s11 + (d02 + d03)*(q1 - s21))/(d03*(-1 + d02*q1*s11))

def calc_q3 (q1,q2,s11,s12,s21,s22):

	return (d02*q2 - d02* q1* q2* s12 - s22)/(-d02 - d03 + d02* q1* s12 + d03* q1* s12 +  d03* q2* s12 + d02* d03* q2* s22)

def calc_K (q): #magnet K-value
	return q/L

def calc_g (q,gamma):  #magnet gradient
	K = q /L
	return K*gamma*mc/299.8


