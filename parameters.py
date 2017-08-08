#FAST quadrupole length
L = 0.167

#Q106 - Q107 
d2 = 0.2 - L

#Q107-Q111
d3 = 1.58 - L

#Total length of the transformer
dT = d2 + d3

#Distances for thin-lens approximation
d02 = d2 + L
d03 = d3 + L
d0T = d02 + d03

mc = 0.511
