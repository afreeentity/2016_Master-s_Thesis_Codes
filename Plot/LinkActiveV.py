import math as mt

def link_active(c,tr):

	Phiz_M=0.0 
	alpha=1.e0 
	beta=1.e5
	Phin_M = c
	eta=1.0

	xp_M = mt.exp(-1*eta*alpha*(Phin_M-Phiz_M))
	gq_M = 1./(1.+beta*xp_M)



	te = False
	if gq_M>tr: te = True
	return te
