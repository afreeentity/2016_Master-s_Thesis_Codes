import matplotlib.mlab as mlab
import matplotlib.pyplot as plt
import numpy as np
import os
import plfit
import powerlaw

#my_path = os.path.abspath(__file__)


theta=[]
with open("ttt") as f:
    for line in f:
        theta.append(float(line))
#print(theta)

plt.xlabel('Degree', fontsize=14)
plt.ylabel('Number of Nodes', fontsize=14)
plt.title('Node Degree Distribution Histogram')
plt.grid(True)
plt.yscale('log')
plt.xscale('log')
plt.hist(theta, 1400, range=[0. , 1400.])
plt.savefig('histogram.eps', facecolor='w', edgecolor='w', format='eps')
plt.show()

a=np.sort(theta)
#x = np.arange(len(a))
print(a)
plt.xlabel('Nodes', fontsize=14)
plt.ylabel('Degrees', fontsize=14)
plt.title('Data_Visualization')
plt.grid(True)
plt.plot(a)
plt.savefig('Plot.eps', facecolor='w', edgecolor='w', format='eps')
plt.show()

thefile = open('test_m.txt', 'w')
for item in a:
  thefile.write("%f\n" % item)
 
pl=plfit.plfit(theta,usefortran=False)
from pylab import *
figure(1)
clf()
pl.plotpdf()
savefig("fit.png")
print "Nodes     : n:%10i degree: mean,std,max: %8.2f,%8.2f,%8.2f degreemin: %8.2f alpha: %8.2f (%8.2f) ntail: %10i p: %5.2f" % (pl.data.shape[0], pl.data.mean(), pl.data.std(), pl.data.max(), pl._xmin, pl._alpha, pl._alphaerr, pl._ngtx, pl._ks_prob)
