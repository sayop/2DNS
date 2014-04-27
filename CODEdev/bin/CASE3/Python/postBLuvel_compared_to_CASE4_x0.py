#!/usr/bin/env python
import os
import sys
import numpy as np
import matplotlib.pyplot as plt

dataDir1 = '../'
dataDir2 = '/home/sayop/data/Devel/GitHub.Clones/2DNS/CODEdev/bin/CASE4'

# setup for DATA file
File1 = 'BoundaryLayer_x0.dat'
File1 = os.path.join(dataDir1,File1)
data1 = np.loadtxt(File1)

File2 = 'BoundaryLayer_x0.dat'
File2 = os.path.join(dataDir2,File2)
data2 = np.loadtxt(File2)

# Plot for RMS log
y1 = data1[:,0]
uvel1 = data1[:,1]

y2 = data2[:,0]
uvel2 = data2[:,1]

MinX = -0.4
MaxX = 1
MinY = 0.0#min(uvel)
MaxY = 0.017#max(uvel)

p = plt.plot(uvel1,y1, 'r-', label='CASE 3 (Adiabatic wall)')
plt.setp(p, linewidth='2.0')
p = plt.plot(uvel2,y2, 'b--', label='CASE 4 (Isothermal wall)')
plt.setp(p, linewidth='2.0')
plt.axis([MinX,MaxX, MinY, MaxY])
plt.xscale('linear')
plt.yscale('linear')
plt.xlabel('Axial velocity [m/s]', fontsize=22)
plt.ylabel('y at x = 0', fontsize=22)
plt.grid(True)
ax = plt.gca()
xlabels = plt.getp(ax, 'xticklabels')
ylabels = plt.getp(ax, 'yticklabels')
plt.setp(xlabels, fontsize=18)
plt.setp(ylabels, fontsize=18)
plt.legend(
          loc='upper right',
          borderpad=0.25,
          handletextpad=0.25,
          borderaxespad=0.25,
          labelspacing=0.0,
          handlelength=2.0,
          numpoints=1)
legendText = plt.gca().get_legend().get_texts()
plt.setp(legendText, fontsize=18)
legend = plt.gca().get_legend()
legend.draw_frame(False)

pltFile = 'BL_uvel_CASE3and4_x0.png'
fig = plt.gcf()
fig.set_size_inches(8,5)
plt.tight_layout()
plt.savefig(pltFile, format='png')
plt.close()

print "%s DONE!!" % (pltFile)

