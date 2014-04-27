#!/usr/bin/env python
import os
import sys
import numpy as np
import matplotlib.pyplot as plt

dataDir1 = '../'
dataDir2 = '/home/sayop/data/Devel/GitHub.Clones/2DNS/CODEdev/bin/CASE4'

# setup for DATA file
File1 = 'BoundaryLayer_x-025.dat'
File1 = os.path.join(dataDir1,File1)
data1 = np.loadtxt(File1)

# setup for DATA file
File2 = 'BoundaryLayer_x-025.dat'
File2 = os.path.join(dataDir2,File2)
data2 = np.loadtxt(File2)

# Plot for RMS log
y1 = data1[:,0]
temp1 = data1[:,2]

y2 = data2[:,0]
temp2 = data2[:,2]

MinX = 300.0
MaxX = 550.0
MinY = 0.0#min(temp)
MaxY = 0.009#max(temp)

p = plt.plot(temp1,y1, 'r-', label='CASE3 (Adiabatic wall)')
plt.setp(p, linewidth='2.0')
p = plt.plot(temp2,y2, 'b--', label='CASE4 (Isothermal wall)')
plt.setp(p, linewidth='2.0')
plt.axis([MinX,MaxX, MinY, MaxY])
plt.xscale('linear')
plt.yscale('linear')
plt.xlabel('Temperature [K]', fontsize=22)
plt.ylabel('y at x = -0.25', fontsize=22)
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

pltFile = 'BL_temp_CASE3and4_x-025.png'
fig = plt.gcf()
fig.set_size_inches(8,5)
plt.tight_layout()
plt.savefig(pltFile, format='png')
plt.close()

print "%s DONE!!" % (pltFile)

