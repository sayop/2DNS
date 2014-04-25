#!/usr/bin/env python
import os
import sys
import numpy as np
import matplotlib.pyplot as plt

dataDir = '../'

# setup for DATA file
File = 'WallData.dat'
File = os.path.join(dataDir,File)
data = np.loadtxt(File)

FileExac = 'ExacMach.dat'
FileExac = os.path.join(dataDir,FileExac)
dataExac = np.loadtxt(FileExac)

# Plot for RMS log
x = data[:,0]
mach = data[:,5]

xExac = dataExac[:,0]
mExac = dataExac[:,1]

MinX = min(x)
MaxX = max(x)
MinY = 1.5#min(mach)
MaxY = 2.5#max(mach)

p = plt.plot(x,mach, 'r-', label='Numerical Solution')
p = plt.plot(xExac,mExac, 'k-', label='Exact Solution')
plt.setp(p, linewidth='1.0')
plt.axis([MinX,MaxX, MinY, MaxY])
plt.xscale('linear')
plt.yscale('linear')
plt.xlabel('x', fontsize=22)
plt.ylabel('Pressure', fontsize=22)
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

pltFile = 'Mach.png'
fig = plt.gcf()
fig.set_size_inches(8,5)
plt.tight_layout()
plt.savefig(pltFile, format='png')
plt.close()

print "%s DONE!!" % (pltFile)

