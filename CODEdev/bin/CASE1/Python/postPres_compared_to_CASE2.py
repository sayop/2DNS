#!/usr/bin/env python
import os
import sys
import numpy as np
import matplotlib.pyplot as plt

dataDir1 = '../'
dataDir2 = '/home/sayop/data/Devel/GitHub.Clones/2DNS/CODEdev/bin/CASE2'

# setup for DATA file
File1 = 'WallData.dat'
File1 = os.path.join(dataDir1,File1)
data1 = np.loadtxt(File1)

File2 = 'WallData.dat'
File2 = os.path.join(dataDir2,File2)
data2 = np.loadtxt(File2)

FileExac = 'ExacP.dat'
FileExac = os.path.join(dataDir1,FileExac)
dataExac = np.loadtxt(FileExac)

# Plot for RMS log
x1 = data1[:,0]
pres1 = data1[:,4]

x2 = data2[:,0]
pres2 = data2[:,4]

xExac = dataExac[:,0]
pExac = dataExac[:,1]

MinX = min(x1)
MaxX = max(x1)
MinY = 0.0#min(pres)
MaxY = 0.4#max(pres)

p = plt.plot(x1,pres1, 'r-', label='CASE1 (1st order accurate)')
plt.setp(p, linewidth='2.0')
p = plt.plot(x2,pres2, 'b--', label='CASE2 (2nd order accurate)')
plt.setp(p, linewidth='2.0')
p = plt.plot(xExac,pExac, 'k-', label='Exact Solution')
plt.setp(p, linewidth='2.0')
plt.axis([MinX,MaxX, MinY, MaxY])
plt.xscale('linear')
plt.yscale('linear')
plt.xlabel('x', fontsize=22)
plt.ylabel('Pressure [Non-Dim]', fontsize=22)
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

pltFile = 'Pressure_CASE1and2.png'
fig = plt.gcf()
fig.set_size_inches(8,5)
plt.tight_layout()
plt.savefig(pltFile, format='png')
plt.close()

print "%s DONE!!" % (pltFile)

