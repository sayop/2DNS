#!/usr/bin/env python
import os
import sys
import numpy as np
import matplotlib.pyplot as plt

dataDir = '../'

# setup for DATA file
File = 'BoundaryLayer.dat'
File = os.path.join(dataDir,File)
data = np.loadtxt(File)

# Plot for RMS log
y = data[:,0]
temp = data[:,2]

MinX = 300.0
MaxX = 550.0
MinY = 0.0#min(temp)
MaxY = 0.02#max(temp)

p = plt.plot(temp,y, 'r-', label='Numerical solution')
plt.setp(p, linewidth='1.0')
plt.axis([MinX,MaxX, MinY, MaxY])
plt.xscale('linear')
plt.yscale('linear')
plt.xlabel('Temperature [K]', fontsize=22)
plt.ylabel('y', fontsize=22)
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

pltFile = 'BL_temp.png'
fig = plt.gcf()
fig.set_size_inches(8,5)
plt.tight_layout()
plt.savefig(pltFile, format='png')
plt.close()

print "%s DONE!!" % (pltFile)

