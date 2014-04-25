#!/usr/bin/env python
import os
import sys
import numpy as np
import matplotlib.pyplot as plt

dataDir = '../'

# setup for RMS log data
rmsFile = 'ErrorLog.dat'
rmsFile = os.path.join(dataDir,rmsFile)
rmsdata = np.loadtxt(rmsFile)

# Plot for RMS log
nIter = rmsdata[:,0]
RMSerr = rmsdata[:,1]

MinX = min(nIter)
MaxX = max(nIter)
MinY = 1.0E-05#min(RMSerr)
MaxY = 1.0#max(RMSerr)

p = plt.plot(nIter,RMSerr, 'r-', label='RMS error')
plt.setp(p, linewidth='1.0')
plt.axis([MinX,MaxX, MinY, MaxY])
plt.xscale('linear')
plt.yscale('log')
plt.xlabel('Number of iteration', fontsize=22)
plt.ylabel('RMS error', fontsize=22)
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

pltFile = 'RMSlog.png'
fig = plt.gcf()
fig.set_size_inches(8,5)
plt.tight_layout()
plt.savefig(pltFile, format='png')
plt.close()

print "%s DONE!!" % (pltFile)

