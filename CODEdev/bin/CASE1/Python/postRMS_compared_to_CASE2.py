#!/usr/bin/env python
import os
import sys
import numpy as np
import matplotlib.pyplot as plt

dataDir1 = '../'
dataDir2 = '/home/sayop/data/Devel/GitHub.Clones/2DNS/CODEdev/bin/CASE2'

# setup for RMS log data
rmsFile1 = 'ErrorLog.dat'
rmsFile1 = os.path.join(dataDir1,rmsFile1)
rmsdata1 = np.loadtxt(rmsFile1)

rmsFile2 = 'ErrorLog.dat'
rmsFile2 = os.path.join(dataDir2,rmsFile2)
rmsdata2 = np.loadtxt(rmsFile2)

# Plot for RMS log
nIter1 = rmsdata1[:,0]
RMSerr1 = rmsdata1[:,1]

nIter2 = rmsdata2[:,0]
RMSerr2 = rmsdata2[:,1]

MinX = min(nIter1)
MaxX = max(nIter1)
MinY = 1.0E-05#min(RMSerr)
MaxY = 1.0#max(RMSerr)

p = plt.plot(nIter1,RMSerr1, 'r-', label='CASE1 (1st order accurate)')
plt.setp(p, linewidth='2.0')
p = plt.plot(nIter2,RMSerr2, 'b--', label='CASE2 (2nd order accurate)')
plt.setp(p, linewidth='2.0')
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

pltFile = 'RMSlog_CASE1and2.png'
fig = plt.gcf()
fig.set_size_inches(8,5)
plt.tight_layout()
plt.savefig(pltFile, format='png')
plt.close()

print "%s DONE!!" % (pltFile)

