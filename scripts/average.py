#!/usr/bin/env python

import sys
import numpy as np
import io

bin_size=50
samples = 5000
if len(sys.argv) != 3:
    print("Usage: %s FILENAME OUTPUT"%sys.argv[0])
    sys.exit(1)

data = np.genfromtxt(sys.argv[1])

# cutoff necessary because i messed up
data[data[:,2]>1e3] = np.NaN
data = data.reshape(-1,samples,3)

data = data.reshape(data.shape[0],-1,bin_size,3)

data = np.nanmean(data,axis=2)

# find out whether there is still relaxation going on
t = np.arange(samples/bin_size)
st = np.std(t)
mt = np.mean(t)
tcorr = (np.nanmean(data[:,:,2]*t,axis=1)-mt*np.nanmean(data[:,:,2],axis=1))/st/np.nanstd(data[:,:,2],axis=1)

VTE = np.nanmean(data,axis=1)
sigE = np.nanstd(data[:,:,2],ddof=1,axis=1)/(samples/bin_size)

# Volume, Temperature, Energy, sigEnergy, tcorrelation
c = np.hstack([VTE,sigE[:,None],tcorr[:,None]])
i = np.lexsort((c[:,0],c[:,1]),0)

np.savetxt(sys.argv[2],c)
