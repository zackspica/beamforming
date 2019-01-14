import glob, os
from pyrocko import model
from obspy import read
import numpy as np

g = glob.glob('/Volumes/MrWhite/DAS_Stanford/DAS_mseed/281/*')
sta = model.load_stations('stations.txt')
bin = np.zeros( (len(sta),len(read(g[0])[0].data)) ) 
for i, s in enumerate(sta):
    tr = read('/Volumes/MrWhite/DAS_Stanford/DAS_mseed/281/DS.%s.281.mseed'%s.station)[0].data
    bin[i]=tr
print (bin)
np.save('DAS.281.npy',bin)


