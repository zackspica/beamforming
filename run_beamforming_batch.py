import os, glob
import numpy as np

for d in xrange(6,25,1):
    os.system('python main_beamforming.py %s'%d)


