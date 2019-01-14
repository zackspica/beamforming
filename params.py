#!/usr/bin/env python

from pyrocko import model

from pylab import *
import os


def get_params():
    # Get Configuration
    params = {}
    params['stations']      = model.load_stations('stations.txt')#all sta
    params['receivers']     = model.load_stations('receivers.txt')#sta to analyze
    params['datdir']        = '/Users/tzompantli/DAS/Stanford/BF/bins/'
    params['outdir']        = 'out/'
    params['nameout']       = 'DS.2018'
    params['df']            = 50.
    params['freq2analyze']  = [1., 3., 5., 8.,10,12, 15, 20]
    params['slowmin']       = 0.05
    params['slowmax']       = 4.55
    params['slowstep']      = 0.05
    params['div_fact']      = 1 # div fact for downsampling
    params['nprep']         = True
    params['nbeam']         = True
    params['plot']          = True
    params['saveplot']      = True#'ttt'
    params['show']          = True
   # params['saveresults']   = True
    params['nsta']          = len(model.load_stations('receivers.txt')) 
    
    params['fftpower']      = 12 #precision
    params['freq_int']      = (0.5,21)

    params['nhours']        = 1 # Input data is cut into chunks of length nhours.
    params['ndays']         = 1 #number of days to compute as a bunch
    params['threshold_std'] = 0.5
    params['onebit']        = True
    params['tempfilter']    = False#True
    params['specwhite']     = True
    params['timenorm']      = False#True

    params['src_param']     = [(90,800000)] #location and vel of source for array resp


    sl                      = arange(params['slowmin'],params['slowmax'],params['slowstep'])
    params['slowness']      = sl.reshape(1,sl.size)
   
    try: os.makedirs(params['outdir'])
    except:pass
    
    return params 


        



if __name__=='__main__':
    params = get_params()
    print (params)

#EDF
