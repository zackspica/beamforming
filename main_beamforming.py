#!/usr/bin/env python
from beamfunctions import *
from params import *
import os, glob, sys

"""
Main code for plane wave beamforming for the vertical component.
""""""""""""""""""""""""""""""""""""""""""""""""""
Based on the original code of Yannik Behr
Last update Zack Spica 15/07/17
"""



def main(p, arg_list):
    """
    main function to calculate the beamforming
    over one binfile.npy
    All params are in params.py
    The julian day for the binfile selection has to be
    parsed through an argument
    """
    print (">> Calculate the beamforming")
    day = '%03d' % int(arg_list[0])
    datdir = p['datdir']
    outdir = p['outdir']
    if datdir.endswith('/'):
        datdir = datdir[0:-1]
    if outdir.endswith('/'):
        outdir = outdir[0:-1]
    try:
        binfile = np.load(glob.glob(datdir + '/*%s*' % day)[0])
    except:
        print ('No file found for day %s.' % day)
        sys.exit()
    
    binfile = get_new_binfile(binfile, p['stations'], p['receivers']) 

    if len(binfile[0]) < p['nsta']:
        print ("Matrix and nsta do not match: %s vs %s" % (file[0], nsta))
        return
    tempname = p['nameout'] + '.' + day
    matfile1 = "%s/%s_%s.mat" % (outdir, 'prep_beam', tempname)
    fseis, freqs, slats, slons, dt = prep_beam(binfile, matfile1, p['receivers'], p['ndays'],
                        nhours=p['nhours'], fsamp=int(p['df']), new=p['nprep'],
                        fact=p['div_fact'], threshold_std=p['threshold_std'], onebit=p['onebit'],
                        tempfilter=p['tempfilter'], specwhite=p['specwhite'], timenorm=p['timenorm'],
                        fftpower=p['fftpower'], freq_int=p['freq_int'])
    zetax, theta, sta_origin_x, sta_origin_y = calc_steer(slats, slons)
    matfile2 = "%s/%s_%s.mat" % (outdir, 'beam', tempname)
    nsources, ntimes, nsub, nfft = fseis.shape
    indices = [argmin(abs(freqs - freq)) for freq in p['freq2analyze']]
    StackedBeams = beamforming(fseis, freqs, p['slowness'], theta, zetax, theta.size, indices,
                               new=p['nbeam'], matfile=matfile2)
    if p['plot']:
        # fout = p['saveplot']
        if p['saveplot']:
            fout = matfile2.replace('.mat', '_vertical')
        polar_plot(StackedBeams, theta, freqs, p['slowness'],  p['freq2analyze'], fout=fout)
        if p['show']:
            show()


def main_response(p, arg_list):
    """
    Calculate the theoretical array response
    """
    print ('>> Calculate the theoretical array response')
    day = '%03d' %int(arg_list[0])
    datdir = p['datdir']
    outdir = p['outdir']
    if datdir.endswith('/'):
        datdir = datdir[0:-1]
    if outdir.endswith('/'):
        outdir = outdir[0:-1]
    try:
        binfile = np.load(glob.glob(datdir + '/*%s*' % day)[0])
    except:
        print ('No file found for day %s.' % day)
        sys.exit()
    
    binfile = get_new_binfile(binfile, p['stations'], p['receivers']) 
    
    if len(binfile[0]) < p['nsta']:
        print ("Matrix and nsta do not match: %s vs %s" % (file[0], nsta))
        return
    tempname = p['nameout'] + '.' + day
    matfile1 = "%s/%s_%s.mat" % (outdir, 'prep_beam', tempname)
    """ gros probleme wiht sample_f"""
    fseis, freqs, slats, slons, dt = prep_beam(binfile, matfile1, p['receivers'], p['ndays'],
                                               nhours=1, fsamp=int(p['df']), new=p['nprep'], fact=p['div_fact'],
                                               threshold_std=p['threshold_std'], onebit=False,
                                               tempfilter=False, specwhite=False, timenorm=False,
                                               fftpower=p['fftpower'], freq_int=p['freq_int'])
    zetax, theta, sta_origin_x, sta_origin_y = calc_steer(slats, slons)

    matfile2 = "%s/%s_%s.mat" % (outdir, 'array_response', tempname)

    nsources, ntimes, nsub, nfft = fseis.shape
    indices = [argmin(abs(freqs - freq)) for freq in p['freq2analyze']]
    
    beam = arr_resp(freqs, p['nsta'], indices, p['slowness'], zetax, theta,
    sta_origin_x, sta_origin_y, new=p['nbeam'], matfile=matfile2, src=True, src_param=p['src_param'])
    
    if p['plot']:
        if p['saveplot']:
            fout = matfile2.replace('.mat', 'resp_vertical')
        polar_plot_resp(beam, theta, freqs, p['slowness'],
                        p['freq2analyze'], polar=True, fout=fout)
        if p['show']:
            show()

if __name__ == '__main__':
    """
    One important thing is that we need only one file (Z) per station
    Use a binfile as for correlation functions, same input (npy)
    usage

    main_beamforming.py arg=#julianday
    """

    p = get_params()
    main(p, sys.argv[1:])
    #main_response(p, sys.argv[1:])
