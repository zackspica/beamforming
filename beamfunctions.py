#!/usr/bin/env python 

"""
Plane wave beamforming for the vertical component.
""""""""""""""""""""""""""""""""""""""""""""""""""
Based on the original code of Yannik Behr
Last update Zack Spica 15/07/17
"""

import os
import sys
import glob
from obspy.core import UTCDateTime
import obspy.signal
from obspy.geodetics import gps2dist_azimuth
from obspy.signal.invsim import cosine_taper
from obspy.signal.filter import bandpass, integer_decimation
import scipy.io as sio
import progressbar as pg
from smooth import smooth
import numpy as np
from params import get_params
from beamplot import *

DEBUG = False


def arr_resp(freqs, nstat, idx, slowness, zetax, theta, sta_origin_x, sta_origin_y,
        new=True, matfile=None, src=True, src_param=[(90, 3000)]):
    """
    Compute the array response.
    """
    StackedBeams = []
    if new:
        beam = zeros((zetax.shape[0], slowness.size, len(freqs)))
        if src_param is not None:
            for couple in src_param:
                theta1, c1 = couple
                zeta_x = -cos(theta1 * pi / 180.)
                zeta_y = -sin(theta1 * pi / 180.)
                zeta_src = zeta_x * sta_origin_x + zeta_y * sta_origin_y
                for ww in idx:
                    FF = freqs[ww]
                    for cc in np.arange(slowness.shape[1]):
                        omega = 2 * pi * FF
                        velocity = 1. / slowness[0][cc] * 1000
                        e_steer = exp(-1j * zetax * omega / velocity).T
                        if src:
                            e_src = exp(1j * zeta_src * omega / c1).T
                            Y = multiply(ones((zetax.shape[1], 1), 'complex128'), atleast_2d(e_src).T)
                            R = dot(Y, conjugate(Y).T)
                        else:
                            R = ones((zetax.shape[1], zetax.shape[1]))
                        beam[:, cc, ww] += 1. / (nstat * nstat) ** 2 * diag(abs(asarray(dot(conjugate(e_steer.T), dot(R, e_steer)))) ** 2)
        if matfile is not None:
            for ind in idx:
                stackedbeam = squeeze(beam[:, :, ind])
                StackedBeams.append(stackedbeam)
                filename = matfile[:-4]+'_%s'%(np.round(freqs[ind],2))+'.mat'
                sio.savemat(filename, {'beam':stackedbeam, 'theta':theta, 'slowness':slowness, 
                    'freqs':freqs, 'f':np.round(freqs[ind],2)})

    else:
        print (' >> Reading %s and passing beamforming response computation...'%matfile)
        for ind in idx:
            filename = matfile[:-4]+'_%s'%(np.round(freqs[ind],2))+'.mat'
            StackedBeams.append(sio.loadmat(filename)['beam'])
    
    return StackedBeams


def calc_steer(slats, slons):
    """
    Compute the steering vector

    :type slats:
    :param slats:
    :type slons:
    :param slons:
    """
    theta = arange(0, 362, 2)
    theta = theta.reshape((theta.size, 1))
    sta_origin_dist = array([])
    sta_origin_bearing = array([])
    meanlat = slats.mean()
    meanlon = slons.mean()
    for lat, lon in zip(slats, slons):
        #get dist between each sta and mean xy of array
        dist, az, baz = gps2dist_azimuth(meanlat, meanlon, lat, lon)
        sta_origin_dist = append(sta_origin_dist, dist)
        sta_origin_bearing = append(sta_origin_bearing, az)
    # ramener dist sur le cadran en x et y (dist*cos(az)) ou dist*sin(az)
    sta_origin_x = sta_origin_dist * cos(sta_origin_bearing * pi / 180.)
    sta_origin_y = sta_origin_dist * sin(sta_origin_bearing * pi / 180.)
    # On le calcule pr plein de theta different car apres on cherche quel est le meilleur.  
    zeta_x = -cos(theta * pi / 180.)
    zeta_y = -sin(theta * pi / 180.)
    # dot product betwen zeta and x: d = -x.cos(theta)-y.sin(theta)   
    zetax = zeta_x * sta_origin_x + zeta_y * sta_origin_y
    return zetax, theta, sta_origin_x, sta_origin_y


def get_new_binfile(binfile, stations_model, receivers_model):
    newbin = np.zeros( (len(receivers_model), len(binfile[0]) ))
    for irec, rec in enumerate(receivers_model):
        for ista, sta in enumerate(stations_model):
            if sta.station == rec.station:
                newbin[irec] = binfile[ista]
    return newbin


def prep_beam(binfile, matfile, stations_model, ndays=1,  nhours=1, fsamp=10., threshold_std=0.5, onebit=False,
              tempfilter=False, specwhite=True, timenorm=False, fact=10, new=False,
              fftpower=7, freq_int=(0.02, 0.4)):
    """
    Prepare the raw data before the beamforming. This comprises bandpass
    filtering, cutting the data into smaller chunks, removing the mean and
    down-weighting strong transient signals from for example earthquakes. One
    can chose between several methods to remove transients: 1-bit normalization,
    time domain normalization, which computes a smoothed traces of the absolute
    amplitude and down-weights the original trace by this smoothed trace, and a
    threshold based method that clips the trace at a predefined factor of the
    traces standard deviation. Note that currently no instrument correction is
    applied so the input traces either have to be already corrected for the
    instrument response, they need to be filtered within a frequency band for
    which all instrument responses are flat or they need to be all from the
    same instruments.

    :param files: Day long vertical component SAC files.
    :param matfile: Name of the file to which the pre-processed traces are
                    written. This saves the time of repeatedly running the
                    pre-processing for beamformer runs with different parameters.
                    The output file is a *.mat file that can also be read with
                    Matlab.
    :param statons_model:   
    :param nhours: Input data is cut into chunks of length of nhours.
    :param fsamp:  Sampling frequency of the input data.
    :param threshold_std: Clipping factor; values greater than
                          treshold_std * std(trace) are set to threshold_std.
    :param onebit: Turn on/off 1-bit normalization
    :param tempfilter: Turn on/off threshold based clipping.
    :param specwhite: Turn on/off spectral whitening (only retain spectral phase
                      and set spectral amplitude to one.
    :param timenomr: Turn on/off time domain normalization.
    :param fact: Decimation factor.
    :param new: If set to false it will try to load all return values from
                matfile. If true it will compute all return values from scratch.
    :param fftpower: Length of data chunks cut before the FFT.
    :param freq_int: Frequency interval for the band-pass filter.

    :return fseis: Pre-processed traces in the frequency domain.
    :return freqs: Frequency array corresponding to fseis
    :return slats: Station latitudes.
    :return slons: Station longitudes.
    :return dt: Sampling interval after decimation.
    :return seissmall: Pre-processed traces in the time domain.
    """
    if new:
        print (' >> Prepare raw data before beamforming...')
        ntimes = int(ndays) * int(round(24 / nhours))
        step = int( nhours * 3600 * fsamp / fact )
        stations = []
        slons = array([])
        slats = array([])
        nfiles = np.shape(binfile)[0]
        seisband = zeros((nfiles, ntimes, step))
        freqs = fftfreq(2 ** fftpower, 1. / (fsamp / fact))
        for i, (t, s) in enumerate(zip(binfile,stations_model)):
            data = t
            print (i, s.station)
            slons = append(slons, s.lon)
            slats = append(slats, s.lat)
            stations.append(s.station)
            data -= data.mean()
            data = bandpass(data, freqmin=freq_int[0], freqmax=freq_int[1], df=fsamp,
                   corners=4, zerophase=True)
            if fact != 1:
                data = integer_decimation(data, fact)
            npts = len(data) 
            df = fsamp / fact  
            dt = 1./df 
            seis0 = zeros(ndays * 24 * 3600 * int(df))
            istart = int(round(((UTCDateTime(0).hour * 60 + UTCDateTime(0).minute) * 60\
                          + UTCDateTime(0).second) * df))
            if timenorm:
                smoothdata = smooth(abs(data), window_len=257, window='flat')
                data /= smoothdata
            if np.isnan(data).any():
                data = zeros(ndays * 24 * 3600 * int(df))
            try:
                seis0[istart:(istart + npts)] = data
            except (ValueError, e):
                print ('Problem with %s'%s.station )
                raise ValueError

            seis0 -= seis0.mean()
            # iterate over nhours
            for j in np.arange(ntimes):
                ilow = j * step
                iup = (j + 1) * step
                seisband[i, j, :] = seis0[ilow:iup]
                seisband[i, j, :] -= seisband[i, j, :].mean()
                if onebit:
                    seisband[i, j, :] = sign(seisband[i, j, :])
            if tempfilter:
                sigmas = seisband[i, :, :].std(axis=1, ddof=1)
                sgm = ma.masked_equal(array(sigmas), 0.).compressed()
                sigma = sqrt(sum(sgm ** 2) / sgm.size)
                threshold = threshold_std * sigma
                seisband[i] = where(abs(seisband[i]) > threshold, threshold * sign(seisband[i]), seisband[i])
                seisband[i] = apply_along_axis(lambda e: e - e.mean(), 1, seisband[i])

        ismall = 2 **fftpower
        ipick = arange(ismall)
        taper = cosine_taper(len(ipick))
        n = ndays * nhours * 3600 * df
        nsub = int(np.floor(n / ismall))  # Number of time pieces -20 mins long each
        seissmall = zeros((nfiles, ntimes, nsub, len(ipick)))
        for ii in np.arange(nfiles):
            for jj in np.arange(ntimes):
                for kk in np.arange(nsub):
                    seissmall[ii, jj, kk, :] = seisband[ii, jj, kk * ismall + ipick] * taper
        fseis = fft(seissmall, n=2**fftpower, axis=3)
        if np.isnan(fseis).any():
            print ("NaN found")
            return
        ind = np.where((freqs > freq_int[0]) & (freqs < freq_int[1]))[0]
        fseis = fseis[:, :, :, ind]
        if specwhite:
            fseis = exp(angle(fseis) * 1j)

        sio.savemat(matfile, {'fseis':fseis, 'slats':slats, 'slons':slons, 'dt':dt, 
                              'freqs':freqs[ind]})
        return fseis, freqs[ind], slats, slons, dt
    else:
        print (' >> Reading %s and passing preprocess...'%matfile)
        a = sio.loadmat(matfile)
        fseis = a['fseis']
        freqs = np.squeeze(a['freqs'])
        slats = np.squeeze(a['slats'])
        slons = np.squeeze(a['slons'])
        dt = a['dt'][0][0]
        return fseis, freqs, slats, slons, dt 


def beamforming(fseis, freqs, slowness, theta, zetax, nsources, idx, new=True, matfile=None):
    """
    Compute the beam in the frequency domain. Main code for BF. 
    """
    StackedBeams = []
    if new:
        print (' >> Compute beamforming...')
        nstat, ntimes, nsub, nfft = fseis.shape
        beam = zeros((nsources, slowness.size, ntimes, nfft))
        if not DEBUG:
            widgets = ['vertical beamforming: ', pg.Percentage(), ' ', pg.Bar('#'),
                       ' ', pg.ETA()]
            pbar = pg.ProgressBar(widgets=widgets, maxval=len(idx) * ntimes * nsub).start()
            count = 0
        for ww in idx:
            FF = freqs[ww]
            for tt in np.arange(ntimes):
                for TT in np.arange(nsub):
                    if not DEBUG:
                        count += 1
                        pbar.update(count)
                    Y = asmatrix(squeeze(fseis[:, tt, TT, ww]))# une matrix avec toutes les sta a une freq donnee
                    YT = Y.T.copy()#transposee
                    R = dot(YT, conjugate(Y))# produit vecotriel de YT * conj(Y)
                    for cc in np.arange(slowness.shape[1]):
                        omega = 2 * pi * FF
                        velocity = 1. / slowness[0][cc] * 1000
                        # this is just the steering vector in the freq domaine tau=r.s
                        e_steer = exp(-1j * zetax * omega / velocity)# phase shift vect
                        e_steerT = e_steer.T.copy()
                        beam[:, cc, tt, ww] += 1. / (nstat * nstat) * \
                        diag(abs(asarray(dot(conjugate(e_steer), dot(R, e_steerT)))) ** 2) / nsub

        if not DEBUG:
            pbar.finish()
        if matfile is not None:
            for ind in idx:
                stackedbeam = squeeze(beam[:, :, :, ind])
                stackedbeam = stackedbeam[:,:,:].mean(axis=2)
                StackedBeams.append(stackedbeam)
                filename = matfile[:-4]+'_%s'%(np.round(freqs[ind],2))+'.mat'
                sio.savemat(filename, {'beam':stackedbeam, 'theta':theta, 'slowness':slowness, 
                                       'freqs':freqs, 'f':np.round(freqs[ind],2) })
    else:
        print (' >> Reading %s and passing beamforming computation...'%matfile)
        for ind in idx:
            filename = matfile[:-4]+'_%s'%(np.round(freqs[ind],2))+'.mat'
            StackedBeams.append(sio.loadmat(filename)['beam'])

    return StackedBeams


def stackbeamfiles(list_beamfiles):
    for iff, f in enumerate(list_beamfiles):
        b = sio.loadmat(f)['beam']
        if iff == 0:
            superstack = b
        else:   
            superstack += b
    superstack /= iff+1

    return superstack

if __name__=='__main__':
    p = get_params()
    import scipy.io as sio
    binfile = np.load(p['datdir']+'daily.NL.2015.330.20sps.Z.npy')
    get_new_binfile(binfile, p['stations'], p['receivers'])


    freq = [0.5]
    for f in freq:
        g = glob.glob('out/beam_DA.2016.*_%s.mat'%f)
    superstack = stackbeamfiles(g)
    a = sio.loadmat(g[0])
    theta = a['theta']
    freqs = a['freqs']
    slowness = a['slowness']
    sio.savemat('superstackbeam_%s.mat'%freq[0], {'beam':superstack, 'theta':theta, 
        'freqs':freqs, 'slowness':slowness, 'f':freq[0]})

    polar_plot_from_file('superstackbeam_%s.mat'%freq[0], title='')
    

#EOF





