#!/usr/bin/env python
from matplotlib.colors import LinearSegmentedColormap
import numpy as np
from matplotlib.colors import Normalize, LogNorm
from matplotlib.colorbar import ColorbarBase
import matplotlib.pyplot as plt
from obspy.imaging.cm import obspy_sequential
from pylab import *
import scipy.io as sio
from params import get_params

"""
Plot function for plane wave beamforming.
""""""""""""""""""""""""""""""""""""""""""""""""""
Based on the original code of Yannik Behr
Last update Zack Spica 15/07/17
"""

def mycolormap():
    """return Nori's colormap"""
    colors = np.loadtxt('color.txt', delimiter=',')
    n_bin = 320  # Discretizes the interpolation into bins
    cmap_name = 'my_list'
    # Create the colormap
    cm = LinearSegmentedColormap.from_list(cmap_name, colors, N=n_bin)
    # Fewer bins will result in "coarser" colomap interpolation
    return cm


def crown(values, thetas=None, radii=None, ax=None, fraction=0.3,cmap=cm.get_cmap('jet_r'), norm=norm, **kwargs):
    """
    Plot extra crown around the polar axis
    e.g., I use it to plot the projection of the maximum and the mean
    """
    values = np.atleast_2d(values)
    if thetas is None:
        thetas = np.linspace(0, 2*np.pi, values.shape[1]).reshape(1, -1)
    if radii is None:
        radii = np.linspace(0, 1, values.shape[0] + 1).reshape(-1, 1)
    if ax is None:
        fig, ax = plt.subplots(1, 1, subplot_kw={'polar':True})
    mesh = ax.pcolormesh(thetas, radii, values, cmap=cmap, norm=norm, **kwargs)
    radrange = radii.ptp()
    ax.set_rlim(radrange * (1 - 1. / fraction), radrange)
    ax.set_axis_off()

    return mesh


def polar_plot(StackedBeams, theta, freqs, slowness, freq2analyze, plotcrown=True, fout=None):
    """
    Plot the beam over a disc:  polar plot
    Function used direclty with the outputs of the main code.
    To plot just a file.mat use polar_plot_file(filename)
    """
    idx = [argmin(abs(freqs - freq)) for freq in freq2analyze]
    theta = theta[:, 0]
    slowness = slowness[0, :]
    for tre, ind in zip(StackedBeams, idx):
        f = np.round(freqs[ind],2)
        print ('>>> Plotting %s Hz' % f)

        tre /= tre.max()
        fig = figure(figsize=(6, 6))
        cax = fig.add_axes([0.85, 0.2, 0.03, 0.5])
        #ax2 = fig.add_axes([0.125, 0.125, 0.65, 0.65], polar=True)
        #ax3 = fig.add_axes([0.1, 0.1, 0.70, 0.70], polar=True)
        ax = fig.add_axes([0.15, 0.15, 0.60, 0.6], polar=True)

        #cmap = obspy_sequential
        #cmap = cm.get_cmap('hot_r')
        cmap = mycolormap()
        phi = (theta[::-1] + 90.) * pi / 180.
        if polar:
            ax.contourf(phi, slowness, tre.T, 100, cmap=cmap, antialiased=True, 
                        linstyles='dotted',zprder=0)
            ax.contour(phi, slowness, tre.T, 100, cmap=cmap, zorder=1)
        #x, y = np.meshgrid(phi, r)
        #projmax = [m.max() for m in tre]
        #projmean = [m.mean() for m in tre]

        #norm = LogNorm(vmin=0.1, vmax=1.)
        #ax.pcolormesh(x, y, tre.T, cmap=cmap, norm=norm)  # ,vmin=0.,vmax=1.)
        #if plotcrown:
        #    crown(projmax, thetas=phi, ax=ax2, fraction=0.05, cmap=cmap, norm=norm)
        #    crown(projmean, thetas=phi, ax=ax3, fraction=0.05, cmap=cmap, norm=norm)
        ax.set_thetagrids([0, 45., 90., 135., 180., 225., 270., 315.],
                          labels=['90', '45', '0', '315', '270', '225', '180', '135'], frac=1.27, fontsize=9)
        ax.set_rgrids([0.5,1.,2, 3], labels=['2','1','0.5','0.3'], color='w',zorder=1000)
        #circle = Circle((0.0, 0.0), 0.2, transform=ax.transData._b, color="k", linestyle='--', fill=False)
        #ax.add_artist(circle)

        ax.grid(True)
        ax.set_title("%s Hz  " % (f), y=1.1)
        #cb = ColorbarBase(cax, cmap=cmap, norm=norm, ticks=[0.2 + 0.001, 1.])
        #cb.ax.set_yticklabels(['<0.25', '1'])
        #cb.set_label('Normalized power')
        if fout is not None:
            fo = '%s_%s.png' %(fout, f)
            savefig(fo, dpi=300, format='png')
        #plt.show()

def polar_plot_resp(StackedBeams, theta, freqs, slowness, freq2analyze, polar=True, fout=None):
    """
    Plot the beam for the array response
    """
    idx = [argmin(abs(freqs - freq)) for freq in freq2analyze]
    theta = theta[:, 0]
    slowness = slowness[0, :]
    for ind, tre in zip(idx, StackedBeams):
        f = np.round(freqs[ind],2)
        print ('>>> Plotting %s Hz' %f)
        fig = figure(figsize=(6, 6))
        if polar:
            ax = fig.add_subplot(1, 1, 1, projection='polar')
        else:
            ax = fig.add_subplot(1, 1, 1)
        
        cmap = mycolormap()#cm.get_cmap('jet')
        
        if polar:
            ax.contourf((theta[::-1] + 90.) * pi / 180., slowness, tre.T,
                    100, cmap=cmap, antialiased=True, linstyles='dotted',zorder=0)
            ax.contour((theta[::-1] + 90.) * pi / 180., slowness, tre.T,
                    100, cmap=cmap, zorder=1)
        else:
            ax.contourf(theta, slowness, tre.T, 100, cmap=cmap, antialiased=True,
                    linstyles='dotted')
            ax.contour(theta, slowness, tre.T, 100, cmap=cmap)

        if polar:
            ax.set_thetagrids([0, 45., 90., 135., 180., 225., 270., 315.],
            labels=['90', '45', '0', '315', '270', '225', '180', '135'])
            #ax.set_rgrids([0.1, 0.2, 0.3, 0.4, 0.5], labels=['0.1', '0.2', '0.3', '0.4', '0.5'], color='r')
            ax.set_rgrids([0.5,1.,2, 3], labels=['2','1','0.5','0.3'], color='w',zorder=1000)
            #ax.set_rmax(0.5)
        else:
            ax.set_xlabel('Azimuth [degrees]')
            ax.set_ylabel('Slowness [s/km]')
        ax.grid(True)
        #ax.set_title("%s %ds" %('Array response at',int(round((1./freqs[ind])))), y=1.075)
        ax.set_title("%s %s Hz" %('Array response at',f), y=1.075)
        if fout is not None:
            fo = '%s_%s.png' %(fout, f)
            savefig(fo, dpi=300, format='png')
        #if show:
        plt.show()



def polar_plot_from_file(filename, title='', polar=True, fout=None):
    """
    Plot the beam for the array response
    """
    a = sio.loadmat(filename)
    tre = a['beam']
    freqs = a['freqs']
    slowness = a['slowness']
    theta = a['theta']
    f = a['f'][0][0]
    theta = theta[:, 0]
    slowness = slowness[0, :]
    print ('>>> Plotting %s Hz' %f)
    fig = figure(figsize=(6, 6))
    if polar:
        ax = fig.add_subplot(1, 1, 1, projection='polar')
    else:
        ax = fig.add_subplot(1, 1, 1)
    cmap = mycolormap()#cm.get_cmap('jet')
    if polar:
        ax.contourf((theta[::-1] + 90.) * pi / 180., slowness, tre.T,
            100, cmap=cmap, antialiased=True,linstyles='dotted',zorder=0)
        ax.contour((theta[::-1] + 90.) * pi / 180., slowness, tre.T,
            100, cmap=cmap,zorder=1)
    else:
        ax.contourf(theta, slowness, tre.T, 100, cmap=cmap, antialiased=True,
               linstyles='dotted')
        ax.contour(theta, slowness, tre.T, 100, cmap=cmap)

    if polar:
        ax.set_thetagrids([0, 45., 90., 135., 180., 225., 270., 315.],
                    labels=['90', '45', '0', '315', '270', '225', '180', '135'])
            #ax.set_rgrids([0.1, 0.2, 0.3, 0.4, 0.5], labels=['0.1', '0.2', '0.3', '0.4', '0.5'], color='r')
        ax.set_rgrids([0.5,1.,2, 3], labels=['2','1','0.5','0.3'], color='w',zorder=1000)
            #ax.set_rmax(0.5)
    else:
        ax.set_xlabel('Azimuth [degrees]')
        ax.set_ylabel('Slowness [s/km]')
    ax.grid(True)
        #ax.set_title("%s %ds" %('Array response at',int(round((1./freqs[ind])))), y=1.075)
    ax.set_title("%s %s Hz" %(title, f), y=1.075)
    if fout is not None:
        fo = '%s_%s.png' %(fout, f)
        savefig(fo, dpi=300, format='png')
    plt.show()

if __name__=='__main__':
    polar_plot_resp_from_file('out/array_response_DA.2016_3.0.mat')

