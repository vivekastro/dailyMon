import numpy as np
from astropy.io import fits
import scipy as sp
import glob
import os
import matplotlib.pyplot as plt
import datetime
from astropy.time import Time
from pyspherematch import *
from scipy.ndimage.filters import gaussian_filter1d,uniform_filter1d

text_font = {'fontname':'Arial', 'size':'6'}

color_list = ['red','blue','green','brown','cyan','magenta','gold','orange','yellow']
def init_plotting():
    plt.rcParams['font.size'] = 15
    plt.rcParams['font.family'] = 'Times New Roman'
    plt.rcParams['axes.labelsize'] = plt.rcParams['font.size']
    plt.rcParams['axes.titlesize'] = 1.5*plt.rcParams['font.size']
    plt.rcParams['legend.fontsize'] = plt.rcParams['font.size']
    plt.rcParams['xtick.labelsize'] = plt.rcParams['font.size']
    plt.rcParams['ytick.labelsize'] = plt.rcParams['font.size']
    plt.rcParams['xtick.major.size'] = 3
    plt.rcParams['xtick.minor.size'] = 3
    plt.rcParams['xtick.major.width'] = 1
    plt.rcParams['xtick.minor.width'] = 1
    plt.rcParams['ytick.major.size'] = 3
    plt.rcParams['ytick.minor.size'] = 3
    plt.rcParams['ytick.major.width'] = 1
    plt.rcParams['ytick.minor.width'] = 1
    plt.rcParams['legend.frameon'] = False
    plt.rcParams['legend.loc'] = 'center left'
    plt.rcParams['axes.linewidth'] = 1

    plt.gca().spines['right'].set_color('none')
    plt.gca().spines['top'].set_color('none')
    plt.gca().xaxis.set_ticks_position('bottom')
    plt.gca().yaxis.set_ticks_position('left')

def currentMJD():
    '''  Returns the current MJD    '''
    return Time(str(datetime.date.today())).mjd


def download_spectra(plate,mjd,fiber,dirname='.'):
        FITS_FILENAME = 'spec-%(plate)04i-%(mjd)05i-%(fiber)04i.fits'
        SDSS_URL = ('https://data.sdss.org/sas/dr14/eboss/spectro/redux/v5_10_0/spectra/%(plate)04i/'
            'spec-%(plate)04i-%(mjd)05i-%(fiber)04i.fits')
        print SDSS_URL % dict(plate=plate,mjd=mjd,fiber=fiber)
        download_url = 'wget '+SDSS_URL % dict(plate=plate,mjd=mjd,fiber=fiber)
        print download_url
        os.system(download_url)
        mv_cmd='mv '+FITS_FILENAME % dict(plate=plate,mjd=mjd,fiber=fiber) + ' '+dirname+'/.'
        print mv_cmd
        os.system(mv_cmd)

def readPlateList(platelistfile='platelist.fits'):
    platelistdir = '.'  # Update with correct path to $BOSS_SPECTRO_REDUX
    platelist= fits.open(os.path.join(platelistdir,platelistfile))[1].data
    doneplatelist = platelist[np.where(platelist['STATUS1D'] == 'Done')[0]]
    mjd = sorted(platelist['MJD'])[-20]
    print 'Working on MJD :',mjd
    plates = doneplatelist[np.where(doneplatelist['MJD'] == mjd)[0]]
    qsoplates = plates[np.where(plates['PROGRAMNAME'] == 'eboss')[0]]['PLATE']
    print 'Plates observed on MJD: {}  are {}'.format(mjd,qsoplates)
    return mjd,qsoplates


def readSpAllfile(plateid,spallfile='spAll-v5_10_7.fits'):
    spAlldir='.'   # Update with correct path to $BOSS_SPECTRO_REDUX
    spAll = fits.open(os.path.join(spAlldir,spallfile))[1].data
    print spAll.columns.names
    platematch = spAll[np.where(spAll['PLATE'] == plateid)[0]]
    balindex = np.where(platematch['EBOSS_TARGET2'] & 2**25)[0]
    print balindex
    return  platematch['FIBERID'][balindex],platematch['RA'][balindex],platematch['DEC'][balindex],platematch['Z'][balindex]

def sphereMatch(ra,dec,tol=0.00001):
    dupplate=[];dupmjd=[];dupfiber=[]
    dr14q = fits.open('DR14Q_v4_4.fits')[1].data
    a,b,ds= spherematch(ra,dec,dr14q['RA'],dr14q['DEC'],tol=tol)
    print a,b,dr14q['RA'][b],ra,ds,dr14q['PLATE'][b]
    for i in range(len(b)):
        dupplate.append(dr14q['PLATE'][b[i]])
        dupmjd.append(dr14q['MJD'][b[i]])
        dupfiber.append(dr14q['FIBERID'][b[i]])
    return dupplate,dupmjd,dupfiber


def plotSpectra(pplate,pmjd,pfiber,splate,smjd,sfiber,z):
    
    # Reading the line lists from speccy
    linelist = np.genfromtxt('linelist_speccy.txt',usecols=(0,1,2),dtype=('|S10',float,'|S5'),names=True)
    bossspectroredux = '.' # Update with correct path to $BOSS_SPECTRO_REDUX
    spPlatefile = os.path.join(bossspectroredux,str(pplate),'spPlate-{0}-{1}.fits'.format(pplate,pmjd))
    spPlate = fits.open(spPlatefile)[0].data
    spheader = fits.open(spPlatefile)[0].header
    pwave = 10**(spheader['CRVAL1']+np.arange(spheader['NAXIS1'])*spheader['CD1_1'])
    pflux = spPlate[pfiber - 1]
    pinvar = fits.open(spPlatefile)[0].data[pfiber - 1]
    plabel = '{0:04d}-{1:5d}-{2:04d}'.format(pplate,pmjd,pfiber)
    fig,ax = plt.subplots(figsize=(20,10))
    init_plotting()
    ax.plot(pwave,gaussian_filter1d(pflux,5),color='black',alpha=0.6,label=plabel)
    ax.plot(pwave,1.0/np.sqrt(pinvar),color='black',alpha=0.1)

    ax.set_xlabel(r'Observed Wavelength $\AA$')
    ax.set_ylabel(r'Flux 10$^{-17}$ ergs s$^{-1}$ cm$^{-2}$ $\AA^{-1}$')
    
    ax.set_ylim(np.median(pflux)-3*np.std(pflux),np.median(pflux)+3*np.std(pflux))
    ax.set_xlim(3500,10400)
    ax2 = ax.twiny()
    ax1Xs = ax.get_xticks()
    ax2Xs = []
    for X in ax1Xs:
        ax2Xs.append(X /(1.0+z))
    ax2Xsm = ["%.1f" % zzz for zzz in ax2Xs]
    ylim = ax.get_ylim()
    xlim = ax.get_xlim()
               
    ax.text(xlim[0]+0.017*(xlim[1]-xlim[0]),ylim[0]+0.95*(ylim[1]-ylim[0]),plabel)
    ax.text(xlim[0]+0.017*(xlim[1]-xlim[0]),ylim[0]+0.9*(ylim[1]-ylim[0]),'Z_PIP:  '+str(z))
    #ax.text(xlim[0]+0.6*(xlim[1]-xlim[0]),ylim[0]+0.75*(ylim[1]-ylim[0]),'SN:'+str(sn2))
    #ax.text(xlim[0]+0.3*(xlim[1]-xlim[0]),ylim[0]+0.8*(ylim[1]-ylim[0]),'RA:'+str(Right_asc)+' DEC:'+str(Decli))
    ax2.set_xticks(ax1Xs)
    ax2.set_xbound(ax.get_xbound())
    ax2.set_xticklabels(ax2Xsm)
    ax2.set_xlabel(r'Rest Wavelength $\AA$')
    obslambda = linelist['lambda']*(1.+z)
    xx = np.where((obslambda > xlim[0]) & (obslambda < xlim[1]))[0]
    plotlambda = obslambda[xx]
    plotname = linelist['Name'][xx]
    plota_e = linelist['a_e'][xx]
    #print plotlambda
    for k in range(len(plotlambda)):
        if plota_e[k].strip() == 'Abs.' :
            ax.axvline(x=plotlambda[k], color='lawngreen', linestyle=':')
            ax.text(plotlambda[k],ylim[0]+0.05*(ylim[1]-ylim[0]),plotname[k],color='Orange',ha='center',rotation=90,**text_font)
        else :
            ax.axvline(x=plotlambda[k], color='lightblue', linestyle=':')
            ax.text(plotlambda[k],ylim[0]+0.15*(ylim[1]-ylim[0]),plotname[k],color='Brown',ha='center',rotation=90,**text_font)

    for k in range(len(splate)):
        download_spectra(splate[k],smjd[k],sfiber[k],'DR14_Spectra')
        dupfile = os.path.join('DR14_Spectra','spec-{0:04d}-{1:5d}-{2:04d}.fits'.format(splate[k],smjd[k],sfiber[k]))
        print dupfile
        slabel = '{0:04d}-{1:5d}-{2:04d}'.format(splate[k],smjd[k],sfiber[k])
        sdata = fits.open(dupfile)[1].data
        sflux = sdata.flux
        swave = 10**sdata.loglam
        serr = 1.0/np.sqrt(sdata.ivar)
        ax.plot(swave,gaussian_filter1d(sflux,5),color=color_list[k],ls='--',alpha=0.9,label=slabel)
        ax.plot(swave,serr,color=color_list[k],ls='--',alpha=0.1)
    ax.legend(loc=1)
    fig.tight_layout()
    destination_dir = os.path.join('plots',str(pmjd))
    dir_cmd = 'mkdir {}'.format(destination_dir)
    if not os.path.exists(destination_dir):
        os.system(dir_cmd)
    savefilename = os.path.join(destination_dir,'Status_{0}_{1}_{2}.jpeg'.format(pplate,pmjd,pfiber))
    fig.savefig(savefilename)
    plt.show()

    
    


currentmjd = int(currentMJD())

mjd,qsoplates = readPlateList()




for plate in qsoplates:
    fiber,ra,dec,z = readSpAllfile(plate)
    for i in range(len(fiber)):
        print plate,mjd,fiber[i],ra[i],dec[i],z[i]
        dupplate,dupmjd,dupfiber=sphereMatch(ra[i],dec[i])
        print 'PLATE-MJD-FIBER: {0}-{1}-{2} has {3} entries in DR14 '.format(plate,mjd,fiber[i],len(dupplate))
        plotSpectra(plate,mjd,fiber[i],dupplate,dupmjd,dupfiber,z[i])
