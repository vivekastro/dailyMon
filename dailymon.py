#!/usr/bin/python

import numpy as np
from astropy.io import fits
import scipy as sp
import glob
import os
import matplotlib.pyplot as plt
from astropy.time import Time
from datetime import date
from pyspherematch import *
from scipy.ndimage.filters import gaussian_filter1d,uniform_filter1d
import calendar

plt.switch_backend('agg')

text_font = {'fontname':'Arial', 'size':'9'}
color_list = ['red','blue','green','brown','cyan','magenta','gold','orange','yellow']

def setRedux():
    try:
        #BossSpectroRedux = '/uufs/chpc.utah.edu/common/home/sdss/ebosswork/eboss/spectro/redux'#os.environ['BOSS_SPECTRO_REDUX']
        BossSpectroRedux = os.environ['BOSS_SPECTRO_REDUX']
    except:
        BossSpectroRedux = '.'
    if BossSpectroRedux == '.':
        spAllfile= sorted(glob.glob('spAll-v5_1*.fits'))[-1]
    else:    
        version = os.path.basename(sorted(glob.glob(os.path.join(BossSpectroRedux,'v5_10*')))[-1])
        BossSpectroRedux = os.path.join(BossSpectroRedux,version)
        spAllfile = 'spAll-'+version+'.fits'
    return BossSpectroRedux,spAllfile

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
    return Time(str(date.today())).mjd

def convertMJD(mjd):
    mydate=Time(mjd,format='mjd')
    mydate.format='fits'
    tdate= date(int(mydate.value[0:4]),int(mydate.value[5:7]),int(mydate.value[8:10]))
    month = calendar.month_name[tdate.month]  
    weekday = calendar.day_name[tdate.weekday()] 
    ret_string = '{0} {1}-{2}-{3} '.format(weekday,int(mydate.value[8:10]),month,int(mydate.value[0:4]))
    return ret_string

def download_spectra(plate, mjd, fiber, dirname='.'):
        FITS_FILENAME = 'spec-%(plate)04i-%(mjd)05i-%(fiber)04i.fits'
        SDSS_URL = ('https://data.sdss.org/sas/dr14/eboss/spectro/redux/v5_10_0/spectra/%(plate)04i/'
            'spec-%(plate)04i-%(mjd)05i-%(fiber)04i.fits')
       # print SDSS_URL % dict(plate=plate,mjd=mjd,fiber=fiber)
        download_url = 'wget -q '+SDSS_URL % dict(plate=plate,mjd=mjd,fiber=fiber)
        #print download_url
        os.system(download_url)
        mv_cmd='mv '+FITS_FILENAME % dict(plate=plate,mjd=mjd,fiber=fiber) + ' '+dirname+'/.'
        #print mv_cmd
        os.system(mv_cmd)

def outPlates(mjd,platelist,status):
    uniqmjd = sorted(np.unique(platelist['MJD']),reverse=True)
    plates = platelist[np.where(platelist['MJD'] == mjd)[0]]
    done = plates[np.where(plates['STATUS1D'] == 'Done')[0]]
    doneplates = plates[np.where(plates['STATUS1D'] == 'Done')[0]]['PLATE']

    if (len(plates) != len(done)) | (len(done) < 1):
        for m in range(len(uniqmjd)):
            mjd = uniqmjd[m+21]
            plates = platelist[np.where(platelist['MJD'] == mjd)[0]]
            done = plates[np.where(plates['STATUS1D'] == 'Done')[0]]
            doneplates = plates[np.where(plates['STATUS1D'] == 'Done')[0]]['PLATE']
            if not (len(plates) != len(done)) | (len(done) < 1):
                break
        for pl in plates:
            status.append('Plate: {} ({})'.format(pl['PLATE'],pl['PROGRAMNAME']))

    else:
        for pl in plates:
            status.append('Plate: {} ({})'.format(pl['PLATE'],pl['PROGRAMNAME']))
    return mjd,doneplates,status

def readPlateList(platelistfile='platelist.fits'):
    platelistdir,dummyname = setRedux()  # Update with correct path to $BOSS_SPECTRO_REDUX
    status=list() #;status.append(' ')
    wait_flag=False ; noplates_flag = False
    platelist= fits.open(os.path.join(platelistdir,platelistfile))[1].data
    mjd = sorted(platelist['MJD'])[-1]
    mjd,doneplates,status = outPlates(mjd,platelist,status)
    #print 'Working on MJD :',mjd
    if  os.path.exists(os.path.join('HTML',str(mjd),'balmonitor_'+str(mjd)+'.html')):
        wait_flag =True
        status=list()
        status.append('Last MJD with BALs = {}. BalMonitor has already seen this MJD'.format(mjd))
    return mjd,doneplates,status,wait_flag


def readSpAllfile(plateid):
    spAlldir,spallfile = setRedux()   # Update with correct path to $BOSS_SPECTRO_REDUX
    spAll = fits.open(os.path.join(spAlldir,spallfile))[1].data
    #print spAll.columns.names
    platematch = spAll[np.where(spAll['PLATE'] == plateid)[0]]
    balindex = np.where(platematch['EBOSS_TARGET2'] & 2**25)[0]
    print balindex
    return  platematch['FIBERID'][balindex],platematch['RA'][balindex],platematch['DEC'][balindex],platematch['Z'][balindex]

def sphereMatch(ra, dec, tol=0.00001):
    dupplate=[];dupmjd=[];dupfiber=[]
    dr14q = fits.open('DR14Q_v4_4.fits')[1].data
    a,b,ds= spherematch(ra,dec,dr14q['RA'],dr14q['DEC'],tol=tol)
    #print a,b,dr14q['RA'][b],ra,ds,dr14q['PLATE'][b]
    for i in range(len(b)):
        dupplate.append(dr14q['PLATE'][b[i]])
        dupmjd.append(dr14q['MJD'][b[i]])
        dupfiber.append(dr14q['FIBERID'][b[i]])
    return dupplate,dupmjd,dupfiber,b[i]


def makehtml(mjd, status, cat):
	images = glob.glob('HTML/{}/balmonitor*.jpeg'.format(mjd))
	html = list()
	
	html.append('<!DOCTYPE html>')
	html.append('<html>')
	html.append('<head>')
	html.append('<link rel="stylesheet" type="text/css" href="../../balmonitor.css">')
	html.append('</head>')
	html.append('<body>')
	html.append(' ')
	html.append('<h1>MJD %d, night of %s</h1>' % (mjd, convertMJD(mjd)))
	html.append('<hr>')
	html.append('<a style="float: left; width: 5%0%;" href= "/Users/vzm83/BALQSO_DailyMonitor/HTML/"+%s+"balmonitor_"+%s+".html">Yesterday: MJD=%d </a>' %(str(mjd-1),str(mjd-mjd-1),mjd-1))
	html.append('<a style="float: right; width: 50%%;text-align: right;"href= "/Users/vzm83/BALQSO_DailyMonitor/HTML/"+%s+"balmonitor_"+%s+".html"> Tomorrow: MJD=%d </a>' %(str(mjd+1),str(mjd+1),mjd+1))
	html.append('<hr>')
	html.append('<h2> Summary </h2>')
	for statusmessage in status:
	        html.append('<h3>%s</h3>' % (statusmessage))
	
	html.append('<hr>')
	html.append(' ')
	for image in images:
		html.append('<div class="row">')
		html.append('<div class="column left" "style="background-color:lightblue;" >')
		html.append('<h3a> %s </br></h3a>' %(image.split('/')[2].split('_')[1]+'-'+str( mjd)+'-'+image.split('/')[2].split('_')[3].split('.')[0]))
		match_string = image.split('/')[2].split('_')[1]+'_'+str( mjd)+'_'+image.split('/')[2].split('_')[3].split('.')[0]
		for outer_key in cat:
			cat_string = '{0}_{1}_{2}'.format(cat[outer_key]['plate'],cat[outer_key]['mjd'],cat[outer_key]['fiber'])
			print match_string,cat_string
			if cat_string == match_string:
				cat1=cat[outer_key]
				html.append('<p>RA : {0:8.4f} DEC : {1:8.4f}</p>'.format(cat1['ra'],cat1['dec']))
				html.append('<p>Z : {0:8.4f} Z_VI : {1:8.4f}</p>'.format(cat1['z'],cat1['zvi']))
				html.append('<p>MI : {0:8.4f} BI_CIV : {1:8.4f}</p>'.format(cat1['mi'],cat1['bi_civ']))
				html.append('<p>u:{0:4.2f} g:{1:4.2f} r:{2:4.2f} i:{3:4.2f} z:{4:4.2f}</p>'.format(cat1['mag'][0],cat1['mag'][1],cat1['mag'][2],cat1['mag'][3],cat1['mag'][4]))
				html.append('<p>N_SPEC : {}</p>'.format(cat1['nspec']))
				html.append('<p>N_SPEC_BOSS : {} N_SPEC_SDSS : {}</p>'.format(cat1['nspec_boss'],cat1['nspec_sdss']))
				html.append('<p>X-ray Flux(0.2-12 KeV) : {0:8.4f} FIRST-Flux : {1:8.4f}</p>'.format(cat1['xray_2to12kevflux'],cat1['first_flux']))
				html.append('<p>NUV : {0:4.2f} FUV : {1:4.2f}</p>'.format(cat1['nuv'],cat1['fuv']))
				html.append('<p>W1:{0:4.2f} W2:{1:4.2f} W3:{2:4.2f} W4:{3:4.2f} </p>'.format(cat1['w1'],cat1['w2'],cat1['w3'],cat1['w4']))
				html.append('</div>')
		html.append('<div class="column right" >')
		html.append('<img  src="%s" width="85%%">' % (os.path.abspath(image)))
		html.append('</div>')
		html.append('</div>')
	html.append('<hr>')
	html.append('</body>')
	html.append('</html>')
	return "\n".join(html)


def plotSpectra(pplate, pmjd, pfiber, splate, smjd, sfiber, z):
    
    # Reading the line lists from speccy
    linelist = np.genfromtxt('linelist_speccy.txt',usecols=(0,1,2),dtype=('|S10',float,'|S5'),names=True)
    bossspectroredux,dumm = setRedux() # Update with correct path to $BOSS_SPECTRO_REDUX
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
        if len(splate) > 0:
            download_spectra(splate[k],smjd[k],sfiber[k],'DR14_Spectra')
            dupfile = os.path.join('DR14_Spectra','spec-{0:04d}-{1:5d}-{2:04d}.fits'.format(splate[k],smjd[k],sfiber[k]))
        #    print dupfile
            slabel = '{0:04d}-{1:5d}-{2:04d}'.format(splate[k],smjd[k],sfiber[k])
            sdata = fits.open(dupfile)[1].data
            sflux = sdata.flux
            swave = 10**sdata.loglam
            serr = 1.0/np.sqrt(sdata.ivar)
            ax.plot(swave,gaussian_filter1d(sflux,5),color=color_list[k],ls='--',alpha=0.9,label=slabel)
            ax.plot(swave,serr,color=color_list[k],ls='--',alpha=0.1)
    ax.legend(loc=1)
    fig.tight_layout()
    destination_dir = os.path.join('HTML',str(pmjd))
    dir_cmd = 'mkdir {}'.format(destination_dir)
    if not os.path.exists(destination_dir):
        os.system(dir_cmd)
    savefilename = os.path.join(destination_dir,'balmonitor_{0}_{1}_{2}.jpeg'.format(pplate,pmjd,pfiber))
    fig.savefig(savefilename)
    #plt.show()


def makeCatalog(nitem, dr14index, plate, mjd, fiber, cat):
    outkey = 'target'+str(nitem)
    dr14Q = fits.open('DR14Q_v4_4.fits')[1].data
    cat[outkey] = {
                   "plate": plate,
                   "mjd" : mjd,
                   "fiber" : fiber,
                   "splate": dr14Q['PLATE'][dr14index],
                   "smjd"  :  dr14Q['MJD'][dr14index],
                   "sfiber" :  dr14Q['FIBERID'][dr14index],
                   "ra"  :  dr14Q['RA'][dr14index],
                   "dec" :  dr14Q['DEC'][dr14index],
                   "z" : dr14Q['Z'][dr14index],
                    "zvi" : dr14Q['Z_VI'][dr14index],
                    "nspec" : dr14Q['N_SPEC'][dr14index],
                    "nspec_boss" : dr14Q['N_SPEC_BOSS'][dr14index],
                    "nspec_sdss" : dr14Q['N_SPEC_SDSS'][dr14index],
                    "bi_civ" :  dr14Q['BI_CIV'][dr14index],
                    "mag":  dr14Q['PSFMAG'][dr14index],
                    "mi" : dr14Q['MI'][dr14index],
                    "xray_2to12kevflux": dr14Q['FLUX_0.2_12.0keV'][dr14index],
                    "first_flux" : dr14Q['FIRST_FLUX'][dr14index],
                    "nuv" : dr14Q['NUV'][dr14index],
                    "fuv" : dr14Q['FUV'][dr14index],
                    "w1" :  dr14Q['W1MAG'][dr14index],
                    "w2" :  dr14Q['W2MAG'][dr14index],
                    "w3" :  dr14Q['W3MAG'][dr14index],
                    "w4" :  dr14Q['W4MAG'][dr14index],
                    }

    return cat


    

def ProgMain():
    currentmjd = int(currentMJD())
    cat = {}
    mjd, qplates, status, wait_flag = readPlateList()
    print 'This',status,qplates,wait_flag
    if   wait_flag :
        html=makehtml(mjd, status,cat)
    else:
        for plate in qplates:
            fiber, ra, dec,z = readSpAllfile(plate)
            print fiber
            if (len(fiber) < 1):
                status.append('No BAL quasars included in Plate:{}'.format(plate))
             #   print status
                html=makehtml(mjd, status, cat)
            else:
                for i in range(len(fiber)):
                    print plate,mjd,fiber[i],ra[i],dec[i],z[i]
                    dupplate,dupmjd,dupfiber,dr14index = sphereMatch(ra[i],dec[i])
                    print 'PLATE-MJD-FIBER: {0}-{1}-{2} has {3} entries in DR14 '.format(plate,mjd,fiber[i],len(dupplate))
                    plotSpectra(plate, mjd, fiber[i], dupplate, dupmjd, dupfiber, z[i])
                    cat = makeCatalog(i, dr14index, plate, mjd, fiber[i], cat)
                    html=makehtml(mjd, status, cat)
                    
    html_dir = os.path.join('HTML',str(mjd))
    html_cmd = 'mkdir {}'.format(html_dir)
    if not os.path.exists(html_dir):
        os.system(html_cmd)
    htmlfile=os.path.join('HTML',str(mjd),'balmonitor_'+str(mjd)+'.html')
    if ((not os.path.isfile(htmlfile)) & (not wait_flag)):    
        fx=open(htmlfile,'w')
        print>>fx,html
        fx.close()
        latesthtmlfile = 'balmonitor_current.html'
        rm_cmd ='rm {}'.format(latesthtmlfile)
        ln_cmd = 'ln -s {} {}'.format(htmlfile,latesthtmlfile)
        latestmjd = 'CurrentMJD'
        lndir_cmd = 'ln -s {} {}'.format(html_dir,latestmjd)
        rmdir_cmd = 'rm {}'.format(latestmjd)
        os.system(rm_cmd)
        os.system(rmdir_cmd)
        os.system(ln_cmd)
        os.system(lndir_cmd)
    return mjd


if __name__ == "__main__":
    mjd = ProgMain()
    file_mtime = os.path.getmtime('balmonitor_current.html')
    
    #cmd = """
    #(
    #echo "From: $USER@`hostname`"
    #echo "To: %s"
    #echo "MIME-Version: 1.0"
    #echo "Content-Type: multipart/mixed;"
    #echo "Subject: Auto BAL Monitor Report %d (%s) "
    #echo ""
    #echo "This is a MIME-encapsulated message"
    #echo ""
    #echo "Content-Type: text/html"
    #echo ""
    #cat  balmonitor_current.html
    #echo ""
    #) | /usr/sbin/sendmail -t
    #""" % ('getkeviv@gmail.com', mjd, convertMJD(mjd), )
    #os.system(cmd)

