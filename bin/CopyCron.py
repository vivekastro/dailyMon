import os
import glob

ucopycmd = 'scp -r u0992250@astro.utah.edu:~/BALQSO_DailyMonitor/CurrentMJD /Users/vzm83/BALQSO_DailyMonitor/UtahCurrentMJD'
print ucopycmd
os.system(ucopycmd)

htmlfile = os.path.basename(glob.glob('/Users/vzm83/BALQSO_DailyMonitor/UtahCurrentMJD/balmonitor*.html')[0])
mjd= htmlfile.split('_')[1].split('.')[0]
cpdir = 'cp -r /Users/vzm83/BALQSO_DailyMonitor/UtahCurrentMJD /Users/vzm83/BALQSO_DailyMonitor/HTML/'+mjd
print cpdir
os.system(cpdir)
latesthtmlfile = 'balmonitor_current.html'
rm_cmd ='rm {}'.format(latesthtmlfile)
ln_cmd = 'ln -s {} {}'.format(os.path.join('/Users/vzm83/BALQSO_DailyMonitor/HTML',mjd,htmlfile),latesthtmlfile)
latestmjd = 'CurrentMJD'
html_dir='/Users/vzm83/BALQSO_DailyMonitor/HTML/'+mjd
lndir_cmd = 'ln -s {} {}'.format(html_dir,latestmjd)
rmdir_cmd = 'rm {}'.format(latestmjd)
os.system(rm_cmd)
os.system(rmdir_cmd)
os.system(ln_cmd)
os.system(lndir_cmd)
