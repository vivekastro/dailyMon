import os
import glob
from time import sleep
run_dailymon = 'ssh u0992250@astro02.astro.utah.edu ./run_dailymon.sh'
os.system(run_dailymon)
print 'Finishing running the script'
sleep(120)
os.system('rm -rf /Users/vzm83/BALQSO_DailyMonitor/UtahCurrentMJD')
ucopycmd = 'scp -r u0992250@astro02.astro.utah.edu:~/BALQSO_DailyMonitor/CurrentMJD /Users/vzm83/BALQSO_DailyMonitor/UtahCurrentMJD'
print ucopycmd
os.system(ucopycmd)
print 'Finished copying the files from utah'
htmlfile = os.path.basename(glob.glob('/Users/vzm83/BALQSO_DailyMonitor/UtahCurrentMJD/balmonitor*.html')[0])
filelines=list()
with open(os.path.join('/Users/vzm83/BALQSO_DailyMonitor/UtahCurrentMJD',htmlfile),'r') as fr:
	for line in fr:
		if 'u0992250' in line:
			line = line.replace('/uufs/astro.utah.edu/common/home/u0992250/BALQSO_DailyMonitor/','')
		filelines.append(line)

with open(os.path.join('/Users/vzm83/BALQSO_DailyMonitor/UtahCurrentMJD',htmlfile),'w') as fw:
	for ll in filelines:
		print >> fw,ll

fw.close()


print 'Finished setting up the correct link'
	

mjd= htmlfile.split('_')[1].split('.')[0]
cpdir = 'cp -r /Users/vzm83/BALQSO_DailyMonitor/UtahCurrentMJD /Users/vzm83/BALQSO_DailyMonitor/HTML/'+mjd
rmdir = 'rm -rf  /Users/vzm83/BALQSO_DailyMonitor/HTML/'+mjd
print cpdir
os.system(rmdir)
os.system(cpdir)
print 'Copied the files from UtahCurrentMJD to HTML/mjd '
latesthtmlfile = '/Users/vzm83/BALQSO_DailyMonitor/balmonitor_current.html'
rm_cmd ='rm {}'.format(latesthtmlfile)
ln_cmd = 'ln -s {} {}'.format(os.path.join('/Users/vzm83/BALQSO_DailyMonitor/HTML',mjd,htmlfile),latesthtmlfile)
latestmjd = '/Users/vzm83/BALQSO_DailyMonitor/CurrentMJD'
html_dir='/Users/vzm83/BALQSO_DailyMonitor/HTML/'+mjd
lndir_cmd = 'ln -s {} {}'.format(html_dir,latestmjd)
rmdir_cmd = 'rm -rf {}'.format(latestmjd)
os.system(rm_cmd)
os.system(rmdir_cmd)
os.system(ln_cmd)
os.system(lndir_cmd)
print 'Finished setting up the correct links'
