# -*- coding: utf-8 -*-
#dmcopy "7435/primary/acisf07435N002_evt2.fits[EVENTS][grade=0,2,3,4,6,status=0,tg_srcid=1,tg_m=-1,1][cols time,tg_lam,tg_m,tg_part]" 7435/evt2_short.fits.gz
#dmcopy "7436/primary/acisf07436N002_evt2.fits[EVENTS][grade=0,2,3,4,6,status=0,tg_srcid=1,tg_m=-1,1][cols time,tg_lam,tg_m,tg_part]" 7436/evt2_short.fits.gz
#dmcopy "7437/primary/acisf07437N002_evt2.fits[EVENTS][grade=0,2,3,4,6,status=0,tg_srcid=1,tg_m=-1,1][cols time,tg_lam,tg_m,tg_part]" 7437/evt2_short.fits.gz
#dmcopy "7438/primary/acisf07438N002_evt2.fits[EVENTS][grade=0,2,3,4,6,status=0,tg_srcid=1,tg_m=-1,1][cols time,tg_lam,tg_m,tg_part]" 7438/evt2_short.fits.gz
basedir = '/data/hguenther/obs/Chandra/'
basedir='/home/moritz/Chandra/TW_Hya/'
basedir='/scratch/hspc47/steh305/XVP/TW_Hya/'
import atpy
import numpy as np
import scipy
import scipy.stats
import scipy.odr

from TWHya_plots import *

d1 = atpy.Table(basedir+'7435/evt2_short.fits.gz', hdu = 1)
d2 = atpy.Table(basedir+'7436/evt2_short.fits.gz', hdu = 1)
d3 = atpy.Table(basedir+'7437/evt2_short.fits.gz', hdu = 1)
d4 = atpy.Table(basedir+'7438/evt2_short.fits.gz', hdu = 1)

heg = d1.where(d1.tg_part == 1)
heg.append(d2.where(d2.tg_part == 1))
heg.append(d3.where(d3.tg_part == 1))
heg.append(d4.where(d4.tg_part == 1))
heg.sort(['time'])

meg = d1.where(d1.tg_part == 2)
meg.append(d2.where(d2.tg_part == 2))
meg.append(d3.where(d3.tg_part == 2))
meg.append(d4.where(d4.tg_part == 2))
meg.add_column('t', meg.time-287906000.)
meg.sort(['t'])

mner = around(meg.tg_lam, 13.45)
mnei = around(meg.tg_lam, 13.55)
mnef = around(meg.tg_lam, 13.70)
mneb = around(meg.tg_lam, 11.545)
mnec = around(meg.tg_lam, 11.00)
mne10 = around(meg.tg_lam, 12.14)

hner = around(heg.tg_lam, 13.45, delta = 0.01)
hnei = around(heg.tg_lam, 13.55, delta = 0.01)
hnef = around(heg.tg_lam, 13.70, delta = 0.01)
hneb = around(heg.tg_lam, 11.54, delta = 0.01)
hnec = around(heg.tg_lam, 11.00, delta = 0.01)
hne10 = around(heg.tg_lam, 12.14)


obs1 = meg.t < 2e5
obs2 = (meg.t > 9e5) & (meg.t < 1.15e6)
obs3 = (meg.t > 1.15e6) & (meg.t < 1.35e6)
obs4 = meg.t > 1.35e6


n10 = np.zeros((16,len(meg)),dtype = np.bool)
for i in range(16): n10[i,i*2000:(i+1)*2000] = True

megp = (meg.tg_m == 1)
megn = (meg.tg_m == -1)
dist = meg.t

ind_meg_ne = {'r':mner, 'i':mnei, 'f': mnef, 'beta': mneb, 'gamma': mnec, 'Ha': mne10}
ind_meg_nen = {'r':mner & megn, 'i':mnei & megn, 'f': mnef & megn, 'beta': mneb & megn, 'Ha': mne10 & megn}

for i in range(len(dist)): dist[i] = np.min(np.abs(meg.time[mner & megp]-meg.time[i]))
plot_rat(mner , mneb, dist<100, (dist>100) & (dist<200),  (dist > 200)&(dist<400), (dist > 400), color = 'y')

plot_rat(mnef , mnei,obs1,obs2,obs3, obs4, color = 'g')

plt.clf()
plot_ratios(ind_meg_ne, obs1, obs2, obs3, obs4)
plt.title('sorted by observation')
plt.savefig('obs.png')

plt.clf()
plot_ratios(ind_meg_ne, *n10)
ylim(0,5)
plt.title('binned to equal MEG count number')
plt.savefig('binned.png')

plt.clf()
plot_ratios(ind_meg_nen, dist<100, (dist>100) & (dist<300),  (dist > 300)&(dist<800), dist > 800)
plt.title('distance to MEG He r photon (<100, 100<t<500, 500<t<1200, >1200)')
plt.savefig('disttoHEG.png')

plt.clf()
plt.hist(dist, bins=20, range=(0,5000))
plt.savefig('hist.png')

plt.clf()
print '(f+i) / ((f+i)/r)'
r,i,f = plot_ratrat(ind_meg_ne['r'], ind_meg_ne['i'], ind_meg_ne['f'], *n10, color ='y')
plt.savefig('fi-vs-fir.png')


#scipy.stats.pearsonr((f+i)/r,f+i)
#scipy.stats.spearmanr((f+i)/r,f+i)

# find flare in counts
#plt.clf()
#hist = plt.hist(meg.t,bins=1000)
## flare 999500-1005000 (element 13740,14498

#scipy.stats.linregress(f/i,(f+i)/r)


def func(B, x):
    ''' Linear function y = m*x + b '''
    return B[0]*x + B[1]


# B is a vector of the parameters.
# x is an array of the current x values.
# x is same format as the x passed to Data or RealData.


x = np.arange(0.3,0.8,0.1)




lya = split_bins(ind_meg_ne['Ha'] ,*n10) 
r = split_bins(ind_meg_ne['r'] ,*n10)
i = split_bins(ind_meg_ne['i'] ,*n10)
f = split_bins(ind_meg_ne['f'] ,*n10)
beta = split_bins(ind_meg_ne['beta'] ,*n10)
gamma = split_bins(ind_meg_ne['gamma'] ,*n10)
#A_eff in cm2 from Nancy's article
flya = lya / 48.5 
fr = r / 22.3
fi = i / 22.7
ff = f / 25.
fbeta = beta / 55.3
fgamma = gamma / 67.4

a=f
b=i
c=f+i
d=r
plt.clf()
lastplot = plt.errorbar(ff/fi,(ff+fi)/fr, xerr = np.abs(error(f,i,25.,22.7)), yerr = np.abs(error2(f,i,r,25.,22.7, 22.3)), fmt = 'o')

linear = scipy.odr.Model(func)
mydata = scipy.odr.RealData(a/b,c/d, sx = gauss_error(a,b), sy = gauss_error(c,d))
myodr = scipy.odr.ODR(mydata, linear, beta0=[1., 2.])
myoutput = myodr.run()
myoutput.pprint()
x = np.arange(0.1,4.,0.1)
plt.plot(x,func(myoutput.beta,x))
myodr = scipy.odr.ODR(mydata, linear, beta0=[-1., 2.])
myoutput = myodr.run()
myoutput.pprint()
x = np.arange(0.1,4.,0.1)
plt.plot(x,func(myoutput.beta,x))

plt.xlabel(r'f/i')
plt.ylabel(r'(f+i)/r')
plt.savefig('rbfir.png')

plt.arrow(1., 2., -0.4, 0., width=0.01, color = 'b')
plt.text(0.7, 1.9, "increasing density", ha="center", size=16)
plt.arrow(1., 2., 0., -0.8, width=0.004, color = 'b')
plt.text(0.9, 1.4, "increasing", ha="center", size=16)
plt.text(0.9, 1.3, "temperature", ha="center", size=16)

plt.savefig('/scratch/hspc47/steh305/Dropbox/my_proposals/Chandra/XVP/firfi.ps')

import scipy.interpolate
import chianti.core as ch

class LineRatioInterpolator():
    def __init__(self, interpolator):
        self.interp = interpolator
    def __call__(self,x):
        x_new = np.minimum(max(self.interp.x),x)
        x_new = np.maximum(min(self.interp.x),x_new)
        return self.interp(x_new)

nH = np.log(rat/rat_0)/(sig_r-sig_beta)
def temp_from_rat_interp(ionstring, nom, denom, temp = 10.**np.arange(5.,8.,.1)):
    ion = ch.ion(ionstring, temperature = 2e6, density = 1e8)
    ion.emiss(temperature = temp, density = 1e8)
    indnom = (ion.Emiss['wvl'] > min(nom) ) & (ion.Emiss['wvl'] < max(nom))
    inddenom = (ion.Emiss['wvl'] > min(denom) ) & (ion.Emiss['wvl'] < max(denom))
    ratio  =  ion.Emiss['emiss'][indnom,:].sum(axis=0) /  ion.Emiss['emiss'][inddenom,:].sum(axis=0)
    ind = ratio.argsort()
    return LineRatioInterpolator(scipy.interpolate.interp1d(ratio[ind], temp[ind]))

def dens_from_rat_interp(ionstring, nom, denom, dens = 10.**np.arange(8.,15.,.1)):
    ion = ch.ion(ionstring, temperature = 2e6, density = dens)
    ion.emiss(temperature = 2e6, density = dens)
    indnom = (ion.Emiss['wvl'] > min(nom) ) & (ion.Emiss['wvl'] < max(nom))
    inddenom = (ion.Emiss['wvl'] > min(denom) ) & (ion.Emiss['wvl'] < max(denom))
    ratio  =  ion.Emiss['emiss'][indnom,:].sum(axis=0) /  ion.Emiss['emiss'][inddenom,:].sum(axis=0)
    ind = ratio.argsort()
    return LineRatioInterpolator(scipy.interpolate.interp1d(ratio[ind], dens[ind]))

def rat_from_temp_interp(ionstring, nom, denom, temp = 10.**np.arange(5.,8.,.1)):
    ion = ch.ion(ionstring, temperature = temp, density = 1e8)
    ion.emiss(temperature = temp, density = 1e8)
    indnom = (ion.Emiss['wvl'] > min(nom) ) & (ion.Emiss['wvl'] < max(nom))
    inddenom = (ion.Emiss['wvl'] > min(denom) ) & (ion.Emiss['wvl'] < max(denom))
    ratio  =  ion.Emiss['emiss'][indnom,:].sum(axis=0) /  ion.Emiss['emiss'][inddenom,:].sum(axis=0)
    ind = temp.argsort()
    return LineRatioInterpolator(scipy.interpolate.interp1d(temp[ind], ratio[ind]))


n_from_f2i = dens_from_rat_interp('ne_9', [13.69, 13.71], [13.54, 13.56])
t_from_fir = temp_from_rat_interp('ne_9', [13.54, 13.71], [13.44, 13.46], temp = 10.**np.arange(5.5,7.,.1))
t_from_betar = temp_from_rat_interp('ne_9', [11.54, 11.56], [13.44, 13.46], temp = 10.**np.arange(5.5,7.,.1))
betar_from_t = rat_from_temp_interp('ne_9', [11.54, 11.56], [13.44, 13.46], temp = 10.**np.arange(5.5,7.,.1))

plt.loglog(n_from_f2i(f/i),t_from_fir((f+i)/r))
plt.xlabel('density [$cm^{-3}$]')
plt.ylabel('temperature [K]')
plt.savefig('tempdens.png')

beta_bg_ind = around(meg.tg_lam, 11.545, delta=0.4)
beta_bg = split_bins(beta_bg_ind ,*n10)
beta_bg = (beta_bg-beta)/(0.4/0.02)

plt.clf()
plt.loglog(t_from_betar(beta/r), t_from_fir((f+i)/r))
plt.clf()
plt.loglog(t_from_fir((f+i)/r),n_from_f2i(f/i), 'bs' )

rat=(beta)/r
rat_0 = betar_from_t(t_from_fir((f+i)/r))

sig_r = 3.07e-22
sig_beta = 2.00e-22

nH = np.log(rat/rat_0)/(sig_r-sig_beta)
T = t_from_fir((f+i)/r)
dens = n_from_f2i(f/i)
for t,d,n in  np.vstack([T,dens,nH]).transpose():
    print "{0:8.1e} {1:8.1e} {2:8.1e}".format(t,d,n)

# get some idea of the errors involved in this
ff/fi,, xerr = np.abs(error(f,i,25.,22.7))
(ff+fi)/fr yerr = np.abs(error2(f,i,r,25.,22.7, 25.))
(beta/55.3)/fr , zerr = np.abs(error(beta,r,55.3,25.))

ratdens = ff/fi
ratdenserr = error(f,i,25.,22.7)[1,:]
rattemp = (ff+fi)/fr
rattemperr = error2(f,i,r,25.,22.7, 25.)[1,:]
ratnh = (beta/55.3)/fr
ratnherr = error(beta,r,55.3,25.)[1,:]

print (ff+fi)/fr
oben = 3
mitte = 4
unten = 14

for sel in [3,4,14]:
    print sel
    T = t_from_fir(rattemp[sel]+rattemperr[sel]*   np.array([0, 0, 0, 0, 1, 1, 1,-1,-1,-1]))
    dens = n_from_f2i(ratdens[sel]+ratdenserr[sel]*np.array([0, 1, 0,-1, 0, 1,-1,0 , 1,-1]))
    nH = np.log((ratnh[sel]+ratnherr[sel]*         np.array([0, 1, 1, 0, 1,-1,-1,-1, 0, 1]))/betar_from_t(T))/(sig_r-sig_beta)
    for t,d,n in  np.vstack([T,dens,nH]).transpose():
        print "{0:8.1e} {1:8.1e} {2:8.1e}".format(t,d,n)