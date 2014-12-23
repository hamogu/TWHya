import glob
import sys

from astropy.io import ascii
from astropy.io import fits

sys.path.append('/data/hguenther/Dropbox/code/python/utils')
import read_IRAF_spec

smartslist = glob.glob('/data/hguenther/TWHya/SMARTSspec/*.ascii')
smartsspec = []
for f in smartslist:
    smartsspec.append(ascii.read(smartslist[0], Reader=ascii.NoHeader, names=['wave','flux','err']))

for spec in smartsspec:
    plt.plot(spec['wave'],spec['flux'])

range = [6553., 6573.]

range = [6536., 6547.]
indices = []
for s in smartsspec:
    indices.append((s.wave > range[0]) & (s.wave < range[1]))

for spec, i in zip(smartsspec, indices):
    plt.plot(spec['wave'][i],spec['flux'][i])

temp = np.correlate(smartsspec[0].flux[indices[0]]-np.mean(smartsspec[0].flux[indices[0]]), smartsspec[1].flux[indices[1]]-np.mean(smartsspec[1].flux[indices[1]]), mode='same')

AAOlist = glob.glob('/data/hguenther/TWHya/AAO/*fits')
AAOspec= []

aao = fits.open(AAOlist[0])
wave = read_IRAF_spec.make_IRAF_wave(aao[0].header)
for i in np.arange(aao[0].data.shape[0]):
    plt.plot(wave[i,:], aao[0].data[i,:])
