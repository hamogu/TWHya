import re

from scipy.io import readsav

def read_pandora_input(filename):
    '''read a Z,TE,NH,VXS structure in PANDORA format'''
    with open(filename) as f:
        content = f.readlines()
    content = [x.strip('\n') for x in content]
    content = ' '.join(content)
    params = content.split(') >')[:4]  # last one is empty string
    out = {}
    for p in params:
        n, v = p.split('( >')
        # replace multiple spaces with single space, othwerwise
        # split return empty values
        # also remove trailing spaces
        v = re.sub(r'\s+', ' ', v.strip())
        out[n.strip()] = np.array(v.split(' '), dtype=np.float)
    return out

sun = read_pandora_input('/data/hguenther/TWHya/eugene/eugene_sun.dat')
model = readsav('/data/guenther/tt_results/large_grid/12_6.sav')
hyd = model['hydrodyn']

# chop of first and last 2 km
# Those have a lot of layers with very little emission measure
# For a simple model it would be a waste of time to include those
ind = (hyd[:,0] > 2e5) & (hyd[:,0] < (hyd[0,0]-2e5))

hyd = hyd[ind,:]
# make coarser grid to save run time
hyd = hyd[::2, :]

# Assembles arrays for writing
out = {}
out['Z'] = hyd[::-1, 0]/1e5 # in km, reverse ordering, reverse direction
out['TE'] = hyd[::-1, 5]  # reverse ordering
out['NH'] = hyd[::-1, 2]  # reverse ordering
out['VXS']= -hyd[::-1,3]/1e5  # in km/s, reverse ordering, reverse direction

# attach solar atmosphere at the bottom
# match by temperature
indsun = sun['TE'] < np.min(hyd[:,5])
sun['Z'] = sun['Z'] - np.min(sun['Z'][indsun]) + np.max(out['Z'])+1 
out['Z'] = np.hstack((out['Z'], sun['Z'][indsun]))
for n in ['TE', 'NH', 'VXS']:
    out[n] = np.hstack((out[n], sun[n][indsun]))

# make the number of layers dividable by 5 for easier writing
for name in out:
    n = out[name].shape[0]//5
    out[name] = out[name][:n*5]


with open('/data/hguenther/TWHya/eugene/model12_6.dat', 'w') as f:
    for section in ['Z', 'TE', 'NH', 'VXS']:
        f.write('{0:5s}             ( >\n'.format(section))
        for i in range(n):
            f.write('{0: 1.8E} {1: 1.8E} {2: 1.8E} {3: 1.8E} {4: 1.8E}\n'.format(*out[section][i*5:(i+1)*5]))
        f.write(') >                \n')
    

### Read and plot the results
# find start and end by hand and put in here - Do not count empty lines!
from astropy.table import Table
civ = Table.read('hgxc4.aaa.23643', format='ascii', data_start=26, data_end=114)
civ = Table.read('hgxc4.aaa.23643.c4', format='ascii', data_start=15)

plt.plot(civ['col2'], civ['col6'])
