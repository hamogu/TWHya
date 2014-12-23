import numpy as np

from astropy.modeling.models import custom_model



'''This class is just 2* Gaussian + constant. 

Can be removed once astropy supports fitting of composite models.
Also (known bug, but not fixed yet), fitting does not work with bounds.
However, without bounds the fit will just put the gauss in the sharpest spike.
Workaround set in eval method
'''

def lineabs_model(x, amplitude=1, mean=1, stddev=1, const=1, amplitude2=-1, mean2=1, stddev2=1):
    if stddev < 0.1: stddev = 0.25
    return const + amplitude * np.exp((-(1 / (2. * stddev**2)) * (x - mean)**2)) - amplitude2 * np.exp((-(1 / (2. * stddev2**2)) * (x - mean2)**2))

def lineabs_deriv(x, amplitude=1, mean=1, stddev=1, const=1, amplitude2=-1, mean2=1, stddev2=1):
    d_amplitude = np.exp((-(1 / (stddev**2)) * (x - mean)**2))
    d_mean = (2 * amplitude *
          np.exp((-(1 / (stddev**2)) * (x - mean)**2)) *
          (x - mean) / (stddev**2))
    d_stddev = (2 * amplitude *
            np.exp((-(1 / (stddev**2)) * (x - mean)**2)) *
            ((x - mean)**2) / (stddev**3))
    d_const = np.zeros_like(x)
    d_amplitude2 = np.exp((-(1 / (stddev2**2)) * (x - mean2)**2))
    d_mean2 = (2 * amplitude2 *
          np.exp((-(1 / (stddev2**2)) * (x - mean2)**2)) *
          (x - mean) / (stddev2**2))
    d_stddev2 = (2 * amplitude2 *
            np.exp((-(1 / (stddev2**2)) * (x - mean2)**2)) *
            ((x - mean2)**2) / (stddev2**3))
    return [d_amplitude, d_mean, d_stddev, d_const, d_amplitude2, d_mean2, d_stddev2]

LineAbsModel = custom_model(lineabs_model, fit_deriv=lineabs_deriv)


'''This class is just Gaussian + constant. 

Can be removed once astropy supports fitting of composite models.
Also (known bug, but not fixed yet), fitting does not work with bounds.
However, without bounds the fit will just put the gauss in the sharpest spike.
Workaround set in eval method.
'''

def line_model(x, amplitude=1, mean=1, stddev=1, const=1):
        if stddev < 0.1: stddev = 0.25
        return const + amplitude * np.exp((-(1 / (2. * stddev**2)) * (x - mean)**2))

def line_deriv(x, amplitude=1, mean=1, stddev=1, const=1):
        d_amplitude = np.exp((-(1 / (stddev**2)) * (x - mean)**2))
        d_mean = (2 * amplitude *
              np.exp((-(1 / (stddev**2)) * (x - mean)**2)) *
              (x - mean) / (stddev**2))
        d_stddev = (2 * amplitude *
                np.exp((-(1 / (stddev**2)) * (x - mean)**2)) *
                ((x - mean)**2) / (stddev**3))
        d_const = np.zeros_like(x)
        return [d_amplitude, d_mean, d_stddev, d_const]


LineModel = custom_model(line_model, fit_deriv=line_deriv)


def FileFormatException(Exception):
    pass

def read(cls, filename):
    '''specific to COS - can be made into a Reader

    Should also work for STIS long slit'''
    data = []

    tab = table.Table.read(filename)
    for c in tab.columns:
        if len(tab[c].shape) == 2:
            data.append(table.Column(data=tab[c].data.flatten(),
                                     unit=tab[c].unit,
                                     name=c))
        else:
            raise FileFormatException('File format does not match COS x1d file standard.')
    flattab = cls(data, meta=tab.meta, dispersion='WAVELENGTH')
    # COS is not an echelle spectrograph, so there never is any overlap
    flattab.sort('WAVELENGTH')
    return flattab


def wave_little_interpol(wavelist):
    '''Make a wavelengths array for merging echelle orders with little interpolation.

    In echelle spectra we often have the situation that neighboring orders overlap
    a little in wavelength space::

        aaaaaaaaaaaa
                 bbbbbbbbbbbbb
                          ccccccccccccc

    When merging those spectra, we want to keep the original wavelength grid where possible.
    This way, we only need to interpolate on a new wavelength grid where different orders
    overlap (here ``ab`` or ``bc``) and can avoid the dangers of flux interpolation in
    those wavelength region where only one order contributes.

    This algorithm has limitations, some are fundamental, some are just due to the 
    implementation and may be removed in future versions:

    - The resulting grid is **not** equally spaced, but the step size should not vary too much.
    - The wavelength arrays need to be sorted in increasing order.
    - There has to be overlap between every order and every order has to have some overlap
      free region in the middle.

    Parameters
    ----------
    wavelist : list of 1-dim ndarrays
        input list of wavelength

    Returns
    -------
    waveout : ndarray
        wavelength array that can be used to co-adding all echelle orders.
    '''
    mins = np.array([min(w) for w in wavelist])
    maxs = np.array([max(w) for w in wavelist])
    
    if np.argsort(mins) != np.arange(len(wavelist)):
        raise ValueError('List of wavelengths must be sorted in increasing order.')
    if np.argsort(mins) != np.arange(len(wavelist)):
        raise ValueError('List of wavelengths must be sorted in increasing order.')
    if not np.all(maxs[:-1] < mins[1:]):
        raise ValueError('Not all orders overlap.')
    if np.any(mins[2:] < maxs[:-2]):
        raise ValueError('No order can be completely overlapped.')


    waveout = []
    waveout.append(wavelist[0][wavelist[0]< mins[1]])
    for i in range(len(wavelist)-1):
        #### overlap region ####
        # No assumptions on how bin edges of different orders match up
        # In overlap region patch in a linear scale with slightly different step.
        dw = maxs[i] - mins[i+1]
        step = 0.5*(np.mean(np.diff(wavelist[i])) + np.mean(np.diff(wavelist[i+1])))
        n_steps = np.int(dw / step + 0.5)
        # overlap start and stop are the last and first "clean" points.
        overlap_start = np.max(waveout[-1])
        overlap_end = np.min(wavelist[i+1][wavelist[i+1] > maxs[i]])
        wave_overlap = np.linspace(overlap_start + step,  overlap_end - step, n_steps-1)
        waveout.append(wave_overlap)

        #### next region without overlap ####
        if i < (len(wavelist) -1):  # normal case
            waveout.append(wavelist[i+1][(wavelist[i+1] > maxs[i]) & (wavelist[i+1]< mins[i+2])])
        else:                       # last array - no more overlap behind that
            waveout.append(wavelist[i+1][(wavelist[i+1] > maxs[i])])

    return waveout
