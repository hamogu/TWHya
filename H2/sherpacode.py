



import os
from itertools import groupby

import numpy as np

import sherpa.astro.ui as ui
# Sherpa currently screws up the logger. This is a workaround:
import logging
import sys
logging.getLogger('sherpa').propagate = 0
sys.tracebacklimit = None

from sherpa.data import Data1D


from spectrum import Spectrum

from H2 import read_Abgrall93

# small helper scripts
def _isfloat(value):
    '''Check if a string can be converted to a float'''
    try:
        float(value)
        return True
    except ValueError:
        return False

def _set_val(model, argstring, val):
    '''Set the value of a Sherpa model parameter

    Parameters
    ----------
    model : sherpa model instance
    argstring : string
        name of the parameter, e.g. 'kT' or 'norm'
    val : string or number
        Value to be set. Can be number or string of the form
        'f 1234' which set the parameter to 1234 and freezes it.
    '''
    par = ui.get_par(model.name+'.'+argstring)
    if isinstance(val, basestring):

        if val[0] == '"':
            val = val[1:-1]
        val = val.strip('" \t').strip("' \t").split()
        par.frozen = False
        for v in val:
            print v
            if v == 'f':
                par.frozen = True
            elif _isfloat(v):
                par.val = float(v)
            else:
                raise ValueError('Not valid format to set a parameter {0}'.format(v))
        
    else:
        par.frozen = False
        par.val = val


def _place_val(line, name, val):
    if name not in line:
        line[name] = val
    else:
        oldval = line[name]
        if oldval[0] == '"':
            oldval = oldval[1:-1]
        oldval = oldval.strip('" \t').strip("' \t").split()
        oldval = oldval[0:-1]
        oldval.append(str(val))
        line[name] = ' '.join(oldval)
        

def _parminmax(res, par):
    if par in res.parnames:
        i = res.parnames.index(par)
        return res.parmins[i], res.parmaxes[i]
    else:
        return 0,0

 
def set_wave(model, wave):
    model.pos.min = 0
    model.pos.max = 1e5
    _set_val(model, 'pos', wave)
    model.pos.min = wave - 0.2
    model.pos.max = wave + 0.2


printfiletype = '.png'

def read_spreadsheet(filename):
    '''Read Gabriel's spreadsheet

    This spreadsheet lists which regions and which lines should be fitted
    together and where the lines are spaced far enough apart to be treated
    separately.
    '''
    regions = []
    Abgrall93 = read_Abgrall93()
    with open(filename) as f:
        header = f.readline()
        filecontent = f.readlines()
    for line in filecontent:
        newline = line.strip().split(',')
        regions.append(Region(newline, Abgrall93=Abgrall93, oldformat=True))
    
    return regions


def read_output_spreadsheet(filename):
    '''Read a fit result spreadsheet as starting point for a new fit

    This can be used to tweak parameters in the fit result in cases where 
    the fit got stuck in a local minimum. 

    Parameters
    ----------
    filename : string

    Returns
    -------
    regions : list of :class:`Region` instances
        

    '''
    Abgrall93 = read_Abgrall93()
    regions = []
    with open(filename) as f:
        header = f.readline()
        filecontent = f.readlines()

    filecontent = [line.strip().split(',') for line in filecontent]

    for name, region in groupby(filecontent, lambda x: x[0]):
        r = list(region)
        regions.append(Region(r, Abgrall93=Abgrall93))

    return regions

def interpret_line_code(code, Abgrall93, oldformat=False):
    '''interpret code for line in Grabriel's table

    Parameters
    ----------
    code : string
        Line name in Gabriel's table. Could be an H_2 line or an
        unidentified line starting with ``?`` followed by an
        approximate wavelength.

    Returns
    -------
    name : string
        some identifying name
    wave : float
        rest wavelength (for H_2 lines) or guessed wavelength (for other lines)
        in Angstroem
    absorption : bool
        True if absorption lines, False for emission lines
    '''

    # genumeric (used to write the table) surrounds all string with "  "
    if code[0] == '"':
        code = code[1:-1]
    code = code.strip()

    if code[0:2] == '?a':
        name = code
        wave = float(code[3:])
        absorption = True
    elif code[0] == '?':
        name = code
        wave = float(code[1:])
        absorption = False
    elif code[0] in 'prPRH':
    #elif (oldformat and (code[0] in 'pr')) or (not oldformat and (code[0:2] == 'H2')):
            
        if oldformat:
            codesplit = code.split()
            name = '{0}({1}) {2}-{3}'.format(code[0].upper(),
                          codesplit[0][1:],codesplit[1], codesplit[2])
        else:
            if code[:3] == 'H2 ':
                name = code[3:]
            else: 
                name = code
        
        wave = Abgrall93[name]
        name = 'H2 '+ name
        absorption = False
    else:
        raise ValueError('code: {0} does not have proper format'.format(code))
    
    return name, wave, absorption
          
        
        
class Datafile(Spectrum):
    def __init__(self, *args, **kwargs):
        super(Datafile, self).__init__(*args, **kwargs)
        self.data = Data1D(1, self.disp, self.flux, self.error)


    def fit_all(self, regions, outfile, conf=False):
        for region in regions:
            if region.FN == self.FN:
                if self.detector != region.detector:
                    self.load_detector(region.detector)
                region.fit(conf=conf)
                region.plot(filename = '_'+self.filename[:-4])
                region.append_to_csv(outfile, conf=conf)


# This should go in a project specific file.
# While the rest of this file is project almost agnostic.
from sherpa import models
from sherpa.astro import models as astromodels

# Make one instance of every type of model that we need
# and set sensible defaults for the numbers.
constbase = models.Const1D('baseconst')
constbase.c0 = 5e-14

linebase = astromodels.Lorentz1D('linebase')
linebase.fwhm = 0.07
linebase.fwhm.min = .04
linebase.fwhm.max = 1.
linebase.ampl.max = 5e-12
linebase.ampl = 2e-13
linebase.ampl.min = 0

abslinebase = linebase.__class__('abslinebase')
copy_pars(linebase, abslinebase)
abslinebase.ampl.min = -2e-12
abslinebase.ampl = -2e-13
abslinebase.ampl.max = 0

H2linebase = linebase.__class__('H2linebase')
copy_pars(linebase, H2linebase)
H2linebase.fwhm.max = 0.09
H2linebase.ampl.max = 1e-12
H2linebase.fwhm.val = 0.05681
H2linebase.fwhm.fixed = True

modelselector = [('const', constbase),
                 ('H2 ', H2linebase),
                 ('\?ab', abslinebase),
                 ('\?', linebase)]

# From here on it's general again
import re


def name2model(name, modelselector):
    '''Instantiate a Sherpa model based on a string model name

    Parameters
    ----------
    name : string
        ``name`` should identify a line.
        This string is also used as the name for a model in the sherpa.ui interface.
    modelselector : list of tuples
        Each tuple should consist of at least two parts:
            1) A regular expression
            2) A Sherpa model **instance** (not class)

    Returns
    -------
    newmodel : sherpa model instance
        The newly created sherpa model.

    Example
    -------
    In this example we want to create models for a constant and Gaussian lines::

        >>> from sherpa import models
        >>> baseline  = models.Gauss1D('baseline')
        >>> baseconst = models.Const1D('baseconst')

    Because the base models are model instances, their parameters
    can be changed and those will be copied to all models created by this
    function::

        >>> baseline.fwhm = 0.1

    So, we now have a list of strings that specify which models we want::

        >>> modelnames = ['consto7', 'O_VII_r', 'O_VII_i, 'O_VII_f']

    All models starting with ``'const'`` shall be constants, all others
    shall be lines. The function tries the list from the first to the last item.
    If ``'const'`` matches, it returns the ``baseconst`` model, if not a
    ``baseline`` model (The regular expression ``'.*'`` matches anything except
    newline)::

        >>> modelselector = [('const', baseconst), ('.*', baseline)]
        >>> o7_modellist = [name2model(s, modelselector) for s in modelnames]
    '''

        
    for m in modelselector:
        if re.match(m[0], name):
            newmodel = m[1].__class__(name)
            copy_pars(m[1], newmodel)
            return newmodel

    raise ValueError('No model definition found for {0}.'.format(name))
        

class Region(object):

    plotpath = '/data/guenther/TWHyaplots'
    printfiletype = '.png'

    def __init__(self, region, Abgrall93 = None):
        '''Read Gabriel's spreadsheet

        This spreadsheet lists which regions and which lines should be fitted
        together and where the lines are spaced far enough apart to be treated
        separately.
        '''
        if Abgrall93 is None:
            Abgrall93=read_Abgrall93()
                          
        name  = region[0][0]
        FN = region[0][1]
        if FN[0] == '"':
             FN = FN[1:-1]

        self.name = name
        self.FN = FN[0]
        self.detector = int(FN[2])
        self.start = float(region[0][2])
        self.stop = float(region[0][3])
        self.const = region[0][5]

        self.H2lines = []
        self.nonH2lines = []
        for line in region:
            print line 
            if len(line) > 0:
                name, wave, absorption = interpret_line_code(line[8], 
                                                         Abgrall93)
                if name[0:2] == 'H2':
                    # ignore the pos value for H2 lines.
                    # will be set from Abgrall line list
                    self.H2lines.append({'name': name, 'wave': wave, 'abs': absorption, 'fwhm': line[11], 'ampl': line[12]})
                else:
                    self.nonH2lines.append({'name': name, 'wave': wave, 'abs': absorption, 'pos': line[10], 'fwhm': line[11], 'ampl': line[12]})


    def set_source(self):
        if self.FN == 'F':
            modelstring = 'empG160M'
        elif self.FN == 'N':
            modelstring = 'tabNUV'
        else:
            raise ValueError('F - FUV/G160M, N - NUV/G285M, not recognized: {0}').format(self.FN)


        modelstring = modelstring + '(const1d.c1'
        for i, line in enumerate(self.H2lines):
            line['source'] = 'lorentz1d.h{0}'.format(i)
            modelstring = modelstring + '+ ' + line['source']
        for i, line in enumerate(self.nonH2lines):
            line['source'] = 'lorentz1d.l{0}'.format(i)
            modelstring = modelstring + '+ ' + line['source']
        modelstring = modelstring + ')'
        print modelstring
        ui.set_source(modelstring)

        # set some reasonable limits for the model parameters
        # to increase the chances of achieving a reasonable fit

        # If there is more than 1 H_2 line, we link the wavelength together
        if len(self.H2lines) > 0:
            wave_base = self.H2lines[0]['wave']
            model_base = ui.get_model_component('h0')
            set_wave(model_base, wave_base)
            model_base.fwhm = 0.07
            model_base.fwhm.min = .04
            model_base.fwhm.max = .09
            model_base.ampl.max = 1e-12
            model_base.ampl = 2e-13
            model_base.ampl.min = 0
            for i, line in enumerate(self.H2lines[1:]):
                modelcomp = ui.get_model_component('h{0}'.format(i+1))
                set_wave(modelcomp, line['wave'])
                modelcomp.pos = model_base.pos + (line['wave'] - wave_base)
                modelcomp.fwhm = 0.07
                modelcomp.fwhm.min = .04
                modelcomp.fwhm.max = .09
                modelcomp.ampl.max = 1e-12
                modelcomp.ampl = 2e-13
                modelcomp.ampl.min = 0

        for i, line in enumerate(self.nonH2lines):
            modelcomp = ui.get_model_component('l{0}'.format(i))
            set_wave(modelcomp, line['wave'])
            modelcomp.fwhm = 0.07
            modelcomp.fwhm.min = .04
            modelcomp.fwhm.max = 1.
            modelcomp.ampl.max = 5e-12
            modelcomp.ampl = 2e-13
            modelcomp.ampl.min = 0
            if line['abs']:
                # The order of these statements is important, because you cannot
                # set a value below the min.
                modelcomp.ampl.min = -2e-12
                modelcomp.ampl = -2e-13
                modelcomp.ampl.max = 0


        # If the input file specified those values, they take precedence
        # over the hard-coded default value      
        for line in (self.H2lines + self.nonH2lines):
            for n in ['pos','fwhm','ampl']:
                if n in line:
                    model = ui.get_model_component(line['source'].split('.')[1])
                    _set_val(model, n, line[n])

        model = ui.get_model_component('c1')
        _set_val(model, 'c0', self.const)

        #ui.get_par('c1.c0').val = self.const



    def fit(self, conf = False):


        ui.ignore(None, None)
        ui.notice(self.start, self.stop)
        self.set_source()
        ui.fit(1)
        if conf:
            ui.conf()
        res = ui.get_fit_results()

        for line in (self.H2lines + self.nonH2lines):
            sourcename = line['source'].split('.')[1]
            print sourcename
            for p in ['pos', 'fwhm', 'ampl']:
                n = '{0}.{1}'.format(sourcename, p)
                _place_val(line, p, ui.get_par(n).val)

        self.const = ui.get_par('c1.c0').val
        self.redchi2 = res.rstat


        if conf:
            res = ui.get_conf_results()
            for line in (self.H2lines + self.nonH2lines):
                sourcename = line['source'].split('.')[1]
                for p in ['pos', 'fwhm', 'ampl']:
                    n = '{0}.{1}'.format(sourcename, p)
                    parmin, parmax = _parminmax(res, n)
                    line[p+'_max'] = parmax
                    line[p+'_min'] = parmin
            # deal with error on const
            parmin, parmax = _parminmax(res, 'c1.c0')
            self.const_min = parmin
            self.const_max = parmax


    def plot(self, filename=''):
        try:
            import chips as ch
            self.plot_chips(filename=filename)
        except ImportError:
            import matplotlib as mpl
            import matplotlib.pyplot as plt
            self.plot(mpl(filename=filename)

    def plot_mpl(self, filename=''):
        raise NotImplementedError

    def plot_chips(self, filename=''):
        ui.plot_fit_delchi(1)
        fig = plt.gcf()
        ch.add_label(0.5, 0.97, self.name + '_'+filename,
                     ["coordsys", ch.FRAME_NORM])

        res = ui.get_fit_results()
        ch.add_label(0.2,0.92, '\\chi^2_{{red}} = {0:3.1f}'.format(res.rstat),
                     ["coordsys", ch.FRAME_NORM])

        for i in range(len(res.parnames)):
            ch.add_label(0.2,0.88-i*0.02,
                         '{0} = {1:8.5g}'.format(res.parnames[i],
                          res.parvals[i]), ["coordsys", ch.FRAME_NORM])

        ch.print_window(os.path.join(self.plotpath, 
                                     self.name + filename + self.printfiletype),
                        {"clobber": True, 'dpi': 300})  




    def append_to_csv(self, filename, conf=False):
        if not os.path.isfile(filename):
            with open(filename, "w") as outfile:
                outfile.write('Region, detector, start, stop, redchi2, const1d, const1d_max, const1d_min, line, wave_in, wave_fit, fwhm, ampl, wave_max, wave_min, fwhm_max, fwhm_min, ampl_max, ampl_min\n')
        with open(filename, "a") as outfile:
            for line in (self.H2lines + self.nonH2lines):
                if conf:
                    outfile.write('{0},{1} {2},{3},{4},{5},{6},{7},{8},{9},{10}, {11}, {12}, {13}, {14}, {15}, {16}, {17}, {18}, {19}\n'.format(self.name, self.FN, self.detector, self.start, self.stop, self.redchi2, self.const, self.const_max, self.const_min, line['name'], line['wave'], line['pos'], line['fwhm'], line['ampl'], line['pos_max'], line['pos_min'], line['fwhm_max'], line['fwhm_min'], line['ampl_max'], line['ampl_min']))

                else:
                    outfile.write('{0},{1} {2},{3},{4},{5},{6},{7},{8},{9},{10},{11},{12},{13}\n'.format(self.name, self.FN, self.detector, self.start, self.stop, self.redchi2, self.const, 'constmin', 'constmax', line['name'], line['wave'], line['pos'], line['fwhm'], line['ampl']))


