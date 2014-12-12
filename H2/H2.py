import re

Abgrall93file = '/data/hguenther/misc/H2/MOLAT/fichiers_BCan94'

def read_Abgrall93(filename=None):
    '''Read the Abgrall 1993 line list of Lyman band transitions

    Parameters
    ----------
    filename : string or None
        If filename is ``None`` use the module level variable 
        ``Abgrall93file``

    Returns
    -------
    out : dictionary
        the keys in the dictionary are the spectroscopic code such as
        ``P(2) 4-5``, the data entry is the wavelength in Ang
    '''
    if filename is None:
        filename = Abgrall93file
    data = np.loadtxt(filename, skiprows=3, dtype=[('vu','i4'), ('Ju','i4'), ('vl','i4'), ('Jl','i4'), ('A','f4'), ('wavenumber','f8')])
    out = {}
    for d in data:
        code = H2numbers2string()
        out[code] = 1./d['wavenumber'] * 1e8
    return out


delta_J_letter2number = {'R': +1, 'P': -1, 'Q': 0}
delta_J_number2letter = dict((v,k) for k, v in delta_J_letter2number.iteritems())

def numbers2string(Ju, Jl, vu, vl):
    '''convert quantum numbers J_up, J_low, v_up, v_low to string like ``3-4 P(2)``
    
    Parameters
    ----------
    J_u : int
        rotational quantum number of upper level
    J_l : int
        rotational quantum number of lower level
    vu : int
        vibrational quantum number of upper level
    vl : int
        vibrational quantum number of lower level

    Returns
    -------
    code : string
        string representation of a ro-vibrational H2 line
    '''
    return '{2}-{3} {0}({1})'.format(delta_J_number2letter[Ju-Jl], Jl, vu, vl)


def string2numbers(code):
    '''convert string like ``3-4 P(2)`` to quantum numbers J_up, J_low, v_up, v_low    
    Parameters
    ----------
    code : string
        string representation of a ro-vibrational H2 line

    Returns
    -------
    J_u : int
        rotational quantum number of upper level
    J_l : int
        rotational quantum number of lower level
    vu : int
        vibrational quantum number of upper level
    vl : int
        vibrational quantum number of lower level
    '''
    H2line = '[H2 ]*(?P<vu>[0-9]+)-(?P<vl>[0-9]+)\s(?P<dJ>[PR])\((?P<Jl>[0-9]+)\)'
    H2_line_order_reversed = '[H2 ]*(?P<dJ>[PR])\((?P<Jl>[0-9]+)\)\s(?P<vu>[0-9]+)-(?P<vl>[0-9]+)'
    m = re.match(H2line, code)
    m1 = re.match(H2_line_order_reversed, code)
    m = m if m is not None else m1
    return int(m.group('Jl')) + delta_J_letter2number[m.group('dJ')], int(m.group('Jl')), int(m.group('vu')), int(m.group('vl'))
