"""
Module with helpers to manage Sherpa model parameters
"""

import os


def copy_pars(oldcomp, newcomp, sametype = True):
    """copy parameters from one component to another
    
    Both components need to be of the same type, e.g. both are gauss1d models
    This routine then copies `val`, `max`, `min`, `frozen` and `link` values.
    
    Example:
    >>> from sherpa.ui import *
    >>> set_model(gauss1d.g1 + gauss1d.g2)
    >>> g1.pos.min = 0.
    >>> copy_pars(g1, g2)
    
    Parameters
    ----------
    oldcomp : Sherpa model instance
        Component with original values. 
        
        .. note:: Instance vs. string name
           This routine expects the model instance objects at input, not just their
           string name (`'g1'`). If you only have their string names, you can do::

               >>> copy_pars(get_model_component('g1'), get_model_component('g2')

    newcomp : Sherpa model instance
        Values of this component will be set
    sametype : bool
        If `False` this passes the check that both model instances have to be of the same
        type. Be careful! Currently, the parameters are copied in the order they appear in
        `oldcomp.pars` and `newcomp.pars`, even if both components have different names!
        
    TBD: replace get_model_component(oldcomp).pars with some way that iterates over names, so that parameters can be copied between two line types, even if pos is once the first and once the second parameter. 
    """
    if sametype:
        if not (type(oldcomp) == type(newcomp)):
            raise TypeError('Old and new model component must be of same type')
    #
    for parold, parnew in zip(oldcomp.pars, newcomp.pars):
        # min cannot be above max.
        # set to -+inf to avoid problems with previously set pars
        setattr(parnew, "min", getattr(parnew, "hard_min"))
        setattr(parnew, "max", getattr(parnew, "hard_max"))
        for elem in ["min", "max", "val", "frozen", "link"]:
            setattr(parnew, elem, getattr(parold, elem))


def set_pars_from_table(table, modelcol):
    '''Set parameters for sherpa model from the values in a table

    This function goes row by row through a table.
    In each row there should be a model


    Parameters
    ----------
    table : masked astropy.table.Table
        
    modelcol : string
        Name of column in ``table`` that holds the model instances.
    '''
    for row in table:
        model = row[modelcol]
        if model is not None:
            for par in model.pars:
                for elem in ["min", "max", "val", "frozen", "link"]:
                    colname = '{0}_{1}'.format(par.name, elem)
                    # Allow fwhm instead of fwhm.par
                    if elem=='val' and not (colname in row):
                        colname = par.name
                    if (colname in row) and not (row[colname].mask):
                        setattr(par, elem, row[colname]



    





