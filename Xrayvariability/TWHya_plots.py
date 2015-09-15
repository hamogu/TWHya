# -*- coding: utf-8 -*-

import numpy as np
import matplotlib.pyplot as plt

def around(x, x0, delta = 0.02):
    '''return an index array of all values in between x0-delta and x0+delta

    :param x: array of values
    :param x0: value at middle of interval
    :keyword delta: half-width of interval 
    '''
    return (x >= (x0 - delta)) & (x <= (x0 + delta))

def split_bins(lam_1, *args):
    bins = np.zeros(len(args))
    for i, arg in enumerate(args):
        bins[i] = (arg & lam_1).sum(dtype = '<f4')
    return bins


def plot_rat(lam_1, lam_2, *args, **kwargs):
    rat = np.zeros(len(args))
    err = np.zeros((2,len(args)))
    for i, arg in enumerate(args):
        x = (arg & lam_1).sum(dtype = '<f4')
        y = (arg & lam_2).sum(dtype = '<f4')
        rat[i] = x/y
        err[:, i] = error(x,y)
        print i, x/y, ' +- ', error(x,y)
    lastplot = plt.errorbar(range(len(args)), rat, yerr = np.abs(err), fmt = '-o', **kwargs)
    plt.xlim(-0.5, len(args)-0.5)
    plt.xlabel('bin')
    plt.ylabel('line ratio')
    return rat, err

def plot_ratrat(lam_r, lam_i, lam_f, *args, **kwargs):
    r = np.zeros(len(args))
    i = np.zeros(len(args))
    f = np.zeros(len(args))
    for j, arg in enumerate(args):
        r[j] = (arg & lam_r).sum(dtype = '<f4')
        i[j] = (arg & lam_i).sum(dtype = '<f4')
        f[j] = (arg & lam_f).sum(dtype = '<f4')
        #print r,i,f
        #print j,': ', str((f/i)/((f+i)/r)), ' +/- ', error(error(f,i),error(f+i,r))
    lastplot = plt.errorbar(f/i,(f+i)/r, xerr = np.abs(error(f,i)), yerr = np.abs(error(f+i, r)), fmt = 'o', **kwargs)
    plt.xlabel('f/i')
    plt.ylabel(r'$\frac{f+i}{r}$')
    return r,i,f

def plot_lineline(lam_x, lam_y, *args, **kwargs):
    x = np.zeros(len(args))
    y = np.zeros(len(args))
    for j, arg in enumerate(args):
        x[j] = (arg & lam_x).sum(dtype = '<f4')
        y[j] = (arg & lam_y).sum(dtype = '<f4')
        #print r,i,f
        #print j,': ', str((f/i)/((f+i)/r)), ' +/- ', error(error(f,i),error(f+i,r))
    lastplot = plt.errorbar(x,y, xerr = np.sqrt(x), yerr = np.sqrt(y), fmt = 'o', **kwargs)
    return x,y

def gauss_error(a, b):
    return np.sqrt(a + a**2./b)/b

def poisson_error(a,b,a_eff, b_eff, n=10000, conf = 0.67):
    '''a[0] - counts
    a[1] - effective area
    a/b
    '''
    x = np.random.poisson(a,n)/a_eff
    y = np.random.poisson(b,n)/b_eff
    rat = 1.*x/y
    rat.sort()
    return rat[np.int((1.-conf)/2.*n)]-(1.*(a/a_eff)/(b/b_eff)), rat[np.int((1.-(1.-conf)/2.)*n)]-(1.*(a/a_eff)/(b/b_eff))

vpoisson_error = np.vectorize(poisson_error)
error = lambda a,b,a_eff, b_eff: np.array(vpoisson_error(a,b,a_eff, b_eff))

def poisson_error2(a,b,c,a_eff, b_eff, c_eff, n=10000, conf = 0.67):
    x = np.random.poisson(a,n)/a_eff
    y = np.random.poisson(b,n)/b_eff
    z = np.random.poisson(c,n)/c_eff
    rat = 1.*(x+y)/z
    rat.sort()
    return rat[np.int((1.-conf)/2.*n)]-(1.*(a/a_eff+b/b_eff)/(c/c_eff)), rat[np.int((1.-(1.-conf)/2.)*n)]-(1.*(a/a_eff+b/b_eff)/(c/c_eff))

vpoisson_error2 = np.vectorize(poisson_error2)
error2 = lambda a,b,c,a_eff, b_eff, c_eff: np.array(vpoisson_error2(a,b,c,a_eff, b_eff, c_eff))



def plot_ratios(lam, *args, **kwargs): 
    if ('r' in lam) & ('beta' in lam):
        print 'He alpha / He beta'
        rbeta = plot_rat(lam['r'], lam['beta'], *args, color ='r', label = r'He$\alpha$/He$\beta$', **kwargs)
    if ('f' in lam) & ('i' in lam):
        print 'f/i'
        fi = plot_rat(lam['f'], lam['i'], *args, color ='b', label = 'f/i', **kwargs)
    if ('r' in lam) & ('i' in lam) & ('f' in lam):
        print '(f+i)/r'
        fir = plot_rat(lam['f'] | lam['i'], lam['r'], *args, color ='g', label = '(f+i)/r', **kwargs)
    if ('Ha' in lam) & ('r' in lam):
        print 'H alpha / He alpha'
        HHe = plot_rat(lam['Ha'], lam['r'], *args, color ='k', label = r'H$\alpha$/He$\alpha$', **kwargs)
    plt.legend()   


