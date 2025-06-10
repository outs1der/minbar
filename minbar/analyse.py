# Classes and routines to analyse X-ray observational data

import numpy as np
from scipy.signal import correlate
from scipy.optimize import curve_fit
import scipy.special
import os, sys
import pickle
import minbar
import astropy.units as u
import matplotlib.pyplot as plt

def idlwhere(bool):
    """
    Routine to replace the IDL where routine, which nicely returns the number of elements, and the complement
    :param bool: boolean array derived from some expression
    Usage:
    a = np.arange(20)
    good, n, bad = idlwhere(a > 15)
    print (a[good], n)
    [16 17 18 19] 4
    """

    i = np.where(bool)[0]
    not_i = np.where(~bool)[0]

    return i, len(i), not_i


def tmatch(t1, t2, exclusive=False, thresh=None):
    '''
    Function to return the indices of array t2 that are the closest to the
    values of t1. Adapted from tmatch.pro

    The exclusive option may not always get the optimal list, since we can
    do replacements where the original best match is set as unmatched (even
    when it might be a better match to one of events in t2 which was
    matched to a different event in t1)
    
    :param thresh: match only if the offset is less than this value
    :param exclusive: when True, only match each entry in t2 once
    
    :return: indices for matches, array of offsets, and array of unmatched elements in t2
    '''

    n = len(t1)
    offset = np.zeros(n)
    xid = [None] * n

    matched = np.zeros(len(t2))
    
    for i in range(n):
        # Loop over each element of t1 and rank the nearest elements of t2
        _imin = np.argsort(np.abs(t1[i]-t2))
        _offset = np.abs(t1[i]-t2[_imin])

        if (not exclusive) & (_offset[0] < (np.inf if thresh is None else thresh)):

            # We have a match, and don't care about exclusivity

            xid[i] = _imin[0]		# pointer to t1
            offset[i] = _offset[0]	# offset
            matched[_imin[0]] += 1	# count of matches in t2

        elif exclusive:

            # work your way down the list of matches

            for j in np.where(_offset < (np.inf if thresh is None else thresh))[0]:

                if (matched[_imin[j]] > 0):
		    # If we've already matched with this event, BUT the
		    # current match is better, we can swap

                    _prev = np.where(xid == _imin[j])[0]
                    assert len(_prev) == 1
                    if (_offset[j] < offset[_prev[0]]):
                        xid[_prev[0]], offset[_prev[0]] = None, 0.0
                        xid[i] = _imin[j]
                        offset[i] = _offset[j]	# offset
                        break

                else:

                    xid[i] = _imin[j]
                    offset[i] = _offset[j]	# offset
                    matched[_imin[j]] += 1
                    break


    unmatched = np.where(matched == 0)[0]
            
    return xid, offset, unmatched


def gtiseg(time, maxgap=None, mingti=None):
    """
    NAME:
      GTISEG

    AUTHOR:
      Craig B. Markwardt, NASA/GSFC Code 662, Greenbelt, MD 20770
      craigm@lheamail.gsfc.nasa.gov
      Python conversion by Duncan Galloway, Monash University
      duncan.galloway@monash.edu

    PURPOSE:
      Convert a list of times to a set of Good Time Intervals (GTIs)

    CALLING SEQUENCE:
      GTI = GTISEG(TIMES, MAXGAP=MAXGAP, MINGTI=MINGTI)

    DESCRIPTION:

      The function GTISEG accepts an array of times and converts
      adjacent data into good time intervals (GTIs).

      Elements of the array are clustered into intervals based on the
      gaps between times.  If the gaps are small enough then the times
      are grouped into a single interval.  If a gap exceeds MAXGAP, then
      an interruption occurs and at least two intervals are formed.
      Thus, the keyword parameter MAXGAP essentially determines how many
      and where the intervals will be formed.

      If the time samples are regularly spaced -- aside from gaps --
      then MAXGAP should be set to a number slightly larger than the
      spacing to prevent roundoff errors.  By default MAXGAP is set to
      the difference between the first and second samples.

      For GTISEG, the samples do not need to be regularly spaced, but
      they *must* be given in ascending order.  Arrays can be sorted
      with the SORT function.  The primary difference between GTISEG and
      MASK2GTI is that MASK2GTI assumes the time samples are regularly
      spaced while GTISEG does not.  Also, MASK2GTI allows intervals to
      be enlarged or shrunk.

      It should be noted that this function is not constrained to
      operation only on time arrays.  It should work on any
      one-dimensional quantity with intervals.

    INPUTS:

      TIME - an array of times in ascending order.


    KEYWORDS:

      MAXGAP - a scalar, the maximum gap between time samples before a
               new interval is created.  Samples with gaps smaller than
               this value are grouped into a single GTI.
               Default: TIME(1) - TIME(0)

      MINGTI - the smallest possible GTI.  Any interval smaller than
               MINGTI is discarded.
               Default: 0   (all intervals are accepted)

      The following keywords are not yet implemented in the Python version

      COUNT - upon return, the number of resulting intervals.  A value
              of zero indicates no good time intervals.

      INDICES - upon return, a 2xCOUNT array of integers which give the
                indices of samples which lie within each interval.  The
                times TIME(INDICES(0,i) : INDICES(1,i)) fall within the
                ith interval.

    RETURNS:

      A new GTI array containing the indices corresponding to the
      enlarged or shrunken intervals (equivalent to the INDICES
      keyword in the IDL version). The array is 2xCOUNT where COUNT
      is the number of resulting intervals.  GTI[i,:] represents the
      indices for the start and stop times of interval number i.
      The intervals are non-overlapping and time-ordered.

      If COUNT is zero then the returned array is a scalar value of
      zero, indicating no good intervals were found.

    SEE ALSO:

      MASK2GTI, GTITRIM, GTIMERGE, GTIWHERE

    MODIFICATION HISTORY:
      Written, CM, 1999-2001
      Documented, CM, Apr 2001
      MINGTI now works as documented, in that segments *equal* to MINGTI
        are now accepted, CM, 30 Oct 2007
      MINGTI now also affects INDICES, CM, 03 Mar 2008
      Handle case of scalar input, CM, 2016-04-15
      Python conversion, dkg, 2019 Oct

     $Id: gtiseg.pro,v 1.8 2016/05/19 16:12:02 cmarkwar Exp $

   -
    Copyright (C) 1999-2001, 2007, 2008, 2016, Craig Markwardt
    This software is provided as is without any warranty whatsoever.
    Permission to use, copy, modify, and distribute modified or
    unmodified copies is granted, provided this copyright and disclaimer
    are included unchanged.
   -
    """

    count = 0

    # No need for usage message, if we have the docstring
    #    if n_params() EQ 0 then begin
    #        message, 'GTI = GTISEG(TIME, [COUNT=,] [INDICES=,] [MAXGAP=,] [MINGTI=])',/info
    #        return, 0
    #    endif
    if maxgap == None:
        maxgap = time[1] - time[0]
    nt = len(time)
    tdiff = [0.0]
    ntseg = 0
    if nt > 1:
        tdiff = time[1:] - time[:-1]
        whdiff = [-1] + list(np.where(tdiff > maxgap)[0]) + [nt - 1]
        ntseg = len(whdiff) - 2
        # print (tdiff*86400., maxgap*86400., whdiff, ntseg) # for debugging
    if (nt <= 1) | (ntseg == 0):
        whdiff = [-1, nt - 1]
    ntseg = ntseg + 1

    mintdiff = min(tdiff) > 0

    # indices = reform(lonarr(2, ntseg), 2, ntseg, /overwrite)
    indices = np.zeros((ntseg, 2), dtype=int)
    # print (ntseg, np.shape(indices)) # more debugging
    # indices(0,*) = whdiff(0:ntseg-1)+1
    # indices(1,*) = whdiff(1:ntseg)
    indices[:, 0] = np.array(whdiff[0:-1]) + 1
    indices[:, 1] = np.array(whdiff[1:])

    #     tgti = reform(make_array(2, ntseg, value=time(0)*0), 2, ntseg, /overwrite)
    #     tgti(0,*) = time(indices(0,*))
    #     tgti(1,*) = time(indices(1,*)) + mintdiff

    if mingti is not None:
        print('** ERROR ** not yet implemented')
        # wh = where(tgti(1,*)-tgti(0,*) >= mingti(0), ntseg)
        # if ntseg GT 0:
        # tgti = reform(tgti(*,wh),2,ntseg)
        # indices = reform(indices(*,wh),2,ntseg)
        # else:
        # tgti = -1L
        # indices = -1L

    count = ntseg

    return indices


def dt_nonzero(time):
    """
    Routine to infer the timestep from an array. Should be robust against duplicate
    values"""

    n = len(time)
    _dt = time[1:] - time[:-1]
    ma = np.ma.masked_equal(_dt, 0.0, copy=False)

    return ma.min()


def fit_burst_times(time=None, entry_obs=None, adjust=[], swt_thresh=0.5*u.hr,
    plot=True, _tunit='hr'):
    """
    Routine to generate a trial model to a set of burst times, and fit to
    get the average recurrence time. Accepts a set of burst times, or a
    set of MINBAR observation IDs from which the bursts are drawn
    
    Adapted from the Analysis of high-state bursts.ipynb notebook
    
    :param time: burst times. Can be astropy column or an array of times (days assumed)
    :param entry_obs: MINBAR observation entry numbers, from which the burst times will be extracted
    :param adjust: list of tuples for adjustments of the default cycle count. E.g. [(3,+1), (5, -2)] will first increment the cycle count for the 3rd and later bursts by 1; and then decrement the 5th and later bursts by 2
    :param swt_thresh: short waiting time threshold, below which burst separations will be ignored for the pursposes of estimating the recurrence time
    :param plot: set to False to skip the diagnostic plot
      
    :return: array of cycle counts, tuple giving best-fit recurrence time and uncertainty, and rms for fit
    """
    
    if (time is None) & (entry_obs is None):
        minbar.logger.error('need to specify burst times or observation IDs')
        return
    elif (time is None) & (entry_obs is not None):
        # when we incorporate this into one of the classes, can replace the below with 
        # self or self.bursts, respectively
        b = minbar.Bursts()
        b.select(entry_obs,'entry_obs')
        b.show(['entry','instr','obsid','bfluen','trec','tdel','alpha'])
        
        time = b['time']
        
    time0 = min(time)
    try:
        unit = time.unit
    except:
        minbar.logger.warning("time has no units, assuming days")
        unit = u.d
        time *= unit

    # print (unit, type(unit))
    tfac = (1.*unit).to(_tunit).value
        
    dt = (time[1:]-time[:-1])
    dt0 = min(dt[dt > swt_thresh.to(unit)])

    # generate the trial cycle count numbers
    
    c = np.round(((time-min(time))/dt0).value)

    # adjust the cycle counts, e.g. where there is an obvious phase wrap, because the
    # initial trial value was too long or too short
    
    for delc in adjust:
        c[delc[0]:] += delc[1]
        
    # do the fit; see
    # https://stackoverflow.com/questions/27283393/linear-regression-slope-error-in-numpy

    coef, cov = np.polyfit(c, time, 1, cov=True)
    # print (coef, time0)
    tdel, e_tdel = coef[0]*tfac*u.Unit(_tunit), np.sqrt(np.diag(cov)[0])*tfac*u.Unit(_tunit)
    minbar.logger.info ('average recurrence time = {:.4f}+/-{:.4f}'.format(tdel.value, e_tdel))
    
    # poly1d_fn is now a function which takes in x and returns an estimate for y

    poly1d_fn = np.poly1d(coef) 
    resid = (time-poly1d_fn(c)*unit).to(_tunit)

    if plot:

        # set up a plot showing the fit and the residuals
        # example here https://matplotlib.org/stable/gallery/text_labels_and_annotations/label_subplots.html#sphx-glr-gallery-text-labels-and-annotations-label-subplots-py

        fig, axs = plt.subplot_mosaic([['a)'], ['a)'], ['b)']],
                                  layout='constrained')

        axs['a)'].plot(c,time, 'yo', c, poly1d_fn(c), '--k') #'--k'=black dashed line, 'yo' = yellow circle marker
        axs['a)'].set_ylabel('Burst time (MJD)')
        axs['b)'].plot(c,resid, 'yo')
        axs['b)'].axhline(0,0,ls='--',c='k')

        axs['b)'].set_xlabel('Burst number')
        axs['b)'].set_ylabel('Residual ({})'.format(_tunit))

        plt.show()
    
    return c, (tdel, e_tdel), np.sqrt(np.mean(resid**2))


def burst_template(x, a):
    """
    Generates a template burst within a time series given by the first
    variable passed, with the specified profile. Includes calculation of
    partial derivatives for use with curvefit. Parameters are:

    | a[0] = baseline count rate
    | a[1] = start time of burst
    | a[2] = rise time of burst
    | a[3] = peak intensity of burst
    | a[4] = exponential decay time

    Copied from jemx/search-new.pro, 23.10.14
    """

    n = len(x)
    dt = dt_nonzero(x)
    f = np.zeros(n)
    rise = max([dt, a[2]])  # Meaningless to have a rise smaller than this
    pder = np.zeros([n, 5])
    pder[:, 0] = 1.
    g = np.where(np.logical_and(x >= a[1], x - a[1] < rise))[0]  # Rise segment
    if len(g) > 0:
        f[g] = a[3] * (x[g] - a[1]) / rise
        #        if n_params() ge 4 then begin
        pder[g, 1] = -a[3] / rise
        pder[g, 2] = -f[g] / rise
        pder[g, 3] = f[g] / a[3]
    #        endif
    #    endif
    g = np.where(x >= a[1] + rise)[0]  # Tail segment
    # print(g)

    if len(g) > 0:
        f[g] = a[3] * np.exp(-(x[g] - (a[1] + rise)) / a[4])

        #        if n_params() ge 4 then begin
        #        pder[g,1:4] = [[f[g]/a[3]],[f[g]/a[2]],[f[g]*(x[g]-a[1])/a[3]**2]]
        #        pder[g,1:4] = np.vstack((f[g]/a[3],f[g]/a[2],f[g]*(x[g]-a[1])/a[3]**2)).T
        pder[g, 1] = f[g] / a[4]
        pder[g, 2] = pder[g, 1]
        pder[g, 3] = f[g] / a[3]
        pder[g, 4] = f[g] * (x[g] - (a[1] + rise)) / a[4] ** 2
    #        endif

    f += a[0]

    # It helps for the fit values if we limit the rise here, but not for the
    # errors

    a[2] = rise

    return f, pder


def burst_template_pow(x, base, start, rise, peak_rate, peak_dur, gamma):
    """
    Generates a model burst within a time series given by the first
    variable passed, with the specified profile. This version replaces
    the exponential decay in the tail with a power law

    Doesn't quite work yet - dkg 2024 Jun

    Copied from burst_template 12.6.2024

    :param x: time bins [s]
    :param base: baseline count rate (background), counts/s or counts/cm^2/s
    :param start: start time of burst, same units as x
    :param rise: rise time of burst, same units as x
    :param peak_rate: peak intensity of burst, counts/s or counts/cm^2/s
    :param peak_dur: duration of peak intensity, same units as x
    :param gamma: power law index for decay portion

    :return: array with burst model evaluated at x values
    """

    n = len(x)
    dt = dt_nonzero(x)
    f = np.zeros(n)
    rise = max([dt, rise])  # Meaningless to have a rise smaller than this
    peak_dur = max([dt, peak_dur])

    g = np.where(np.logical_and(x >= start, x < start + rise))[0]  # Rise segment
    if len(g) > 0:
        f[g] = peak_rate * (x[g] - start) / rise

    g = np.where(np.logical_and(x >= start + rise, x <= start + rise + peak_dur))[0]  # Peak segment
    if len(g) > 0:
        f[g] = peak_rate

    g = np.where(x > start + rise + peak_dur)[0]  # Tail segment
    if len(g) > 0:
        # with the convention below, gamma is positive, for consistency with edt
        f[g] = peak_rate * (dt/(x[g] - (start + rise + peak_dur)))**gamma

    return f+base #, pder


def burst_template_flat(x, base, start, rise, peak_rate, peak_dur, edt):
    """
    Generates a model burst within a time series given by the first
    variable passed, with the specified profile. Rather than the sharp peak of
    burst_template, this profile has a peak that lasts for at least one bin,
    likely better suited to realistic bursts and particularly PRE bursts

    Copied from burst_template 12.6.2024

    :param x: time bins [s]
    :param base: baseline count rate (background), counts/s or counts/cm^2/s
    :param start: start time of burst, same units as x
    :param rise: rise time of burst, same units as x
    :param peak_rate: peak intensity of burst, counts/s or counts/cm^2/s
    :param peak_dur: duration of peak intensity, same units as x
    :param edt: exponential decay time, same units as x

    :return: array with burst model evaluated at x values
    """

    n = len(x)
    dt = dt_nonzero(x)
    f = np.zeros(n)
    rise = max([dt, rise])  # Meaningless to have a rise smaller than this
    peak_dur = max([dt, peak_dur])

    g = np.where(np.logical_and(x >= start, x < start + rise))[0]  # Rise segment
    if len(g) > 0:
        f[g] = peak_rate * (x[g] - start) / rise

    g = np.where(np.logical_and(x >= start + rise, x < start + rise + peak_dur))[0]  # Peak segment
    if len(g) > 0:
        f[g] = peak_rate

    g = np.where(x >= start + rise + peak_dur)[0]  # Tail segment
    if len(g) > 0:
        f[g] = peak_rate * np.exp(-(x[g] - (start + rise + peak_dur)) / edt)

    return f + base  # , pder


def findburst(_time, _rate, _error=None, tunit=None,
    rise=None, edt=20., refine=False, fit=False, fit_func=burst_template_flat,
    # sig=sig, param=param, error=sigma, resid=resid,
    plot=False, # xrange=xrange,
    min_burst_sep=3., sigthresh=5., dbg=False):
    """
    $Id: findburst.pro,v 1.1 2011/06/29 12:54:00 duncan Exp duncan $

    This routine is designed to find one or more bursts in noisy data.
    Extracted from search.pro v1.3, 29/6/2011

    This version moved from minbar/jemx, and simplified 23.10.2014

    :param _time: array with burst lightcurve time bins
    :param _rate: array with burst rate [counts/sec, or counts/sec/cm^2]
    :param _error: array with error on burst rate
    :param tunit: string giving time units, 's' or 'd'
    :param rise: (optional) template burst rise time [sec]
    :param edt: (optional) template burst exponentional decay timescale [sec]
    :param refine: boolean, set to True to 'refine' the burst time based on a fit (not yet implemented)
    :param fit: boolean, set to True to also fit the burst profile
    :param fit_func: function to fit the burst profiles
    :param min_burst_sep: minimum burst separation [min]
    :param sigthresh: significance threshold for burst detection (no. of sigmas)
    :param dbg: boolean, set to True to display some diagnostic information

    :return: <t_cand> list of candidate burst times [same units as _time], and (if fit=True) the fit parameters and sigmas

    | sig estimated significance of burst detection
    | param list of burst parameters: peak rate [count/s], base rate [count/s],
    |         fitted start time [MJD], rise [s], fitted peak rate [count/s],
    |         exponential decay time [s]
    | sigma error on fit parameters (from curvefit)

    Example usage:

    | import minbar
    | from astropy.io import fits
    |
    | file = '../xmm/data/1RXSJ180408.9-342058/0741620101/2to7good.fits'
    | hdul = fits.open(file)
    | data = hdul[1].data
    | test = minbar.findburst(data['TIME'], data['RATE'], data['ERROR'])
    | print(test)
    | [5.42058957e+08 5.42067368e+08 5.42075296e+08 5.42083081e+08 5.42090903e+08]

    Routines called:

    | burst_template_flat
    | gtiseg
    | dt_nozero
    """

    # Some parameters passed to the IDL routine via common statements
#  common findburst,sigthresh,min_burst_sep
#  common jemx,mjdref    ; MJDREF not currently used

    pre_burst = 32.0     # interval prior to the burst to include in the fit

    if rise == None:
        rise=edt/4.     # Default burst rise [sec]

    if refine: # or fit:
        minbar.logger.error("refining/fitting of template not yet fully implemented")
        return None

# Make copies of the input arrays to avoid unintentionally modifying them

    time = _time
    rate = _rate
    n = len(time)
    dt = dt_nonzero(time)         # Nominal dt

# The error is required for the fitting; we should have an option where it
# can be estimated from the lightcurve, to make the code easier to use

    if _error is None:
        minbar.logger.error("errors are required")
        return None
    else:
        error=_error

# Also search for non-finite burst rate measurements

    fin = np.where(np.isfinite(rate) & np.isfinite(error))[0]
    if len(fin) < n:
        minbar.logger.warning("non-finite rate measurements, excising from analysis")
        time = time[fin]
        rate = rate[fin]
        error = error[fin]

# Try to establish what the time unit here is, and set the scaling factors
# and timescales appropriately

#  sd=864d2 & md=144d1
    if tunit is None:
        if dt >= 0.1:
            tunit='s'
        else:
            tunit='d'
        minbar.logger.warning("setting tunit='{}', set explicitly if incorrect".format(tunit))
    if tunit == 's':
        sd=1.0
    elif tunit == 'd':
        sd=86400.0
    else:
        minbar.logger.error('unrecognised time unit')
        return
    md=sd/60.

# Define the fitting window here

    pb=pre_burst/sd            # interval prior to the burst to include in the
                          #   fitting window
    fwin = np.array([-pb,max([min_burst_sep/md,edt/sd])-pb])
    fwin_ind = np.where((time-time[0] > fwin[0]) & (time-time[0] < fwin[1]))[0]
    if dbg:
        print("{}: defined fitting window as {} relative to the burst time".format(
            'findburst', fwin*sd))

# Some screening here. This gap selection used for the April and June 2011
# searches [of JEM-X data?], but subsequently Jerome suggested it was too restrictive...
# although in actual fact it only results in one exception

    g = np.where((dt*sd >= 16.) & (rate > 4000.))[0]
    if len(g) > 0:
        minbar.logger.error ("exception in lightcurve")
        return None

# First define the burst template

    m = np.median(rate)
    p = max(rate)-m
#  t0=0.5*(max(time)+min(time)) ; Trial burst is in the middle of the segment
    # t0=time[int(n/2)]
    # trial, pder = burst_template(time,[m,t0,rise/sd,p,edt/sd])
    t0 = time[0]+pb
    # trial, pder = burst_template(time[fwin_ind],[m,t0,rise/sd,p,edt/sd])
    # set peak duration to the same value as the rise
    trial = burst_template_flat(time[fwin_ind], m, t0, rise/sd, p, rise/sd, edt/sd)
    # we store the template values as the guess for curve_fit
    # TODO modify these guesses etc. based on the specific form of the fit function
    p0 = [m, t0, rise/sd, p, rise/sd, edt/sd]
    if dbg & (not fit):
        plt.plot(time[fwin_ind],trial)

# Now cross-correlate the template, sliding the window over the entire lightcurve

    c = correlate(rate, trial, mode='same')
    # plt.plot(time,c)

# Estimate the significance here by shuffling the rate
# measurements and determining the statistics of the resulting correlation
# Now uses more efficient Fisher-Yates shuffle

#    seed=systime(/sec)
    o = np.arange(n)
    for i in range(n):
        r = np.random.choice(n)
        tmp = o[i]
        o[i] = o[r]
        o[r] = tmp
    rshuf = rate[o]
    ctest = correlate(rshuf, trial, mode='same')
    ctest_m = np.mean(ctest)
    ctest_s = np.std(ctest)

# Find all candidates for bursts, from peaks of the cross-correlation
# function, and use gtiseg to identify each contiguous region exceeding
# the threshold
# The cexc array gives the indices of the cross-correlation array that are
# in excess of the significance limit.

#  cmax=max(c,imax)
#  t_cand=time[imax]
#  sig=(cmax-ctest_m)/ctest_s
    csig=(c-ctest_m)/ctest_s
    cexc=np.where(csig > sigthresh)[0]
    if len(cexc) == 0:

        # No bursts
        return None

    # Define an array giving the start and end indices for each region
    # The 5% padding below avoids issues with rounding errors

    # print (1.05*min_burst_sep/md,min_burst_sep)
    t_cand_ind=gtiseg(time[cexc], maxgap=1.05*min_burst_sep/md)#,indices=t_cand_ind)
    t_cand_ind=cexc[t_cand_ind]
    # print (t_cand_ind,np.shape(t_cand_ind))
    # print (len(time),len(c))

    # And finally assemble the list of start times, best estimates as the 
    # time at which the cross-correlation is maximum, less the window
    # offset and the rise time

    t_cand = [time[x[0]+np.argmax(c[x[0]:x[1]])] for x in t_cand_ind]+fwin[0]-rise/sd
    if not fit:
        # return time[t_cand_ind[:,0]]-fwin[0]
        return t_cand
    else:

        # Now set up the result arrays, based on the number of contiguous regions
        # exceeding the search threshold, identified in the previous step

        ncand=len(t_cand_ind)
        param = np.zeros((ncand,8)) # fit parameters plus significance & maximum rate
        result = -1     # originally the fits for each burst; not currently used
        sigma = np.zeros((ncand,7))
        resid = rate.copy()
        for i in range(ncand):

            # Calculate peak intensity in the window over which the significance
            # exceeds the threshold
            # Initial estimate of t_cand is probably sufficiently good

            ipmax = np.argmax(rate[t_cand_ind[i,0]:t_cand_ind[i,1]])+t_cand_ind[i,0]
            # p=max(rate[t_cand_ind[0,i]:t_cand_ind[1,i]],ipmax)
            # ipmax += t_cand_ind[0,i]      # note ipmax is a pointer into time, not win
            # pe=error[ipmax]
            peak_rate, peak_rate_error = rate[ipmax], error[ipmax] # store the maximum rate and error
            # t_cand[i]=time[ipmax]
            sig = max(csig[t_cand_ind[i,0]:t_cand_ind[i,1]])

            # Now determine the window over which to do the fitting. This has a
            # minimum size of fwin

            win = np.where((time-t_cand[i] >= min([fwin[0],time[t_cand_ind[i,0]]-t_cand[i]]))
                & (time-t_cand[i] < max([fwin[1],time[t_cand_ind[i,1]]-t_cand[i]])))[0]

# Start the fitting loop

            # assemble the initial guess
            # _param=[m,t_cand[i],rise/sd,p,rise/sd,edt/sd]
            p0[1] = t_cand[i]

            # if len(win) > len(_param):
            if len(win) > 6:

                # try to do a fit, now just within a window surrounding the burst
                # originally this code used IDL's curvefit, replaced with scipy.optimize's curve_fit, see
                # https://docs.scipy.org/doc/scipy/reference/generated/scipy.optimize.curve_fit.html

                # check for zero errors, and fill with the median, if any are present
                _error = error[win]
                zerr = _error <= 0.
                if np.any(zerr):
                    _error[zerr] = np.median(_error[~zerr])

                # previously we'd set the initial parameters, but now we guide the fit with the
                # bounds parameter
                # _param=[_param[0],t_cand[i],_param[2:]]
                # _result=curvefit(time[win],rate[win],weights[win],_param,_sigma,
                #     function_name='burst_template',iter=1000,status=_s)
                # burst_template_flat parameters are base, start, rise, peak_rate, peak_dur, edt
                _param, pcov, idict, mesg, _s = curve_fit(fit_func, time[win], rate[win], sigma=_error,
                    absolute_sigma=True, full_output=True, p0=p0,
                    bounds=([0,   t_cand[i]-edt/sd, 0,      0.9*(p-m),                 0,      0],
                            [2*m, t_cand[i]+edt/sd, edt/sd, peak_rate+3*peak_rate_error, edt/sd, 100/sd]))
                _sigma = np.sqrt(np.diag(pcov))

# if dbg then begin
#   plot,time[win],rate[win]
#   oplot,time[win],_result,color='00ff00'x
#   a=get_kbrd(1)
# endif

            else:
                _s=-1    # don't have enough points to fit

# Try to refine the detection criteria here by redoing the search with a
# decay value closer to the fit value
# This will also require refining the significance threshold exceedance array,
# AND omitting the segments to which a burst has already been fitted

            if refine & (ncand == 1) & (_s in (1,2,3,4)) & (_param[4] < 1000./sd):
                trial2 = fit_func(time,[_param[0],t0,_param[2:5]])


                c2 = correlate(rate, trial2, mode='same')
                ctest2 = correlate(rshuf, trial2, mode='same')
                cexc2 = np.where((c2-mean(ctest2))/np.std(ctest2) > sigthresh)[0]

#      if cexc2[0] gt -1 then begin
                if len(cexc2) == 1:
                    if (t_cand[ncand-1] >= time[cexc2]) & (t_cand[ncand-1] < time[cexc2+1]):
                        cexc2 = []
                elif len(cexc2) >= 2:
                    intvl=gtiseg(cexc2,maxgap=1.01)#,indices=intvl)
                    for i in range(len(t_cand))-1:
                        if len(cexc2) > 0:
                            for j in range(len(intvl)):
                                if ((t_cand[i] > time[intvl[0,j]]) &
                                    (t_cand[i] <= time[intvl[1,j]-1])):
                                    cexc2 = cexc2[np.where(cexc2 < intvl[0,j] | cexc2 >= intvl[1,j])[0]]
                cmax2=max(c2[cexc2])
                imax2 = np.argmax(c2[cexc2])

                sig2=(cmax2-mean(ctest2))/np.std(ctest2)

# If this works, we re-fit with the better parameter values and use the
# revised template for subsequent fitting

                if sig2 > sig[ncand-1]:
# if dbg then print,'refined ',sig[ncand-1],' -> ',sig2
                    sig[ncand-1]=sig2
                    t_cand[ncand-1]=time[cexc2[imax2]]
                    c=c2
                    cexc=cexc2
                    ncexc=ncexc2
                    trial=trial2
                    _param[1]=t_cand[ncand-1]
                    ctest_m=mean(ctest2)
                    ctest_s=np.std(ctest2)
#;        result=curvefit(time,rate,weights,_param,_sigma, $
#;          function_name='burst_template',iter=1000)
                    if len(win) > len(_param):
                        _result=curvefit(time[win],rate[win],weights[win],_param,_sigma,
                                         function_name='fit_func',iter=1000,status=_s)
# if dbg then begin
#   plot,time[win],rate[win]
#   oplot,time[win],_result,color='00ff00'x
#   a=get_kbrd(1)
# endif
                else:
                    _s=-1      # don't have enough points to fit

# If the fit has failed, or some of the uncertainties are non-finite, then
# the results can't be trusted; so blank them out

            if (_s not in (1,2,3,4)) or (np.any(~np.isfinite(_sigma))) or (_param[4] < 0.0):
                _param=np.zeros(6)
                _param[1]=t_cand[ncand-1]
                _sigma=np.zeros(6)

            # Assemble parameter, result list

            param[i,:]=np.concatenate(([peak_rate],       _param*[1.,1.,sd,1.,sd,sd], [sig]))
            sigma[i,:]=np.concatenate(([peak_rate_error], _sigma*[1.,1.,sd,1.,sd,sd]))

            if dbg:
                # plot the burst data versus the model
                # TODO might want to save this plot as a matter of course for each burst
                plt.plot(time[win],rate[win],'b-')
                plt.plot(time[win], fit_func(time[win], *_param), 'g--')
#   if ncexc gt 0 then oplot,[time[cexc]],[rate[cexc]],color='00ff00'x
#   a=get_kbrd(1)
# endif

# if dbg then begin
#   plot,time,rate
#   if ncexc gt 0 gt -1 then $
#     oplot,[time[cexc]],[rate[cexc]],color='00ff00'x
#   a=get_kbrd(1)
# endif

# if dbg then stop

    return t_cand, param, sigma


def generate_bins(x_b, x_o, nperbin=10, log=True, x_lim=None, x_add=None):
    """
    Function to generate bin boundaries for calc_alpha, to maintain roughly the
    same number of events per bin, with the option of adding a few by hand
    :param x_b: array of values that are to be distributed evenly between bins
    :param x_o: auxiliary array used to determine overall range of bins
    :param nperbin: number of values of x_b to include per bin
    :param log: determines how to calculate the bin boundaries
    :param x_lim: limits on the maximum and minimum bin
    :param x_add: additional bin boundaries to include
    """

    # fn = 'generate_bins'
    fn = sys._getframe().f_code.co_name

    # Depending on whether you intend to plot as log or linear spacing, this
    # option will correctly determine the bin boundaries

    _log = 0
    if log:
        _log = 1

    # Calculate the number of bins

    n = len(x_b)
    si = np.argsort(x_b)
    nbins = int(n / nperbin) + 1

    # If there's only 2 bins and the second one has less than half of the first
    # just divide them evenly instead

    if nbins == 2 and n - nperbin <= nperbin / 2:
        _nperbin = n / 2
    else:
        # Even out the population in the bins
        _nperbin = int(n / nbins)
        if ((n % _nperbin) < _nperbin / 2) and ((n % (_nperbin + 1)) > (n % _nperbin)):
            _nperbin += 1

    if x_add is not None:
        if not isinstance(x_add, list):
            x_add = [x_add]
        n_add = len(x_add)

    if nbins > 1:
        print("{}: Dividing up into {} bins with approx. {} bursts each".format(fn, nbins, nperbin))
    nbins = nbins + 1  # +gap

    # Define the first bin

    bins = np.array(min(np.concatenate((x_b, x_o))))
    inext_old = 0  # index for counting bursts

    if x_lim is not None:
        if not isinstance(x_lim, list):
            x_lim = [x_lim]
        if x_lim[0] > bins:
            inext_old = len(np.where(x_b < x_lim[0])[0])
            print('** WARNING ** {} events below hard bin limit'.format(inext_old))
            bins = np.array(x_lim[0])

    # print (min(x_b),min(x_o), bins)

    # Find the largest gap between observations (bursts?). This block returns
    # _img and img;
    # nflux[_si[_img]] = obs flux at the lower edge of the biggest gap,
    # nflux2[si[img]] = biggest burst flux that is below the gap

    lastgap = 0
    #  if gap then begin
    # ;    if log then intvl=nflux2[si[1:n-1]]/nflux2[si[0:n-2]] else $
    # ;      intvl=nflux2[si[1:n-1]]-nflux2[si[0:n-2]]
    #    _n=n_elements(nflux) & _si=sort(nflux)
    #    if log then intvl=nflux[_si[1:_n-1]]/nflux[_si[0:_n-2]] else $
    #      intvl=nflux[_si[1:_n-1]]-nflux[_si[0:_n-2]]
    #    _img=where(intvl eq max(intvl))	; Biggest gap
    #    _img=_img[0]
    #    img=max(where(nflux_bursts[si] le nflux[_si[_img[0]]]))
    #    img=img[0]
    #  endif

    # Now loop to create the bin boundaries

    while inext_old < n - 1:

        # Set the upper limit of this bin

        inext = min([n - 1, inext_old + nperbin])

        if (inext == n - 1) and (lastgap == 0):

            # Make the last bin with a maximum 1% bigger than the maximum either of the
            # bursts or observations

            #      bins=[bins,1.01*max([nflux,nflux2])]

            # Here's an alternative approach, which may (unfortunately) leave an extra
            # (empty) bin if there are observations at gammas exceeding the maximum burst
            # (i.e. max(nflux)>max(nflux_bursts):

            bins = np.append(bins, 1.01 * max(np.concatenate((x_b, x_o))))
            # print('last bin', bins[-2], '-', bins[-1], max(x_b), max(x_o))
            if (x_lim is not None):
                if len(x_lim) > 1:
                    if (bins[-1] > x_lim[1]):
                        # Add an extra bin here corresponding to the limit
                        bins = np.insert(bins, -1, x_lim[1])

        else:

            # Here we can optionally check for gaps
            # If the biggest gap is within the current bin, try to correct by putting
            # a bin to cover the gap, which will then be empty

            #      if gap then begin
            #        if debug then print,fn,': ',i,lastgap
            # ;        if nflux2[si[img]] ge bins[i] and $
            # ;          nflux2[si[img+1]] le nflux2[si[inext]] and $
            # ;          n_elements(where(nflux2 ge bins[i] and lastgap eq 0 and $
            # ;          nflux2 le nflux2[si[img+1]])) gt nperbin/2 then begin
            # ;;            bins=[bins,(1-log)*0.5*(nflux2[si[img]]+nflux2[si[img+1]]) + $
            # ;;                       log*sqrt(nflux2[si[img]]*nflux2[si[img+1]])] $
            # ;            bins=[bins,nflux2[si[img]]+.01*(nflux2[si[img]]+nflux2[si[img+1]])]
            #        if nflux[_si[_img]] ge bins[i] and $
            # ;          nflux[_si[_img+1]] le nflux2[si[inext]] and $
            #          nflux[_si[_img]] le nflux_bursts[si[inext]] and $
            #          n_elements(where(nflux_bursts ge bins[i] and lastgap eq 0 and $
            #          nflux_bursts le nflux[_si[_img+1]])) gt nperbin/2 then begin
            #            bins=[bins,nflux[_si[_img]]+ $
            #              .01*(nflux[_si[_img]]+nflux[_si[_img+1]])]
            #            lastgap=1
            #            if debug then print,fn,': gap edge bin ',i
            #        endif else if lastgap eq 1 then begin
            # ;          bins=[bins,nflux_bursts[si[img+1]]]
            #          bins=[bins,nflux[_si[_img+1]]]
            #          lastgap=0
            #          if debug then print,fn,': gap other edge bin ',i
            #        endif else $
            #          bins=[bins,(1-log)*0.5*(nflux_bursts[si[inext-1]]+nflux_bursts[si[inext]]) + $
            #                   log*sqrt(nflux_bursts[si[inext-1]]*nflux_bursts[si[inext]])]
            #      endif $	; end block for gap checking
            #
            #      else $

            # Set the upper limit for the next bin

            bins = np.append(bins, (1 - _log) * 0.5 * (x_b[si[inext - 1]] + x_b[si[inext]]) +
                             _log * np.sqrt(x_b[si[inext - 1]] * x_b[si[inext]]))

            # Add any custom bins here

            ninrange = 0
            if x_add is not None:
                inrange, ninrange, dummy = idlwhere((np.array(x_add) > bins[-2]) & (np.array(x_add) < bins[-1]))
                if ninrange > 0:
                    # Add a custom bin here
                    bins[-1] = x_add[inrange[0]]
                    inext = inext_old + len(np.where((x_b >= bins[-2]) & (x_b < bins[-1]))[0])

            # If you've set a hard lower limit, you need to make sure the 2nd (and subsequent)
            # bins are at larger values
            # Don't get this mixed up with the additional bins though

            if (x_add is None or ninrange <= 0) and x_lim is not None:
                assert bins[-1] > x_lim[0]
                if len(x_lim) > 1:
                    # Also impose the upper limit
                    if bins[-1] > x_lim[1]:
                        bins[-1] = min([bins[-1],x_lim[1]])
                        inext = inext_old + len(np.where((x_b >= bins[-2]) & (x_b < bins[-1]))[0])
                        # print ('imposing bin upper limit', inext, inext_old)

        inext_old = inext

    return bins


def calc_alpha(src, bc=None, nperbin=16, limit=None, custom=None,
               filter=None, exclude_short=True, gap=False, log=False, s_z=False, conf=0.68,
               verify=False, debug=False, verbose=True, write=False, clobber=False):
    """
    This function calculates mean burst rates for a given source as well as
    alphas etc. binned as a function of source flux

    Copied from burst/minbar/analysis/calc_alpha.pro and converted to Python

    The original routine returned a structure like this, where the entries in
    square brackets are only when you run on multiple sources:
     [ SRC        STRING Array[19] ]	Source list (excludes "deleted" sources)
      BINS       FLOAT  Array[20]		Boundary values for gamma/S_Z bins
      PBINS      FLOAT  Array[19]		Midpoint of each bin
      PBINE      FLOAT  Array[19]		Half-width of each bin
      N          INT    Array[19]		Number of bursts in each bin
     [ N_SRC      INT    Array[19, 19] ]	Number of bursts per bin, per source
     [ DUR_SRC    FLOAT  Array[19, 19] ]	Duration per source per bin (indexes
                                          in that order!)
      AVRATE     FLOAT     0.120416	Mean rate (/hr) for all bins
      AVRATE_ERR FLOAT     0.00562668	Uncertainty on mean rate
      RATE       FLOAT  Array[19]		Burst rate (/hr) per bin
      RATE_ERR   FLOAT  Array[19]		Uncertainty on rate
      ALPHA      FLOAT  Array[19]		Alpha per bin
      ALPHA_ERR  FLOAT  Array[19]		Uncertainty on alpha
     [ BADBURSTS  LONG   Array[5] ]	Bursts omitted because no gamma/S_Z val
     [ BADOBS     LONG   Array[328] ]	Obs omitted because no gamma/S_Z val

    Example usage:
    res = minbar.calc_alpha('4U 1636-536', nperbin=50, limit=0.03, bc=1.51)
    minbar.binplot(res['bins'], res['rate'], res['rate_err'],xrange=[0.01,0.5],yrange = [0.04,0.8])

    :param src: source name to bin bursts
    :param bc: Bolometric correction for alpha-calculation
    :param nperbin:  Number of bursts per rate/alpha bin
    :param limit: 1 or 2 element list giving the hard limits for binning
    :param gap: Avoid big gaps within bins (not yet implemented)
    :param log: Bin boundary is geometric mean rather than arithmetic mean
    :param s_z: Plot vs. position on color-color diagram instead of flux
    :param exclude_short: Exclude the short bursts (default: yes)
    :param filter: Allows you to filter the observations, e.g. on some other parameter
    :param conf: Confidence level for uncertainties
    :param verify: Check for completeness of the data
    :param debug: Provide debugging information
    :param write: Save the analysis results to a (pickle) file
    :param clobber: Force overwrite of the save file
    """

    fn = sys._getframe().f_code.co_name

    eta = 1e-6
    q = 0.5 * (1. - conf)   # confidence interval for errors
    conf_lim = 0.95            # confidence value for limits
    tdel_thresh = 1800. / 3600.  # [hr] threshold for short-recurrence time bursts

    if s_z:
        print('** WARNING ** binning on S_Z is not tested')
        gap = False  # Gap won't work if we're using colors
        log = False  # Better to do this in linear space

    _filter=''
    if filter is not None:
        print('''
** WARNING ** filtering on gamma produces inexplicably variable results,
              and should probably be avoided''')
        _filter = filter + ','

    # Set bolometric correction for alpha-calculation

    if bc is not None:
        if not isinstance(bc,list):
            bc=[bc]
        if abs(bc[0]-minbar.BOL_CORR[src][1]) > eta:
            minbar.logger.warning('overriding default MINBAR bolometric correction for {}'.format(src))
        _bc = bc[0]
        if len(bc) > 1:
            _bce = bc[1]
        else:
            _bce = 0.0
    else:
        _bc, _bce = minbar.BOL_CORR[src][1:]

    gflag = True  # Flag to tell if we have fedd defined for each source. If
    #   True, then we bin using gammas; if False, we use fluxes

    # ----------------------------------------------------------------------------
    # Load the source table 
    s = minbar.Sources()
    s.name_like(src)

    # Find the F_Edd value for the source
    # this bit is redundant since we're only working with one source at at
    # time (and we're working with gamma, when we're normalising)

    # fedd = None
    # if s.local_files:
    #     if debug:
    #         print('{}: extracting Eddington luminosity via get_FEdd'.format(fn))
    #     test = s['F_Edd']
    #     if test > 0.:
    #         fedd = test

    # Check that all sources have an F_Edd value, and adjust if not

    # if (fedd is None):
    #     print("{}: ** WARNING ** no F_Edd value for {}".format(fn, src))
    #     # return None

    # Get the anisotropy factor

    assert len(s['type']) == 1
    if 'D' in (s['type'].values[0]):
        aniso_fac = minbar.ANISOTROPY['dipper'][1] / minbar.ANISOTROPY['dipper'][0]
    else:
        aniso_fac = minbar.ANISOTROPY['non-dipper'][1] / minbar.ANISOTROPY['non-dipper'][0]

    s.clear()

    # ----------------------------------------------------------------------------
    # Generate the list of bursts.

    # print("{}: generating the list of bursts matching criteria {}".format(fn, criteria + filter))
    if verbose:
        print("{}: generating the list of bursts".format(fn))

    # Build up the burst list

    m = minbar.Bursts()

    m.name_like(src)
    m.unique()  # Eliminate the multiple bursts
    lburst = m['entry']

    # Exclude short recurrence-time bursts

    if exclude_short:
        _tdel = m['tdel'].data
        good = np.where((_tdel == 0.0) | (_tdel > tdel_thresh))[0]
        if len(lburst) > len(good):
            if verbose:
                print("{}: omitting {} of {} short del-t bursts".format(fn, len(lburst) - len(good), len(lburst)))
            if len(good) > 0:
                lburst = lburst[good]
            else:
                minbar.logger.error('calc_alpha: exclude_short left no bursts')
                lburst = []
                return None

    n = len(lburst)

    # Exclude bursts here where there are multiple active sources in the (RXTE)
    # field

    #  _o=obs_id(lburst,/nodup)
    #  dbopen,!minbar_root+'/minbar-obs'
    #  dbext,_o,'flag',_flag
    #  conf=where((_flag and 2) gt 0,nconf,compl=_ok)
    #  if nconf gt 0 then begin
    #    print,nconf,n, $
    #      format='("** WARNING ** omitting ",i4," of ",i4," bursts in confused fields")'
    #    lburst=lburst[_ok]
    #  endif

    n = len(lburst)

    # Extract the values on which to bin. Normally this is gamma, but if the
    # S_Z flag is set, it will be s_z. (The variable is always called gamma, though)

    # dbext,lburst,'name,gamma,fluen',name,gamma,fluen
    name = m[lburst]['name'].data
    gamma = m[lburst]['gamma'].data
    fluen = m[lburst]['bfluen'].data
    fluen_err = m[lburst]['e_bfluen'].data

    if s_z:
        gamma = m[lburst]['s_z'].data
    #  name=strtrim(name,2)

    # Now loop over the sources and check that there are bursts from this source
    # This section is for where we are calculating the rate averaged over more
    # than one source

    #  for i=0,nsrc-1 do begin
    # ;    _l=dbfind(criteria+'cat_num>0,type=1,flag<15,name='+src[i],/silent)
    # ;    if _l[0] gt 0 then begin
    # ;      if lburst[0] eq -1 then lburst=_l else lburst=[lburst,_l]
    # ;    endif else begin
    #    tmp=where(name eq src[i])
    #    if tmp[0] eq -1 then $
    #
    # ; If we find no bursts for this source, we should remove it from the list!
    # ; Or should we? I'm not sure. I think No.
    #
    #      print,fn,src[i], $
    #        format='(a,": ** WARNING ** no bursts found for ",a," with these criteria")'
    # ;    endelse
    #  endfor

    # Check the gamma values for the bursts

    bad, nbad, good = idlwhere(gamma <= 0.)  # ,nbad,compl=nz)

    if verify:
        if (nbad > 0) and verbose:
            print("{}: ** WARNING ** {}/{} burst gamma-values not set".format(fn, len(bad), n))

        # Recalculate the gamma-values from minbar-obs, to check how up-to-date the
        # values from MINBAR are

        oid = m['entry_obs'].data

        missing = np.where(oid <= 0)[0]  # ,nmissing)
        if len(missing) > 0:
            if verbose:
                print("{}: ** WARNING ** {} bursts without matching observation in minbar-obs".format(fn, len(missing)))
                print(lburst[missing])
            if (max(gamma[missing]) > 0.) and verbose:
                print(fn + ': ** ERROR ** burst with missing obs has non-zero gamma value (how?)')

    if len(good) == 0:
        # There are no good gamma values, so we have to fall back to the flux
        if verbose:
            minbar.logger.warning('no good gamma values, so falling back on the flux')
        gflag = False
        gamma = m[lburst]['pflux'].data * _bc * aniso_fac
        bad, nbad, good = idlwhere(gamma <= 0.)  # ,nbad,compl=nz)

    # Extract all the burst parameters. The important parameters here are
    #   lburst -> burst IDs, including those with no gamma values
    #
    # next three arrays are for the subset with non-zero gamma
    #   name -> source name
    #   gamma -> becomes binning parameter nflux_bursts for bursts (flux/gamma/s_z)
    #   fluen -> integrated energy in burst (normalized if gflag=1)
    #   n = final number of "good" bursts (i.e. with gamma[good]=nflux_bursts>0)

    min_gamma = min(gamma[good])  # To calculate lower limit for obs

    #  min_gamma=0.01	; special here (for 1636-536?)

    # Here we limit to the number of bursts for which we have gamma values

    nflux_bursts = gamma[good]  # Normalized flux for bursts, where applicable
    fluen = fluen[good]
    name = name[good]
    n = len(good)

    # Also clean up the various pointers

    badbursts = lburst[bad]
    lburst = lburst[good]

    # ----------------------------------------------------------------------------
    # Get the list of observations, as well as an F_Edd value for each obs
    # appropriate for the source. This is so we can calculate gamma

    obs_limit = [min_gamma * 0.9, 10.]

    criteria = 'flag<1,'  # This eliminates observations with multiple active
    #   sources in the FOV, observations with no good
    #   times, and observations with no Standard 2 data
    # see http://burst.sci.monash.edu/wiki/index.php?n=MINBAR.Minbar-obsColumn#flag
    criteria = criteria + 'sig>3,fluxe>0,'

    if verbose:
        print("{}: generating the list of observations matching criteria {}".format(fn, criteria))

    o = minbar.Observations()

    o.name_like(src)
    o.good()
    lobs = o['entry']

    # Now extract the relevant columns

    #  dbext,lobs,'name,gamma,flux,s_z,exp',_name,gamma,flux,_s_z,dur
    _name = o[lobs]['name'].data
    gamma = o[lobs]['gamma'].data
    flux = o[lobs]['flux'].data
    flux_err = o[lobs]['e_flux'].data
    _s_z = o[lobs]['s_z'].data
    dur = o[lobs]['exp'].data

    # Extract all the observation parameters. The important parameters here are
    #   flux -> binning param nflux for observations (flux/ofedd,s_z or just flux)
    #        -> also (appropriately rescaled) flux value, nflux3; e.g. for s_z
    #             may be different from nflux!
    #   dur -> duration of each observation
    #   oid -> corresponding observation id (may have two columns!)
    #   badbursts -> list of IDs for those bursts without a gamma value
    # Here we extract our binning parameter, as well as the (appropriately
    # normalized) source flux

    #  if gflag then nflux3=flux/ofedd else nflux3=flux
    if s_z:
        nflux = _s_z
    else:
        if gflag:
            nflux = gamma
        else:
            nflux = flux

    bad, nbad, nz = idlwhere(nflux <= 0.)
    if verify and (not s_z):
        if (nbad > 0) and verbose:
            print("{}: ** WARNING ** {}/{} obs flux/gamma-values not set".format(fn, nbad, len(lobs)))

    # Now that we have a maximal set of fluxes, select the "good" observations here

    gfl, nobs, oexcl = idlwhere((nflux > obs_limit[0])
                                & (nflux < obs_limit[1]))  # Good fluxes, was tmp

    if (len(oexcl) > 0) and verbose:
        print("{}: ** WARNING ** excluded {} observations based on gamma".format(fn, len(oexcl)))
    # if debug and (limit_low) and (not s_z):
    #     print('** WARNING ** omitting observations with gamma/s_z < {}'.format(min_gamma*0.9))
    if nobs == 0:
        minbar.logger.error('{}: ** WARNING ** can''t calculate rates, no valid fluxes'.format(fn))
        return None

    # This section may be partially redundant

    badobs = []
    if len(nflux) > nobs:
        #    bad=where(nflux le 0.0)
        if verbose:
            if s_z:
                print("{}: Omitting {} of {} observations for binning (on S_Z)".format(fn, len(oexcl), len(nflux)))
            else:
                print("{}: Omitting {} of {} observations for binning (on flux/gamma)".format(fn, len(oexcl), len(nflux)))
        badobs = lobs[oexcl]

    # Now reduce the arrays to the good values

    nflux = nflux[gfl]
    dur = dur[gfl]
    _name = _name[gfl]
    lobs = lobs[gfl]

    # ----------------------------------------------------------------------------
    # Now we actually do the binning.

    bins = generate_bins(nflux_bursts, nflux, nperbin=nperbin, x_lim=limit, x_add=custom)
    nbins = len(bins)
    nb = np.zeros(nbins)    # Number of bursts in each bin
    nb_lo = np.zeros(nbins) # 1-sigma lower and upper limits
    nb_hi = np.zeros(nbins)
    bin_id = {}             # Keep track of which bursts and observations are in each bin
    bin_oid = {}
    yall = np.zeros(nbins)  #
    rateall = np.zeros(nbins)       # Average burst rate in each bin
    rateall_lo = np.zeros(nbins)    # 1-sigma lower and upper limits on rate
    rateall_hi = np.zeros(nbins)
    alpha = np.zeros(nbins)         # Average alpha for each bin
    alpha_err = np.zeros(nbins)
    alpha_err_lo = np.zeros(nbins)  # 1-sigma lower and upper limits on alpha
    alpha_err_hi = np.zeros(nbins)

    # Assign the bursts, and observations, to each bin

    ind_b = np.searchsorted(bins, nflux_bursts, side='right')

    ind_o = np.searchsorted(bins, nflux, side='right')

    # And loop over the bins to calculate the rate etc.

    for i in range(nbins):

        if debug:
            print(i, bins[i + 1])

        # Find all the bursts that fall into this bin. If we are using more than
        # one source, split up the results per source

        gb, ngb, dummy = idlwhere(ind_b == i + 1)
        if ngb > 0:
            nb[i] = ngb

            # Accumulate the burst IDs in this bin

            bin_id[i] = lburst[gb]

            # Calculate the errors; see e.g. https://gist.github.com/JohannesBuchner/ea0269f567e2d1f062b0304b575884fc

            nb_lo[i] = scipy.special.gammaincinv(nb[i] + 1, q)
            nb_hi[i] = scipy.special.gammaincinv(nb[i] + 1, 1. - q)

            # What if we have a few zero fluences? This should work fine for either
            # gflag=1 or 0, since in the former case we're using U_b

            bad, nbad, nz = idlwhere(fluen[gb] <= 0.0)
            if nbad > 0:
                mfluen = np.mean(fluen[gb[nz]])
                fluen[gb[bad]] = np.zeros(len(bad)) + mfluen
                if verbose:
                    print("{}: ** WARNING ** one or more zero fluences for bin {}...".format(
                        fn, i))
                    print("{}: replaced {} of {} fluence values with bin mean of {}e-9 ergs/cm^2".format(
                        fn, len(bad), len(gb), mfluen))
                # print(nb[i], fluen[gb])

        else:

            # Always have an upper limit, even when nb[i] is zero

            nb_hi[i] = scipy.special.gammaincinv(nb[i] + 1, conf_lim)

        # Find all the observations that fall into this bin

        g, ng, dummy = idlwhere(ind_o == i + 1)
        if ng > 0:
            bin_oid[i] = lobs[g]

            # Calculate the rate and error, and alpha and error

            yall[i] = sum(dur[g] / 1000.)  # Duration in ks
            rateall[i] = (nb[i]) / (yall[i] / 3.6)
            # Gaussian version
            #      ratealle[i]=max([1.,sqrt(float(nb[i]))])/(yall[i]/3.6)
            # Define instead asymmetric errors
            # print (k,scipy.special.gammaincinv(k+1,q),scipy.special.gammaincinv(k+1,1.-q))

            rateall_lo[i] = rateall[i] - nb_lo[i] / (yall[i] / 3.6)
            rateall_hi[i] = nb_hi[i] / (yall[i] / 3.6) - rateall[i]

            # print (nb[i], ng, rateall[i], rateall_lo[i], rateall_hi[i])

            # Calculate the alpha value and error

            if nb[i] > 0:
                flux_int = sum(dur[g] * flux[g]) * _bc
                fluen_sum = sum(fluen[gb]) * 1e3
                alpha[i] = aniso_fac * flux_int / fluen_sum
                # alpha_err[i] = alpha[i]/sqrt(float(nb[i]))
                alpha_err[i] = aniso_fac * np.sqrt(sum((dur[g] * _bc / fluen_sum) ** 2 * flux_err[g] ** 2)
                                       + sum((flux_int / fluen_sum ** 2) ** 2 * (fluen_err[gb[nz]]*1e3) ** 2))

                # Factor in the bolometric correction uncertainty if present

                if _bce > 0.:
                    alpha_err[i] = np.sqrt( alpha_err[i]**2 + alpha[i]**2*(_bce/_bc)**2 )

                # Factor in the Poisson uncertainties here
                alpha_err_lo[i] = np.sqrt( alpha_err[i]**2 + alpha[i]**2*((nb[i]-nb_lo[i])/nb[i])**2 )
                alpha_err_hi[i] = np.sqrt( alpha_err[i]**2 + alpha[i]**2*((nb_hi[i]-nb[i])/nb[i])**2 )

    # return everything

    result =  {'src': src, 'minbar_version': minbar.VERSION, 'date': minbar.DATE,
               'conf': conf, 'conf_lim': conf_lim, 'bc': (_bc, _bce), 'aniso': aniso_fac, 'nperbin': nperbin,
               'bins': bins, 'n': nb, 'rate': rateall, 'rate_err': (rateall_lo, rateall_hi),
               'alpha': alpha, 'alpha_err': (alpha_err_lo, alpha_err_hi),
               'id': lburst, 'id_obs': lobs, 'badbursts': badbursts, 'badobs': badobs,
               'bin_id': bin_id, 'bin_oid': bin_oid }

    # Optionally save the results to a pickle file. You can read this in as follows:
    # import pickle
    # result = pickle.load( open( "4U 1608-522_rate.p", "rb" ) )
    if write:
        file = src+'_rate.p'
        if os.path.isfile(file) and not clobber:
            minbar.logger.warning('savefile exists, and clobber=False; skipping save')
        else:
            if verbose:
                print('{}: writing rate, alpha values to file {}'.format(fn, file))
            pickle.dump(result, open(file, "wb"))

    return result

def binplot(bins, rate, rate_err, panel=None, pbins=None,
            log=True, xrange=[0.009, 0.99], yrange=[0.05, 2],
            color='green', noxlab=False, bc=1.0, ylabel='Burst rate (hr$^{-1}$)'):
    '''This is the replacement for the IDL version of binplot, part of
    analysis/burst_rate.pro
    To get the right behaviour, we need to use matplotlib's step routine,
    with option "post", and duplicate the last (non-zero) value of the
    rate array
    There are also a limited number of options for controlling the plot
    appearance, including the option to omit the x-labels for stacked plots
    (noxlab = True)'''

    if panel is None:
        panel = plt.plot()

    # Bolometric correction is handled in the gamma-calculation stage, now

    #    assert (r['bc'] == 1)

    #    y=r['rate'][0] # burst rate array
    y = rate
    y = y[np.nonzero(y)] # non-zero rates for plotting
    y = np.append(y, y[-1])  # augment by one for step plot
    nbins = len(y)

    # need to modify the length of the rate_err tuples here to omit the last
    # bin

    asym_errors = False
    if isinstance(rate_err, tuple):
        asym_errors = True
        rate_lo, rate_hi = rate_err
        # Define a combined rate "array-like" of shape(2,N): Separate - and + values for each bar.
        # First row contains the lower errors, the second row contains the upper errors.
        _rate_err = (rate_lo[0:nbins - 1], rate_hi[0:nbins - 1])
        # print (np.shape(_rate_err))
        # We take one element off the arrays below to match pbins
        lim = np.logical_and(rate <= 0.0, rate_hi > 0.0)[:-1]
        rate_hi = rate_hi[:-1]
        # print(lim)
    else:
        _rate_err = rate_err[0:nbins - 1]
        lim = []

    #    print (y)
    nbins = len(y)

    # Generate the bin midpoints for plotting

    if pbins is None:
        if log:
            pbins = np.sqrt(bins[1:] * bins[:-1])
        else:
            pbins = 0.5 * (bins[1:] + bins[:-1])

    #    bins=r['bins'][0][0:nbins]*bc # bins for step plot

    #    pbins=r['pbins'][0][0:nbins-1]*bc # where to plot symbols (with errors)
    #    rate_err=[r['rate_err_lo'][0][0:nbins-1]*bc,r['rate_err_hi'][0][0:nbins-1]*bc]

    #    print (nbins)
    #    y=np.append(r['rate'][0],r['rate'][0][nbins-1])
    #    print (r['rate'][0],y)

    #    plt.step(r['pbins'][0],r['rate'][0][0:nbins])
    #    print (len(bins),len(y))
    plt.step(bins[0:nbins], y, where='post', color=color)

    # Plot the points with errorbars, if pbins is defined

    try:
        plt.errorbar(pbins[0:nbins - 1], rate[0:nbins - 1], yerr=_rate_err,
                     fmt='o', color=color)
        # Plot the limits

        if np.any(lim):
            plt.errorbar(pbins[lim], rate_hi[lim],
                         yerr=0.3 * rate_hi[lim], uplims=True, fmt='o', color=color)
            tmp = np.where(lim)[0][0]
            #        print (pbins[lim],rate_hi[lim],tmp,type(tmp))

            plt.plot(bins[tmp:tmp + 2], np.array([rate_hi[tmp], rate_hi[tmp]]),
                     '-', color=color)

            plt.plot(pbins[lim], rate_hi[lim], '.')

    except:
        pass

    plt.xlim(xrange)
    plt.ylim(yrange)

    # plt.xlabel('Accretion rate $\dot{m}/\dot{m}_{\rm Edd}$')
    if noxlab == False:
        plt.xlabel('Accretion rate $\dot{m}/\dot{m}_\mathrm{Edd}$')

    if log:
        plt.xscale('log')
        plt.yscale('log')
    plt.ylabel(ylabel)

def reduce(obs, to_level=None, clean=False, clobber=False, dryrun=False):
    '''
    Routine to reduce & analyse a given observation. Will check for different components
    (lightcurve, spectrum, analysis results) and if they are not there, regenerate them.
    This behaviour is moderated by the to_level option, with the descriptions below (see
    also the Observation class):

    | level = 0 # no analysis products (known to be) present
    | level = 1 # raw analysis/reduction steps have been completed but no lightcurves or spectra yet created
    | level = 2 # lightcurve(s) are available
    | level = 3 # burst search has been completed and any burst data has been generated
    | level = 4 # time-resolved spectroscopy is available
    | level = 5 # persistent spectra are available
    | level = 9 # defined from a MINBAR observation entry, may not be any raw data available

    Calls analyse_persistent, and the burst search; if any bursts are found, will also
    call analyse_burst

    :param obs: target observation to analyse (Observation object)
    :param to_level: target observation level to reach
    :param clean: boolean, set to True to redo all steps, or if False, will only do incremental analysis
    :param clobber: boolean, set to True to overwrite previous analyses files
    :param dryrun: if True, won't create any files, will just write what's going to happen

    :return:
    '''

    from datetime import date

    header_default = '''MINBAR DR1+ data file
Created {} with the pyminbar package, see https://github.com/outs1der/minbar
'''.format(date.today())
    burst_file = 'bursts.lis'
    burst_dir = 'burst{}'
    burst_win = (-50, 206) # [s] time window to extract burst lightcurve over
    burst_lc_file = 'burst.lc'

    # dict listing all the files created, so we can clean them later

    files_created = {1: (), 2: (),
        3: (burst_file, burst_dir.format('*'), '/'.join((burst_dir.format('*'),burst_lc_file))),
        4: (),
        5: () }

    assert type(obs) == minbar.Observation
    levels = (0, 1, 2, 3, 4, 5, 9)

    if to_level is not None:
        if obs.level >= to_level:
            minbar.logger.warning('analysis already complete to level {}, nothing to do'.format(obs.level))
            return True
        if to_level not in levels:
            minbar.logger.error('no such processing level {}'.format(to_level))
            return False

    path = obs.get_path()
    minbar.logger.info('Reducing {} observation #{} in directory {}...'.format(obs.instr.name, obs.obsid, path))
    for level in levels[(1-int(clean))*(obs.level+1):levels.index(to_level)+1]:
        # loop over the levels here
        minbar.logger.info('Level {} processing...'.format(level))

        # For Python 3.10 you could use the match-case statement; see e.g.
        # https://stackoverflow.com/questions/11479816/what-is-the-python-equivalent-for-a-case-switch-statement
        if level == 1:
            # perform raw analysis/reduction steps have been completed but no lightcurves or spectra yet created
            pass
        elif level == 2:
            # extract lightcurve(s)
            pass
        elif level == 3:
            #################################################################
            # run burst search and generate subdirectories etc.
            #
            # runs the burst file and creates burst.lis in the observation directory,
            # and then for each detected burst a subdirectory burst[n] with a lightcurve file
            # around the burst time

            if obs.time is None:
                obs.get_lc()

            bursts, param, sigma = findburst(obs.mjd.value, obs.rate.value, obs.error.value, tunit='d', fit=True)

            time_unit = '{} {}'.format(obs.mjd.format.upper(), obs.mjd.scale.upper())
            # verify step here
            if len(bursts) == 0:
                logger.info('no bursts detected, nothing more to do')
            else:
                print ('Got burst(s) at {}, triggering data extraction'.format(bursts))
                # useful diagnostic plot to show burst candidates, could have a loop here that would
                # allow the user to verify the events
                fig = obs.plot(show=False)
                plt.plot(bursts, np.full(len(bursts), max(obs.rate)*1.05),'|r')
                plt.show()

                # write the burst file
                obs.burst_file = burst_file
                _path_and_file = '/'.join((path,obs.burst_file))
                obs.burst_cand = bursts # can't use bursts attribute, as that's just for MINBAR bursts
                if os.path.isfile(_path_and_file) & (not clobber):
                    minbar.logger.warning('{} exists and clobber is set to False, skipping burst file write'.format(obs.burst_file))
                    minbar.logger.info('Level {} processing interrupted, need to reconcile burst search results and existing file'.format(level))
                    break
                else:
                    f = open(_path_and_file, 'w')
                    np.savetxt(f, bursts, fmt='%11.5f',
                               header=header_default+'''File {}

Output of burst search on {} observation {}

Created with minbar.reduce v{}

Burst times are in {}'''.format(obs.burst_file, obs.instr.name, obs.obsid, minbar.__version__, time_unit))
                    f.close()

                for i, burst in enumerate(bursts):
                    minbar.logger.info('processing burst #{}, {} {}'.format(i+1, obs.mjd.format.upper(), burst))
                    _burst = burst_dir.format(i+1)
                    _burst_path = '/'.join((path,_burst))
                    # os.mkdir(_burst_path)
                    sel = (obs.mjd.value > burst+burst_win[0]/86400) & (obs.mjd.value < burst+burst_win[1]/86400)
                    _path_and_file = '/'.join((_burst_path,burst_lc_file))

                    # assemble the extraction command
                    # this is supposed to be generic so that it would work with any instrument
                    _extract = ' '.join((obs.filter_events['invoke'], obs.filter_events['outfile'])).format(
                        obs.obsid,
                        obs.filter_events['timesel'].format(min(obs.time[sel].value), max(obs.time[sel].value)),
                        '/'.join((_burst_path,'burst.evt')))
                    breakpoint()
                    # this will run once you import subprocess, but throws an error:
                    # ERROR: Device not configured
                    # Task niextract - events 0.0 terminating with status 6
                    _result = subprocess.run(';'.join((obs.reduce_init,'cd {}'.format(obs.get_path()),_extract)), shell=True, executable='/bin/csh')
                    breakpoint()
                    if os.path.exists(_path_and_file) & (not clobber):
                        minbar.logger.warning('{} exists and clobber is set to False, skipping burst lightcurve file write'.format(burst_lc_file))
                        minbar.logger.info('Level {} processing interrupted, need to reconcile burst search results and existing file'.format(level))
                        break
                    else:
                        # original way of generating the burst lightcurve, extracted from the main lightcurve
                        f = open('/'.join((_burst_path,burst_lc_file)), 'w')
                        np.savetxt(f, np.column_stack([(obs.mjd.value[sel]-burst)*86400, obs.rate.value[sel], obs.error.value[sel]]),
                                   fmt='%9.5f', header=header_default+'''File {}

Extracted from {} observation {}, all-observation lightcurve
  {}
  
Rate is NOT background subtracted, assumed effective area is {}

Columns:
  Time (s relative to {} {:.5f}), rate & error (per {})'''.format(
                            burst_lc_file, obs.instr.name, obs.obsid, obs.file, obs.instr.effarea, time_unit, burst, (1./obs.rate[0]).unit))
                        f.close()
        elif level == 4:
            # extract time-resolved spectroscopy
            pass
        elif level == 5:
            # extract persistent spectra
            pass
        elif level == 9:
            # write to the database?
            pass

        minbar.logger.info('... level {} processing complete.'.format(level))

    minbar.logger.info('...reduce step complete')
    return True # assuming no errors


def analyse_burst(self, obs, burst):
    """
    Function to analyse each burst in a single
    observation (not yet implemented)

    :param obs: Observation in which the burst occurred
    :param burst: [str] directory with burst data

    :return:
    """

    pass


def analyse_persistent(self, obs):
    """
    Function to analyse the persistent emission (lightcurve and spectrum) for a single
    observation (not yet implemented)

    :param obs: Observation in which the burst occurred

    :return:
    """

    pass

