# Classes and routines to analyse X-ray observational data

import numpy as np
from scipy.signal import correlate

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


def burst_template(x, a):
    """
    Generates a template burst within a time series given by the first
    variable passed, with the specified profile. Includes calculation of
    partial derivatives for use with curvefit. Parameters are:
    a[0] = baseline count rate
    a[1] = start time of burst
    a[2] = rise time of burst
    a[3] = peak intensity of burst
    a[4] = exponential decay time

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


def findburst(_time, _rate, _error=None, tunit=None,
    rise=None, edt=20., refine=False, fit=False,
    # sig=sig, param=param, error=sigma, resid=resid,
    plot=False, # xrange=xrange,
    min_burst_sep=3., sigthresh=5., mjdref=0.0, dbg=False):
    """
    $Id: findburst.pro,v 1.1 2011/06/29 12:54:00 duncan Exp duncan $

    This routine is designed to find one or more bursts in noisy data.
    Extracted from search.pro v1.3, 29/6/2011

    This version moved from minbar/jemx, and simplified 23.10.2014

    Input parameters:
      time [days]
      rate & error [counts/sec]
      rise [sec] (optional) burst rise time
      edt [sec] (optional) burst timescale
    Outputs:
      <t_cand> list of candidate burst times [MJD]
      sig estimated significance of burst detection
      param list of burst parameters: peak rate [count/s], base rate [count/s],
              fitted start time [MJD], rise [s], fitted peak rate [count/s],
              exponential decay time [s]
      sigma error on fit parameters (from curvefit)
      resid residual lightcurve, after all detected bursts have been subtracted
              out

    Example usage:

    import minbar
    from astropy.io import fits

    file = '../xmm/data/1RXSJ180408.9-342058/0741620101/2to7good.fits'
    hdul = fits.open(file)
    data = hdul[1].data
    test = minbar.findburst(data['TIME'], data['RATE'], data['ERROR'])
    print(test)
    [5.42058957e+08 5.42067368e+08 5.42075296e+08 5.42083081e+08 5.42090903e+08]

    Routines called:
      burst_template
      gtiseg
      curvefit
      dt_nozero
    """

    # Some parameters passed to the IDL routine via common statements
#  common findburst,sigthresh,min_burst_sep
#  common jemx,mjdref    ; MJDREF not currently used

    sn='findburst'

    if rise == None:
        rise=edt/4.     # Default burst rise [sec]

    if refine:
        print("{}: ** ERROR ** refining of template not yet implemented".format(sn))
        return None

# Make copies of the input arrays to avoid unintentionally modifying them

    time = _time
    rate = _rate
    n = len(time)
    dt = dt_nonzero(time)         # Nominal dt

# The error is required for the fitting; we should have an option where it
# can be estimated from the lightcurve, to make the code easier to use

    if _error is None:
        print("{}: ** WARNING ** errors are required".format(sn))
        return None
    else:
        error=_error

# Also search for non-finite burst rate measurements

    fin = np.where(np.isfinite(rate) & np.isfinite(error))[0]
    if len(fin) < n:
        print("{}: ** WARNING ** non-finite rate measurements, excising from analysis".format(sn))
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
        print("{}: ** WARNING ** setting tunit='{}', set explicitly if incorrect".format(sn,tunit))
    if tunit == 's':
        sd=1.0
    elif tunit == 'd':
        sd=86400.0
    else:
        stop,'** ERROR ** unrecognised time unit'
    md=sd/60.

# Define the fitting window here

    pb=32.0/sd            # interval prior to the burst to include in the
                          #   fitting window
    fwin = np.array([-pb,max([min_burst_sep/md,edt/sd])-pb])
    fwin_ind = np.where((time-time[0] > fwin[0]) & (time-time[0] < fwin[1]))[0]
    if dbg:
        print("{}: defined fitting window as {} relative to the burst time".format(
            sn, fwin*sd))

# Some screening here. This gap selection used for the April and June 2011
# searches, but subsequently Jerome suggested it was too restrictive...
# although in actual fact it only results in one exception

    g = np.where((dt*sd >= 16.) & (rate > 4000.))[0]
    if len(g) > 0:
        print ("{}: ** WARNING ** exception in lightcurve".format(sn))
        return None

# First define the burst template

    m=np.median(rate)
    p=max(rate)-m
#  t0=0.5*(max(time)+min(time)) ; Trial burst is in the middle of the segment
    # t0=time[int(n/2)]
    # trial, pder = burst_template(time,[m,t0,rise/sd,p,edt/sd])
    t0 = time[0]+pb
    trial, pder = burst_template(time[fwin_ind],[m,t0,rise/sd,p,edt/sd])
    if dbg:
        plt.plot(time[fwin_ind],trial)

# Now cross-correlate the template, sliding the window over the entire lightcurve

    c=correlate(rate, trial, mode='same')
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
        return None

# Define an array giving the start and end indices for each region
# The 5% padding below avoids issues with rounding errors

    # print (1.05*min_burst_sep/md,min_burst_sep)
    t_cand_ind=gtiseg(time[cexc], maxgap=1.05*min_burst_sep/md)#,indices=t_cand_ind)
    t_cand_ind=cexc[t_cand_ind]
    # print (t_cand_ind,np.shape(t_cand_ind))

    if not fit:
        return time[t_cand_ind[:,0]]
    else:

        # Now set up the result arrays, based on the number of contiguous regions
        # exceeding the search threshold, identified in the previous step

        ncand=len(t_cand_ind)
        t_cand = np.zeros(ncand)
        param = np.zeros((ncand,6))
        sig = np.zeros(ncand)
        result = -1     # originally the fits for each burst; not currently used
        sigma = param.copy()
        weights = 1./error^2
        resid = rate.copy()
        for i in range(ncand):

# Calculate peak intensity in the window over which the significance
# exceeds the threshold

            p=max(rate[t_cand_ind[0,i]:t_cand_ind[1,i]],ipmax)
            ipmax += t_cand_ind[0,i]      # note ipmax is a pointer into time, not win
            pe=error[ipmax]
            t_cand[i]=time[ipmax]
            sig[i]=max(csig[t_cand_ind[0,i]:t_cand_ind[1,i]])

# Now determine the window over which to do the fitting. This has a
# minimum size of fwin

            win = np.where(time-t_cand[i] >= min([fwin[0],time[t_cand_ind[0,i]]-t_cand[i]])
              & time-t_cand[i] < max([fwin[1],time[t_cand_ind[1,i]]-t_cand[i]]))[0]

# Start the fitting loop

            _param=[m,t_cand[i],rise/sd,p,edt/sd]

            if len(win) > len(_param):

# try to do a fit, now just within a window surrounding the burst

                _param=[_param[0],t_cand[i],_param[2:4]]
                _result=curvefit(time[win],rate[win],weights[win],_param,_sigma,
                    function_name='burst_template',iter=1000,status=_s)

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

            if refine & (ncand == 1) & (_s == 0) & (_param[4] < 1000./sd):
                trial2 = burst_template(time,[_param[0],t0,_param[2:5]])


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
                                         function_name='burst_template',iter=1000,status=_s)
# if dbg then begin
#   plot,time[win],rate[win]
#   oplot,time[win],_result,color='00ff00'x
#   a=get_kbrd(1)
# endif
                else:
                    _s=-1      # don't have enough points to fit

# If the fit has failed, or some of the uncertainties are non-finite, then
# the results can't be trusted; so blank them out

        if (_s > 0) or (min(finite(_sigma)) == 0) or (_param[4] < 0.0):
            _param=[0.0,t_cand[ncand-1],fltarr(3)]
            _sigma=fltarr(5)

# Assemble parameter, result list

        sigma[:,i]=[pe,_sigma*[1.,1.,sd,1.,sd]]
# if dbg then begin
#   plot,time,rate
#   if ncexc gt 0 then oplot,[time[cexc]],[rate[cexc]],color='00ff00'x
#   a=get_kbrd(1)
# endif

# Subtract off the model from the residual

        burst_template,time[win],[0.0,_param[1:4]],model
        resid[win]=resid[win]-model


# if dbg then begin
#   plot,time,rate
#   if ncexc gt 0 gt -1 then $
#     oplot,[time[cexc]],[rate[cexc]],color='00ff00'x
#   a=get_kbrd(1)
# endif

# if dbg then stop

        if param[0] == -1:
            param=p
                        # If we don't fit, we only return the burst peak

    return t_cand
