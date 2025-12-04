# folded out of the Bursts class so that it's available as a stand-alone
# tool

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib import gridspec
import matplotlib.ticker as mticker

import minbar

def burstplot(_bdata, param='flux', show=True, **kwargs):
    """
    General-purpose routine to plot burst data. Would like to be able
    to call this in a number of ways, both with a burst ID from
    MINBAR, but also with a pandas table (as read in with
    get_burst_data, for example). And do a bunch of different plots,
    including the three-panel "in 't Zand" plot (see
	    e.g. Fig 3, `in 't Zand et al. 2012, A&A 547, A47 <https://ui.adsabs.harvard.edu/abs/2012A%26A...547A..47I>`_), as well as the
    three-panel "HR-diagram" version

    TODO add some annotation identifying the burst, somewhere...

    Usage:

    | burstplot(burst,param='rad',xlim=[-5,25])
    | burstplot(burst,param=['flux','rad','kT','chisq'])

    :param bdata: pandas table of input data, including for bursts not in MINBAR
    :param param: parameter or list of parameters to plot; special 'hr' to plot the three-panel HR-diagram version
    :param show: show the figure immediately if ``True``, otherwise withhold it (e.g. to add further annotation etc.)
    """

    def plot_param(bdata, ax, param='flux', ylabel=None, color=None):
        """
        This routine is called by burstplot to display each panel of data
        It plots binned data as steps, with errors

        :param bdata: burst data to plot
        :param ax: axis to plot on
        :param param: parameter to plot; has to be present in bdata
        :param ylabel: dict of y-labels, param names as key
        :param color: dict of colors, param names as key
        """

        # boolean to determine whether the parameter has an error, but doesnt' tell you 
        # what the _name_ of that error parameter is

        has_error = np.all([x in bdata for x in [param+'_min',param+'_max']])\
            | (param+'err' in bdata) | (param+'_err' in bdata) \
            | ('e_'+param in bdata) 
            # | (param == 'r') redundant, for MINBAR data

        # filter on good data
        # _gd = bdata.flux*bdata.fluxerr > 0
        _gd = bdata.flux > 0
        if 'fluxerr' in bdata:
            _gd = bdata.flux*bdata.fluxerr > 0

        # add an extra value copy here to plot that last step
        ax.step(np.append(bdata.time[_gd].values,
            bdata.time[_gd][-1:].values+bdata.dt[_gd][-1:].values),
            np.append(bdata[param][_gd].values,bdata[param][_gd][-1:].values),
            where='post',color=color[param])

        if has_error:
            # plot errors
            yerr = None
            if param == 'r':
                # special convention here for the radius
                yerr = bdata['re'][_gd]
            else:
                # want to accommodate more than one label on the errors
                try:
                    yerr = np.stack((bdata[param][_gd]-bdata[param+'_min'][_gd],
                        bdata[param+'_max'][_gd]-bdata[param][_gd]))
                    # yerr = np.stack((bdata[param][_gd]-bdata[param+'_min'][_gd],
                    #     bdata[param+'_max'][_gd]-bdata[param][_gd]))
                except:
                    # we don't have <param>_min, <param>_max so have to
                    # use whatever's present
                    if param+'err' in bdata:
                        yerr = bdata[param+'err'][_gd]
                    elif param+'_err' in bdata:
                        yerr = bdata[param+'_err'][_gd]
                    elif 'e_'+param in bdata:
                        yerr = bdata['e_'+param][_gd]

            assert yerr is not None
            ax.errorbar(bdata.time[_gd]+bdata.dt[_gd]/2., bdata[param][_gd], yerr,
                         fmt='none',ecolor=color[param])
        ax.set_ylabel(ylabel[param])

    # Set the label names and colours here. To be passed also to
    # plot_param
    # Might need to set up some custom labels for the different
    # conventions of the SAX and RXTE data
    ylabel = {'r': 'Count rate [s$^{-1}$]',
              'flux': 'Flux [$10^{-9} \mathrm{erg\,cm^{-2}\,s^{-1}}$]',
              'kT': 'kT [keV]',
              'rad': 'Blackbody normalisation\n[$(R_{\mathrm{km}}/d_{10\ \mathrm{kpc}})^2$]',
              'chisq': 'Fit $\chi^2/n_{\mathrm{DOF}}$'}
    color = {'r': 'k', 'flux': 'k', 'kT': 'r', 'rad': 'b', 'chisq': 'g'}

    if type(_bdata) != pd.DataFrame:
        # we can also work on a concord ObservedBurst, but need to have
        # that package available; see https://github.com/outs1der/concord
        try:
            import concord as cd
        except:
            minbar.logger.error('required input not pandas DataFrame and concord not available')
            return None

        if type(_bdata) != cd.ObservedBurst:
            minbar.logger.error('required input not pandas DataFrame or concord ObservedBurst')
            return None

        # convert ObservedBurst to pandas table
        # new names for error component of flux, below - 2025 Dec
        bdata = pd.DataFrame({'time': _bdata.time, 'dt': _bdata.dt,
            'flux': _bdata.flux/minbar.FLUX_U, 'flux_err': _bdata.e_flux/minbar.FLUX_U})
        for _param in ylabel:
            for key in (_param, 'e_'+_param, 'E_'+_param,
                        # problem with _min, _max params coming from concord.ObservedBurst
                        # _param+'_min', _param+'_max',
                _param+'err', _param+'_err'):
                if hasattr(_bdata, key) & (key not in bdata.columns):
                    # print ('hasattr: {}'.format(key))
                    # if assigning new columns like this, don't want units
                    if hasattr(getattr(_bdata, key), 'unit'):
                        bdata[key] = getattr(_bdata, key).value
                    else:
                        bdata[key] = getattr(_bdata, key)
            
    else:
        bdata = _bdata.copy()

    xlabel='Time [s]'

    # check that the passed labels match one of the above

    test_param = param
    if type(param) == str:
        test_param = [param]
    for _param in test_param:
        if (_param not in ylabel.keys()) & (_param != 'hr'):
            print ('** ERROR ** plot param/type {} not recognized; allowed choices are:'.format(_param))
            for k in ylabel.keys():
                print ('  {}: {}'.format(k, ylabel[k]))
            print ('  hr: show 3-panel plot with HR-style diagram')
            return None

    fig = plt.figure()

    # Use GridSpec to constrain the layout, for maximum flexibility; see
    # https://matplotlib.org/stable/tutorials/intermediate/gridspec.html

    if param == 'hr':
        # this is the special three-panel plot with flux, blackbody
        # radius, and the H-R diagram on the right, with temperature
        # vs. flux

        gs = gridspec.GridSpec(2, 2)

        # kT - flux plot
        ax0 = fig.add_subplot(gs[:,1])
        kT_err = np.stack((bdata['kT']-bdata['kT_min'],bdata['kT_max']-bdata['kT']))
        flux_err = np.stack((bdata['flux']-bdata['flux_min'],bdata['flux_max']-bdata['flux']))

        ax0.errorbar(bdata.kT, bdata.flux, flux_err, kT_err)
        ax0.set_yscale('log')
        ax0.set_xscale('log')
        ax0.invert_xaxis()
        # ax0.ticklabel_format(useOffset=False, style='plain')
        ax0.xaxis.set_minor_formatter(mticker.ScalarFormatter())
        # ax0.ticklabel_format(style='plain', axis='x')
        ax0.set_xlabel(ylabel['kT'])
        ax0.yaxis.tick_right()
        ax0.yaxis.set_ticks_position('both')
        ax0.yaxis.set_label_position('right')
        ax0.set_ylabel(ylabel['flux'])

        ax1 = fig.add_subplot(gs[0,0])
        plot_param(bdata, ax1, 'flux', ylabel, color)

        ax2 = fig.add_subplot(gs[-1,0], sharex=ax1)
        plot_param(bdata, ax2, 'rad', ylabel, color)
        ax2.set_xlabel(xlabel)
    else:
        # generic plot
        if type(param) != list:
            param = [param]

        gs = gridspec.GridSpec(len(param), 1)

        for i, _param in enumerate(reversed(param)):

            if i == 0:
                ax0 = plt.subplot(gs[len(param)-i-1])
                plot_param(bdata, ax0, _param, ylabel, color)
                this_ax = ax0

            else:
                axi = plt.subplot(gs[len(param)-i-1], sharex = ax0)
                plot_param(bdata, axi, _param, ylabel, color)
                plt.setp(axi.get_xticklabels(), visible=False)
                this_ax = axi

            # This won't work if you have chisq as the first parameter
            if _param == 'chisq':
                this_ax.axhline(1, color='grey', linestyle='--')

        ax0.set_xlabel(xlabel)
        plt.subplots_adjust(hspace=.0)

        # interpret kwargs here
        if 'xlim' in kwargs:
            plt.xlim(kwargs['xlim'])
        # print (kwargs)

    plt.tight_layout()

    if show:
        plt.show()

    return fig



