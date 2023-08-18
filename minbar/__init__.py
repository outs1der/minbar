"""
Python class to easily extract information from MINBAR.
Provides the burst catalog (Bursts), observation catalog (Observations),
and source catalog (Sources). See the doc strings of those classes
for example usage.

The data files required for this package can be downloaded from the
bridges (figshare) repository at https://doi.org/10.26180/5e4a697d9b8b6
The files minbar.txt, minbar-obs.txt and minbar_sources.fits
should be placed in the 'data' subdirectory of this package.

The table and attribute descriptions, and the data analysis procedures,
are all described in the accompanying paper:
The Multi-INstrument Burst ARchive (MINBAR), by D.K. Galloway et al. (the
MINBAR Collaboration) 2020, ApJS 249, 32; available at
https://iopscience.iop.org/article/10.3847/1538-4365/ab9f2e

(c) 2020, Duncan Galloway duncan.galloway@monash.edu & Laurens Keek,
  laurens@xrb.space

Updated for MINBAR DR1, 2020, Duncan Galloway, duncan.galloway@monash.edu
Updated for MINBAR v0.9, 2017, Laurens Keek, laurens.keek@nasa.gov
"""

__author__ = """Laurens Keek and Duncan Galloway"""
__email__ = 'duncan.galloway@monash.edu'
__version__ = '1.15.0'

from .idldatabase import IDLDatabase
from .analyse import *
import numpy as np
import pandas as pd
import os, re
from astropy.io import fits, ascii
import astropy.units as u
from astropy.time import Time
from datetime import datetime
import logging
import sys

import matplotlib.pyplot as plt
from matplotlib import gridspec
import matplotlib.ticker as mticker

# kpc = 3.086e21 # cm
kpc = u.kpc.to('cm')*u.cm # cm

# Record the version and current date for analysis timestamps

VERSION = 1.0
DATE = datetime.now()

# Local paths for MINBAR data; you need to update this when installing
# on a new system
# This is NOT the directory where the table data files are found
# LOCAL_DATA flag is now determined dynamically as part of minbar.__init__

MINBAR_ROOT = '/Users/Shared/burst/minbar'
# LOCAL_DATA = True
MINBAR_URL = 'https://burst.sci.monash.edu/'

# Bytearr entries are for consistency between the IDL and ASCII
# versions of the data

MINBAR_INSTR_LABEL = {'PCA': 'XP', 'WFC': 'SW', 'JEM-X': 'IJ'}
MINBAR_INSTR_PATH = {'XP': 'xte', 'SW': 'wfc', 'IJ': 'jemx',
                     'PCA': 'xte', 'WFC': 'wfc', 'JEM-X': 'jemx',
                     b'XP': 'xte', b'SW': 'wfc', b'IJ': 'jemx'}

# derived effective areas used to renormalise PCA and JEM-X lightcurves;
# see the MINBAR paper (Galloway et al. 2020), Table 4

PCA_EFFAREA = 1400*u.cm**2
JEMX_EFFAREA = 64*u.cm**2
JEMX_EFFAREA_BURSTS = 100*u.cm**2

# Units for selected atributes in the MINBAR tables
# This is redundant for the MRT tables as those have their own units

FLUX_U = 1e-9*u.erg/u.cm**2/u.s
CFLUX_U = u.ct/u.cm**2/u.s
FLUEN_U = 1e-6*u.erg/u.cm**2
UNITS = { 'time': u.d, 'tstart': u.d, 'tstop': u.d, 'exp': u.s,
          'angle': u.arcmin, 'rise': u.s, 'tau': u.s, 'taue': u.s, 'e_tau': u.s,
          'dur': u.s, 'dure': u.s, 'e_dur': u.s, 'edt': u.s, 'edte': u.s, 'e_edt': u.s,
          'tdel': u.hr, 'trec': u.hr,
          'perflx': FLUX_U, 'perflxe': FLUX_U, 'e_perflx': FLUX_U,
          'flux': FLUX_U, 'fluxe': FLUX_U, 'e_flux': FLUX_U,
          'pflux': CFLUX_U, 'pfluxe': CFLUX_U, 'e_pflux': CFLUX_U,
          'count': CFLUX_U, 'counte': CFLUX_U, 'e_count': CFLUX_U,
          'fluen': u.ct/u.cm**2, 'fluene': u.ct/u.cm**2, 'e_fluen': u.ct/u.cm**2,
          'bpflux': FLUX_U, 'bpfluxe': FLUX_U, 'e_bpflux': FLUX_U,
          # Spectral model parameters
          'kT': u.keV, 'kTe': u.keV, 'e_kT': u.keV,
          'T_0': u.keV, 'T_0e': u.keV, 'e_T_0': u.keV,
          'kT_e': u.keV, 'kT_ee': u.keV, 'e_kT_e': u.keV,
          'line': u.keV, 'linee': u.keV, 'e_line': u.keV,
          'sigma': u.keV, 'sigmae': u.keV, 'e_sigma': u.keV,
          'gnorm': u.ct/u.cm**2/u.s, 'gnorme': u.ct/u.cm**2/u.s, 'e_gnorm': u.ct/u.cm**2/u.s,
          'rad': u.km/(10.*u.kpc), 'rade': u.km/(10.*u.kpc), 'e_rad': u.km/(10.*u.kpc),
          'bbnorm': (u.km/(10.*u.kpc))**2, 'bbnorme': (u.km/(10.*u.kpc))**2, 'e_bbnorm': (u.km/(10.*u.kpc))**2,
          'plnorm': u.ct/u.keV/u.cm**2/u.s, 'plnorme': u.ct/u.keV/u.cm**2/u.s, 'e_plnorm': u.ct/u.keV/u.cm**2/u.s,
          'bfluen': FLUEN_U, 'bfluene': FLUEN_U, 'e_bfluen': FLUEN_U }

# Bolometric corrections adopted for different sources, based on Table 9 from
# the paper (Galloway et al. 2020)
# Each entry gives the # of measurements, the mean, and the standard deviation
# Two special values are defined for "burster" and "pulsar" without their
# own measurements.

BOL_CORR = {'4U 0513-40':         (1, 1.47, 0.02),
            'EXO 0748-676':       (4, 1.6, 0.3),
            '4U 0836-429':        (1, 1.82, 0.02),
            '4U 1254-69':         (3, 1.30, 0.15),
            '4U 1323-62':         (3, 1.65, 0.05),
            'Cir X-1':            (4, 1.12, 0.06),
            '4U 1608-522' :       (4, 1.59, 0.13),
            '4U 1636-536' :      (44, 1.51, 0.12),
            'XTE J1701-462':      (2, 1.44, 0.07),
            'MXB 1658-298':      (17, 1.32, 0.05),
            '4U 1702-429':       (11, 1.4, 0.3),
            '4U 1705-44':         (7, 1.51, 0.15),
            'XTE J1709-267':      (3, 1.45, 0.05),
            'XTE J1710-281':      (1, 1.42, 0.13),
            'IGR J17191-2821':    (1, 1.36, 0.04),
            'XTE J1723-376':      (1, 1.05, 0.02),
            '4U 1728-34':        (43, 1.40, 0.15),
            'MXB 1730-335':      (33, 1.30, 0.05),
            'KS 1731-260':        (6, 1.62, 0.13),
            '4U 1735-444':       (10, 1.37, 0.12),
            'XTE J1739-285':      (1, 1.30, 0.06),
            'SAX J1747.0-2853':   (1, 1.93, 0.06),
            'IGR J17473-2721':    (3, 1.6, 0.5),
            'SLX 1744-300':       (4, 1.45, 0.14),
            'GX 3+1':             (2, 1.458, 0.008),
            'IGR J17480-2446':   (34, 1.21, 0.02),
            'EXO 1745-248':       (7, 1.8, 0.3),
            'SAX J1748.9-2021':  (15, 1.43, 0.08),
            '4U 1746-37':         (5, 1.33, 0.07),
            'SAX J1750.8-2900':   (2, 1.338, 0.008),
            'GRS 1747-312':       (1, 1.34, 0.04),
            'SAX J1806.5-2215':   (1, 1.30, 0.05),
            'GX 17+2':            (9, 1.35, 0.10),
            '4U 1820-303':        (2, 1.45, 0.17),
            'GS 1826-24':        (13, 1.66, 0.11),
            'Ser X-1':           (15, 1.45, 0.08),
            'Aql X-1':            (7, 1.65, 0.10),
            'XB 1916-053':        (1, 1.37, 0.09),
            'XTE J2123-058':      (2, 1.35, 0.06),
            'Cyg X-2':           (54, 1.41, 0.05),
            'SAX J1808.4-3658':   (1, 2.14, 0.03),
            'XTE J1814-338':      (1, 1.86, 0.03),
            'burster':           (-1, 1.42, 0.17), # Mean value for all bursters
            'pulsar':            (-1, 2.00, 0.2) } # Mean value for both pulsars

# Standard anisotropy factors for non-dippers and dippers; the tuple below
# gives \xi_b, \xi_p for each class of source (see sec. 5.4 of Galloway et al. 2020)

ANISOTROPY = {'non-dipper': (0.898, 0.809),
              'dipper': (1.639, 7.27)}

# List of ultra compacts based on In 't Zand (2007)
# Includes all candidates.  Updated as of MINBAR source list v2.6
# Can also generate using the Sources object, with method .type('C')
UCXBS = ['4U 0513-40', '4U 0614+09', '2S 0918-549', '4U 1246-588',
           '4U 1543-624', 'IGR J17062-6143', '4U 1705-32',
           'XTE J1709-267', 'SAX J1712.6-3739', 'RX J1718.4-4029',
           'IGR J17254-3257', '4U 1722-30', '4U 1728-34', 'SLX 1735-269',
           'SLX 1737-282', 'IGR J17464-2811', 'SLX 1744-299',
           'XMMU J181227.8-181234', '4U 1812-12', '4U 1820-303',
           'XB 1832-330', '4U 1850-086', 'XB 1905+000', 'XB 1916-053',
           'M15 X-2']

# Set up for publication-quality plots

plt.rcParams.update({
    "text.usetex": True,
    "font.family": "Times"
})


def create_logger():
    """
    Create a logger instance where messages are sent
    See https://docs.python.org/3/library/logging.html
    """
    logger = logging.getLogger(__name__)
    if not logger.handlers: # Check if created before, otherwise a reload will add handlers
        logger.setLevel(logging.INFO)
        handler = logging.StreamHandler()
        formatter = logging.Formatter('%(levelname)s:%(name)s:%(message)s')
        handler.setFormatter(formatter)
        logger.addHandler(handler)
    return logger

logger = create_logger()

def mjd_to_ss(time):
    """
    Converts MJD UTC (as used for the MINBAR burst start times) to RXTE
    Spacecraft seconds, as used for the time-resolved spectroscopic files derived
    from MIT "packet" data
    NOT for use with the FITS data, which use a slightly different convention
    Spacecraft Clock Seconds (SCCS), or Mission Elapsed Time (MET), will represent
    true elapsed seconds since January 1, 1994, at 0h0m0s UTC, which corresponds to
    MJD = 49353.0 (UTC)
    see https://heasarc.gsfc.nasa.gov/docs/xte/abc/time_tutorial.html
    Conversion below is chosen for consistency with the IDL tconv routine, as well
    as the xtime page, with the "Apply Clock Offset Correction(s) for RXTE and Swift"
    option; see https://heasarc.gsfc.nasa.gov/cgi-bin/Tools/xTime/xTime.pl
    :param time:
    :return:
    """

    _time = Time(time, scale='utc', format='mjd')

    # Omitted the 3.378 s offset to line up the WFC and PCA time-resolved
    # spectral files
    return (_time.tt.mjd - 49353.000696574074) * 86400. #- 3.378431

def verify_path(source, path, source_path, verbose=True):
    """
    Generic routine to verify the data path, and try to match up the source names with the
    directories in the tree
    Should be incorporated into the Instrument class
    :param source:
    :param path:
    :param source_path:
    :return:
    """

    has_dir = [False] * len(source)

    if path is None: # flag for nonexistent source data path
        return has_dir, 0, len(source)

    nonmatched = 0
    with os.scandir('/'.join([MINBAR_ROOT, path, 'data'])) as entries:
        for entry in entries:
            # print (entry, entry.name)
            if os.path.isdir(entry) and not (entry.name in source_path):
                if verbose:
                    logger.warning("possible non-compliant source path {}".format(entry.name))
                nonmatched += 1
            else:
                # self.has_dir[self.source_path.index(entry.name)] = True
                _i = np.where(source_path == entry.name)[0]
                # print (entry.name,_i)
                assert len(_i) <= 1
                # print (entry.name, _i)
                if len(_i) == 1:
                    has_dir[_i[0]] = True

    nmissing = len(np.where(np.logical_not(has_dir))[0])
    if (nonmatched > 0) & verbose:  # or self.nmissing > 0:
        print("""
Possible inconsistencies in the directory tree; 
  {} directories without source matches
  {} sources without data directories
You should check your source to directory name mapping""".format(nonmatched, nmissing))

    return has_dir, nonmatched, nmissing


class Minbar(IDLDatabase):
    """
    This class is provided to access EITHER the IDL database or the information via the MRT tables
    Gradually moving the methods from the Burst class here, so they can also be used by the
    Observations class
    Behaviour is controlled by the selection attribute, which indicates which of the table entries
    are selected for display or data extraction
    Ordering is persistent and is controlled by the order attribute, by default time ordering
    """

    def __init__(self, filename=None, IDL=False):

        if filename==None:
            # If no filename specified, then choose the bursts
            filename = self.get_default_path('minbar')

        if not IDL:
            filename += '.txt'

        if IDL:
            # Read in the IDL version of the database
            IDLDatabase.__init__(self, filename)

            self.fix_labels()

        else:
            # read in the MRT version, and fake the rest of the attributes; need to return
            # a recarray, as described here:
            #   https://docs.scipy.org/doc/numpy/reference/generated/numpy.recarray.html
            # attributes to generate:
            #   field_labels - dictionary of names, descriptions
            #   field_names - just the field names
            #   field_format - not used

            data = ascii.read(filename)
            # astropy seems to treat bytes as str, so can't convert strings back to byte arrays
            # for compatibility with the IDL version
            self.records = data

            self.field_names = data.colnames
            info = data.info(out=None)
            self.field_labels = dict(zip(self.field_names, info['description']))

        # Keep the flag so we know what kind of data we're using
        self.IDL = IDL

        # Generate the list of names
        self.names = self.get_names()

        # Define an index; see https://docs.astropy.org/en/stable/table/indexing.html
        self.records.add_index('entry')

        # Also want to determine if we have local data
        # By default, MINBAR includes data from RXTE/PCA, BeppoSAX/WFC, and INTEGRAL/JEM-X
        # so define those instruments here. We pass the names to avoid having to create
        # a Sources object each time (noting that this limits the source name list to those
        # sources with bursts in MINBAR)

        self.instruments = {'XP': Instrument('PCA', source_name=self.names),
                            'SW': Instrument('WFC', source_name=self.names),
                            'IJ': Instrument('JEM-X', source_name=self.names)}

        # local_data is set to true if *any* of the source paths are found
        self.local_data = os.path.isdir(MINBAR_ROOT)
        if self.local_data:
            # the path exists, so check if there's any data there
            self.local_data = False
            for key in self.instruments.keys():
                self.local_data = self.local_data | np.any(self.instruments[key].has_dir)

        # set the default attributes for displaying via the show() method
        self.attributes_default = ['entry','name','obsid','instr','sflag']
        if self.entryname == 'burst':
            self.attributes_default += ['time','rexp']
        else:
            self.attributes_default += ['tstart','tstop']


    def show(self, attributes=None, all=False):
        """
        Display the object in a user-friendly way
        :param attributes:
        :param all:
        :return:
        """

        if attributes is None:
            attributes = self.attributes_default
        else:
            # check all the attributes are in the field names
            in_attr_list = True
            for attr in attributes:
                in_attr_list = in_attr_list & (attr in self.field_names)
                # print (attr, in_attr_list)
            if not in_attr_list:
                logger.error("attribute not present in table")
                return
        print (self)
        if all:
            self.records[self.selection][attributes].pprint_all()
        else:
            self.records[self.selection][attributes][self.order].pprint()

    def fix_labels(self):
        """
        Fix the whitespace in the labels when reading in the IDL version.
        """
        pat = re.compile('\s+')
        for k in self.field_labels:
            self.field_labels[k] = re.sub(pat, ' ', self.field_labels[k])

    def get_names(self):
        """
        Returns a list of all source names in the archive. Ordered by right ascension.
        RA is not a part of the MRT tables, so can only do this with the IDL version.
        Should really replace with a read of the FITS table, which is the definitive
        version
        """
        names = np.unique(self.records.field('name').data)

        if self.IDL:
            ra = np.array([self.records.field('ra')[self.records.field('name') == name][0] for name in names])
            ind = np.argsort(ra)

            return names[ind]

        return names

    def get_default_path(self, filename):
        """
        Return the default path of the minbar data files with prefix
        filename
        """
        return os.path.join(os.path.dirname(__file__), 'data', filename)

    def __len__(self):
        """
        Return the number of entries in the current selection.
        """
        return len(self.ind)

    def get_type(self):
        """
        Return an index array selecting the specified burst type (self.type).
        """
        if self.type == None:
            return np.ones(len(self.records), bool)
        else:
            return self.records.field('type') == self.type


    def clear(self):
        """
        Clear the selection. If self.type is not None, only bursts of the given
        type are selected.
        """
        self.name = ''
        self.selection = self.get_type()
        self.order_field = self.timefield
        self.order = np.argsort(self.records[self.selection].field(self.order_field).value)
        self.ind = np.where(self.selection)[0][self.order]


    def select(self, value, attribute='name', exclude=False):
        """
        Set the selection to the source with given name.
        Previously this would reset the previous selection, but now it respects any prior filter
        Allows multiple items and uses implied or; also implements exclusion
        """
        if attribute not in self.field_names:
            logger.error('attribute {} not present in table'.format(attribute))
            return None

        # Handle multiple value items here
        if np.shape(value) != ():
            if len(value) == 0:
                logger.error("can't filter {} on empty list".format(attribute))
                return None
            _selection = (self.records.field(attribute) == value[0])
            # print (value[0], len(_selection))
            for _value in value[1:]:
                _selection = _selection | (self.records.field(attribute) == _value)
                # print (_value, len(_selection))
        else:
            if value == '':
                logger.error("can't filter {} on empty string".format(attribute))
            # Match attribute including existing selection, and burst type (for backwards compatibility)
            _selection = self.records.field(attribute) == value

        # Don't return a null selection
        # Might instead use attr_like in case that's what was meant
        if ~np.any(_selection) | (exclude & np.all(_selection)):
            logger.warning('criteria left no {}s, skipping'.format(self.entryname))
            return None

        if exclude:
            n_select = len(np.where(self.selection)[0])
            self.selection = self.selection & ~_selection & self.get_type()
            action, n_action = 'excluded', n_select-len(np.where(self.selection)[0])
        else:
            self.selection = self.selection & _selection & self.get_type()
            action, n_action = 'selected', len(np.where(self.selection)[0])

        # Retain the time order for the index

        self.order = np.argsort(self.records[self.selection].field(self.order_field).value)
        self.ind = np.where(self.selection)[0][self.order]
        if attribute == 'name':
            self.name = value
            logger.info('{} {} {}s from {}'.format(action, n_action, self.entryname, self.name))
        else:
            logger.info('{} {} {}s with {}={}'.format(action, n_action, self.entryname, attribute, value))

        # Return self so we can "cascade" selections

        return self

    def sort(self, attribute='time', ascending=True, descending=False):
        """
        Basic sort functionality; will modify the default ordering of the selection,
        until another sort parameter is chosen, or a call to clear() is made
        (which replaces the default sort field as time)

        Usage:
        b.sort('tdel',descending=True) # sort the selected bursts in descending order
                                       #   of the time since the last burst

        :param attribute: (single) table attribute to sort on
        :param ascending: set to False to sort descending instead
        :param descending: set to True to sort descending instead
        :return: sorted object
        """

        if attribute not in self.field_names:
            logger.error('attribute {} not present in table, skipping sort'.format(attribute))
            return self

        self.order_field = attribute
        self.order = np.argsort(self.records[self.selection].field(self.order_field).value)
        if (not ascending) or descending:
            self.order = self.order[::-1]
        self.ind = np.where(self.selection)[0][self.order]

        return self


    def obsid(self, obsid):
        """
        Select observations with a particular obsid
        :param obsid:
        :return:
        """

        self.select(obsid, attribute='obsid')

    def name_like(self, name):
        """
        Like select, but get the name from get_name_like().
        If multiple names are found, the first is selected
        and the rest are reported.
        """
        names = self.get_name_like(name)
        if not names:
            logger.info('No matching source')
        else:
            self.select(names[0])
            if len(names) > 1:
                logger.info('{} more matching sources: {}'.format(len(names) - 1, ', '.join(names[1:])))

        # Return self so we can "cascade" selections

        return self


    def attr_like(self, substring, attribute='name'):
        """
        Multi-purpose function to find attributes matching a particular string
        Output can then be used as input to select
        :param value:
        :param attribute:
        :return:
        """
        to_search = self[attribute]  # includes the selection
        if attribute == 'name':
            to_search = self.names
        elif attribute == 'instr':
            # This is not very efficient for PCA observations, as it will return every distinct
            # instrument code...
            alias = {'pca': 'XP', 'wfc': 'SW', 'jemx': 'IJ'}
            if substring in alias:
                substring = alias[substring]

        selection = []
        for i in to_search:
            if i.find(substring) > -1:
                selection.append(i)

        return list(set(selection))

    def get_name_like(self, name):
        """
        Return a list of sources that have 'name' in their archive
        identifier.
        """
        return self.attr_like(name, 'name')


    def _pad_name(self, name):
        """
        Source names in records are padded with spaces. This routine pads with
        given name with spaces for comparison to the records.
        """
        size = self.records['name'].dtype.itemsize
        return name.ljust(size)


    def get(self, field):
        """
        Get the field as a numpy array while applying the
        current (source) selection and time ordering the
        result.
        Modified to allow getting data from an attribute, if it is
        not available in the standard database.
        Mainly for use by __getitem__
        :param field: field name or array of field names
        :return: astropy MaskedColumn giving the data subject to the
                 current selection
        """
        fields_ok = True
        if np.shape(field) != ():
            for _field in field:
                fields_ok = fields_ok & (_field in self.field_names)

        else:
            fields_ok = field in self.field_names

        if fields_ok:
            # The field method restricts returns to a single attribute
            # return self.get_records().field(field)
            return self.get_records()[field]
        else:
            # Really need to check that the attributes are all present here
            return getattr(self, field)[self.ind]


    def get_records(self):
        """
        Get the time ordered records that are currently
        selected.
        """
        return self.records[self.selection][self.order]


    def __getitem__(self, field):
        """
        Originally this routine was just a shorthand for get(), but has been
        expanded to to also return the item with the entry number; return selected
        items from the current selection; and return
        Note that the ID selections below ignore the data selection, and also
        return astropy Tables (rather than Bursts or Observations objects) so
        you unfortunately can't use the show() method after selecting

        Usage:
        b['time'] - return all the times for the current selection
        b[2094] - return the entry for ID 2094 (ignores current selection)
        b[[2094,2095]] - return the entry for IDs 2094, 2095
        b[['time','tau','bpflux']] - return all the attributes for the current selection
        b[b['time'] > 54000.] - return only those items with the selected properties
        """

        if (np.shape(field) != ()) & (np.shape(field) != (1,)):
            # Array argument
            if type(field[0]) == str:
                return self.get(field)
            elif (type(field[0]) == bool) or (type(field[0]) == bool) or (type(field) == np.ma.core.MaskedArray) or (type(field[0]) == np.bool_):
                # Boolean argument
                return self.get_records()[field]
            elif np.issubdtype(type(field[0]), np.integer):
                # integer arrays, interpret as IDs
                return self.records.loc[field]
            else:
                logger.error("can't interpret argument of type {}, skipping".format(type(field)))
                return None
        else:
            # Scalar argument here
            if type(field) == str:
                return self.get(field)

            # The code below requires that entry is set up as the index
            # return self.records[self.records['entry'] == field]
            result = self.records.loc[field]

            if (type(self) == Observations):
                if (result['nburst'] > 0):
                    # add the burst details to the record, via the meta flag
                    result.meta['bursts'] = self.bursts[self.bursts['entry_obs'] == field]
                else:
                    result.meta['bursts'] = None

            return result

    def instr_like(self, instrument):
        """
        Return an array that selects all entries where the instrument name
        begins with the given instrument string. E.g., 'XP' for PCA. For convenience,
        the following aliases are provided: 'pca', 'wfc', 'jemx'.
        """
        alias = {'pca': 'XP', 'wfc': 'SW', 'jemx': 'IJ'}
        if instrument.lower() in alias:
            instrument = alias[instrument.lower()]

        # This is a crappy way to try to make things work with both strings and byte arrays
        try:
            return np.char.array(self['instr']).startswith(instrument)
        except:
            return np.char.array(self['instr']).startswith(instrument.encode('ascii'))


    def instr(self, instrument, exclude=False):
        """
        More general-purpose routine to allow selections (and exclusions) on instr
        :param instrument:
        :param exclude:
        :return:
        """

        alias = {'pca': 'XP', 'wfc': 'SW', 'jemx': 'IJ'}

        if np.shape(instrument) != ():
            _instrument = [alias[x.lower()] if x.lower() in alias else x for x in instrument]
            instr_all = []
            for i in _instrument:
                instr_all += self.attr_like(i, 'instr')
            self.select(instr_all, 'instr', exclude=exclude)
            return self
        elif instrument.lower() in alias:
            _instrument = alias[instrument.lower()]
        else:
            _instrument = instrument

        instr_all = self.attr_like(_instrument, 'instr')
        if len(instr_all) > 0:
            self.select(instr_all, 'instr', exclude=exclude)

        return self

    def instr_exclude(self, instrument):
        """
        Wrapper for instr
        :param instrument:
        :return:
        """
        self.instr(instrument, exclude=True)

    def select_all(self, names):
        """
        Select multiple sources with given names. Now redundant
        """
        self.select(names, 'name')

    def exclude(self, name):
        """
        Removes source with given name from the current selection. Replaced by select
        """

        self.select(name, exclude=True)

    def exclude_like(self, name):
        """
        Like exclude, but get the name from get_name_like().
        If multiple names are found, the first is excluded
        and the rest are reported.
        """
        names = self.get_name_like(name)
        if not names:
            logger.info('No matching source')
        else:
            self.exclude(names[0])

            if len(names) > 1:
                logger.info('{} more matching sources: {}'.format(len(names) - 1, ', '.join(names[1:])))

    def exclude_flag(self, flags):
        """
        Removes entries with flags matching one or more labels from the current selection.
        """
        selection = [re.search('['+flags+']',x) is None for x in self.records['sflag']]
        self.selection = np.logical_and(self.selection, selection)
        self.order = np.argsort(self.records[self.selection].field(self.order_field).value)
        if len(np.where(self.selection)[0]) < len(self.ind):
            logger.info('excluded {} {}s by excluding flag(s) {}'.format(len(self.ind)-len(np.where(self.selection)[0]),
                    self.entryname, flags))
        self.ind = np.where(self.selection)[0][self.order]


class Bursts(Minbar):
    """
    Read the MINBAR IDL database and give access to its contents.
    
    Example usage:
    import minbar
    mb = minbar.Bursts()
    mb.name_like('1636') # Select a source using part of its name
    print mb.field_labels # See which fields are available
    time = mb['time'] # Get a field as a numpy array (automatically time-ordered)
    flux = mb['pflux']*1e-9 # Flux in erg/s/cm2
    mb.create_distance_correction() # Include distance information from Sources()
    distance = mb['dist']
    luminosity = flux*mb['distcor'] # Luminosity in erg/s
    
    mb.name_like('1826') # Replace selection by another source
    mb.select_all(['GS 1826-24', '4U 1636-536']) # Select multiple sources; requires exact names
    mb.clear() # Clear the selection so all sources are included
    mb.exclude_like('1636') # Exclude source from selection
    mb.exclude_like('1826') # Now two sources are excluded
    """
    
    timefield = 'time' # The field used for determining time order
    entryname = 'burst'
    
    def __init__(self, filename=None, type=1, IDL=False):
        """
        Create a new Bursts instance using the data from the minbar database.
        
        filename: path to the database files, excluding their extension.
                  By default the minbar database in the directory of this
                  script is used.
        type: burst type. The default, 1, selects all vetted Type I bursts. Setting it
              to None means no type is selected.
        """

        # IDLDatabase.__init__(self, filename)
        Minbar.__init__(self, filename, IDL=IDL)

        if IDL:
            self.type = type
        else:
            self.type = None

        # Among other things, this call sets the selection, order, order_field and ind arrays
        self.clear()


    def get_burst_data(self, id):
        """
        Retrieve time-resolved spectroscopy for this burst
        Replacement for the get_burst_data IDL function
        :param id:
        :return:
        Example usage:
        data1 = b.get_burst_data(2380)
        data2 = b.get_burst_data(1600)
        """

        instr = self[id]['instr']
        if instr[0:2] == 'IJ':
            logger.error('Time-resolved spectroscopy not available for JEM-X bursts')
            return None

        wfcdataroot = 'wfcspec_vs2'

        # First we set the _file variable, with the location of the data
        # This can be the name of a file stored locally, or a URL pointing
        # at the MINBAR website
        if self.local_data:
            # Try to get the file locally

            # Get the source path
            # Can probably package this into the Minbar class for wider use

            _match = self.instruments[instr[0:2]].source_name == self[id]['name']
            if not np.any(_match):
                logger.error('No source directory for this burst')
                return None

            # Different naming conventions for the local files, sadly
            if instr[0:2] == 'SW':
                # Get WFC data
                # File naming convention is a bit messed up in the directory Jean provided
                # - EXO 0748-676 files are EXO0748-678_*
                # - +'s in source names are replaced by p
                _path = '/'.join([MINBAR_ROOT, MINBAR_INSTR_PATH[instr[0:2]], wfcdataroot])

                file_pref = self.instruments[instr[0:2]].source_path[_match][0]
                if file_pref == '4U2129+12':
                    file_pref='M15'
                elif file_pref == 'EXO0748-676':
                    file_pref='EXO0748-678'
                elif file_pref == '4U1246-588':
                    file_pref='A1246-588'
                # this is so fucking painful
                file_pref = str(np.char.replace(file_pref, '+', 'p'))

                _file = '/'.join([_path, file_pref ]) \
                        + '_{}w{}_{:02d}.spec'.format( self[id]['obsid'].zfill(5), instr[2:], self[id]['bnum'])


            elif instr[0:2] == 'XP':
                # Get PCA data
                # this is a bit of a challenge with the historical arrangement of the burst data under the
                # alternate file heirarchy, but we can overcome this with some judicious softlinks
                _path = '/'.join([MINBAR_ROOT, MINBAR_INSTR_PATH[instr[0:2]], 'data',
                                  self.instruments[instr[0:2]].source_path[_match][0],
                                  self[id]['obsid'], 'burst{}'.format(self[id]['bnum']) ])
                _file = '/'.join([_path, 'analysis/bbfit_kabs.log'])

            else:
                logger.error('local files not available for instrument {}'.format(instr))
                return None

            # Check that local files exist
            if not os.path.isfile(_file):
                logger.error('time-resolved spectroscopy file {} not found!'.format(_file))
                return None

        else:
            # Try to get the data remotely
            # URLs for the time-resolved spectroscopic data look like
            # https://burst.sci.monash.edu/minbar/data/trs/0001_trs.dat

            # logger.error('Remote data retrieval not yet implemented')

            _file = MINBAR_URL+'minbar/data/trs/{0:04d}_trs.dat'.format(id)

            # Should also check here that the remote file exists...

        # Now that we've defined the file location, read in the data
        # The format is differrent for different files
        if instr[0:2] == 'SW':
            # Now read in the data...format is
            # (1) MJD interval
            # (2) kT of black body model[keV]
            # (3) 1 sigma error in kT[keV]
            # (4) black body radius R[km for d=10 kpc]
            # (5) 1 sigma error in R[km for d=10 kpc]
            # (6) 3 - 25 keV flux[erg / cm ^ 2 / s]
            # (7) 1 sigma error on 3 - 25 keV flux[erg / cm ^ 2 / s]
            # (8) Unabsorbed bolometric flux of black body[erg / cm ^ 2 / s]
            # (9) Error in unabsorbed bolometric flux, based on delta(chi2) = 2.7[erg / cm ^ 2 / s]
            # (10) chi ^ 2 - red

            _data = pd.read_csv(_file,
                                names=['trange', 'kT', 'kT_err', 'rad', 'rad_err', 'flux_3_25', 'flux_3_25_err',
                                       'flux', 'fluxerr', 'chisq'], sep='\s+')
            nspec = len(_data)
            time = np.zeros(nspec)
            dt = np.zeros(nspec)
            for j in range(nspec):
                tmp = _data['trange'][j].split('-')
                time[j] = float(tmp[0])
                dt[j] = (float(tmp[1]) - time[j]) * 86400.
            _data['dt'] = dt
            # This calculation assumes UTC for the time ranges
            _data['time'] = (time - self[id]['time']) * 86400.

            # These flux errors are FAR too big
            _data['flux'] *= 1e9
            _data['fluxerr'] *= 1e9
            _data['flux_min'] = _data['flux'] - _data['fluxerr']
            _data['flux_max'] = _data['flux'] + _data['fluxerr']
            logger.warning('BeppoSAX/WFC flux is 2-10 keV only!')
            lz = np.where(_data['flux_min'] < 0.)[0]
            if len(lz) > 0:
                _data['flux_min'][lz] = 0.1
            _data['kT_min'] = _data['kT'] - _data['kT_err']
            _data['kT_max'] = _data['kT'] + _data['kT_err']
            _data['rad_min'] = _data['rad'] - _data['rad_err']
            _data['rad_max'] = _data['rad'] + _data['rad_err']

            # Need to be aware of different conventions here for the files; SAX files
            # have radius in units of km/10kpc, while RXTE is (km/10kpc)^2

            _data['rad'] = _data['rad'] ** 2
            _data['rad_min'] = _data['rad_min'] ** 2
            _data['rad_max'] = _data['rad_max'] ** 2

        elif instr[0:2] == 'XP':

            _data = pd.read_csv(_file, comment='#',
                                names=['time', 'r', 're', 'dt', 'nH', 'nH_min', 'nH_max', 'kT', 'kT_min', 'kT_max',
                                       'rad', 'rad_min', 'rad_max', 'chisq', 'rawflux', 'flux', 'flux_min', 'flux_max'],
                                sep='\s+')

            # _data = _path # for testing
            _data['fluxerr'] = 0.5 * (_data['flux_max'] - _data['flux_min'])

            # Need to adjust time from SS to seconds post star time

            _data['time'] -= mjd_to_ss(self[id]['time'])

            # The dt value includes the corrections for deadtime. For consistency with the BeppoSAX data, we
            # also want a "well-behaved" dt array that corresponds to the difference between the time bins

            _data['exp'] = _data['dt']
            _dt = _data['time'].values[1:] - _data['time'].values[:-1]
            _dt = np.append(_dt, _dt[-1])
            assert _data.exp.values[-1] < _dt[-1] # check if our guess is wrong
            _data['dt'] = _dt

        # We'd like to have some information about the data, so add that here
        # This is apparently not the best way to do this necessarily, but we're left
        # with little choice, as the DataFrame is the only thing returned:
        # https://stackoverflow.com/questions/14688306/adding-meta-information-metadata-to-pandas-dataframe

        _data.attrs.update({'file': _file})

        return _data


    def burstplot(self, entry=None, param='flux', bdata=None, show=True, 
        **kwargs):
        """
	General-purpose routine to plot burst data. Would like to be able
	to call this in a number of ways, both with a burst ID from
	MINBAR, but also with a pandas table (as read in with
	get_burst_data, for example). And do a bunch of different plots,
        including the three-panel "in 't Zand" plot (see
	e.g. Fig 3, in 't Zand et al. 2012, A&A 547, A47), as well as the
        three-panel "HR-diagram" version

        TODO add some annotation identifying the burst, somewhere...

        Usage:
        burstplot(burst,param='rad',xlim=[-5,25])
        burstplot(burst,param=['flux','rad','kT','chisq'])

        :param entry: MINBAR burst entry # to plot
        :param param: parameter or list of parameters to plot; special
          'hr' to plot the three-panel HR-diagram version
        :param bdata: alternative table of input data, e.g. for bursts not
          in MINBAR
        :param show: show the figure immediately if True, otherwise
          withhold it (e.g. to add further annotation etc.)
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
    
            has_error = np.all([x in bdata for x in [param+'_min',param+'_max']]) | (param == 'r')
    
            # filter on good data
            _gd = bdata.flux*bdata.fluxerr > 0

            # add an extra value copy here to plot that last step
            ax.step(np.append(bdata.time[_gd].values, 
                bdata.time[_gd][-1:].values+bdata.dt[_gd][-1:].values), 
                np.append(bdata[param][_gd].values,bdata[param][_gd][-1:].values),
                where='post',color=color[param])
            if has_error:
                if param == 'r':
                    yerr = bdata['re'][_gd]
                else:
                    yerr = np.stack((bdata[param][_gd]-bdata[param+'_min'][_gd],
                        bdata[param+'_max'][_gd]-bdata[param][_gd]))
                ax.errorbar(bdata.time[_gd]+bdata.dt[_gd]/2., bdata[param][_gd], yerr,
                             fmt='none',ecolor=color[param])
            ax.set_ylabel(ylabel[param])

        xlabel='Time [s]'

        # Set the label names and colours here. To be passed also to
        # plot_param
        # Might need to set up some custom labels for the different
        # conventions of the SAX and RXTE data
        ylabel = {'r': 'Count rate [s$^{-1}$]',
                  'flux': 'Flux [$10^{-9} \mathrm{erg\,cm^{-2}\,s^{-1}}$]',
                  'kT': 'kT [keV]', 
                  'rad': 'Blackbody normalisation [$(R_{\mathrm{km}}/d_{10\ \mathrm{kpc}})^2$]',
                  'chisq': 'Fit $\chi^2/n_{\mathrm{DOF}}$'}
        color = {'r': 'k', 'flux': 'k', 'kT': 'r', 'rad': 'b', 'chisq': 'g'}
    
        # Get the data here

        if (bdata is None) & (entry is None):
            logger.error('please specify either the burst data or MINBAR ID')
            return
        elif (bdata is None) & (entry is not None):
            bdata = self.get_burst_data(entry)

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
                assert _param in bdata

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

        if show:
            plt.show()

        return fig


    def get_lc(self, id, pre=16., post=None):
        """
        Preliminary routine to return the lightcurve corresponding to a burst
        from the lightcurve for the host observation

        Later this should probably be incorporated into a Burst object or similar
        Usage:
        b = minbar.Bursts()
        mjd, rate, error = b.get_lc(2257)
        :param id: Burst ID to retrieve
        :param pre: pre-burst interval (in seconds) to include
        :param post: post-burst interval (in seconds) to include
        :return:
        """

        mjd, rate, error = None, None, None

        if not (id in self['entry']):
            logger.error('not a valid burst ID')

        t0 = self[id]['time']
        if self.local_data:
            # Try to get the file locally. At the moment this relies on the observation lightcurves,
            # which are in different units; so is not really consistent with the version that gets
            # the lightcurves from the online repo
            logger.warning('getting burst data segment from observation files, careful with units')
            oid = self[id]['entry_obs']
            if post is None:
                post = 5. * self[id]['dur']

            o = minbar.Observations()
            obs = minbar.Observation(obs_entry=o[oid])
            lc = obs.get_lc()
            sel_burst = np.where((obs.mjd.mjd > t0 - (pre * u.s).to('d').value) &
                                 (obs.mjd.mjd < t0 + (post * u.s).to('d').value))[0]

            mjd, rate, error = obs.mjd.mjd[sel_burst], obs.rate[sel_burst], obs.error[sel_burst]

        else:
            # Try to get the data remotely
            # URLs for the data look like
            # https://burst.sci.monash.edu/wiki/uploads/MINBAR/bursts/0001_lc.csv

            # logger.error('Remote data retrieval not yet implemented')

            url = MINBAR_URL+'wiki/uploads/MINBAR/bursts/{0:04d}_lc.csv'.format(id)

            data = pd.read_csv(url)
            if type(data) == pd.DataFrame:
                mjd = (data['time'].values*u.s).to('d') + t0*u.d
                rate = data['flux']*CFLUX_U # bit of a misnomer, units of counts/cm^2/s
                error = data['error']*CFLUX_U
            else:
                logger.error('Some problem with the remote data file')

        return mjd, rate, error

    def PRE(self, rexp_thresh=1.629, marginal=False):
        """
        Select only photospheric radius-expansion (PRE) bursts, according
        to the standard threshold
        Additional marginal flag will optionally include those "marginal"
        events
        :param rexp_thresh:
        :return:
        """

        marginal_incl = ('excluding','including')[int(marginal)]
        selection = self.selection \
            & (self.records.field('rexp') >= rexp_thresh) \
            & ((self.records.field('rexp') < 3.0) | marginal)
        if ~np.any(selection):
            logger.warning('no PRE bursts in current selection')
        else:
            logger.info('selected {} bursts with PRE (rexp > {}), {} marginal cases'.format(
                len(np.where(selection)[0]), rexp_thresh, marginal_incl))
            self.selection = selection

        self.order = np.argsort(self.records[self.selection].field(self.order_field).value)
        self.ind = np.where(self.selection)[0][self.order]

    def unique(self):
        """
        Removes multiply-observed bursts, with a priority for XTE
        """
        instr_label = [x[0:2] for x in self['instr']]
        set_instr_label = set(instr_label)
        if len(set_instr_label.difference({'IJ','SW','XP'})) > 0:
            logger.error('unique filtering may not function correctly with additional instruments')

        selection = np.logical_or(self.records.field('mult') == 1,
                                  np.char.array(self.records['instr']).startswith('XP'))
        self.selection = np.logical_and(self.selection, selection)

        self.order = np.argsort(self.records[self.selection].field(self.order_field).value)
        n_excluded = len(self.ind) - len(np.where(self.selection)[0])
        if n_excluded > 0:
            logger.info('selected {} unique {}s by excluding {}'.format(len(np.where(self.selection)[0]),
                    self.entryname, n_excluded))
        self.ind = np.where(self.selection)[0][self.order]

        # Return self so we can "cascade" selections

        return self

    def __str__(self):
        """
        Return a nice string.
        """
        return "Multi-INstrument Burst ARchive (MINBAR) ({} {}s from {} sources)".format(
            len(self.records), self.entryname, len(self.names) )


    def has_error(self, field):
        """
        Return if there exists a field with error
        data on field.
        """
        error_name = self.get_error_name(field)
        return error_name in self.field_names


    def is_error(self, field):
        """
        Return whether given field represents the
        uncertainty in another.
        """
        if field[-1] == 'e':
            return field[:-1] in self.field_names
        else:
            return False


    def get_error(self, field):
        """
        Return an array representing the uncertainty in
        field. has_error() can be used to check whether
        errors are available. If not available, an array
        of zeros will be returned.
        """
        if self.has_error(field):
            return self.get(self.get_error_name(field))
        else:
            return np.zeros(len(self.get(field)))


    def get_error_name(self, field):
        """
        Return the name of the field containing the
        error for field.
        """
        return field + 'e'


    def create_distance_correction(self):
        """
        Create an array of distance corrections, 4 pi d^2, for each burst.
        This can be used to easily convert flux to luminosity, assuming
        isotropic X-ray emission and neglecting the bolometric correction
        
	Returns correction factor (cm^2), error, distance (kpc, from the
	Source table), error. Furthermore, these arrays are stored and can
	be accessed as self['fieldname'], with fieldnames distcor,
        distcore, dist, diste.
        """
        s = Sources()
        dist = np.zeros(len(self.records))
        diste = np.zeros_like(dist)
        cor = np.zeros_like(dist)
        core = np.zeros_like(dist)
        
        names = self.records.field('name')
        for name in self.names:
            s.name_like(name.strip(), verbose=False)
            if s.selection!=None:
                ind = names==name
                dist[ind] = s['dist']
                diste[ind] = s['diste']
            s.clear()
        
        cor = 4*np.pi*(dist*kpc)**2
        ind = dist>0.0 # Prevent division by 0 in next line
        core[ind] = cor[ind]*2.0*diste[ind]/dist[ind]
        
        self.dist = dist
        self.diste = diste
        self.distcor = cor
        self.distcore = core
        
        return cor, core, dist, diste
    
    def clean(self, ucxbs=True, high_mdot=True, gx354=True, rapidburster=True, clear=True, fishy=True):
        """
        Exclude (groups of) sources, such as UCXBs, sources with high mdot (Cyg X-2
        and GX 17+2), GX 254-0 (aka 4U 1728-34, frequent burster and possibly a UCXB)
        and the Rapid Burster.
        """
        if clear:
            self.clear()
        
        if rapidburster:
            self.exclude('MXB 1730-335') # Rapid burster
        
        if gx354:
            self.exclude('4U 1728-34') # GX 354-0, close to rapid burster and possible UCXB
        
        if ucxbs:
            # Exclude (candidate) UCXBs
            for name in UCXBS:
                self.exclude(name)
        
        if high_mdot:
            # Exclude systems that accrete near the Eddington limit
            self.exclude('GX 17+2')
            self.exclude('Cyg X-2') # type-II bursts?
        
        if fishy:
            self.exclude('4U 1746-37') # Possibly 2 sources
            self.exclude('EXO 1745-248') # Possibly Type II


class Observations(Minbar):
    """
    Load MINBAR database of observations.
    """
    
    timefield = 'tstart' # The field used for determining time order
    entryname = 'observation'


    def __init__(self, filename=None, type=None, IDL=False, bursts=True, verbose=True):
        """
        Load the database of observations.
        """
        if filename==None:
            filename = self.get_default_path('minbar-obs')

        if verbose:
            logger.info('loading observations, please wait...')

        Minbar.__init__(self, filename, IDL=IDL)

        if bursts:
            # by default also import the bursts
            # I think filename is only used for the IDL option, but try to fix that here
            bfilename = filename
            if bfilename is not None:
                bfilename = bfilename[:-4] # drop the -obs part
            self.bursts = Bursts(filename=bfilename, IDL=IDL)
        else:
            self.bursts = None

        self.type = type

        # Among other things, this call sets the selection, order, order_field and ind arrays
        self.clear()

    def good(self):
        """
        Filter for only the "good" observations, excluding bad flags and non detections
        :return:
        """
        self.exclude_flag('bcdefg')
        selection = (self.records['flux'] > 3.*self.records['e_flux']) & (self.records['sig'] >= 3.)
        self.selection = np.logical_and(self.selection, selection)
        self.order = np.argsort(self.records[self.selection].field(self.order_field).value)
        if len(np.where(self.selection)[0]) < len(self.ind):
            logger.info('Restricted to {} "good" {}s by also excluding nondetections'.format(
                len(np.where(self.selection)[0]), self.entryname))
        self.ind = np.where(self.selection)[0][self.order]

    def __str__(self):
        """
        Return a nice string.
        """
        return "Multi-INstrument oBservation ARchive (MINBAR) ({} {}s from {} sources)".format(
            len(self.records), self.entryname, len(self.names))

    def plot(self, bursts=True, entry=None, lightcurve=None, **kwargs):
        """
        Simple plotting interface for MINBAR observations, useful for
        producing (for example) summary plots of source behaviour over the
        MINBAR sample

        This routine works off the fluxes in the current selection, which
        could include multiple sources (no current method to plot for >1
        sources)

        :param bursts: boolean flag to toggle plotting of bursts
        :param entry: alternative selection method, by just passing an
          array of the observation ID numbers
        :param lightcurve: boolean to trigger plotting of the high-time
          resolution lightcurves. By default, this will plot the
          lightcurves if there are less than 10 observations or they cover
          10d or less

        Example usage:
        import minbar
        o = minbar.Observations()
        o.select('4U 1636-536')
        o.select(0,'flux',exclude=True) # don't plot zero fluxes
        o.plot()

        # or alternatively, for individual observations
        o.plot(entry=[14637,14639])
        """

        if (entry is not None):
            sel_old = self.selection # store to preserve the selection
            self.clear()
            self.select(entry, 'entry')

        # by default, we adopt the flux range for the table data; if we're
        # reading in lightcurves later on this may be modified
        yrange = (min(self['flux']-self['e_flux']), 
                max(self['flux']+self['e_flux']))

        if (len(set(self['name'])) > 1):
            logger.warning('multiple sources in current selection, can\'t plot')
            return
        src = self['name'][0]

        fig = plt.figure()

        extent = max(self['tstop'])-min(self['tstart'])
        if lightcurve is None:
            lightcurve = self.local_data & ((len(self['entry']) <= 10) | (extent <= 10.))

        instr = np.array([x[:2] for x in self['instr']])
        label = {'IJ': '{\it INTEGRAL}/JEM-X', 'SW': '{\it BeppoSAX}/WFC',
            'XP': '{\it RXTE}/PCA'}
        colors = {'IJ': 'C0', 'SW': 'C1', 'XP': 'C2'}

        for _instr in set(instr):
            _s = instr == _instr
            if lightcurve:
                # plot individual high-time resolution lightcurves
                # not yet tested
                _label = label[_instr]
                for _id in self['entry'][_s]:
                    _obs = Observation(self[_id])
                    _obs.plot(fig, show=False, show_bursts=False, 
                        label=_label, color=colors[_instr])
                    # have to calculate the yrange progressively here with
                    # each lightcurve
                    yrange = (min([yrange[0], min(_obs.rate.value)]),
                        max([yrange[1], max(_obs.rate.value)]))
                    _label = None
            else:
                # just plot the averaged fluxes over the entire observation
                plt.errorbar(0.5*(self['tstart'][_s]+self['tstop'][_s]),
                    self['flux'][_s], yerr = self['e_flux'][_s],
                    xerr = 0.5*(self['tstop'][_s]-self['tstart'][_s]),
                    fmt='.', label=label[_instr], color=colors[_instr])

        # need a reverse lookup for the bursts; this is pretty slow, but I
        # can't think of a better way to do it
        # Actually this is already implemented in select!

        if bursts:
            self.bursts.clear()
            self.bursts.select(self['entry'], 'entry_obs')

            nburst = len(self.bursts['time'])
            if nburst > 0:
                burst_pos = yrange[1]*1.05 - 0.05*yrange[0]
                plt.plot(self.bursts['time'],
                    np.full(len(self.bursts['time']), burst_pos),'|r',
                    label='type-I bursts')
            else:
                logger.info('no bursts to show in current selection')
            self.bursts.clear()

        plt.xlabel('Time [MJD]')
        plt.ylabel('Flux [3-25 keV, $10^{-9}\, \mathrm{erg\,cm^{-2}\,s^{-1}}$]')
        if (nburst > 0) | (len(set(instr)) > 1):
            plt.legend()

        plt.show()

        if entry is not None:
            # restore the original selection
            self.selection = sel_old

        return fig


class Observation:
    """
    This object is intended to allow all the possible actions you might have on an
    observation. You can create it from a minbar entry, or given an instrument, source name and obs ID

    Example usage:
    import minbar
    o = minbar.Observations()
    obs = minbar.Observation(o[14637])
    """


    def __init__(self, obs_entry=None, instr=None, name=None, obsid=None):
        """
        Create an observation instance, either from a MINBAR obs entry, or by-hand
        Ideally this object should make available every parameter in the MINBAR observation table
        It'd be nice to be able to create this just given the obs ID (for example), but
        then we'd need to keep a copy of the Observations object, which seems wasteful

        :param obs_entry:
        :param instr: 
        :param source: 
        :param obsid: 
        """

        if obs_entry is not None:
            # logger.warning('initialisation from obs entry not yet completely implemented')

            # create an instrument here
            label = [key for key, value in MINBAR_INSTR_LABEL.items() if value == obs_entry['instr'][:2]][0]
            # print (label)

            # Don't set these parameters yet, as they'll be defined outside this block
            instr = Instrument(label, obs_entry['instr'])
            name = obs_entry['name']
            obsid = obs_entry['obsid']

            # self.tstart = obs_entry['tstart']
            # self.tstop = obs_entry['tstop']
            # Copy all the columns to the new object. NOTE this will includ
            # the instrument label, which will be overwritten with the 
            # instrument object, below; but we keep the label by specifying it
            # as a configuration/camera above
            for col in obs_entry.columns:
                if col in UNITS.keys():
                    setattr(self, col, obs_entry[col]*UNITS[col])
                else:
                    setattr(self, col, obs_entry[col])

            # add any bursts present
            if self.nburst > 0:
                # self.bursts =
                self.bursts = obs_entry.meta['bursts']

        else:
            # this observation isn't in MINBAR, so it has no entry
            self.entry = None

        # Potentially need to check here that all of the passed parameters are set
        self.instr = instr
        self.name = name
        self.obsid = obsid

        # Define parameters for the lightcurve; later this might be a class
        # The lightcurve might also not be available, so don't force it to be read in now

        self.time = None
        self.rate = None
        self.error = None


    def __str__(self):
        """
        Return a nice string.
        """

        output = "MINBAR observation of {}\nInstrument: {}\nObsID: {}".format(
            self.name, self.instr.name, self.obsid)
        if hasattr(self,'tstart') & hasattr(self,'tstop'):
            output += "\nTime range: MJD {}-{}".format(self.tstart.value, self.tstop)
        _path = self.get_path()
        if _path is not None:
            output += "\nData path: {}".format(_path)

        return output


    def plot(self, figure=None, show=True, show_bursts=True, sym_burst='^r',
            **kwargs):
        """
        Plot the observation lightcurve, reading it in first if need be.
        Keyword arguments are passed on to the plot command (see the example
        below).
        :param figure: existing figure to add to, if multiple observations are
          to be plotted on the same axis (for example)
        :param show: display the figure immediately or not (latter case for 
          multiple observations to be plotted together)
        :return: plot object

        Example usage, showing a combined plot of two subsequent observations
        from 4U 1254-69:

        obs1=mb.Observation(o[17920])
        obs2=mb.Observation(o[17921], color='C0')
        fig = obs1.plot(show=False)
        obs2.plot(fig)

        """
        ylabel = 'Rate (count s$^{-1}$ cm$^{-2}$)'

        if self.time is None:
            self.get_lc()

        if figure is None:
            figure = plt.figure()

        if self.time is None:
            logger.info('Showing schematic plot for flux')
            ylabel = 'Flux (3-25 keV, $10^{-9}\ \rm{erg\,cm^{-2}\,s^{-1}$)'
            plt.errorbar(0.5*(self.tstart+self.tstop), self.flux, 
                xerr=0.5*(self.tstop-self.tstart), yerr=self.e_flux, **kwargs)
            rate_max = (self.flux+self.e_flux).value
            rate_range = rate_max
            plt.ylim((0,rate_max+rate_range*0.1))
            tscale = ''

        else:
            # plt.plot(lc['TIME'],lc['RATE'])
            # plot can't work with "raw" time units
            # the duplication of .mjd is not a typo below! Also the scale method
            # is for Time objects, which mjd should be defined as
            plt.plot(self.mjd.mjd, self.rate, **kwargs)

            rate_max = np.nanmax(self.rate)
            rate_range = rate_max-np.nanmin(self.rate)
            tscale = self.mjd.scale.upper()

        plt.xlabel('Time (MJD '+tscale+')')
        # plt.ylabel(ylabel)

        if hasattr(self, 'bursts') & show_bursts:
            # also plot the bursts
            plt.plot(self.bursts['time'], 
                rate_max+0.05*rate_range * np.full(len(self.bursts), 1), sym_burst)
        if show:
            plt.show()

        return figure


    def get_path(self, split=False):
        """
        Return the path for MINBAR observations, assuming you have them stored locally
        :param entry:
        :return:
        """

        instr = self.instr.label
        if not self.instr.local_data:
            logger.warning('no local data is present, check your MINBAR_ROOT')
            return None

        _match = np.where(self.instr.source_name == self.name)[0]
        # print (_match)
        assert len(_match) == 1
        # if len(_match) > 1:
            # logger.warning("multiple source name matches for path")
        if self.instr.label == 'SW':
            # WFC doesn't have individual directories for sources
            path = '/'.join([MINBAR_ROOT, self.instr.path, 'data'])
            if os.path.isdir(path):
                return path
        elif self.instr.has_dir[_match[0]]:
            # this only implies that there is a source directory for this
            # object, NOT that the obs directory exists
            path = '/'.join([MINBAR_ROOT, self.instr.path, 'data',
                     self.instr.source_path[_match[0]],
                     self.obsid])
            if os.path.isdir(path):
                # the "split" option is only relevant for PCA/XTE data
                if (not split) or (instr != 'XP'):
                    return path
                else:
                    # PCA data are sometimes split over multiple paths, so
                    # try to find those here
                    i, path_arr = 0, [path]
                    while (os.path.isdir(path+str(i))): 
                        path_arr.append(path+str(i))
                        i +=1
                    # if no additional directories are found, this routine will
                    # still return a (single-element) list, but I think that's
                    # OK
                    return path_arr
            else:
                logger.warning('directory {} not found locally'.format(path))
                return None
  

        return None


    def get_lc(self):
        """
        Return the lightcurve for a particular observation; this is a replacement for the IDL routine get_lc.pro
        This routine also populates the time, mjd_tt, mjd, rate, and error attributes for the observation
        :param entry:
        :return:
        """

        def pca_time_to_mjd_tt(time, header):
            """
            convert raw times to MJD (TT) here; see https://heasarc.gsfc.nasa.gov/docs/xte/abc/time_tutorial.html
            can check results using https://heasarc.gsfc.nasa.gov/cgi-bin/Tools/xTime
            """

            return Time( time.to('d') + (header['MJDREFI'] + header['MJDREFF'])*u.d, format='mjd', scale='tt') # TT

        def read_fits_lc(file, effarea = 1.*u.cm**2):
            """
            Utility routine to read in the important bits of a lightcurve

            :param file: name of file to read in
            :param effarea: effective area adopted to convert count rate to
              counts/cm^2/s (standard for MINBAR lightcurves)

            :return: time, rate, error, timesys, timeunit, mjd, gti
            """

            # print (path+'/'+filename)
            lcfile = fits.open(file)

            # For XTE files, convention is to have the first extension RATE and a second extension STDGTI
            header = lcfile[0].header
            if 'INSTRUME' not in header:
                # For WFC, TIMESYS etc. are in the first extension header
                header = lcfile[1].header
            instrument = header['INSTRUME']
            lc = lcfile[1].data

            gti_ext = None
            if len(lcfile) > 2:
                gti_ext = lcfile[2].data

            lcfile.close()

            # Can clean up the table here if necessary; i.e. create a Lightcurve
            # object (not yet defined), adopt uniform time scale etc.
            # Using astropy to keep track of the time scale and units, read from the
            # header file; see https://docs.astropy.org/en/stable/time

            timesys = header['TIMESYS']
            timeunit = header['TIMEUNIT']
            time = lc['TIME']*u.Unit(timeunit)

            rate = lc['RATE']/u.s/effarea
            error = lc['ERROR']/u.s/effarea

            # print (header['INSTRUME'])
            if instrument == 'PCA':

                mjd_tt = pca_time_to_mjd_tt(time, header)
                mjd = mjd_tt.utc

                # not sure if GTI arrays exist for other instruments
                gti = np.array([])
                for _gti in gti_ext:
                    gti = np.append(gti, _gti)#, axis=1)

                return time, rate, error, timesys, timeunit, mjd, pca_time_to_mjd_tt(gti*u.Unit(timeunit), header).utc
            else:
                if 'MJDREF' in header:
                    time += header['MJDREF']*u.d
                mjd = Time( time, format='mjd', scale=timesys.lower())
                if timesys == 'TT':
                    mjd = mjd.utc

                return time, rate, error, timesys, timeunit, mjd, None

        path = self.get_path()
        if path is None:
            # indicates no local data present, so skip the read
            return None

        add_files = []
        if self.instr.label == 'XP':
            # The XTE path uses a function to get the lightcurve name
            filename = self.instr.lightcurve(self.obsid)
            logger.warning('PCA lightcurves are standard products which may not show all bursts')

            # there is also the possibility of additional paths, where the
            # observation data is split over
            add_paths = self.get_path(split=True)[1:]
            for _path in add_paths:
                add_files.append(self.instr.lightcurve(self.obsid+_path[-1:]))
            # print (add_paths, add_files)
        elif self.instr.label == 'SW':
            # WFC also uses a (different) function to get the name
            # also have to translate to the name set used in the filenames
            __name = self.instr.source_path[np.where(self.instr.source_name == self.name)[0]][0]
            filename = self.instr.lightcurve(__name, self.obsid, self.instr.instr[2:])
        elif self.instr.label == 'IJ':
            # JEM-X can have multiple lightcurves, but that's not implemented
            # yet. Instead we pick one of them for the JEM-X1 or -X2 here
            filename = self.instr.lightcurve[int(self.instr.instr[2:])-1]
        else:
            # logger.warning("other instruments not yet implemented")
            filename = self.instr.lightcurve

        self.file = '/'.join((path, filename))
        self.time, self.rate, self.error, self.timesys, self.timeunit, self.mjd, self.gti = read_fits_lc(self.file, self.instr.effarea)

        for i, file in enumerate(add_files):
            # here we read in any additional files and append them to the already-created arrays
            _time, _rate, _error, _timesys, _timeunit, _mjd, _gti = read_fits_lc('/'.join((add_paths[i], file)), self.instr.effarea)
            self.file = np.append(self.file, '/'.join((add_paths[i], file)))
            # The ordering here will not necessarily be in time
            self.time = _time.insert(0, self.time)
            self.rate = _rate.insert(0, self.rate)
            self.error = _error.insert(0, self.error)
            self.mjd = _mjd.insert(0, self.mjd)
            assert _timeunit == self.timeunit
            assert _timesys == self.timesys
            if _gti is not None:
                self.gti = _gti.insert(0, self.gti)
        if self.gti is not None:
            self.gti = np.reshape(self.gti.sort(), (-1,2))
        i = self.time.argsort()
        self.time = self.time[i]
        self.rate = self.rate[i]
        self.error = self.error[i]
        self.mjd = self.mjd[i]

        return #lc


class Sources:
    """
    Contains all information on the sources present in MINBAR, via
    the file minbar_sources.fits.
    
    Example:
    s = Sources()
    print s.field_names # Show available data fields
    ra = s['ra_obj'] # Right ascension for all sources
    s.name_like('1636')
    ra = s['ra_obj'] # Right ascension for selected source only
    s.clear() # Clear selection
    """
    
    def __init__(self, filename=None, X=0.0, Gaia=True):
        """
        Load source list from FITS file.
        """
        if filename==None:
            filename = self.get_default_path()
        self._f = fits.open(filename)
        self.header = self._f[1].header
        self.field_names = [i.lower() for i in self._f[1].data.dtype.names]
        self.field_labels = self._get_field_labels()
        self._fits_names = list(self.field_names) # Keep track of which fields are in fits file
        self.clear()
        # Now has a bit more information about the distances
        # self.dist, self.diste = self.get_distances(X=X, Gaia=Gaia)
        self.dist, self.diste, self.diste_lo, self.dist_method = self.get_distances(X=X, Gaia=Gaia)
        self.field_names += ['dist', 'diste','diste_lo','dist_method']
        self.field_labels['dist'] = 'Distance (kpc)'
        self.field_labels['diste'] = 'Error on distance (kpc)'
        self.field_labels['diste_lo'] = 'Lower error on distance (kpc, where defined)'
        self.field_labels['diste_method'] = 'Distance method'
        self.X = X
        self.Gaia = Gaia

        # Extract the version information

        version_string = [x for x in self.header['history'] if re.search('minbar_sources.fits', x)][0]
        self.version = version_string[26:]

        # Get the Eddington flux information. This step is kind of a hack since
        # the fluxes should really be included in the source FITS file.

        F_Edd = self.get_F_Edd()
        if F_Edd is None:
            self.local_files = False
        else:
            self.local_files = True
            self.F_Edd, self.F_Edd_err = F_Edd


    def get_default_path(self):
        """
        Return the default path of the source list
        """
        return os.path.join(os.path.dirname(__file__), 'data', 'minbar_sources.fits')


    def _get_field_labels(self):
        """
        Get the field labels from the comment fields in the fits header
        """
        columns = [field[5:] for field in self.header if field.startswith('TTYPE')]
        field_labels = {}
        for column in columns:
            label = self.header.get('TCOMM'+column, self.header['TTYPE'+column])
            unit = self.header.get('TUNIT'+column, '')
            if unit:
                label = '{} ({})'.format(label, unit)
            field_labels[self.header['TTYPE'+column].lower()] = label
        return field_labels


    def __len__(self):
        """
        Return the number of entries
        """
        return self._f[1].data.shape[0]


    def get(self, field, all=False):
        """
        Return field with given name.
        all: whether to return all values or only the selected source
        """
        if field.lower() in self._fits_names:
            data = self._f[1].data[field]
        else: # If not in the fits file, see if it is an attribute
            data = getattr(self, field)
        
        if all or self.selection is None:
            return data
        else:
            return data[self.selection]


    def __getitem__(self, field):
        """
        Return field with given name. See self.get().
        """
        return self.get(field)

    def type(self, types):
        """
        Select only objects matching a particular type code
        :param types:
        :return:
        """
        type_names = {'atoll': 'A', 'atolls': 'A', 'ultracompact': 'C', 'UCXB': 'C', 'dipper': 'D',
                      'dippers': 'D', 'eclipsing': 'E', 'globular': 'G', 'cluster': 'G',
                      'intermittent': 'I', 'microquasar': 'M', 'oscillation': 'O', 'pulsar': 'P',
                      'pulsars': 'P', 'radio': 'R', 'superburst': 'S', 'transient': 'T',
                      'transients': 'T', 'burster': 'B', 'bursters': 'B'}
        if types in type_names:
            types=type_names[types]

        src_types = self._f[1].data['type']
        sel = np.array([types[:1] in x for x in src_types])
        for type in types[1:]:
            sel = (sel & [type in x for x in src_types])

        self.selection = sel

        return len(np.where(sel)[0])

    def get_name_like(self, name):
        """
        Return a list of source indices that have 'name' in their name or name_2 fields.
        Case insensitive.
        """
        name = name.lower()
        selection = []
        for i, (name1, name2) in enumerate(zip(self._f[1].data['name'], self._f[1].data['name_2'])):
            if name1.lower().find(name)>-1:
                selection.append(i)
            elif name2.lower().find(name)>-1:
                selection.append(i)
        return np.array(selection)


    def name_like(self, name, verbose=True):
        """
        Select the source with given name. Uses first result from self.get_name_like()
        """
        ind = self.get_name_like(name)
        if len(ind)>0:
            self.selection = ind[0]
            if verbose:
                logger.info('Selected source {}'.format(self['name']))
                if len(ind)>1:
                    logger.info('{} more matching sources: {}'.format(len(ind) - 1, ', '.join(self.get('name', True)[ind[1:]])))
        else:
            logger.info('No matching source')


    def clear(self):
        """
        Clear current source selection
        """
        self.selection = None


    def get_F_Edd(self):
        """
        Read in information from Eddington fluxes file
        :return:
        """

        self.EddingtonFluxes = MINBAR_ROOT + '/EddingtonFluxes.dat'
        if not os.path.isfile(self.EddingtonFluxes):
            logger.warning('can\'t read file {}, Eddington flux values unavailabile'.format(self.EddingtonFluxes))
            return None

        d = pd.read_csv(self.EddingtonFluxes, header=None, engine='python',
                        names=["BursterNo", "SourceName", "RXTEname", "F_Edd", "F_Edd_Err", "Chisqdf",
                               "Timestamp", "PCA", "SWFC", "LBC", "flag0", "flag1", "flag2", "flag3", "flag4", "flag5",
                               "flag6", "flag7", "flag8"], sep=', ')
        d.set_index('SourceName', inplace=True)

        # Map Eddington fluxes onto the name array, as for the distances

        F_Edd = np.zeros_like(self['ra_obj'])
        F_Edd_Err = np.zeros_like(F_Edd)
        for i, name in enumerate(self['name']):
            if name in d.index:
                F_Edd[i] = d.loc[name]['F_Edd']
                F_Edd_Err[i] = d.loc[name]['F_Edd_Err']
                # print(name, F_Edd[i], F_Edd_Err[i])

        return F_Edd, F_Edd_Err


    def get_distances(self, X=0.0, Gaia=True):
        """
        Create an array of distances for all sources in self.name
        
        Distances taken from Table 8 of the MINBAR paper, which includes
        distances derived from the Eddington flux measured for MINBAR
        bursts, as well as distances measured from Gaia and other sources
        """

        idist = 3 # index for distances to adopt
        if (X != 0.7) & (X != 0.0):
            logger.error('distances not defined for other than X=(0.0,0.7)')
            return None
        elif X == 0.7:
            idist = 2

        # The table below defines a "distance tuple" for each source; the components are:
        # 1. the average peak luminosity and 1-sigma error (in 1E-9 erg/cm^2/s);
        # 2. the "measured" distance in kpc, with two-sided error (upper and lower), and a string
        #    indicating the source;
        # 3. the distance in kpc inferred from the peak Eddington flux for X=0.7, and 1-sigma error;
        # 4. the distance in kpc inferred from the peak Eddington flux for X=0.0, and 1-sigma error
        # The flags to this routine are used to populate the (simpler) array distances

        table8 = {  '4U 0513-40': ((14.4, 6.7), (10.32, 0.2, -0.24, '1'), (8.5, 1.5), (11.1, 1.9)),
                    '4U 0614+09': ((266.0, 6.0), (3.27, 2.42, -1.3, 'G'), (1.99, 0.02), (2.59, 0.03)),
                    'EXO 0748-676': ((46.5, 4.3), (0.0, 0.0, -0.0, ''), (4.7, 0.2), (6.2, 0.3)),
                    '4U 0836-429': ((0.0, 0.0), (3.18, 2.25, -1.4, 'G'), (0.0, 6.9), (0.0, 9.0)),
                    '2S 0918-549': ((119.1, 14.4), (5.77, 2.77, -1.6, 'G'), (3.0, 0.2), (3.9, 0.2)),
                    '4U 1246-588': ((120.3, 11.9), (2.03, 2.37, -1.17, 'G'), (3.0, 0.1), (3.8, 0.2)),
                    '4U 1254-69': ((0.0, 0.0), (3.18, 3.16, -1.33, 'G'), (0.0, 6.0), (0.0, 7.9)),
                    'SAX J1324.5-6313': ((0.0, 0.0), (0.0, 0.0, -0.0, ''), (0.0, 4.7), (0.0, 6.1)),
                    '4U 1323-62': ((0.0, 0.0), (0.0, 0.0, -0.0, ''), (0.0, 5.2), (0.0, 6.8)),
                    'Cir X-1': ((0.0, 0.0), (6.17, 2.86, -1.96, 'G'), (0.0, 13.6), (0.0, 17.7)),
                    '4U 1608-522': ((169.0, 41.2), (5.8, 1.8, -2.0, '2'), (2.5, 0.3), (3.2, 0.3)),
                    '4U 1636-536': ((72.5, 18.8), (4.42, 3.08, -1.63, 'G'), (3.8, 0.4), (5.0, 0.5)),
                    'XTE J1701-462': ((43.4, 1.4), (0.0, 0.0, -0.0, ''), (4.9, 0.1), (6.4, 0.1)),
                    'MXB 1658-298': ((17.0, 15.9), (0.0, 0.0, -0.0, ''), (7.9, 2.2), (10.2, 2.9)),
                    '4U 1702-429': ((88.7, 45.0), (0.0, 0.0, -0.0, ''), (3.4, 0.6), (4.5, 0.8)),
                    '4U 1705-32': ((0.0, 0.0), (0.0, 0.0, -0.0, ''), (0.0, 4.4), (0.0, 5.7)),
                    '4U 1705-44': ((41.3, 17.5), (0.0, 0.0, -0.0, ''), (5.0, 0.8), (6.6, 1.1)),
                    'XTE J1709-267': ((0.0, 0.0), (0.0, 0.0, -0.0, ''), (0.0, 2.7), (0.0, 3.6)),
                    'XTE J1710-281': ((7.1, 1.8), (0.0, 0.0, -0.0, ''), (12.2, 1.3), (15.9, 1.7)),
                    'SAX J1712.6-3739': ((76.0, 46.9), (0.0, 0.0, -0.0, ''), (3.7, 0.8), (4.8, 1.0)),
                    '2S 1711-339': ((0.0, 0.0), (0.0, 0.0, -0.0, ''), (0.0, 4.9), (0.0, 6.4)),
                    'RX J1718.4-4029': ((47.2, 6.2), (0.0, 0.0, -0.0, ''), (4.7, 0.3), (6.1, 0.4)),
                    'IGR J17191-2821': ((0.0, 0.0), (0.0, 0.0, -0.0, ''), (0.0, 5.9), (0.0, 7.7)),
                    'XTE J1723-376': ((0.0, 0.0), (0.0, 0.0, -0.0, ''), (0.0, 3.7), (0.0, 4.8)),
                    '4U 1722-30': ((61.8, 16.6), (7.4, 0.5, -0.5, '3,4,5'), (4.1, 0.5), (5.4, 0.6)),
                    '4U 1728-34': ((94.0, 35.9), (0.0, 0.0, -0.0, ''), (3.3, 0.5), (4.4, 0.7)),
                    'MXB 1730-335': ((28.0, 7.0), (7.87, 0.56, -0.5, '5'), (6.1, 0.6), (8.0, 0.8)),
                    'KS 1731-260': ((50.5, 20.4), (0.0, 0.0, -0.0, ''), (4.6, 0.7), (5.9, 0.9)),
                    'SLX 1735-269': ((52.9, 21.0), (0.0, 0.0, -0.0, ''), (4.5, 0.7), (5.8, 0.9)),
                    '4U 1735-444': ((34.2, 22.0), (5.65, 3.62, -2.14, 'G'), (5.5, 1.2), (7.2, 1.6)),
                    'XTE J1739-285': ((0.0, 0.0), (4.06, 4.25, -2.44, 'G'), (0.0, 6.1), (0.0, 7.9)),
                    'SLX 1737-282': ((68.1, 12.4), (0.0, 0.0, -0.0, ''), (3.9, 0.3), (5.1, 0.4)),
                    'KS 1741-293': ((0.0, 0.0), (0.0, 0.0, -0.0, ''), (0.0, 4.2), (0.0, 5.4)),
                    'GRS 1741.9-2853': ((35.3, 9.8), (0.0, 0.0, -0.0, ''), (5.5, 0.6), (7.1, 0.8)),
                    '1A 1742-294': ((37.7, 1.3), (0.0, 0.0, -0.0, ''), (5.3, 0.1), (6.9, 0.1)),
                    'SAX J1747.0-2853': ((52.7, 31.4), (0.0, 0.0, -0.0, ''), (4.5, 0.9), (5.8, 1.2)),
                    'IGR J17473-2721': ((113.6, 11.7), (0.0, 0.0, -0.0, ''), (3.0, 0.1), (4.0, 0.2)),
                    'SLX 1744-300': ((13.7, 3.2), (0.0, 0.0, -0.0, ''), (8.7, 0.9), (11.4, 1.1)),
                    'GX 3+1': ((53.3, 15.2), (0.0, 0.0, -0.0, ''), (4.4, 0.5), (5.8, 0.7)),
                    'IGR J17480-2446': ((36.1, 9.0), (6.9, 0.5, -0.5, '3,4,6'), (5.4, 0.6), (7.0, 0.7)),
                    '1A 1744-361': ((0.0, 0.0), (0.0, 0.0, -0.0, ''), (0.0, 7.0), (0.0, 9.1)),
                    'SAX J1748.9-2021': ((38.0, 6.1), (8.4, 1.5, -1.3, '3,4'), (5.3, 0.4), (6.9, 0.5)),
                    'EXO 1745-248': ((63.4, 10.9), (6.9, 0.5, -0.5, '3,4,6'), (4.1, 0.3), (5.3, 0.4)),
                    'IGR J17498-2921': ((51.6, 1.6), (0.0, 0.0, -0.0, ''), (4.5, 0.1), (5.9, 0.1)),
                    '4U 1746-37': ((5.4, 0.8), (0.0, 0.0, -0.0, ''), (13.9, 1.0), (18.2, 1.3)),
                    'SAX J1750.8-2900': ((54.3, 6.1), (0.0, 0.0, -0.0, ''), (4.4, 0.2), (5.7, 0.3)),
                    'GRS 1747-312': ((13.4, 7.3), (6.7, 0.5, -0.5, '3,4,6'), (8.8, 1.7), (11.5, 2.2)),
                    'IGR J17511-3057': ((0.0, 0.0), (0.0, 0.0, -0.0, ''), (0.0, 4.1), (0.0, 5.4)),
                    'SAX J1752.3-3138': ((0.0, 0.0), (0.0, 0.0, -0.0, ''), (0.0, 6.8), (0.0, 8.9)),
                    'SAX J1753.5-2349': ((0.0, 0.0), (0.0, 0.0, -0.0, ''), (0.0, 4.4), (0.0, 5.7)),
                    'IGR J17597-2201': ((15.7, 0.8), (0.0, 0.0, -0.0, ''), (8.2, 0.2), (10.7, 0.3)),
                    'SAX J1806.5-2215': ((0.0, 0.0), (0.0, 0.0, -0.0, ''), (0.0, 4.8), (0.0, 6.2)),
                    '2S 1803-245': ((0.0, 0.0), (0.0, 0.0, -0.0, ''), (0.0, 4.2), (0.0, 5.4)),
                    'SAX J1808.4-3658': ((230.2, 26.3), (0.0, 0.0, -0.0, ''), (2.1, 0.1), (2.8, 0.1)),
                    'XTE J1810-189': ((54.2, 1.8), (0.0, 0.0, -0.0, ''), (4.4, 0.1), (5.7, 0.1)),
                    'SAX J1810.8-2609': ((111.3, 7.2), (0.0, 0.0, -0.0, ''), (3.1, 0.1), (4.0, 0.1)),
                    'XMMU J181227.8-181234': ((2.4, 0.3), (14.0, 2.0, -2.0, '7'), (20.9, 1.2), (27.2, 1.5)),
                    'XTE J1814-338': ((0.0, 0.0), (0.0, 0.0, -0.0, ''), (0.0, 8.6), (0.0, 11.3)),
                    '4U 1812-12': ((203.1, 40.1), (0.0, 0.0, -0.0, ''), (2.3, 0.2), (3.0, 0.3)),
                    'GX 17+2': ((14.6, 5.0), (0.0, 0.0, -0.0, ''), (8.5, 1.2), (11.1, 1.5)),
                    'SAX J1818.7+1424': ((0.0, 0.0), (0.0, 0.0, -0.0, ''), (0.0, 5.9), (0.0, 7.7)),
                    '4U 1820-303': ((60.5, 22.6), (7.6, 0.4, -0.4, '3,4,8'), (4.2, 0.6), (5.4, 0.8)),
                    'SAX J1828.5-1037': ((0.0, 0.0), (0.0, 0.0, -0.0, ''), (0.0, 5.8), (0.0, 7.6)),
                    'GS 1826-24': ((40.0, 3.0), (0.0, 0.0, -0.0, ''), (5.1, 0.2), (6.7, 0.2)),
                    'XB 1832-330': ((33.8, 4.5), (9.6, 0.4, -0.4, '3,4,9'), (5.6, 0.3), (7.3, 0.4)),
                    'Ser X-1': ((29.4, 13.8), (4.31, 2.54, -1.61, 'G'), (6.0, 1.0), (7.8, 1.4)),
                    'HETE J1900.1-2455': ((123.9, 8.6), (0.0, 0.0, -0.0, ''), (2.9, 0.1), (3.8, 0.1)),
                    'Aql X-1': ((103.3, 19.6), (2.97, 2.64, -1.32, 'G'), (3.2, 0.3), (4.2, 0.3)),
                    'XB 1916-053': ((30.6, 3.5), (0.0, 0.0, -0.0, ''), (5.8, 0.3), (7.6, 0.4)),
                    'XTE J2123-058': ((0.0, 0.0), (0.0, 0.0, -0.0, ''), (0.0, 12.3), (0.0, 16.0)),
                    'M15 X-2': ((40.8, 7.2), (10.38, 0.15, -0.15, '1'), (5.1, 0.4), (6.6, 0.5)),
                    'Cyg X-2': ((13.1, 2.3), (6.95, 1.16, -0.91, 'G'), (8.9, 0.7), (11.6, 0.9)),
                    'SAX J2224.9+5421': ((0.0, 0.0), (0.0, 0.0, -0.0, ''), (0.0, 6.0), (0.0, 7.9)) }

        # also define a string which determines how the distance was derived
        distances = { x: (table8[x][1] if (table8[x][1][0] > 0.0 and Gaia) 
            else table8[x][idist]) for x in table8 }

        # logger.warning('distances are outdated, use with caution')
        addl = ''
        if Gaia:
            addl = ', using Gaia/cluster distances where available'

        logger.info('adopted distances for X={}{}'.format( X, addl))

        dist = np.zeros_like(self['ra_obj'])
        diste_hi = np.zeros_like(dist)
        diste_lo = np.zeros_like(dist)
        dist_method = np.empty(len(dist), dtype="<U10")
        for i, name in enumerate(self['name']):
            if name in distances:
                dist[i] = distances[name][0]
                if len(distances[name]) == 2:
                    diste_hi[i] = distances[name][1]
                    diste_lo[i] = distances[name][1]
                    dist_method[i] = 'B'
                else:
                    # two-sided distances are converted to single-sided, for now
                    # diste[i] = 0.5*(distances[name][1]-distances[name][2])
                    diste_hi[i] = distances[name][1]
                    diste_lo[i] = -distances[name][2]
                    dist_method[i] = distances[name][3]

        return dist, diste_hi, diste_lo, dist_method

    def __str__(self):
        """
        Display some information about this object
        The "local files" option refers to auxiliary information including the
        list of Eddington fluxes, which is not part of the online dataset
        :return:
        """

        sel_str = ""
        if self.selection is not None:
            # Generate a string to give information about the current selection
            n_sel = len(np.where(self.selection)[0])
            if n_sel == 1:
                sel_str = "{} selected ({}".format(n_sel, self['name'])
            else:
                sel_str = "{} selected ({}".format(n_sel, ", ".join(self['name'][:min([3,n_sel])]) )
            if n_sel > 4:
                sel_str += " ... "
            elif n_sel == 4:
                sel_str += ", "
            if n_sel > 3:
                sel_str += self['name'][-1]
            sel_str += ")\n"

        # the "Local files" option here refers only to the F_Edd file, which is not part of the
        # public distribution
        available = ['unavailable','present']

        return """Multi-INstrument Burst ARchive (MINBAR) source table v{}
{} sources ({} with bursts in MINBAR)
{}
Local files: {}""".format( self.version, len(self),
                           # this version will include the selection
                           # len(np.where(self['nburst'] > 0)[0]),
                           len(np.where(self._f[1].data['nburst'] > 0)[0]),
                           sel_str,
                           available[int(self.local_files)])


class Instrument:
    """
    Here's a generic instrument class which can be adapted and/or duplicated for different
    instruments. This class is kept pretty lean to avoid having to replicate lots of code
    Defines the properties of an instrument with data that we're going to analyse and add to MINBAR
    :param name: name of the instrument
    :param source_name: list of source names, to avoid having to load a Sources object
    Example
    import minbar
    jemx = minbar.Instrument('JEM-X', 'jemx', 'IJ')
    """

    def __init__(self, name, camera=None, path=None, label=None,
                 lightcurve=['lc1_3-25keV_1.0s.fits','lc2_3-25keV_1.0s.fits'],
                 source_name=None,
                 spectrum=None,
                 verbose=False,
                 effarea=None):

        self.name = name
        if name in MINBAR_INSTR_LABEL.keys():
            # These are the known MINBAR instruments
            self.label = MINBAR_INSTR_LABEL[name]
            self.instr = camera # copy of the instr row
            if path is None:
                path = MINBAR_INSTR_PATH[self.label]

        else:
            # all other instruments
            if (label is None) or (path is None):
                logger.error('for "new" instruments need to specify path and label')
                return None
            logger.warning('new instrument {} may not be fully implemented'.format(name))
            self.label = label
            if effarea is None:
                self.effarea = 100*u.cm**2
                logger.warning('effective area not supplied, assuming {}'.format(self.effarea))

        if os.path.isdir('/'.join([MINBAR_ROOT, path])):
            self.path = path
        else:
            # sys.exit(1)
            # logger.error('need a valid path for the data files')
            # avoid printing the warning, will be annoying for multiple instruments
            # logger.warning('MINBAR_ROOT path for raw data does not exist')
            # return None
            self.path = None

        # Set the list of names. If not provided we'll create a Sources object and read it
        # (as a chararray) from there; otherwise you can pass it as a parameter

        if source_name is None:
            self.source = Sources()
            self.source_name = self.source['name']
        else:
            self.source_name = source_name

        # Define the paths corresponding to each source here
        # Default is just the source name with spaces removed; this won't work for RXTE
        # (without some judicious softlinks)

        # print (type(self.source_name), type(self.source_name[0]), self.source_name[0])
        self.source_path = np.char.replace(self.source_name, " ", "")
        if self.label == 'XP':
            # for RXTE/PCA, most sources just omit the prefix (but not the GX sources)
            # self.source_path = np.array([x.split()[1] for x in self.source_name])
            # The standard product lightcurve also includes the obsid, so define the lightcurve
            # parameter as a function here
            lightcurve = self.pca_lightcurve_filename
            self.effarea = PCA_EFFAREA
            self.effarea_bursts = PCA_EFFAREA

        if self.label == 'IJ':
            # for JEM-X, the convention is to also replace '+' with 'p'
            # self.source_path = self.source_path.replace("+","p")
            self.source_path = np.char.replace(self.source_path, "+","p")
            self.effarea = JEMX_EFFAREA
            self.effarea_bursts = JEMX_EFFAREA_BURSTS

        if self.label == 'SW':

            # BeppoSAX/WFC
	    # list of WFC source file names here. There are many more names
	    # than the bursters
            self.wfc_lightcurve_source_name = ['1RXHJ173523.7-354013',
             '1RXSJ180408.9-342058', 'AXJ1745.6-2901', 'AXJ1754.2-2754',
             'HETEJ1900.12455', 'IGRJ17062-6143', 'IGRJ17191-2821',
             'IGRJ17254-3257', 'IGRJ17380-3749', 'IGRJ17464-2811',
             'IGRJ17498-2921', 'IGRJ17511-3057', 'ks1724-35', 'M28',
             'MAXIJ1647-227', 'RXJ1718.4-4029', 'RXSJ170854-3219',
             'SAXJ1324.5-6313', 'SAXJ1603.9-7753', 'SAXJ1712.6-3739',
             'SAXJ1747.0-2853', 'SAXJ1748.9-2021', 'SAXJ1750.8-2900',
             'SAXJ1752.4-3138', 'SAXJ1753.5-2349', 'SAXJ1806.5-2215',
             'SAXJ1808.4-3658', 'SAXJ1810.8-2609', 'SAXJ1818.6-1703',
             'SAXJ1818.7+1424', 'SAXJ1828.5-1037', 'SAXJ2103.5+4545',
             'SAX J2224.9+5421', 'SAXJ2224.9+5421', 'SWIFTJ1749.4-2807',
             'SWIFTJ185003.2-005627', 'SWIFTJ1922.7-1716', 'Terzan5',
             'UWCrb', 'velapulsa', 'velax1', 'X0025+641', 'X0042+327',
             'X0050-727', 'X0052-739', 'X0053+604', 'X0103-720', 'X0103-762',
             'X0104-720', 'X0114+650', 'X0115+634', 'X0115-737', 'X0116-737',
             'X0142+614', 'X0143+611', 'X0240+611', 'X0352+309', 'X0422+328',
             'X0449-055', 'X0512-401', 'X0521-720', 'X0522-696', 'X0527-328',
             'X0531+219', 'X0531-661', 'X0532-664', 'X0535-668', 'X0535-692',
             'X0538-641', 'X0540-693', 'X0540-697', 'X0543-682', 'X0544-665',
             'X0547-711', 'X0614+091', 'X0656-072', 'X0726-260', 'X0746-532',
             'X0748-676', 'X0755-609', 'X0812-311', 'X0818-52', 'X0821-42',
             'X0834-430', 'X0836-429', 'X0918-549', 'X0921-630', 'X1008-570',
             'X1010-584', 'X1011-447', 'X1024-575', 'X1028-567', 'X1037-564',
             'X1037-592', 'X1048-596', 'X1059-771', 'X1118-616', 'X1119-603',
             'X1124-684', 'X1137-651', 'X1145-616', 'X1145-619', 'X1148-665',
             'X1223-624', 'X1235-751', 'X1239-599', 'X1244-603', 'X1246-588',
             'X1248-411', 'X1249-289', 'X1249-637', 'X1251-568', 'X1254-690',
             'X1258-613', 'X1259-600', 'X1322-427', 'X1323-618', 'X1344-326',
             'X1344-603', 'X1354-644', 'X1414+250', 'X1417-624', 'X1455-314',
             'X1516-569', 'X1524-617', 'X1538-522', 'x1543-475', 'X1543-624',
             'X1544-475', 'X1550-552', 'X1550-564', 'X1553-542', 'X1556-605',
             'X1608-522', 'X1617-155', 'X1624-375', 'X1624-490', 'X1627-673',
             'X1630-472', 'X1632-497', 'X1633-170', 'X1636-536', 'X1642-455',
             'X1644+500', 'X1652+390', 'X1655-400', 'X1656+354', 'X1657-415',
             'X1658-298', 'X1659-487', 'X1700-377', 'X1702-363', 'X1702-429',
             'X1704+240', 'X1705-440', 'X1706-266', 'X1708-233', 'X1708-407',
             'X1711-339', 'X1715-321', 'X1716-249', 'X1722-363', 'X1724-307',
             'X1727-214', 'X1728-169', 'X1728-247', 'X1728-337', 'X1728-338',
             'X1730-220', 'X1730-311', 'X1730-333', 'X1731-260', 'X1732-300',
             'X1732-304', 'X1734-275', 'X1734-292', 'X1735-269', 'X1735-28',
             'X1735-444', 'X1736-297', 'X1736-310', 'X1737-132', 'X1737-282',
             'X1739-277', 'X1741-286', 'X1741-289', 'X1741-293', 'X1741-297',
             'X1741-322', 'X1742-289', 'X1742-290', 'X1742-294', 'X1742-295',
             'X1742-326', 'X1743-287', 'X1743-288', 'X1743-290', 'X1743-322',
             'X1744-265', 'X1744-28', 'X1744-299', 'X1744-300', 'X1744-361',
             'X1745-203', 'X1745-24', 'X1745-248', 'X1746-324', 'X1746-370',
             'X1747-214', 'X1747-299', 'X1747-312', 'X1747-341', 'X1749-285',
             'X1750-248', 'X1750-27', 'X1754-32', 'X1755-338', 'X1758-205',
             'X1758-25', 'X1758-250', 'X1758-258', 'X1803-245', 'X1805-204',
             'X1806-202', 'X1807-10', 'X1811-171', 'X1812-12', 'X1813-140',
             'X1814+498', 'X1820-303', 'X1822-000', 'X1822-371', 'X1826-235',
             'X1832-330', 'X1833-103', 'X1837+049', 'X1839-06', 'X1843+00',
             'X1845-024', 'X1845-03', 'X1846-031', 'X1850-087', 'X1851-034',
             'X1851-312', 'X1854+052', 'X1905+000', 'X1907+097', 'X1908+005',
             'X1908+075', 'X1909+048', 'X1915+105', 'X1916-053', 'X1918-32',
             'X1942+274', 'X1947+300', 'X1949+32', 'x1953+319', 'X1956+350',
             'X1957+115', 'X2017-01', 'X2023+338', 'X2030+3', 'X2030+375',
             'X2030+407', 'X2057+3', 'X2058+42', 'X2100+41', 'X2127+119',
             'X2127+1191', 'X2129+470', 'X2138+567', 'X2140+433',
             'X2142+380', 'X2159+499', 'X2206+542', 'X2215-086', 'X2250+165',
             'X2252-034', 'X2259+587', 'X2323+585', 'XMMJ174457-2850.3',
             'XTEJ1701-407', 'XTEJ1701-462', 'XTEJ1709-267', 'XTEJ1710-281',
             'XTEJ1723-056', 'XTEJ1723-376', 'XTEJ1739-285', 'XTEJ1747-274',
             'XTEJ1759-221', 'XTEJ1810-189', 'XTEJ1814-338', 'XTEJ2123-058']

            self.nmatched = 0
            for i, name in enumerate(self.source_path):

                if not (name in self.wfc_lightcurve_source_name):
                    # if we can't find a match, we might need to modify a bit
                    alt = name.replace('4U','X')

                    if (alt in self.wfc_lightcurve_source_name):
                        self.source_path[i] = alt
                        self.nmatched += 1
                    else:
                        # 2nd last try, in case file version has an extra digit
                        imatch = [j for j, x in enumerate(self.wfc_lightcurve_source_name) if alt == x[:len(alt)]]
                        if imatch == []:
                            # absolutely last try, in case source version has
                            # an extra digit
                            imatch = [j for j, x in enumerate(self.wfc_lightcurve_source_name) if alt.ljust(len(x),'0') == x]
                        if imatch != []:
                            self.source_path[i] = self.wfc_lightcurve_source_name[imatch[0]]
                            self.nmatched += 1
                else:
                    self.nmatched += 1
            # above code matches 51/115 sources
            if self.nmatched < len(self.source_path):
                logger.warning('only got matched source path names for {}/{} sources'.format(self.nmatched, len(self.source_path)))
            lightcurve = self.wfc_lightcurve_filename
            # Lightcurves are already normalised
            self.effarea = 1.*u.cm**2
            self.effarea_bursts = 1.*u.cm**2 # I think this is correct! - dkg

        # Check that the directories in the data path all correspond to source_paths
        # You don't want to miss observations in some directory because of inconsistencies
        # with the directory names

        self.has_dir, self.nonmatched, self.nmissing = verify_path(
            self.source_name, self.path, self.source_path, verbose=verbose)
        if self.label == 'SW':
            # special here for WFC as we don't have source subdirectories
            self.has_dir = os.path.isdir('/'.join([MINBAR_ROOT, path]))

        # as for MINBAR

        self.local_data = np.any(self.has_dir)

        # Define the lightcurve and spectral files

        assert lightcurve is not None
        self.lightcurve = lightcurve
        self.spectrum = spectrum


    def __str__(self):
        """
        Method to display information about this object
        TODO: add PCU decode function to show active PCUs
        :return:
        """

        return """
MINBAR instrument definition

Name: {} ({})
Camera/detector flag: {}
Data path: {}/{}/data
Lightcurve(s): {}
Spectra: {}""".format(self.name, self.label, self.instr[2:], MINBAR_ROOT, self.path, self.lightcurve, self.spectrum)


    def filename_with_obsid(self, template, obsid, exclude=None):
        """
        Function to return a filename (possibly including path)
        incorporating the obsservation ID, as used for RXTE and NuSTAR
        Also have the option of excluding a character from the obsid string
        See pca_lightcurve_filename for example usage
        :return:
        """

        if not ('{}' in template):
            logger.warning("template doesn't seem to include obsid placeholder")

        if exclude:
            obsid = obsid.replace(exclude, "")

        return template.format(obsid)


    def pca_lightcurve_filename(self, obsid):
        """
        Function to return the lightcurve name for PCA observations
        :return:
        """

        # return 'stdprod/xp{}_n1.lc.gz'.format(obsid.replace("-", ""))
        return self.filename_with_obsid('stdprod/xp{}_n1.lc.gz', obsid, "-")


    def wfc_lightcurve_filename(self, name, obsid, camera):
        """
        Function to return the lightcurve name for WFC observations
        This convention is for the 2013 data
        :return:
        """

        return '{}_{}w{}_e1_31.lcv'.format(name, obsid, camera)


    def analyse_persistent(self, src, obsid):
        """
        Function to analyse the persistent emission (lightcurve and spectrum) for a single
        observation
        :param src:
        :param obsid:
        :return:
        """

    def analyse_burst(self, bursts):
        """
        Function to analyse the persistent emission (lightcurve and spectrum) for a single
        observation
        :param src:
        :param obsid:
        :return:
        """
