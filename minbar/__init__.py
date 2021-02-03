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
MINBAR_URL = 'https://burst.sci.monash.edu/wiki/uploads/MINBAR/'

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

    return (_time.tt.mjd - 49353.000696574074) * 86400. - 3.378431

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
        :return:
        """

        if attributes is None:
            attributes = self.attributes_default
        else:
            # check all the attributes are in the field names
            in_attr_list = True
            for attr in attributes:
                in_attr_list = in_attr_list & (attr in self.field_names)
                print (attr, in_attr_list)
            if not in_attr_list:
                logger.error("attribute not present in table")
                return
        print (self)
        if all:
            self.records[self.selection][attributes].pprint_all()
        else:
            self.records[self.selection][attributes].pprint()

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
            return np.ones(len(self.records), np.bool)
        else:
            return self.records.field('type') == self.type

    def clear(self):
        """
        Clear the selection. If self.type is not None, only bursts of the given
        type are selected.
        """
        self.name = ''
        self.selection = self.get_type()
        self.time_order = np.argsort(self.records[self.selection].field(self.timefield))
        self.ind = np.where(self.selection)[0][self.time_order]

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

        self.time_order = np.argsort(self.records[self.selection].field(self.timefield))
        self.ind = np.where(self.selection)[0][self.time_order]
        if attribute == 'name':
            self.name = value
            logger.info('{} {} {}s from {}'.format(action, n_action, self.entryname, self.name))
        else:
            logger.info('{} {} {}s with {}={}'.format(action, n_action, self.entryname, attribute, value))

        # Return self so we can "cascade" selections

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
        return self.records[self.selection][self.time_order]


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
            elif type(field[0]) == bool or type(field[0]) == np.bool_:
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
        self.time_order = np.argsort(self.records[self.selection].field(self.timefield))
        if len(np.where(self.selection)[0]) < len(self.ind):
            logger.info('excluded {} {}s by excluding flag(s) {}'.format(len(self.ind)-len(np.where(self.selection)[0]),
                    self.entryname, flags))
        self.ind = np.where(self.selection)[0][self.time_order]


class Bursts(Minbar):
    """
    Read the MINBAR IDL database and give access to its
    contents.
    
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
        self.clear()

    def get_burst_data(self, id):
        """
        Retrieve time-resolved spectroscopy for this burst
        Replacement for the get_burst_data IDL function
        :param id:
        :return:
        Example usage:
        data = b.get_burst_data(1918)
        """

        wfcdataroot = 'wfcspec_vs2'

        if self.local_data:
            # Try to get the file locally
            instr = self[id]['instr']
            if instr[0:2] == 'IJ':
                logger.error('Time-resolved spectroscopy not available for JEM-X bursts')
                return None

            # Get the source path
            # Can probably package this into the Minbar class for wider use

            _match = self.instruments[instr[0:2]].source_name == self[id]['name']
            if not np.any(_match):
                logger.error('No source directory for this burst')
                return None

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

                if not os.path.isfile(_file):
                    logger.error('time-resolved spectroscopy file not found!')
                else:
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
                        names=['trange','kT','kT_err','rad','rad_err','flux_3_25','flux_3_25_err',
                               'flux','fluxerr','chisq'], sep='\s+')
                    nspec=len(_data)
                    time=np.zeros(nspec)
                    dt=np.zeros(nspec)
                    for j in range(nspec):
                        tmp=_data['trange'][j].split('-')
                        time[j]=float(tmp[0])
                        dt[j]=(float(tmp[1])-time[j])*86400.
                    _data['dt'] = dt
                    # This calculation assumes UTC for the time ranges
                    _data['time'] = (time-self[id]['time'])*86400.

                    # These flux errors are FAR too big
                    _data['flux'] *= 1e9
                    _data['fluxerr'] *= 1e9
                    _data['flux_min'] = _data['flux']-_data['fluxerr']
                    _data['flux_max'] = _data['flux']+_data['fluxerr']
                    logger.warning('BeppoSAX/WFC flux is 2-10 keV only!')
                    lz=np.where(_data['flux_min'] < 0.)[0]
                    if len(lz) > 0:
                        _data['flux_min'][lz]=0.1
                    _data['kT_min'] = _data['kT']-_data['kT_err']
                    _data['kT_max'] = _data['kT']+_data['kT_err']
                    _data['rad_min'] = _data['rad']-_data['rad_err']
                    _data['rad_max'] = _data['rad']+_data['rad_err']

                    # Need to be aware of different conventions here for the files; SAX files
                    # have radius in units of km/10kpc, while RXTE is (km/10kpc)^2

                    _data['rad'] = _data['rad']**2
                    _data['rad_min'] = _data['rad_min']**2
                    _data['rad_max'] = _data['rad_max']**2

            elif instr[0:2] == 'XP':
                # Get PCA data
                # this is a bit of a challenge with the historical arrangement of the burst data under the
                # alternate file heirarchy, but we can overcome this with some judicious softlinks
                _path = '/'.join([MINBAR_ROOT, MINBAR_INSTR_PATH[instr[0:2]], 'data',
                                  self.instruments[instr[0:2]].source_path[_match][0],
                                  self[id]['obsid'], 'burst{}'.format(self[id]['bnum']) ])
                _file = '/'.join([_path, 'analysis/bbfit_kabs.log'])

                _data = pd.read_csv(_file, comment='#',
                    names=['time', 'r', 're', 'dt', 'nH', 'nH_min', 'nH_max', 'kT', 'kT_min', 'kT_max',
                        'rad', 'rad_min', 'rad_max', 'chisq', 'rawflux', 'flux', 'flux_min', 'flux_max'], sep='\s+')

                # _data = _path # for testing
                _data['fluxerr'] = 0.5*(_data['flux_max']-_data['flux_min'])

                # Need to adjust time from SS to seconds post star time

                _data['time'] -= mjd_to_ss(self[id]['time'])

            else:
                logger.error('time-resolved spectroscopy not yet implemented for this instrument')

        else:
            # Try to get the data remotely
            logger.error('Remote data retrieval not yet implemented')

        return _data

    def get_lc(self, id, pre=16., post=None):
        """
        Preliminary routine to return the lightcurve corresponding to a burst
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

            url = MINBAR_URL+'bursts/{0:04d}_lc.csv'.format(id)

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

        self.time_order = np.argsort(self.records[self.selection].field(self.timefield))
        self.ind = np.where(self.selection)[0][self.time_order]

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

        self.time_order = np.argsort(self.records[self.selection].field(self.timefield))
        n_excluded = len(self.ind) - len(np.where(self.selection)[0])
        if n_excluded > 0:
            logger.info('selected {} unique {}s by excluding {}'.format(len(np.where(self.selection)[0]),
                    self.entryname, n_excluded))
        self.ind = np.where(self.selection)[0][self.time_order]

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
        This can be used to easily convert flux to luminosity.
        
        Returns correction factor (cm^2), error, distance (kpc), error. Furthermore,
        these arrays are stored and can be accessed as self['fieldname'], with fieldnames
        distcor, distcore, dist, diste.
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
        self.clear()

    def good(self):
        """
        Filter for only the "good" observations, excluding bad flags and non detections
        :return:
        """
        self.exclude_flag('bcdefg')
        selection = (self.records['flux'] > 3.*self.records['e_flux']) & (self.records['sig'] >= 3.)
        self.selection = np.logical_and(self.selection, selection)
        self.time_order = np.argsort(self.records[self.selection].field(self.timefield))
        if len(np.where(self.selection)[0]) < len(self.ind):
            logger.info('Restricted to {} "good" {}s by also excluding nondetections'.format(
                len(np.where(self.selection)[0]), self.entryname))
        self.ind = np.where(self.selection)[0][self.time_order]

    def __str__(self):
        """
        Return a nice string.
        """
        return "Multi-INstrument oBservation ARchive (MINBAR) ({} {}s from {} sources)".format(
            len(self.records), self.entryname, len(self.names))


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
            instr = Instrument(label)
            name = obs_entry['name']
            obsid = obs_entry['obsid']

            # self.tstart = obs_entry['tstart']
            # self.tstop = obs_entry['tstop']
            # Copy all the columns to the new object
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
            output += "\nTime range: {}-{}".format(self.tstart, self.tstop)
        _path = self.get_path()
        if _path is not None:
            output += "\nData path: {}".format(_path)

        return output


    def plot(self):
        """
        Plot the lightcurve, reading it in first if need be
        :return:
        """
        if self.time is None:
            lc = self.get_lc()

        # plt.plot(lc['TIME'],lc['RATE'])
        # plot can't work with "raw" time units
        plt.plot(self.mjd.mjd, self.rate)
        plt.xlabel('Time (MJD '+self.mjd.scale.upper()+')')
        plt.ylabel('Rate (count s$^{-1}$ cm$^{-2}$)')
        plt.show()


    def get_path(self):
        """
        Return the path for MINBAR observations, assuming you have them stored locally
        :param entry:
        :return:
        """

        instr = self.instr.label
        if self.instr.local_data:
            _match = np.where(self.instr.source_name == self.name)[0]
            # print (_match)
            assert len(_match) == 1
            # if len(_match) > 1:
                # logger.warning("multiple source name matches for path")
            if self.instr.has_dir[_match[0]]:
                return '/'.join([MINBAR_ROOT, self.instr.path, 'data',
                         self.instr.source_path[_match[0]],
                         self.obsid])

        return None


    def get_lc(self):
        """
        Return the lightcurve for a particular observation; this is a replacement for the IDL routine get_lc.pro
        This routine also populates the time, mjd_tt, mjd, rate, and error attributes for the observation
        :param entry:
        :return:
        """

        path = self.get_path()
        if self.instr.label == 'XP':
            filename = self.instr.lightcurve(self.obsid)
        else:
            # logger.warning("other instruments not yet implemented")
            filename = self.instr.lightcurve

        # print (path+'/'+filename)
        lcfile = fits.open(path+'/'+filename)

        # For XTE files, convention is to have the first extension RATE and a second extension STDGTI
        header = lcfile[0].header
        lc = lcfile[1].data

        lcfile.close()

        # Can clean up the table here if necessary; i.e. create a Lightcurve
        # object (not yet defined), adopt uniform time scale etc.
        # Using astropy to keep track of the time scale and units, read from the
        # header file; see https://docs.astropy.org/en/stable/time

        self.timesys = header['TIMESYS']
        self.timeunit = header['TIMEUNIT']
        self.time = lc['TIME']*u.Unit(self.timeunit)

        effarea = 1.*u.cm**2    # dummy value
        if self.instr.name == 'PCA':
            # convert raw times to MJD (TT) here; see https://heasarc.gsfc.nasa.gov/docs/xte/abc/time_tutorial.html
            # can check results using https://heasarc.gsfc.nasa.gov/cgi-bin/Tools/xTime
            self.mjd_tt = Time( self.time.to('d') + (header['MJDREFI'] + header['MJDREFF'])*u.d, format='mjd', scale='tt') # TT
            self.mjd = self.mjd_tt.utc

        self.rate = lc['RATE']/u.s/self.instr.effarea
        self.error = lc['ERROR']/u.s/self.instr.effarea

        return lc


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

    def __init__(self, name, path=None, label=None,
                 lightcurve=['lc1_3-25keV_1.0s.fits','lc2_3-25keV_1.0s.fits'],
                 source_name=None,
                 spectrum=None,
                 verbose=False):

        self.name = name
        if name in MINBAR_INSTR_LABEL.keys():
            # These are the known MINBAR instruments
            self.label = MINBAR_INSTR_LABEL[name]
            if path is None:
                path = MINBAR_INSTR_PATH[self.label]

        else:
            # all other instruments
            if (label is None) or (path is None):
                logger.error('for "new" instruments need to specify path and label')
                return None
            logger.warning('new instrument {} may not be fully implemented'.format(name))
            self.label = label

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
            # Lightcurves are already normalised
            self.effarea = 1.*u.cm**2

        # Check that the directories in the data path all correspond to source_paths
        # You don't want to miss observations in some directory because of inconsistencies
        # with the directory names

        self.has_dir, self.nonmatched, self.nmissing = verify_path(
            self.source_name, self.path, self.source_path, verbose=verbose)

        # as for MINBAR

        self.local_data = np.any(self.has_dir)

        # Define the lightcurve and spectral files

        assert lightcurve is not None
        self.lightcurve = lightcurve
        self.spectrum = spectrum


    def __str__(self):
        """
        Method to display information about this object
        :return:
        """

        return """
MINBAR instrument definition

Name: {} ({})
Data path: {}
Lightcurve(s): {}
Spectra: {}""".format(self.name, self.label, self.path, self.lightcurve, self.spectrum)


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
