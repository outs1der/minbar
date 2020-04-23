"""
Python class to easily extract information from MINBAR.
Provides the burst catalog (Bursts), observation catalog (Observations),
and source catalog (Sources). See the doc strings of those classes
for example usage.

The data files should be placed in the 'data' subdirectory of this
package.

(c) 2007, Laurens Keek, l.keek@sron.nl
Updated for MINBAR DR1, 2020, Duncan Galloway, duncan.galloway@monash.edu
Updated for MINBAR v0.9, 2017, Laurens Keek, laurens.keek@nasa.gov
"""

from .idldatabase import IDLDatabase
from .analyse import *
import numpy as np
import os, re
from astropy.io import fits
from astropy.io import ascii
import astropy.units as u
from astropy.time import Time
import logging
import sys

import matplotlib.pyplot as plt

# kpc = 3.086e21 # cm
kpc = u.kpc.to('cm') # cm

# Local paths for MINBAR data
# Bytearr entries are for consistency between the IDL and ASCII
# versions of the data

MINBAR_ROOT = '/home/burst/minbar'
MINBAR_INSTR_LABEL = {'PCA': 'XP', 'WFC': 'SW', 'JEM-X': 'IJ'}
MINBAR_INSTR_PATH = {'XP': 'xte', 'SW': 'wfc', 'IJ': 'jemx',
                     'PCA': 'xte', 'WFC': 'wfc', 'JEM-X': 'jemx',
                     b'XP': 'xte', b'SW': 'wfc', b'IJ': 'jemx'}

# List of ultra compacts from In 't Zand (2007) (1850-087 and 1905+000 no bursts in MINBAR)
# Includes all candidates. Should I add 1728?
# Should generate this list dynamically based on the type code --- dkg
UCXBS = ['4U 0513-40',
         '4U 0614+09',
         '2S 0918-549',
         '4U 1246-588',
         '4U 1705-32', # = 1RXS J170854.4-321957
         'SAX J1712.6-3739',
         'RX J1718.4-4029',
         '4U 1722-30',
         'SLX 1735-269',
         'SLX 1737-282',
         'SLX 1744-300',
         '4U 1812-12',
         '4U 1820-303',
         'XB 1832-330',
         '4U 1850-086',
         '1905+000', # no bursts in MINBAR
         'XB 1916-053',
         '4U 2129+12',
         ]

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

    def fix_labels(self):
        """
        Fix the whitespace in the labels when reading in the IDL version.
        """
        pat = re.compile('\s+')
        for k in self.field_labels:
            self.field_labels[k] = re.sub(pat, ' ', self.field_labels[k])

    def get_names(self):
        """
        Get a list of all source names in the archive. Ordered by right ascension.
        RA is not a part of the MRT tables, so can only do this with the IDL version.
        Should really replace with a read of the FITS table, which is the definitive
        version
        """
        names = np.unique(self.records.field('name'))

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


    def select(self, name):
        """
        Set the selection to the source with given name.
        """
        self.name = name
        self.selection = (self.records.field('name') == name) & self.get_type()  # Match name and burst type
        self.time_order = np.argsort(self.records[self.selection].field(self.timefield))
        self.ind = np.where(self.selection)[0][self.time_order]
        logger.info('Selected {} {}s from {}'.format(len(np.where(self.selection)[0]), self.entryname, self.name))


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


    def get_name_like(self, name):
        """
        Return a list of sources that have 'name' in their archive
        identifier.
        """
        selection = []
        for i in self.names:
            if i.find(name) > -1:
                selection.append(i)
        return selection


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
        if field in self.field_names:
            return self.get_records().field(field)
        else:
            return getattr(self, field)[self.ind]


    def get_records(self):
        """
        Get the time ordered records that are currently
        selected.
        """
        return self.records[self.selection][self.time_order]


    def __getitem__(self, field):
        """
        Shorthand for get().
        Expand the usage here to also return the item with the entry number
        """

        if type(field) == str:
            return self.get(field)

        return self.records[self.records['entry'] == field]


    def instr_like(self, instrument):
        """
        Return an array that selects all entries where the instrument name
        begins with the given instrument string. E.g., 'XP' for PCA. For convenience,
        the following aliases are provided: 'pca', 'wfc', 'jemx'.
        """
        alias = {'pca': 'XP', 'wfc': 'SW', 'jemx': 'IJ'}
        if instrument in alias:
            instrument = alias[instrument]

        # This is a crappy way to try to make things work with both strings and byte arrays
        try:
            return np.char.array(self['instr']).startswith(instrument)
        except:
            return np.char.array(self['instr']).startswith(instrument.encode('ascii'))


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

    
    def select_all(self, names):
        """
        Select multiple sources with given names.
        """
        selection = np.zeros_like(self.selection)
        for name in names:
            selection = np.logical_or(selection, self.records.field('name') == name)
        self.selection = selection&self.get_type() # Only bursts of specified type
        self.time_order = np.argsort(self.records[self.selection].field(self.timefield))
        self.ind = np.where(self.selection)[0][self.time_order]
        logger.info('Selected {} {}s'.format(len(np.where(self.selection)[0]), self.entryname))
    
    def exclude(self, name):
        """
        Removes source with given name from the current selection.
        """
        selection = np.logical_not(self.records.field('name') == name)
        self.selection = np.logical_and(self.selection, selection)
        self.time_order = np.argsort(self.records[self.selection].field(self.timefield))
        self.ind = np.where(self.selection)[0][self.time_order]
        logger.info('Selected {} {}s by excluding {}'.format(len(np.where(self.selection)[0]), self.entryname, name))
    
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
            
            if len(names)>1:
                logger.info('{} more matching sources: {}'.format(len(names) - 1, ', '.join(names[1:])))
    
    def __str__(self):
        """
        Return a nice string.
        """
        return "Multi-INstrument Burst ARchive (MINBAR) ({} {}s from {} sources)".format(len(self.records), self.entryname, len(self.names))
    
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
            s.name_like(name.strip())
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


    def __init__(self, filename=None, type=None, IDL=False):
        """
        Load the database of observations.
        """
        if filename==None:
            filename = self.get_default_path('minbar-obs')

        Minbar.__init__(self, filename, IDL=IDL)

        self.type = type
        self.clear()


    def __str__(self):
        """
        Return a nice string.
        """
        return "Multi-INstrument oBservation ARchive (MINBAR) ({} {}s from {} sources)".format(len(self.records), self.entryname, len(self.names))


class Observation:
    """
    This object is intended to allow all the possible actions you might have on an
    observation. You can create it from a minbar entry, or given an instrument, source name and obs ID
    """


    def __init__(self, obs_entry=None, instr=None, name=None, obsid=None):
        """
        Create an observation instance, either from a MINBAR obs entry, or by-hand
        Ideally this object should make available every parameter in the MINBAR observation table
        :param obs_entry: 
        :param instr: 
        :param source: 
        :param obsid: 
        """

        if obs_entry is not None:
            logger.warning('initialisation from obs entry not yet completely implemented')
            # Really need to create an instrument here
            label = [key for key, value in MINBAR_INSTR_LABEL.items() if value == obs_entry['instr'][0][:2]][0]
            # print (label)
            instr = Instrument(label)
            name = obs_entry['name'][0]
            obsid = obs_entry['obsid'][0]

        # Potentially need to check here that all of the passed parameters are set
        self.instr = instr
        self.name = name
        self.obsid = obsid

        # Define parameters for the lightcurve; later this might be a class
        # The lightcurve might also not be available, so don't force it to be read in now

        self.time = None
        self.rate = None
        self.error = None


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
        plt.ylabel('Rate (count s$^{-1}$)')
        plt.show()


    def get_path(self):
        """
        Return the path for MINBAR observations, assuming you have them stored locally
        :param entry:
        :return:
        """

        instr = self.instr.label
        return '/'.join([MINBAR_ROOT, self.instr.path, 'data',
                         self.instr.source_path[self.instr.source_name == self.name][0],
                         self.obsid])


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
        if self.instr.name == 'PCA':
            # convert raw times to MJD (TT) here; see https://heasarc.gsfc.nasa.gov/docs/xte/abc/time_tutorial.html
            # can check results using https://heasarc.gsfc.nasa.gov/cgi-bin/Tools/xTime
            self.mjd_tt = Time( self.time.to('d') + (header['MJDREFI'] + header['MJDREFF'])*u.d, format='mjd', scale='tt') # TT
            self.mjd = self.mjd_tt.utc

        self.rate = lc['RATE']/u.s
        self.error = lc['ERROR']/u.s

        return lc


class Sources:
    """
    Contain all information on the sources present in Minbar in
    the file minbar_sources.fits.
    
    Example:
    s = Sources()
    print s.field_names # Show available data fields
    ra = s['ra_obj'] # Right ascension for all sources
    s.name_like('1636')
    ra = s['ra_obj'] # Right ascension for selected source only
    s.clear() # Clear selection
    """
    
    def __init__(self, filename=None):
        """
        Load source list from filename.
        """
        if filename==None:
            filename = self.get_default_path()
        self._f = fits.open(filename)
        self.field_names = [i.lower() for i in self._f[1].data.dtype.names]
        self.field_labels = self._get_field_labels()
        self._fits_names = list(self.field_names) # Keep track of which fields are in fits file
        self.clear()
        self.dist, self.diste = self.get_distances()
        self.field_names += ['dist', 'diste']
        self.field_labels['dist'] = 'Distance (kpc)'
        self.field_labels['diste'] = 'Error on distance (kpc)'
    
    def get_default_path(self):
        """
        Return the default path of the source list
        """
        return os.path.join(os.path.dirname(__file__), 'data', 'minbar_sources.fits')
    
    def _get_field_labels(self):
        """
        Get the field labels from the comment fields in the fits header
        """
        header = self._f[1].header
        columns = [field[5:] for field in header if field.startswith('TTYPE')]
        field_labels = {}
        for column in columns:
            label = header.get('TCOMM'+column, header['TTYPE'+column])
            unit = header.get('TUNIT'+column, '')
            if unit:
                label = '{} ({})'.format(label, unit)
            field_labels[header['TTYPE'+column].lower()] = label
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
        
        if all or self.selection==None:
            return data
        else:
            return data[self.selection]


    def __getitem__(self, field):
        """
        Return field with given name. See self.get().
        """
        return self.get(field)


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
    
    def name_like(self, name):
        """
        Select the source with given name. Uses first result from self.get_name_like()
        """
        ind = self.get_name_like(name)
        if len(ind)>0:
            self.selection = ind[0]
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
    
    def get_distances(self):
        """
        Create an array of distances for all sources in self.name
        
        Distances taken from Kuulkers et al. (2003) for Globular Cluster
        sources, otherwise from the Liu et al. (200?) lmxb catalog. In
        individual cases, the reference is indicated in a comment. If
        uncertainty is unknown, it defaults to 0.

        TODO: check for missing entries; add distance reference, similar to NH.
        """
        distances = {'4U 0513-40': (12.1, 0.3),
                     '4U 0614+09': (3.0, 0),
                     'EXO 0748-676': (8.3, 0),
                     '4U 0836-429': (11, 0),
                     '2S 0918-549': (4.0, 0),
                     '4U 1246-588': (13.0, 0),
                     '4U 1254-69': (0, 0),
                     '4U 1323-62': (11.0, 0),
                     'SAX J1324.5-6313': (6.2, 0),
                     '4U 1608-522': (5.16, 0.72),
                     '4U 1636-536': (5.95, 0.12), # Galloway et al. 2008 solar PRE
                     'MXB 1658-298': (12, 0),
                     '4U 1702-429': (5.5, 0),
                     '4U 1705-32': (13.0, 2.0),
                     '4U 1705-44': (7.5, 1.1),
                     'IGR J17062-6143': (7.3, 0.5), # Keek et al. 2017
                     '4U 1708-40': (0.0, 0),
                     'XTE J1709-267': (8.8, 0),
                     'XTE J1710-281': (14.0, 2.0),
                     '2S 1711-339': (7.5, 0),
                     'SAX J1712.6-3739': (7.0, 0),
                     '1H 1715-321': (0.0, 0),
                     'RX J1718.4-4029': (6.5, 0),
                     'XTE J1723-376': (13.0, 0),
                     '4U 1722-30': (9.5, 2.5),
                     '4U 1728-34': (5.2, 0),
                     'MXB 1730-335': (8.8, 3.3),
                     'KS 1731-260': (7.0, 0),
                     'SLX 1732-304': (4.3, 2.3),
                     'SLX 1735-269': (7.3, 0),
                     '4U 1735-444': (8.0, 1.0),
                     'SLX 1737-282': (6.5, 1.5),
                     'XTE J1739-285': (10.6, 0),
                     'KS 1741-293': (0.0, 0),
                     'GRS 1741.9-2853': (7.5, 2.5),
                     '1A 1742-294': (10.0, 0),
                     '1A 1744-361': (0.0, 0),
                     'SLX 1744-300': (6.7, 0),
                     '1A 1744-361': (11.0, 0),
                     'EXO 1745-248': (5.5, 0.9), # 2007A&A...470.1043O
                     '4U 1746-37': (11.0, 0.9),
                     'IGR J17473-2721': (0, 0),
                     'SAX J1747.0-2853': (8, 0), # 2006ApJS..165..173M
                     'GRS 1747-312': (9.5, 3.3),
                     'SAX J1748.9-2021': (8.4, 1.5),
                     'SAX J1750.8-2900': (5.0, 0),
                     'SAX J1752.3-3138': (9.0, 0),
                     'SAX J1753.5-2349': (6.0, 0),
                     'IGR J17597-2201': (5.0, 0),
                     '2S 1803-245': (0, 0),
                     'SAX J1806.5-2215': (8.0, 0),
                     'SAX J1808.4-3658': (2.5, 0),
                     'SAX J1810.8-2609': (4.9, 0),
                     '4U 1812-12': (4.1, 0),
                     'XTE J1814-338': (8.0, 1.6),
                     'SAX J1818.7+1424': (9.4, 0),
                     '4U 1820-303': (7.6, 0.4),
                     'AX J1824.5-2451': (5.6, 0),
                     'GS 1826-24': (6.0, 0),
                     'SAX J1828.5-1037': (6.2, 0),
                     'XB 1832-330': (9.6, 0.4),
                     'Ser X-1': (8.4, 0),
                     '4U 1850-086': (8.2, 0.6),
                     'HETE J1900.1-2455': (4.7, 0.6),
                     'XB 1905+000': (0.0, 0),
                     'XB 1916-053': (8.9, 0),
                     'XTE J2123-058': (8.5, 2.5),
                     '4U 2129+12': (10.3, 0.4),
                     'SAX J2224.9+5421': (7.1, 0),
                     'Aql X-1': (5.0, 0.0),
                     'Cir X-1': (0.0, 0),
                     'Cyg X-2': (11.6, 0),
                     'GX 17+2': (9.8, 0),
                     'GX 3+1': (6.5, 0),
                     'GX 13+1': (7.0, 0), # Christian & Swank 1997
                     }
        
        dist = np.zeros_like(self['ra_obj'])
        diste = np.zeros_like(dist)
        for i, name in enumerate(self['name']):
            if name in distances:
                dist[i] = distances[name][0]
                diste[i] = distances[name][1]
        return dist, diste


class Instrument:
    """
    Here's a generic instrument class which can be adapted and/or duplicated for different
    instruments. This class is kept pretty lean to avoid having to replicate lots of code
    Defines the properties of an instrument with data that we're going to analyse and add to MINBAR
    Example
    import minbar
    jemx = minbar.Instrument('JEM-X', 'jemx', 'IJ')
    """

    def __init__(self, name, path=None, label=None,
                 lightcurve=['lc1_3-25keV_1.0s.fits','lc2_3-25keV_1.0s.fits'],
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
            logger.error('need a valid path for the data files')
            return None

        # Define the paths corresponding to each source here
        # Default is just the source name with spaces removed; this won't work for RXTE

        self.source = Sources()
        self.source_name = self.source['name']
        if self.label == 'XP':
            # for RXTE/PCA, most sources just omit the prefix
            self.source_path = np.array([x.split()[1] for x in self.source_name])
            # The standard product lightcurve also includes the obsid, so define the lightcurve
            # parameter as a function here
            lightcurve = self.pca_lightcurve_filename
        else:
            # more commonly we just remove the spaces
            self.source_path = self.source_name.replace(" ", "")

        if self.label == 'IJ':
            # for JEM-X, the convention is to also replace '+' with 'p'
            self.source_path = self.source_path.replace("+","p")

        # Define the lightcurve and spectral files

        assert lightcurve is not None
        self.lightcurve = lightcurve
        self.spectrum = spectrum

        # Check that the directories in the data path all correspond to source_paths
        # You don't want to miss observations in some directory because of inconsistencies
        # with the directory names

        self.has_dir, self.nonmatched, self.nmissing = verify_path(
            self.source, self.path, self.source_path, verbose=verbose)


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


    def pca_lightcurve_filename(self, obsid):
        """
        Function to return the lightcurve name for PCA observations
        :return:
        """

        return 'stdprod/xp{}_n1.lc.gz'.format(obsid.replace("-", ""))


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
