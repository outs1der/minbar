"""
Python class to easily extract information from MINBAR.
Provides the burst catalog (Bursts), observation catalog (Observations),
and source catalog (Sources). See the doc strings of those classes
for example usage.

The data files should be placed in the 'data' subdirectory of this
package.

(c) 2007, Laurens Keek, l.keek@sron.nl
Updated for MINBAR v0.9, 2017, Laurens Keek, laurens.keek@nasa.gov
"""

from idldatabase import IDLDatabase
import numpy as n
import os, re
from astropy.io import fits
import logging

logger = logging.getLogger('minbar')
logger.setLevel(logging.INFO)

kpc = 3.086e21 # cm

# List of ultra compacts from In 't Zand (2007) (1850-087 and 1905+000 no bursts in MINBAR)
# Includes all candidates. Should I add 1728?
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

class Bursts(IDLDatabase):
    """
    Read the MINBAR IDL database and give access to its
    contents.
    
    Example usage:
    import minbar
    mb = minbar.Bursts()
    mb.select_like('1636') # Select a source using part of its name
    print mb.field_labels # See which fields are available
    time = mb['time'] # Get a field as a numpy array (automatically time-ordered)
    flux = mb['pflux']*1e-9 # Flux in erg/s/cm2
    mb.create_distance_correction() # Include distance information from Sources()
    distance = mb['dist']
    luminosity = flux*mb['distcor'] # Luminosity in erg/s
    
    mb.select_like('1826') # Replace selection by another source
    mb.select_all(['GS 1826-24', '4U 1636-536']) # Select multiple sources; requires exact names
    mb.clear() # Clear the selection so all sources are included
    mb.exclude_like('1636') # Exclude source from selection
    mb.exclude_like('1826') # Now two sources are excluded
    """
    
    timefield = 'time' # The field used for determining time order
    entryname = 'burst'
    
    def __init__(self, filename=None, type=1):
        """
        Create a new Bursts instance using the data from the minbar database.
        
        filename: path to the database files, excluding their extension.
                  By default the minbar database in the directory of this
                  script is used.
        type: burst type. The default, 1, selects all vetted Type I bursts. Setting it
              to None means no type is selected.
        """
        if filename==None:
            filename = self.get_default_path('minbar')
        IDLDatabase.__init__(self, filename)
        
        self.fix_labels()
        
        self.names = self.get_names()
        self.type = type
        self.clear()
    
    def get_default_path(self, filename):
        """
        Return the default path of the minbar data files with prefix
        filename
        """
        return os.path.join(os.path.dirname(__file__), 'data', filename)
    
    def fix_labels(self):
        """
        Fix the whitespace in the labels.
        """
        pat = re.compile('\s+')
        for k in self.field_labels:
            self.field_labels[k] = re.sub(pat, ' ', self.field_labels[k])
    
    def get_names(self):
        """
        Get a list of all source names in the archive. Ordered by right
        ascension.
        """
        names = n.unique(self.records.field('name'))
        
        ra = n.array([self.records.field('ra')[self.records.field('name')==name][0] for name in names])
        ind = n.argsort(ra)
        
        return names[ind]
    
    def get_name_like(self, name):
	"""
	Return a list of sources that have 'name' in their archive
	identifier.
	"""
	selection = []
	for i in self.names:
		if i.find(name)>-1:
			selection.append(i)
	return selection
    
    def get_type(self):
        """
        Return an index array selecting the specified burst type (self.type).
        """
        if self.type==None:
            return n.ones(len(self.records), n.bool)
        else:
            return self.records.field('type')==self.type
    
    def clear(self):
        """
        Clear the selection. If self.type is not None, only bursts of the given
        type are selected.
        """
        self.name = ''
        self.selection = self.get_type()
        self.time_order = n.argsort(self.records[self.selection].field(self.timefield))
        self.ind = n.where(self.selection)[0][self.time_order]

    def _pad_name(self, name):
        """
        Source names in records are padded with spaces. This routine pads with
        given name with spaces for comparison to the records.
        """
        size = self.records['name'].dtype.itemsize
        return name.ljust(size)
    
    def select(self, name):
        """
        Set the selection to the source with given name.
        """
        self.name = name
        name = self._pad_name(name)
        self.selection = (self.records.field('name') == name)&self.get_type() # Match name and burst type
        self.time_order = n.argsort(self.records[self.selection].field(self.timefield))
        self.ind = n.where(self.selection)[0][self.time_order]
        logger.info('Selected {} {}s from {}'.format(len(n.where(self.selection)[0]), self.entryname, name))
    
    def select_like(self, name):
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
            if len(names)>1:
                logger.info('{} more matching sources: {}'.format(len(names) - 1, ', '.join(names[1:])))
    
    def select_all(self, names):
        """
        Select multiple sources with given names.
        """
        selection = n.zeros_like(self.selection)
        for name in names:
            name = self._pad_name(name)
            selection = n.logical_or(selection, self.records.field('name') == name)
        self.selection = selection&self.get_type() # Only bursts of specified type
        self.time_order = n.argsort(self.records[self.selection].field(self.timefield))
        self.ind = n.where(self.selection)[0][self.time_order]
        logger.info('Selected {} {}s'.format(len(n.where(self.selection)[0]), self.entryname))
    
    def exclude(self, name):
        """
        Removes source with given name from the current selection.
        """
        selection = n.logical_not(self.records.field('name') == name)
        self.selection = n.logical_and(self.selection, selection)
        self.time_order = n.argsort(self.records[self.selection].field(self.timefield))
        self.ind = n.where(self.selection)[0][self.time_order]
        logger.info('Selected {} {}s by excluding {}'.format(len(n.where(self.selection)[0]), self.entryname, name))
    
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
    
    def __len__(self):
        """
        Return the number of entries in the current selection.
        """
        return len(self.ind)
    
    def __getitem__(self, field):
        """
        Shorthand for get().
        """
        return self.get(field)
    
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
            return n.zeros(len(self.get(field)))
    
    def get_error_name(self, field):
        """
        Return the name of the field containing the
        error for field.
        """
        return field + 'e'
    
    def instr_like(self, instrument):
        """
        Return an array that selects all entries where the instrument name
        begins with the given instrument string. E.g., 'XP' for PCA. For convenience,
        the following aliases are provided: 'pca', 'wfc', 'jemx'.
        """
        alias = {'pca': 'XP', 'wfc': 'SW', 'jemx': 'IJ'}
        if instrument in alias:
            instrument = alias[instrument]
        
        return n.char.array(self['instr']).startswith(instrument)
    
    def create_distance_correction(self):
        """
        Create an array of distance corrections, 4 pi d^2, for each burst.
        This can be used to easily convert flux to luminosity.
        
        Returns correction factor (cm^2), error, distance (kpc), error. Furthermore,
        these arrays are stored and can be accessed as self['fieldname'], with fieldnames
        distcor, distcore, dist, diste.
        """
        s = Sources()
        dist = n.zeros(len(self.records))
        diste = n.zeros_like(dist)
        cor = n.zeros_like(dist)
        core = n.zeros_like(dist)
        
        names = self.records.field('name')
        for name in self.names:
            s.select_like(name.strip())
            if s.selection!=None:
                ind = names==name
                dist[ind] = s['dist']
                diste[ind] = s['diste']
            s.clear()
        
        cor = 4*n.pi*(dist*kpc)**2
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

class Minbar(Bursts):
    """
    Bursts used to be called Minbar. This class is defined for backwards compatibility.
    """
            
class Observations(Bursts):
    """
    Load MINBAR database of observations.
    """
    
    timefield = 'tstart' # The field used for determining time order
    entryname = 'observation'
    
    def __init__(self, filename=None, type=None):
        """
        Load the database of observations.
        """
        if filename==None:
            filename = self.get_default_path('minbar-obs')
        Minbar.__init__(self, filename, type)

    def __str__(self):
        """
        Return a nice string.
        """
        return "Multi-INstrument oBservation ARchive (MINBAR) ({} {}s from {} sources)".format(len(self.records), self.entryname, len(self.names))

class Sources:
    """
    Contain all information on the sources present in Minbar in
    the file minbar_sources.fits.
    
    Example:
    s = Sources()
    print s.field_names # Show available data fields
    ra = s['ra_obj'] # Right ascension for all sources
    s.select_like('1636')
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
	return n.array(selection)
    
    def select_like(self, name):
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
        
        dist = n.zeros_like(self['ra_obj'])
        diste = n.zeros_like(dist)
        for i, name in enumerate(self['name']):
            if name in distances:
                dist[i] = distances[name][0]
                diste[i] = distances[name][1]
        return dist, diste
