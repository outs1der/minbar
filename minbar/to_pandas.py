"""
Load Minbar databases as pandas DataFrames

Example:
mb = load_bursts()
mo = load_observations()
ms = load_sources()
"""

import os
import numpy as n
import pandas as pd
import minbar
from astropy.table import Table

class MinbarFrame(pd.DataFrame):
    """
    A DataFrame with Minbar functionality
    
    Example:
    df = minbar2pandas()
    aql = df.get_source('aql') # Select source with partial name
    pca = df.get_instrument('pca') # Select instrument 'pca', 'wfc', or 'jemx'
    """
    
    def get_source(self, name):
        """
        Return a MinbarFrame with only given source
        name: source name (case insensitive). If a partial name is 
              given, the first match is used.
        """
        if not hasattr(self, '_source_names'):
            self._source_names = self['name'].cat.categories.astype('str').str.lower()
        names = self._source_names
        name = name.lower()
        ind = names.str.contains(name)
        ind = n.where(ind)[0]
        if len(ind)==0:
            minbar.logger.info('No matching source')
        else:
            selected_name = self['name'].cat.categories[ind[0]]
            df = self[self.name==selected_name]
            minbar.logger.info('Selected {} entries from {}'.format(len(df), selected_name.strip()))
            if len(ind)>1:
                other_names = [i.strip() for i in self['name'].cat.categories[ind[1:]]]
                minbar.logger.info('{} more matching sources: {}'.format(len(other_names), ', '.join(other_names)))
            return df
    
    def get_instrument(self, name):
        """
        Return a MinbarFrame with only given instrument
        name: instrument name. This may be only the start to match all
              configurations. E.g, 'XP' selects all XTE/PCA configurations.
              The following aliases are provided: 'pca', 'wfc', 'jemx'
        """
        name = {'pca': 'XP', 'wfc': 'SW', 'jemx': 'IJ'}.get(name, name) # Aliases
        ind = self.instr.str.decode('ascii').str.startswith(name)
        return self[ind]

def _minbar2pandas(bursts=True):
    """
    Return minbar as a pandas' DataFrame
    bursts: whether to load bursts or observations instead
    """
    if bursts:
        mb = minbar.Bursts()
    else:
        mb = minbar.Observations()
    data = mb.records
    data = data.byteswap().newbyteorder('N') # Switch to native byteorder for pandas
    df = MinbarFrame.from_records(data)
    df['name'] = df['name'].astype('category')
    return df

def load_bursts():
    """
    Load the Minbar burst database as a pandas DataFrame
    """
    return _minbar2pandas()

def load_observations():
    """
    Load the Minbar observation database as a pandas DataFrame
    """
    return _minbar2pandas(False)

def load_sources():
    """
    Load the Minbar Sources database as a DataFrame. Columns with multi-d 
    fields are not supported by pandas.
    Note: field 'NAME' is renamed to 'name', to match the field in the other
    Minbar databases.
    """
    table = Table.read(os.path.join(os.path.dirname(minbar.__file__), 'data/minbar_sources.fits'))
    # Remove multi-d columns, which are not supported by pandas
    table.remove_columns([name for name in table.colnames if len(table[name].shape)>1])
    df = table.to_pandas()
    df = df.rename(columns={'NAME': 'name'})
    return df
