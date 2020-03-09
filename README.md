# README #

### python-minbar ###

* Python package for analysing observational data of thermonuclear (type-I) X-ray bursts
* Includes code to update and extend the MINBAR sample; see https://arxiv.org/abs/2003.00685 
and the online data table repository at https://doi.org/10.26180/5e4a697d9b8b6
* Also provides facilities for reading MINBAR IDL database files (not released publicly)
* Requires Python 2 or 3, NumPy, and Astropy
* Compatible with MINBAR data release 1

### Getting Started ###

* Place the `minbar` directory in your python path, e.g., the `site-packages` directory.
* Copy the MINBAR database files (normally this would be the MRT tables) into the `data` directory.
* Start python and `import minbar`. Check `__init__.py` for instructions on usage.
* All data are presented as numpy arrays, which allows for easy selection of subsets of the data and operations on such selections. See [this introduction](https://docs.scipy.org/doc/numpy-dev/user/quickstart.html) to get started.

### Example Usage ###

The examples below originally utilised the IDL database files, and may not yet be implemented with the MRT tables.

```
#!python

import minbar
mb = minbar.Bursts() # Load the burst database
mb.select_like('1636') # Select a source using part of its name
print (mb.field_labels) # See which fields are available
time = mb['time'] # Get a field as a numpy array (automatically time-ordered)
flux = mb['pflux']*1e-9 # Flux in erg/s/cm2
mb.create_distance_correction() # Include distance information from Sources()
distance = mb['dist']
luminosity = flux*mb['distcor'] # Luminosity in erg/s
pca = mb.instr_like('pca') # Get index array for bursts observed with PCA
pca_luminosity = luminosity[pca] # Luminosity of PCA bursts

mb.select_like('1826') # Replace selection by another source
mb.select_all(['GS 1826-24', '4U 1636-536']) # Select multiple sources; requires exact names
mb.clear() # Clear the selection so all sources are included
mb.exclude_like('1636') # Exclude source from selection
mb.exclude_like('1826') # Now two sources are excluded

mo = minbar.Observations() # Load the observation database
mo.select_like('1636') # Same source selection options as for burst database
time = mo['tstart'] # And fields are accessed in the same way
print mo.field_labels # See which fields are available

ms = minbar.Sources() # Load the source database
print ms.field_labels # Show available data fields
ra = ms['ra_obj'] # Right ascension for all sources
ms.select_like('1636') # Select a source using part of its name
ra = ms['ra_obj'] # Right ascension for selected source only
ms.clear() # Clear selection
# Note that Sources() does not have select_all() or exclude_like()
```

### pandas ###

[pandas](https://pandas.pydata.org/) is a data analysis library with many
features. The Minbar databases can be loaded as pandas DataFrames:

```
#!python

import minbar.to_pandas
mb = minbar.to_pandas.load_bursts() # Create a DataFrame with burst data
mo = minbar.to_pandas.load_observations()
ms = minbar.to_pandas.load_sources()
df = mb.get_source('1636') # Select a source using part of its name
df = mb.get_instrument('pca') # Select an instrument ('pca', 'wfc', or 'jemx')
```