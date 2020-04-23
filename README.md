# README #

### python-minbar ###

* Python package for analysing observational data of thermonuclear (type-I) X-ray bursts
* Provides facilities for reading MINBAR ASCII and IDL database files (latter are not available publicly); see https://arxiv.org/abs/2003.00685 
and the online data table repository at https://doi.org/10.26180/5e4a697d9b8b6
* Includes code, under development, to extend the MINBAR sample
* Requires Python 2 or 3, NumPy, and Astropy
* Compatible with MINBAR data release 1

### Getting Started ###

* Place the `minbar` directory in your python path, e.g., the `site-packages` directory.
* Copy the MINBAR database files (normally this would be the MRT tables) into the `data` directory, or softlink to the locations
* Start python and `import minbar`. Check `__init__.py` for instructions on usage.
* All data are presented as numpy arrays, which allows for easy selection of subsets of the data and operations on such selections. See [this introduction](https://docs.scipy.org/doc/numpy-dev/user/quickstart.html) to get started.

### Example Usage ###

The examples below originally utilised the IDL database files, and may not yet be implemented with the MRT tables.

```
#!python

import minbar
mb = minbar.Bursts() # Load the burst database
mb.name_like('1636') # Select a source using part of its name
print (mb.field_labels) # See which fields are available
time = mb['time'] # Get a field as a numpy array (automatically time-ordered)
flux = mb['pflux']*1e-9 # Flux in erg/s/cm2
mb.create_distance_correction() # Include distance information from Sources()
distance = mb['dist']
luminosity = flux*mb['distcor'] # Luminosity in erg/s
pca = mb.instr_like('pca') # Get index array for bursts observed with PCA
pca_luminosity = luminosity[pca] # Luminosity of PCA bursts

mb.name_like('1826') # Replace selection by another source
mb.select_all(['GS 1826-24', '4U 1636-536']) # Select multiple sources; requires exact names
mb.clear() # Clear the selection so all sources are included
mb.exclude_like('1636') # Exclude source from selection
mb.exclude_like('1826') # Now two sources are excluded

mo = minbar.Observations() # Load the observation database
mo.name_like('1636') # Same source selection options as for burst database
time = mo['tstart'] # And fields are accessed in the same way
print (mo.field_labels) # See which fields are available

ms = minbar.Sources() # Load the source database
print (ms.field_labels) # Show available data fields
ra = ms['ra_obj'] # Right ascension for all sources
ms.name_like('1636') # Select a source using part of its name
ra = ms['ra_obj'] # Right ascension for selected source only
ms.clear() # Clear selection
```

Note that Sources() does not have select_all() or exclude_like()

Below are some basic examples to analyse some (new?) X-ray data and search for bursts

```
import minbar

xte = minbar.Instrument('PCA') # Create an instrument definition
obs = minbar.Observation(None, xte, '4U 1636-536', '10088-01-07-02')

obs.plot()
print (obs.mjd, obs.rate)
```

Can also define a new instrument for analysis of data from instruments not originally part of MINBAR

```
xmm = minbar.Instrument('XMM-Newton', 'xmm', 'XN', '2to7good.fits')
obs = minbar.Observation(None, xmm, '1RXS J180408.9-342058', '0741620101')
lc =obs.get_lc()

import matplotlib.pyplot as plt
plt.plot(lc['TIME'], lc['RATE'])
plt.show()

test = minbar.findburst(lc['TIME'], lc['RATE'], lc['ERROR'])
print(test)
[5.42058957e+08 5.42067368e+08 5.42075296e+08 5.42083081e+08
 5.42090903e+08]
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