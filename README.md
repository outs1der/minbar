# README #

### python-minbar ###

* Python package for analysing observational data of thermonuclear (type-I) X-ray bursts observed by _RXTE_/PCA, _BeppoSAX_/WFC and _INTEGRAL_/JEM-X
* Provides facilities for reading MINBAR ASCII data tables, available via the online data table repository at https://doi.org/10.26180/5e4a697d9b8b6; see https://arxiv.org/abs/2003.00685 
* Requires Python 3, NumPy, and Astropy (originally adapted from Python 2)
* Compatible with MINBAR data release 1
* Includes code, under development, to extend the MINBAR sample
* Can also read in IDL database files (not  publicly available), from which the ASCII tables are derived

### Getting Started ###

* Place the `minbar` directory in your python path, e.g., the `site-packages` directory.
* Copy the MINBAR table files `minbar.txt` and `minbar-obs.txt` into the `data` directory, or softlink to the locations
* Start python and `import minbar`. 
* Base classes are `Bursts`, `Observations`, and `Sources`; creating instances of each one will read the data from the corresponding table file
* Try the suggested example commands below; see the tutorial jupyter notebook; and also check `__init__.py` for instructions on usage

### Example Usage ###

The examples below originally utilised the IDL database files, and may not yet be implemented with the MRT tables.

```
#!python

import minbar
mb = minbar.Bursts() # Load the burst database
mb.name_like('1636') # Select a source using part of its name
print (mb.field_labels.keys()) # See which fields are available
mb.show() # List the selected bursts

time = mb['time'] # Get a field as a numpy array (automatically time-ordered)
mb[time > 54000.].show() # show all the bursts after the specified time
flux = mb['bpflux'] # Flux in 1e-9 erg/s/cm2
sub = mb[['time','bpflux']] # extract a subset of the columns, for the given selection
mb.create_distance_correction() # Include distance information from Sources()
luminosity = (flux*mb['distcor']).to('erg s-1') # Isotropic peak luminosity in erg/s
pca = mb.instr_like('pca') # Get index array for bursts observed with PCA
pca_luminosity = luminosity[pca] # Luminosity of PCA bursts

mb.clear() # Reset the selection
mb.select_all(['GS 1826-24', '4U 1636-536']) # Select multiple sources; requires exact names
mb.clear() # Clear the selection so all sources are included
mb.exclude_like('1636') # Exclude source from selection
mb.exclude_like('1826') # Now two sources are excluded

mo = minbar.Observations() # Load the observation database
mo.name_like('1636') # Same source selection options as for burst database
time = mo['tstart'] # And fields are accessed in the same way
print (mo.field_labels.keys()) # See which fields are available

ms = minbar.Sources() # Load the source database
print (ms.field_labels.keys()) # Show available data fields
ra = ms['ra_obj'] # Right ascension for all sources
ms.name_like('1636') # Select a source using part of its name
ra = ms['ra_obj'] # Right ascension for selected source only
ms.clear() # Clear selection
```

Note that Sources() methods do not include `select_all()` or `exclude_like()`

Below are some basic examples to analyse some (new?) X-ray data and search for bursts
(under development)

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

