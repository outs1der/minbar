# README #

### python-minbar ###

Python package for analysing observational data of thermonuclear (type-I) X-ray bursts observed by _RXTE_/PCA, _BeppoSAX_/WFC and _INTEGRAL_/JEM-X

Provides facilities for reading MINBAR ASCII data tables, available via the online data table repository at https://doi.org/10.26180/5e4a697d9b8b6; see https://arxiv.org/abs/2003.00685 

* Requires Python 3, NumPy, and Astropy (originally adapted from Python 2)
* Compatible with MINBAR data release 1
* Includes code, under development, to extend the MINBAR sample
* Can also read in IDL database files (not  publicly available), from which the ASCII tables are derived

Full documentation can be found at https://burst.sci.monash.edu/minbar/pydoc/index.html

### Getting Started ###

* Clone/download the repository to your local machine
* Place the `minbar` directory in your python path, e.g., the `site-packages` directory.
* Alternatively, install the package using `python3 -m pip install ./` from the distribution directory
* Copy the MINBAR table files `minbar.txt`, `minbar-obs.txt`, and `minbar_sources.fits` (rename the downloaded `minbar_sources_vXX.fits` to `minbar_sources.fits`) into the data directory, or softlink them to the respective locations.
* Start python and `import minbar`. 
* Base classes are `Bursts`, `Observations`, and `Sources`; creating instances of each one will read the data from the corresponding table file

Try the suggested example commands below; see the [online documentation](https://burst.sci.monash.edu/minbar/pydoc) and the tutorial jupyter notebook; and also check `__init__.py` for instructions on usage

### Example Usage ###

The examples below illustrate creating examples of the available `Bursts`, `Observations`, and `Sources` classes, and performing selections and data extractions from the tables.

```
#!python

import minbar
mb = minbar.Bursts() # Load the burst database
mb.name_like('1636') # Select a source using part of its name
print (mb.field_labels.keys()) # See which fields are available
mb.show() # List the selected bursts

time = mb['time'] # Get a field as a numpy array (automatically time-ordered)
id = mb[time > 54000.]['entry'] # extract ID #s for all the bursts after the specified time
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

