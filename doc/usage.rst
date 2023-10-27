=====
Usage
=====

Here we provide basic usage; more examples are provided in the tutorial
jupyter notebook


1. Working with bursts
----------------------

Here we initialise a :class:`minbar.Bursts` object, and select all the bursts from 4U 1636-536

.. code-block:: python

    import minbar
    mb = minbar.Bursts() # Load the burst database
    mb.name_like('1636') # Select a source using part of its name
    print (mb.field_labels.keys()) # See which fields are available
    mb.show() # List the selected bursts

The selection made with the :py:meth:`minbar.Minbar.name_like` and equivalent commands is *persistent*, and any subsequent query (e.g. extracting one of the table columns) will be restricted to the same set of events. You can reset to the full sample with :func:`minbar.Minbar.clear()`

Let's explore that below

.. code-block:: python

    time = mb['time'] # Get a field as a numpy array (automatically time-ordered)
    id = mb[time > 54000.]['entry'] # extract ID #s for all the bursts after the specified time
    flux = mb['bpflux'] # Flux in 1e-9 erg/s/cm2
    sub = mb[['time','bpflux']] # extract a subset of the columns, for the given selection
    mb.create_distance_correction() # Include distance information from Sources()
    luminosity = (flux*mb['distcor']).to('erg s-1') # Isotropic peak luminosity in erg/s
    pca = mb.instr_like('pca') # Get index array for bursts observed with PCA
    pca_luminosity = luminosity[pca] # Luminosity of PCA bursts

Now some slightly more complex includes and excludes

.. code-block:: python

    mb.clear() # Reset the selection
    mb.select_all(['GS 1826-24', '4U 1636-536']) # Select multiple sources; requires exact names
    mb.clear() # Clear the selection so all sources are included
    mb.exclude_like('1636') # Exclude source from selection
    mb.exclude_like('1826') # Now two sources are excluded

2. Working with observations
----------------------------

.. code-block:: python

    mo = minbar.Observations() # Load the observation database
    mo.name_like('1636') # Same source selection options as for burst database
    time = mo['tstart'] # And fields are accessed in the same way
    print (mo.field_labels.keys()) # See which fields are available

3. Working with sources
----------------------------

.. code-block:: python

    ms = minbar.Sources() # Load the source database
    print (ms.field_labels.keys()) # Show available data fields
    ra = ms['ra_obj'] # Right ascension for all sources
    ms.name_like('1636') # Select a source using part of its name
    ra = ms['ra_obj'] # Right ascension for selected source only
    ms.clear() # Clear selection
