=====
Usage
=====

The base classes to access bursts, burst observations, and burst sources are :class:`minbar.Bursts`, :class:`minbar.Observations`, and :class:`minbar.Sources`, respectively. Creating instances of each one will read the data from the corresponding table file.

Here we provide basic usage; more examples are provided in the tutorial jupyter notebook


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

Let's explore that below with some slightly more complex includes and excludes

.. code-block:: python

    mb.clear() # Reset the selection
    mb.select_all(['GS 1826-24', '4U 1636-536']) # Select multiple sources; requires exact names
    mb.clear() # Clear the selection so all sources are included
    mb.exclude_like('1636') # Exclude source from selection
    mb.exclude_like('1826') # Now two sources are excluded

Finally, here's some analysis which involves working with the fluxes, and
estimating peak luminosity from the bolometric peak flux (table attribute
``bpflux``)

.. code-block:: python

    time = mb['time'] # Get a field as a numpy array (automatically time-ordered)
    id = mb[time > 54000.]['entry'] # extract ID #s for all the bursts after the specified time
    flux = mb['bpflux'] # Flux in 1e-9 erg/s/cm2
    sub = mb[['time','bpflux']] # extract a subset of the columns, for the given selection
    mb.create_distance_correction() # Include distance information from Sources()
    luminosity = (flux*mb['distcor']).to('erg s-1') # Isotropic peak luminosity in erg/s
    pca = mb.instr_like('pca') # Get index array for bursts observed with PCA
    pca_luminosity = luminosity[pca] # Luminosity of PCA bursts

A complete list of the :class:`minbar.Bursts` `table attributes <https://burst.sci.monash.edu/static/html/columns.html#bursts>`_, along with those for the :class:`Observations` and :class:`Sources` tables, can be found at the alternate web interface for the MINBAR data, `http://burst.sci.monash.edu <http://burst.sci.monash.edu>`_

2. Working with observations
----------------------------

As both the :class:`minbar.Bursts` and :class:`minbar.Observations`
classes are built on the underlying :class:`minbar.Minbar` class, many of
the methods are common to both classes. 

Below whe show an example of selecting all the observations from 4U
1636-536 and extracting the start times.

.. code-block:: python

    mo = minbar.Observations() # Load the observation database
    mo.name_like('1636') # Same source selection options as for burst database
    time = mo['tstart'] # And fields are accessed in the same way
    print (mo.field_labels.keys()) # See which fields are available

3. Working with sources
----------------------------

The :class:`minbar.Sources` class is a little different as it is really just a wrapper for the underlying FITS table. Still, some of the methods as for the other two classes are available, including :meth:`minbar.Sources.name_like`.

You can also select sources by type, e.g. ``C`` for ultracompact, or ``S``
for sources that have shown a superburst, or combinations of the two

.. code-block:: python

    ms = minbar.Sources() # Load the source database
    print (ms.field_labels.keys()) # Show available data fields
    ra = ms['ra_obj'] # Right ascension for all sources
    ms.name_like('1636') # Select a source using part of its name
    ra = ms['ra_obj'] # Right ascension for selected source only
    ms.clear() # Clear selection
    ms.type('SC') # Select all ultracompacts that have shown a superburst
    ms['name'] # ... and show their names
