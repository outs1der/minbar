.. highlight:: shell

============
Installation
============

* Clone/download the repository to your local machine
* Place the ``minbar`` directory in your python path, e.g., the ``site-packages`` directory.
* Alternatively, install the package using ``python3 -m pip install .`` from the distribution directory
* Copy the MINBAR table files ``minbar.txt``, ``minbar-obs.txt``, and ``minbar_sources.fits`` (rename the downloaded ``minbar_sources_vXX.fits`` to ``minbar_sources.fits``) into the ``data`` directory, or softlink them to the respective locations.
* Start python and ``import minbar``. 
* Base classes are :class:`minbar.Bursts`, :class:`minbar.Observations`, and :class:`minbar.Sources`; creating instances of each one will read the data from the corresponding table file
* Try the suggested example commands in the next section; see the tutorial jupyter notebook; and also check ``__init__.py`` for additional examples and full parameter lists
