.. highlight:: shell

============
Installation
============

* Clone/download the repository to your local machine
* Copy the MINBAR table files ``minbar.txt``, ``minbar-obs.txt``, and ``minbar_sources.fits`` (rename the downloaded ``minbar_sources_vXX.fits`` to ``minbar_sources.fits``) from the `data repository <https://doi.org/10.26180/5e4a697d9b8b6>`_ into the ``data`` directory, or create softlinks in that directory to the respective locations for each file.
* Start python and ``import minbar``. 

With the instructions above you can use the ``minbar`` package from the
directory in which you installed it. If you want to make it available from
any directory, you should

* add the path to your ``PYTHONPATH`` environment variable; or,
* install the package using ``python3 -m pip install .`` from the distribution directory.

**NOTE** that for the second option above, you will also need to copy the table files to a ``data`` subdirectory in the installation directory; you can find where this is by running (after installation)

.. code-block:: python

    import minbar
    print (minbar.__file__)

Which should return a path that looks something like ``.../miniforge3/envs/minbar/lib/python3.9/site-packages/minbar/__init__.py`` (the bit to the left of ``site_packages`` will depend on your specific installation). Create the ``data`` directory as a subdirectory of ``site-packages/minbar``, and put a copy of (or softlink to) the data files there.

Try the suggested example commands in the next section; see the tutorial jupyter notebook; and also check ``__init__.py`` for additional examples and full parameter lists
