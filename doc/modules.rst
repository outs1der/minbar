=====================
Classes and functions
=====================

:class:`minbar.Minbar` is the base class for :class:`minbar.Bursts` and :class:`minbar.Observations`, so the functions below are common to both those objects

:class:`minbar.Sources` is the other main class that is derived from the tables accompanying the paper, but it's architecture is a bit different, and is based on the FITS table.

.. autoclass:: minbar.Minbar
   :members: show, get_default_path, __len__, get_type, clear, select, sort, obsid, name_like, attr_like, get_name_like, _pad_name, get, get_records, __getitem__, instr_like, instr, instr_exclude, select_all, exclude, exclude_like, exclude_flag

.. autoclass:: minbar.Bursts
   :members: __init__, get_burst_data, burstplot, get_lc, PRE, unique, __str__, has_error, is_error, get_error, get_error_name, create_distance_correction

.. autoclass:: minbar.Observations
   :members: __init__, good, __str__, plot

.. autoclass:: minbar.Sources
   :members: __init__, get_default_path, _get_field_labels, get, __getitem__, type, get_name_like, name_like, clear, get_F_Edd, get_distances, __str__

