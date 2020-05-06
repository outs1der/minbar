ananew.f
========

- contains fortran code for fitting the light curves
- can be compiled with gfortran -ffixed-line-length-none
- is badly commented (sorry)

doall
=====

- c-shell script to fit all light curves
- the data (ascii light curves) need to be in DATA/WFC, DATA/RXTE
  and DATA/JEMX. These directories now contain only the files for
  4U 1722-30. The other ~7000 files were already sent
- calls cshell script dofitit

dofitit
=======

- cshell script that sets parameters, calls ananew and makes plots
  through idl

runidl
=====

- simple non-interactive call to idl

plotfit.pro
===========

- idl plotting routine

CAT
===

- catalog of all burst sources

burstid.txt
===========

- identifies RXTE burstid's with MJD


