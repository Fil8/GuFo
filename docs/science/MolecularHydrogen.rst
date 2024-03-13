.. GaNGiaLF documentation file, created by
   sphinx-quickstart on Mon Feb 18 15:04:26 2019.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.
 
.. _Molecular Hydrogen:

****************
Molecular Hydrogen
****************
 
.. toctree::
   :maxdepth: 2
 


######################
CO (1-0) emission 
######################


++++++++++++++++++++
Flux & Flux Density
++++++++++++++++++++


--------------------------------------------------




+++++++++++++++++++++++++++++++++++
CO-to-H2 conversion
+++++++++++++++++++++++++++++++++++


--------------------------------------------------




+++++++++++++++++++++++++++++++++++
H2 Surface Brightness
+++++++++++++++++++++++++++++++++++


------------------------------------------------


-------------------------------------------------

:code:`hiPlay.nhiMap` converts an input **flux density map** into an HI-column density map. The input is in units of **Jy/beam km/s**, or in units of   **Jy/beam m/s**. The spectral unit also needs to be specified by the user, along with the redshift and the path of the map. Beam size is read from the header of the input.

:code:`hiPlay.nhiMap` converts an input **flux cube** into an HI-column density units (:math:`cm^{-2}`). The input datacube has units of **Jy/beam**. The channel width is read from the header, while its units (either m/s or km/s) need to be specified by the user.

Surface Brightness limits and maps
----------------------------------
:code:`hiPlay.surfBrightHI` converts a column density into HI surface brightness values. While the module :code:`hiPlay.surfBrightHIMap` can be used to convert a column density map into a surface brightness map.


++++
Mass
++++



Mass detection limit
--------------------


- :code:`hiPlay.hiMass` computes the total HI mass from a flux density, channel widht and beam size specified by the user.

- :code:`hiPlay.hiMassFromMom0` computes the total HI mass seen in a moment-0 map enclosed within a threshold cutoff (the input map is in Jy/beam km/s). Redshift and luminosity distance are specified by the user.

- :code:`hiPlay.hiMassFromSofiaTable` computes the total HI mass from the output table of the :code:`SoFiA` source finder. Columns **z**, **dL** and **f_sum** (total flux) are used. The spectral unit (either m/s or km/s) needs to be specified by the user.

**References & links**


################################
CO higher excitation states 
################################



+++++++++++++++
CO line ratios
+++++++++++++++



--------------------

+++++++++++++++
Optical depth 
+++++++++++++++

--------------------


++++
Mass
++++


--------------------
		
**Example**:

**Further information**:
:code:`[modelCube]`: The model must be given as a datacube at the same **spectral** and **spatial** resolution as :code:`[general][inputCube]`. Pixel size and datacube sizes must also be the same. 

