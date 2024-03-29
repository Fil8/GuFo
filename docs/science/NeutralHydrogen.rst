.. GaNGiaLF documentation file, created by
   sphinx-quickstart on Mon Feb 18 15:04:26 2019.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

.. _Neutral Hydrogen:

****************
Neutral Hydrogen
****************

.. toctree::
   :maxdepth: 2

   hi_emission

###########
HI emission
###########

+++++++++++++++++++
Flux & Flux Density
+++++++++++++++++++

The observed flux density

For a source at redshift z, with a rest frame total HI luminosity of L, its observed **flux S** will be:

.. math::

    S_\nu  & = \int_{\Omega_{\rm src}}I_{\nu}(\theta,\phi)d\Omega \\
           & = \frac{1}{\Omega_{\bm}}\int_{\Omega_{\rm src}}I_{\nu}(\theta,\phi)\times P_n(\theta,\phi) d\Omega \\
           & = \frac{1}{\Omega_{\bm}}\int_{\Omega_{\rm src}}I^{\rm obs}_{\nu}(\theta,\phi)d\Omega

To recover the correct flux for an extended region, :math:`\int_{\Omega_{\rm src}I^{obs}_\nu(\theta,\phi)d\Omega)` is measured by summing values of the pixels in the source region (in Jy per beam) multiplied by the pixel area, from which :math:`S_\nu` can then be recovered by dividing by the area of the telescope beam.

:math:`S_\nu` can be measured by simply converting the specific intensity values of the image to Jy per pixel, and then summing over the region of interest.

Determine total flux of a source from moment-0 map
--------------------------------------------------

:code:`[hiPlay.totalFlux]`: measures the total flux in a moment-0 map within a cutoff [Jy beam:math:`{-1}` km s:math:`{-1}`] given by the user.
The beam and pixel sizes (expressed in degrees) are read from the header of the moment-0 map.

.. math::

    S_\nu & = \frac{1}{\Omega_{\bm}}\int_{\Omega_{\rm src}}I^{\rm obs}_{\nu}(\theta,\phi)d\Omega \\
           & = \sum^{\rm pix}_{S{\rm pix}>{\rm cutoff}}{I^{\rm obs}_{\nu}(pix)} \times {\frac{1}{\Omega_{\bm (pixel)}}

which in terms of the beam angular major and minor axis, and pixel size `pixSize` is:

.. math::

    S_\nu & = \frac{1}{\Omega_{\bm}}\int_{\Omega_{\rm src}I^{\rm obs}_{\nu}(\theta,\phi)d\Omega} \\
           & = \sum^{\rm pix}_{S{\rm pix}>{\rm cutoff}}{I^{\rm obs}_{\nu}(pix)} \times {\frac{1}{2 \pi (ab)/(2.35482^2)}}\times {\rm pixSize}^2

+++++++++++++++++++++++++++++++++++
Column density & Surface Brightness
+++++++++++++++++++++++++++++++++++

The HI column density gives the number of atoms per unit area along the line of sight through an astronomical object. The column density :math:`N_{\rm HI}` for a flux `S` measured over a solid angel is given by:

.. math::

    \frac{N_{\rm HI}}{cm^{-2}} = 2.64\times 10^{20}(1+z)^4\big(\frac{S}{\rm JyHz}\big)\big(\frac{\Omega_{\rm bm}}{\rm arcsec}^2\big)^{-1}

which in terms of the beam angular major and minor axis :math:`a,b` is:


.. math::

    \frac{N_{\rm HI}}{cm^{-2}}  = 2.33\times 10^{20}(1+z)^4\big(\frac{S}{\rm JyHz}\big)\big(\frac{ab}{\rm arcsec}^2\big)^{-1}

Hence, the **column density** in terms of observed **frame velocity integrated flux** is given by:

.. math::

    \frac{N_{\rm HI}}{cm^{-2}}  = 1.10\times 10^{24}(1+z)^2\big(\frac{S}{\rm Jy\,km s^{-1}}\big)\big(\frac{ab}{\rm arcsec}^2\big)^{-1}

Knowing the mass of a single hydrogen atom, :math:`m_{\rm H} = 1.6735575 \times 10^-24 g`, converting it to solar masses along with transforming :math:`cm^{-2}` into :math:`pc^{-2}`, one determines the the expression for the HI surface brightness as:

 .. math::

    \frac{\Sigma_{\rm HI}}{M_\odot \,\,pc^{-2}}  &= 9.01\times 10^{-21}\big\frac{N_{\rm HI}}{cm^{-2}}\big) \\
                                                              &= 1.00\times 10^4(1+z)^2\times\big(\frac{S^{\rm V_{obs}}{\rm Jy \,\, km\, s^{-1}}\big)\big(
                                                              \frac{\Omega_{\rm bm}}{\rm arcsec})^{-1}

Column density detection limit of an observation
------------------------------------------------

Given the r.m.s noise measured in a datacube the column density detection limit for the HI emission is computed assuming a flux threshold and the minimum number of continous channels within which this flux is measured. Knowing the beam size (:math:`ab`) and channel width (:math:`\Delta\,\,v`) the :math:`3\sigma` column density detection limit over one channel is given by: 

 .. math::

    \frac{3\sigma N_{\rm HI}}{cm^{-2}}  = 2.33\times 10^{24}(1+z)^2\big(\frac{3\sigma {\rm rms}}{\rm Jy}\cdot\frac{\Delta\,\,v}{km s^{-1}}\big)\big(\frac{ab}{\rm arcsec}^2\big)^{-1}

The column density detection limit can be estimated using :code:`hiPlay.nhi`. The sigma threshold and noise level (in Jy) are given as input, as well as the channel-width (in km/s) and the beam size (in :code:`astropy.units` arcseconds).

Convert from flux density to column density units
-------------------------------------------------

:code:`hiPlay.nhiMap` converts an input **flux density map** into an HI-column density map. The input is in units of **Jy/beam km/s**, or in units of   **Jy/beam m/s**. The spectral unit also needs to be specified by the user, along with the redshift and the path of the map. Beam size is read from the header of the input.

:code:`hiPlay.nhiMap` converts an input **flux cube** into an HI-column density units (:math:`cm^{-2}`). The input datacube has units of **Jy/beam**. The channel width is read from the header, while its units (either m/s or km/s) need to be specified by the user.

Surface Brightness limits and maps
----------------------------------
:code:`hiPlay.surfBrightHI` converts a column density into HI surface brightness values. While the module :code:`hiPlay.surfBrightHIMap` can be used to convert a column density map into a surface brightness map.


++++
Mass
++++
Taking 3/4 of HI atoms to be in the upper hyperfine state, with a spontaneous emission rate of $A_{\rm HI}$, an emitted photon energy of $h\nu_{rm HI}}$, and an HI source with luminosity :math:`L` to be optically thin, the number of HI atoms, :math:`\mathcal{N}_{\rm HI}` is given by:

.. math::

    \mathcal{N}_{\rm HI} & =\frac{L}{\frac{3}{4}h\nu_{\rm HI}A_{\rm HI}} \\
                          & =\frac{16 D_L^2 S}{3 h\nu_{HI}A_{HI}}

hence:

.. math ::

    \big(\mathcal{N}_{\rm HI}{h^-2})=5.91\times10^{58}\big(\frac{D_L}{h^{-1} Mpc}\big)^2\big(\frac{S}{\rm Jy\,Hz}\big)

where :math:`h` is the Hubble Constant (:math:`\frac{H_0}{100 km\,s^{-1}\,Mpc^{-1}}`).

The total HI mass is given by: :math:`M_{\rm HI}=\mathcal{N}_{\rm HI}m_H`, therefore:

\big({M_{\rm HI}}{h^{-2}M_\odot}) &=49.7\big(\frac{D_L}{h^{-1}{\rm Mpc}}\big)^2\big(\frac{S}{\rm Jy\,Hz}\big) \\
&=\frac{2.35\times1-^5}{(1+z)^2}\big(\frac{D_L}{h^{-1}{\rm Mpc}}\big)^2\big(\frac{S^{\rm V_{\rm obs}}}{\rm Jy\,km/s}\big)


Mass detection limit
--------------------
Given the r.m.s noise measured in a datacube the HI mass emission detection limit is computed assuming a flux threshold and the minimum number of continous channels within which this flux is measured. Knowing the beam size (:math:`ab`) and channel width (:math:`\Delta\,\,v`) the :math:`3\sigma` mass detection in one channel is given by: 

\big({M_{\rm HI}}{h^{-2}M_\odot}) = \frac{2.35\times1-^5}{(1+z)^2}\big(\frac{D_L}{h^{-1}{\rm Mpc}}\big)^2 \big(\frac{3\sigma {\rm rms}}{\rm Jy}\cdot\frac{\Delta\,\,v}{km s^{-1}}\big).


- :code:`hiPlay.hiMass` computes the total HI mass from a flux density, channel width and beam size specified by the user.

- :code:`hiPlay.hiMassFromMom0` computes the total HI mass seen in a moment-0 map enclosed within a threshold cutoff (the input map is in Jy/beam km/s). Redshift and luminosity distance are specified by the user.

- :code:`hiPlay.hiMassFromSofiaTable` computes the total HI mass from the output table of the :code:`SoFiA` source finder. Columns **z**, **dL** and **f_sum** (total flux) are used. The spectral unit (either m/s or km/s) needs to be specified by the user.

**References & links**
This section on HI emission science faithfully follows the information shown in:

`Meyer, M. et al. 2017 <https://www.cambridge.org/core/journals/publications-of-the-astronomical-society-of-australia/article/tracing-hi-beyond-the-local-universe/AFDDDD5C391967752B3229E65559F824#>`_

`HI fidelity calculator <https://hifi.icrar.org>`_

Absorption
##########



++++
Flux
++++

++++++++++++++++++++++++++++++
Optical depth & column density
++++++++++++++++++++++++++++++

++++
Mass
++++



**Example**:

**Further information**:
:code:`[modelCube]`: The model must be given as a datacube at the same **spectral** and **spatial** resolution as :code:`[general][inputCube]`. Pixel size and datacube sizes must also be the same. 
