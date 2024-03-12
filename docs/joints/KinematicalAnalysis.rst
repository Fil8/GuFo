.. GaNGiaLF documentation file, created by
   sphinx-quickstart on Mon Feb 18 15:04:26 2019.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.
 
.. _KinematicalAnalysis:

********************
Kinematical analysis
********************
 
.. toctree::
   :maxdepth: 2
 
This set of tools can be used to analyse some kinematical properties of the observed gas. The analysed line emission must be specified in :code:`Name` which is assumed to be the base-string for all needed columns of the `inputTable`. :code:`TableFromMoms` and the default spectral line fitting routines of GanGIaLF generate the `inputTable` in the right format.

Select-rotation
###############

If a model of the kinematics is available, GanGIaLF can look for all lines of sight of the spectral line cube that match with the model. 

All l.o.s. where the full width at zero intensity of the observed/fitted line matches within :code:`rotationPercent` (:math:`\%`) with the corresponding emission line of the model are stored with 1 in the RotID column of the :code:`Ancels` extension of the *outputTable*. 


Schema of configuration file parameters within :code:`selectRotation`::

 kinematicalAnalysis: 
	Name: 'name of gas in input table'
  	selectRotation:
  		enable: True
		rotationID: True / False
		modelCube: 'path-to-model-cube'
		rotationPercent: 90.

**Note**:
Typically, in the kinematical anaylsys of spectral line observations tilted-ring models the gas are generated to describe the regular rotation of the gas. This function can be used to identify which all lines of sight are or not compatible with regularly rotating kinematics. The `Kinematical-plot`_ can be used to see all-in-one plot the result of this analysis.

**Requirements**:
:code:`[modelCube]`: The model must be given as a datacube at the same **spectral** and **spatial** resolution as :code:`[general][inputCube]`. Pixel size and datacube sizes must also be the same. 	

Kinematical-plot
################

Once the kinematical properties of a spectral line in the field of view of the observations have been measured and stored in the :code:`Ancels` extension of the output table, is possible to draw the **kinematical plot** (**k-plot**): dispersion of the line (i.e. :math:`\sigma`) vs. centroid (i.e. :math:`v_{\rm los}-v_{\rm sys}`, the line shift with respect to the systemic velocity).

Schema of configuration file parameters within :code:`kinematicalAnalysis` ::

 kinematicalAnalysis: 
		Name: 'name of gas in input table'
		computeSigmaCentroid: True / False	
		CCAanalysis: 
			enable: True / False
			sigmaInCCA: 1.5
		kPlot:
			enable: True / False
			theoreticalCCA: True / False
			ellipsisKind: Pencil / Ensemble
			plotRotation: True / False
			legendLabels: ['legend label col1','legend label col2','legend label col3' ]
		kPlotMoment:
			enable: True / False
			plotRotation: True
			CCAOverlay: True
			rotOutCCA: False
			plotElse: False


In the k-plot is easy to understand the kinematical differences between lines of sight that are rotating (i.e. have been found compatible with a rotating model, see `Select-rotation`_) with lines of sight that are compatible with the predictions of CCA (see :ref:`CCA`). The distribution of these points in the f.o.v. can be plotted specifying the option :code:`kPlotMoment`.
 
* estimate parameters of the k-plot
	:code:`computeSigmaCentroid: True` the velocity dispersion and the centroid of the given line are computed and stored in the output table. 
	*Note* This function is useful if these information have not been stored in the table during the fit.

* estimate compatibility with CCA 
	:code:`CCAanalsysis: True` for each point of the k-plot compatibility within **sigmaInCCA**sigma is computed (Fig. 10 of Gaspari et al. 2018). Results are stored in column CCA of the *ancels_BF* extension of the :code:`outputTable`. 

* draw the k-plot
	:code:`[kPlot][enable]: True` plots the k-plot taking the velocity dispersion and centroid from the :code:`ancels_BF` extension of the input table. All succcessfully characterised independent pixels are plotted. 
	If multi-fit component has been run, and columns called in 'dispIntr-Name-g1'/'centroid-Name-g1' 'dispIntr-Name-g2/centroid-Name-g2' are present these are plotted first. Then the k-plot of the total integrated line (dispIntr-Name-BF/centroid-Name-BF) is plotted. Labels of legend are specified in :code:`legendLabels`. 
	*Options*

	1. :code:`theoreticalCCA: True`: plots the 1,2 sigma ellipsis predicted by CCA for the kinematical properties of the condensing gas. When high spectral resolution observations are available, :code:ellipsisKind: 'pencil'` should be used. Otherwise, :code:`ellipsisKind:`.

	2. :code:`plot Rotation: True``: plots in blue the lines that have been found compatible with rotation, i.e. column *RotMod* is present in the extension *ancels_BF* of the input table.

	3. :code:`CCAOverlay: True`: plots the :code:`sigmaInCCA` ellipsis to highlight the lines of sight compatible with CCA.

* distribution in the f.o.v. of the k-plot's loci
	:code:`[kPlotMoment][enable]: True** plots the line distribution color-coded based on the analysis from the k-plot. Plots either all rotating lines of sight or the CCA-compatible lines of sight, and combinations of both.
	*Options*
	
	1. :code:`plot Rotation: True`: plots points compatible with rotation in blue (see above).	
	
	2. :code:`CCAOverlay: True`: plots points that are compatible with CCA in green (see above).
	
	3. :code:`rotOutCCA: True`: plots points that are compatible with rotation but not with CCA in grey.
	
	4. :code:`plotElse: True`: plots points not compatible either with rotation and CCA in grey.	

Multi-region k-plot
###################

GanGIaLF can draw the k-plot of multiple regions simultaneously. A figure showing the location of these regions in the field of view can be generated :code:`[multiRegionMoment][enable]: True` providing an moment map for the header (:code:`inputMoment`. Colors and labels of regions are kept the same as in the k-plot.

Schema of configuration file parameters within :code:`kinematicalAnalysis`::

 kinematicalAnalysis:
 	Name: 'name of gas in input table'
	kPlotMultiReg:
		enable: True
		tableNames:['table1.fits','table2.fits','table3.fits']
		colors: ['colorMap 1','colorMap 2','colorMap 3']
		regionLabels: ['label Table1','label Table2','labele Table3']
		multiRegionMoment:
			enable: True / False
			inputMoment: 'path-to-input-moment-map'
			contOver: True
			contName: 'path-to-continuum-image'
			contLevels: [continuum-contour-values (list of float)]
		
**Requirements** 

	* :code:`sigmaCentroid` option must be run on the total table, then subtables must be saved externally. 
	* Multi-region looks for 'sigma-'+Name and 'centroid+'Name columns in the :code:`ext=Ancels` of :code:`tableNames`.

**Options**

	* :code:`[contOver]:True` contours of a second image (typically a continuum image, but also a moment map) can be overlaid.

**Needed updates**

	* :code:`[kPlotMultiReg][Name]` becomes :code:`[kinematicalAnalysis][Name]`

Cold Chaotic Accretion Aanalysis tools
######################################

GanGIaLF can compute kinematical and physical parameters that can be directly compared to cold chaotic
acccretion criteria. Specifically, the C-ratio (Gaspari et al. 2017) and the Taylor number (Gaspari et al. 2018) which both measure if thermal instabilities can generate a condensation cascade in the gas and consequent infall towards the SMBH. 

Taylor Number
*************

Abba::

	CCAcriteria:
	    taylorNumber:
	      enable: False
	      maxRadius: 5.8 #kpc
	      inTables: ['IonRotAllSmall.fits','HIRotAllSmall.fits','CORotAllSmall.fits']
	      yColumns: ['dispIntr_NII6583','sigma_HI','sigma_CO']
	      gasColors: ['dodgerblue','darkorange','seagreen']
	      gasNames: ['[NII]6583 inside disk','HI inside disk','CO inside disk']
	      vRot: 340.
		  ...

C-Ratio
*******

ABBBa::

	CCAcriteria:    
    	cRatio:
	      enable: False
	      maxRadius: 5.
	      inTables: ['IonRotAllSmall.fits','IonOutAllSmall.fits','IonCCAAllSmall.fits']
	      tEddy: ['teddyAlt','teddyAlt','teddyAlt']
	      rColumns: ['r','r','r']
	      tCool: ['tCool','tCool','tCool']
	      gasColors: ['dodgerblue','firebrick','seagreen']
	      gasNames: ['[NII]6583 in rotation','[NII]6583 outflow','[NII]6583 CCA candidate']
		  ...
