## GaNGiaLF

Gas N-Gaussians Line Fitter is a framework of tools for the analysis of spectral observations in the radio and optical wavelengths

### Installation

- Clone this repository 

`git clone https://github.com/Fil8/GaNGiaLF.git`

- start and activate virtual-environment (or conda environment) outside of GaNGiaLF directory

```
virtualenv Genv
source Genv/bin/activate
pip install pip wheel setuptools -U
```

- install GaNGiaLF

```
pip install <absolute path to GaNGiaLF folder>
export PYTHONPATH='' # Ensure that you use venv Python
```

### Help
``` 
	*************		---   GuFo   ---	*************

	 +---------+		    Help		+---------+

usage: gufo [-h] [-v] [-gd GENERATE_DEFAULT] [-stS] [-vor] [-gp] [-gMp] [-gPl]
            [-gPs] [-mm] [-mPl] [-lr] [-res] [-bptPl] [-bptMp] [-cD] [-af]
            [-clT] [-bin BINID] [-dF] [-c CFGFILE]

gufo: tools to fit gaussian lines in cubes

version 1.0.0

install path /Users/maccagni/programs/gufo/bin

Filippo Maccagni <filippo.maccagni@gmial.com>

optional arguments:
  -h, --help            Print help message and exit
  -v, --version         show program's version number and exit
  -gd GENERATE_DEFAULT, --generate_default GENERATE_DEFAULT
                        Generate a default configuration file (YAML format)
  -stS, --starSub       make stat subtracted datacube from PPXF output
  -vor, --vorPlay       voronoi binning of line subtracted datacube
  -gp, --gPlay          fit gaussian components for each Voronoi bin
  -gMp, --gMpPlay       multi-process fit gaussian components for each Voronoi
                        bin
  -gPl, --gPlot         loop to plot fit results for each Voronoi bin
  -gPs, --gPlotSingle   plot fit results from a single bin
  -mm, --moments        compute moment maps from fit results of each line
  -mPl, --momPlot       plot moment maps from fit results of each line
  -lr, --lineRatios     compute line ratios of OIII/Hbeta, NII/Halpha,
                        SII/Halpha
  -res, --resPlot       compute and plot residuals of fits
  -bptPl, --bptPlots    compute and draw BPT plots
  -bptMp, --bptMaps     tool to draw BPT maps
  -cD, --cDist          compute and plot eta-parameter
  -af, --ancillaryInfo  compute sigma, w80 and centroid of fitted line
  -clT, --clTbl         clean fit table from everything except fit results
  -bin BINID, --binID BINID
                        bin to plot
  -dF, --doFit          fit single bin
  -c CFGFILE, --cfgFile CFGFILE
                        input .fits file

Run one of the following commands:

        gufo			-gd <configFileName>:	generate config file in working directory
        gufo			-c  <configFileName>:	run GuFo with modules enabled from config file
        gufo -stS		-c  <configFileName>:	make line cube from PPXF outputs
        gufo -gp		-c  <configFileName>:	run fit with parameters from config file and lineList file (outdated)
        gufo -gpMp		-c  <configFileName>:	run fit on multiprocess with parameters from config file and lineList file
        gufo -gPl		-c  <configFileName>:	plot fit results for each bin
        gufo -gPs -bin <binNum>	-c  <configFileName>:	plot fit results for binNum
        gufo -gPs -bin <binNum>	-c  <configFileName> -dF :	do fit for binNum with specs from configFile
        gufo -mom		-c  <configFileName>:	make moment maps of fitted lines
        gufo -momPl		-c  <configFileName>:	 plot moment maps from fits files
        gufo -res		-c  <configFileName>:	make residual maps for each fitted line
        gufo -lr		-c  <configFileName>:	estimate lineRatios, make lineRatio maps and BPT plots of fitted lines
        gufo -clT		-c  <configFileName>:	clean table leaving only fit results
        gufo -vor		-c  <configFileName>:	voronoi binning of output of PPXF (star subtracted datacube) 
        gufo -stS		-c  <configFileName>:	make stellar subtracted datacube (output of PPXF is a table)
        gufo -af		-c   <configFileName>:	compute sigma, centroid, w80 for each voronoi bin
            

	*************		--- GuFo End ---	*************
```
