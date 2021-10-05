# List of tools

This is the list of all routines available within the GaNGiaLF framework. Each routine is activated by with the keyword `enable: True` in the configuration file. Some quick usage routines can be called directly from the command line (run `gufo -h` for further information to run GaNGiaLF from the command line). GaNGiaLF can be imported and run using `jupyter notebook`, with `from scavengers import *` (see `Tutorial`_ for more information).


```
GaNGiaLF' scavengers
|    
└─── General
|    |    Store inputs from general block of configuration file
|    |    Create output directory structure
|    |    Optional: Generate default configuration file ()
|
└─── Optical spectral observations -- data analysys    
|    |
|    └─── Stellar subtraction (`starSub`)
|    |    lorem ipsum
|    |
|    └─── Voronoi Binning (`vorBin`)
|    |    lorem ipsum
|    |
|    └─── Multi-Gaussian Fitter (`gPlayMp`)
|    |    lorem ipsum
|    |    
|    └─── Plot multi-gaussian results (`gPlot`)
|    |    lorem ipsum
|
└─── 2D kinematical analysis
|    |    
|    └─── Moment maps (momPlay.momplay())
|    |    lorem ipsum
|    |
|    └─── Plot Moment maps (momPlot.MOMPlot())
|    |    lorem ipsum
|    |
|    └─── Kinematical analysis (ancelsPlot.ancelsplot())
|    |    lorem ipsum
|    
└─── 3D kinematical analysis
|    |    lorem ipsum
|    
└─── Physical parameters diagnostics 
|    |    lorem ipsum
|    |
|    └─── Optical spectral observations    
|    |    |
|    |    └─── Line Ratios (linePlay.lineplay())
|    |    |    lorem ipsum
|    |    |    Dust extinction from Balmer Decrement
|    |    |    Electron density from [SII] doublet
|    |    |         
|    |    └─── Plot BPTs and line ratios diagnostics (linePlay.lineplay())
|    |    |    lorem ipsum
|    |
|    └─── HI observations    
|    |    |    
|    |    └─── HI emission ()
|    |    |    Column density (linePlay.lineplay())
|    |    |    Surface brightness (linePlay.lineplay())
|    |    |    HI mass
|    |    |
|    |    └─── HI absorption ()
|    |    |    lorem ipsum    
|    |
|    └─── CO observations    
|    |    Mass (linePlay.lineplay())
|    |    .... 
|
└─── Datacube specific operations   
|    |    lorem ipsum
|    |    rebin
|    |    regrid
|    |    convolve (you wish)    
|    |
|    └─── header manipulations (headPlay.headplay())
|    |    unit conversions
|
└─── GaNGiaLF output table specific operations   
|    |    lorem ipsum
|    |    rebin
|    |    make table from moment maps      
|
└─── GaNGiaLF generic operations (cvPlay.cvplay())   
     unit conversions
```
