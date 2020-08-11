## List of tools

This is the list of all routines available within the GaNGiaLF framework. Each routine is activated by with the keyword `enable: True` in the configuration file, or directly from the command line (run `gufo -h` for further information to run GaNGialF from the command line).


```
scavengers
|    
└─── Generate default configuration file ()
|    
|  
└─── Stellar subtraction (`starSub`)
|    |    rfi_base_000m.fits
|    |    rfi_base_XXXm.fits
|    |    ... for each time chunk of width XXX  
|    |        % of RFI or flags par baseline and frequency channel    
│
└─── Voronoi Binning (`vorBin`)
|
|
|
└─── Multi-Gaussian Fitter (`gPlayMp`)
|    |    msname_root_rfi_full.fits
|    |    ... table with % of rfi flagged and noise per channel for all, short, long baselines 
|    |
|    |    msname_root_rfi_full_spwbin.fits
|    |    ... table with % of rfi flagged and noise per channel for all, short, long baselines 
|    |        binned over channels of width given by the spwbin parameter    
|    |
|    └───time_chunks
|        |    msname_root_rfi_XXXm.fits
|        |    ... for each time chunk of width XXX  
|        |        table with % of rfi flagged and noise per channel for all, short, long baselines 
|        |
|        |    msname_root_rfi_XXXm_spwbin.fits
|        |    ... for each time chunk of width XXX  
|        |        table with % of rfi flagged and noise for all, short, long baselines
|        |        binned over channels of width given by the spwbin parameter
|        |
|        |    msname_root_flags_XXXm.fits
|        |    ... for each time chunk of width XXX  
|        |        table with % of flags and consequent noise per channel for all, short, long baselines 
|        |
|        |    msname_root_rfi_XXXm_spwbin.fits
|        |    ... for each time chunk of width XXX  
|        |        table with % of flags and consequent noise for all, short, long baselines
|        |        binned over channels of width given by the spwbin parameter
|
└─── Plot multi-gaussian results (`gPlot`)
     │    rfi_base_full.png
     |    ... % of rfi par baseline and channel full observation
     |
     |    noise_full_sl_rfi.png
     │    ... expected noise par channel for all,short,long baselines because of RFI clipped
     |
     |    noise_factor_full_sl_rfi.png
     │    ... factor of noise increase over over channel for all,short,long baselines
     |
     │    flags_base_full.png
     |    ... % of flags par baseline and channel
     |
     |    noise_full_sl_flags.png
     │    ... expected noise par channel for all,short,long baselines because of flags
     |
     |    noise_factor_full_sl_flags.png
     │    ... factor of noise increase over over channel for all,short,long baselines because of flags
     │    
     └─── altaz
     |    │    AltAz_rfiXXXX-XXXXMHz.png
     |    │    ... Altitude vs. Azimuth plot with colorcoded the % of RFI flagged per time chunk 
     |    |        1 image for each spectral window of width spwbin in which the observation has been divided 
     |    |      
     |    │    AltAz_flagsXXXX-XXXXMHz.png
     |    │    ... Altitude vs. Azimuth plot with colorcoded the % of flags identified in each time chunk 
     |    |        1 image for each spectral window of width spwbin in which the observation has been divided
     |    
     └─── movies
     |    |
     |    | 
     |    
     └─── time_chunks
          |
          └─── 1D
          |    |    flags_XXXm_sl_rfi.png
          |    |    ... % of rfi par channel for all,short,long baselines because of RFI clipped
          |    |        for each time chunk of width 0XX minutes.                
          |    |   
          |    |    noise_XXXm_sl_rfi.png
          |    |    ... expected noise par channel for all,short,long baselines because of RFI clipped
          |    |
          |    |    noise_factor_XXXm_sl_rfi.png
          |    |    ... factor of noise increase over over channel for all,short,long baselines
          |    |
          |    |    flags_XXXm_sl_flags.png
          |    |    ... % of rfi par channel for all,short,long baselines because of flags
          |    |        for each time chunk of width 0XX minutes.                
          |    |   
          |    |    noise_XXXm_sl_flags.png
          |    |    ... expected noise par channel for all,short,long baselines because of flags
          |    |
          |    |    noise_factor_XXXm_sl_flags.png
          |    |    ... factor of noise increase over over channel for all,short,long baselines because of flags
          |    |
          └─── 2D
               |    rfi_base_XXXm.png
               |    ... % of rfi par baseline and channel
               |        for each time chunk of width 0XX minutes.
               |
               |    flags_base_XXXm.png
               |    ... % of flags par baseline and channel
               |        for each time chunk of width 0XX minutes.
```
