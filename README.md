## Radiobs

tools for radioastronomy data analysis

### Installation

Add `<path to radiobs>` to `PYTHONPATH`
  
  - in your `.cshrc`
  `setenv PYTHONPATH <path_to_radiobs>:${PYTHONPATH}`
  - in your `.bashrc`
  `export PYTHONPATH=<path_to_radiobs>:$PYTHONPATH`
  
### Usage
 In `python` environment type `import radiobs`
 
 find which classes are available.
 
 `radiobs.__dict__.keys()`
 
 We choose to determine the expected frequency of the HI line at a given redshift, we need  `hiline` contained in the class `hi`. Let's load the class:
 
 `hi = radiobs.hi()`
 
  use the module:
 `freq = hi.hiline(z=0.1)`

 
