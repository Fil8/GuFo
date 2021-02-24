.. GaNGiaLF documentation file, created by
   sphinx-quickstart on Mon Feb 18 15:04:26 2019.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.
 
.. _JupyterTutorial:

******************************
GaNGiaLF with jupyter notebook
******************************

.. toctree::
   :maxdepth: 2
 
User settings are given to GaNGiaLF through a configuration file consisting of a
sequence of blocks, each corresponding to the run of a GaNGiaLF joint, i.e. a set of
function to generate a specific astrophysical observational analysis. 
The `Manual`_ lists the available GanGiaLF joints.

A default configuration file can be generated in the working directory running ``gufo -gd``.

The ``general`` block of parameters must always appear in the configuration file. It stores the paths
and names of inputs and outputs.

The parameters of each block are arranged in a nested structure following the YAML syntax rules 
(see https://yaml.readthedocs.io). As an example, a block of the config file may look like::

  joint_name:
    enable: true
    parameter_1: value_1
    parameter_2:
      parameter_2_1: value_2_1
      parameter_2_2: value_2_2
    parameter_3: value_3
    ...

* Each joint is activated setting ``enable: True``. 
* More that one joint can be enabled in a single run. In this 
  case GaNGiaLF will follow the logical order stored in ``gufo.py`` (see `AutomatedGaNGiaLF`_ for more info on automated routines). 

**jupyter notebook usage**
Parameters of a single GaNGiaLF function can  be specified through the configuration file calling 
``gufo.util.loadCfg()``.
