# OpenMM SDM Plugin

A plugin to implement the Single-Decoupling alchemical protocol in OpenMM

## Contributors

Emilio Gallicchio <egallicchio@brooklyn.cuny.edu>  
Rajat Pal <rajatfor2014@gmail.com>  
Baofeng Zhang <BZhang@brooklyn.cuny.edu>  


## License

This software is released under the LGPL license. See LICENSE.

## Credits

This software is maintained by the Gallicchio's laboratory at Department of Chemistry of Brooklyn College of CUNY. Development and maintenance of this software is supported in part from a grant from the National Science Foundation ([CAREER 1750511](https://www.nsf.gov/awardsearch/showAward?AWD_ID=1750511&HistoricalAwards=false)).

The plugin interface is based on the [openmmexampleplugin](https://github.com/peastman/openmmexampleplugin) by Peter Eastman.

## Requirements

OpenMM 7. Known to work with the OpenMM 7.2.2 and OpenMM 7.3.1

Platforms supported: Reference and OpenCL (CUDA will be supported soon).

Our workflow to set up SDM calculations with this plugin [openmm_sdm_workflow](https://github.com/egallicc/openmm_sdm_workflow) assumes input from structure and parameter files in desmond file format. This requires the latest [desmond file reader](https://github.com/egallicc/openmm/blob/master/wrappers/python/simtk/openmm/app/desmonddmsfile.py). Copy the [desmond file reader](https://github.com/egallicc/openmm/blob/master/wrappers/python/simtk/openmm/app/desmonddmsfile.py) in your OpenMM source tree and rebuild it.

## Installation

Locate the OpenMM installation directory, otherwise it will default to `/usr/local/openmm`. See above regarding patching OpenMM's desmond file reader.

Download the package from github:

```
git clone https://github.com/rajatkrpal/openmm_sdm_plugin.git
```


Build and install the plugin with cmake. Assuming a unix system:

```
mkdir build_openmm_sdm_plugin
cd build_openmm_sdm_plugin
ccmake -i ../openmm_sdm_plugin
```

Hit `c` (configure) until all variables are correctly set, then `g` to generate the makefiles. `OPENMM_DIR` should point to an existing OpenMM installation. `CMAKE_INSTALL_PREFIX` normally is the same as `OPENMM_DIR`. The SDM plugin requires the python API. You need `python` and `swig` to install it.

Once the configuration is done do:

```
make
make install
make PythonInstall
```

The last two steps may need superuser access depending on the installation target. It is recommended to to build the plugin under a `conda` environment to install the python modules without superuser access.

## Test


```

cd example
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:<openmm_dir>/lib:<openmm_dir>/lib/plugins
python test.py

`<openmm_dir>` is the OpenMM installation directory. Again, the last step is best accomplished under the same `virtualenv` environment used to build the python modules.

