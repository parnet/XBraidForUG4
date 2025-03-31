# XBraidForUG4Lib

## Installation
use [ughub](https://github.com/UG4/ughub) to install this plugin for [ug4](https://github.com/UG4).
depending on the ughub version it might be necessary to manually initialize submodules so that all used libraries are fetched from github
to do so open the terminal and change the directory to the plugin-directory

```bash
cd $UGRootDirectory/plugins/XBraidForUG4/
git submodule init
git submodule update
```

follow then the ughub-readme to activate and build this plugin
make sure that the option PARALLEL=ON is activated
```bash
cmake -DPARALLEL=ON .
cmake -DXBraidForUG4=ON .
make -j2
```

the external library [XBraid](https://github.com/XBraid/xbraid) which is used within this plugin should be compiled while ug4 compiles.
it could be necessary to set the compiler (for instance by setting CC, CXX, MPICC, MPICXX environment variable) manually because an additional, external make is called during the cmake and this may not inherit the compiler from ug4-toolchains.


## Usage
Examples for the usage of this plugin are given by ParallelConvectionDiffusion or ParallelPoroelasticity

## Notes for Development
When using this plugin the SpaceTimeCommunicator has to be set so that routines of UG4 and those of XBraid are working properly.
By doing this the Process Communication Layer will be overwritten, and therefore it is not possible to use -outproc for every proc, but only for those that are in the same spatial processor group.
Other feature like the profiling (or profiling pcl) might not work as expected.

It is possible to define variables in BraidVectorStruct.h that are respected by the compiler to define the behavior of this plugin.
For example to
* get execution times for some operations (clone, init, step, assemble operator, assemble rhs, solve, residual).
* write out every result as a vector which can be used in matlab
* write out a script in matlab notation
* get timestamps when receiving and sending is executed