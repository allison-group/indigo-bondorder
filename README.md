# indigo-bondorder

Bondorder and formal charge determination for molecules.

## Requirements
- CMake version >= 3.5 (An older version will probably work, you will just have to modify the minimum required version in CMakeLists.txt)
- C++14 compilent compiler. Tested with AppleClang >= 9.0, Clang >= 3.9, GCC >= 5.5
- [optional] Boost. The required libraries are provided from boost version 1.66.0
- [optional] Python >= 3.4, for building the python bindings
- [optional] Java runtime, for using the QuickBB method to generate tree decompositions. Not really worthwhile.
- [optional] doxygen for building the (sparse) documentation

## Installation
Tested on macOS and Linux. If you're on Windows, you're on your own. We have no means of testing on a Windows system.

mkdir build && cd build

cmake .. && make

make install

Examples are provided in the example directory.

