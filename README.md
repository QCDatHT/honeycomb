# tw3evolution


### Installation and throubleshooting 
0. The windows platform are not supported in the current version.
1. Please first bootstrap the build using 
   ```sh
   ./configure.sh <opt: compiler>
   ```
   where the `<opt>`
   optional argument can be used to select a specific compiler (or compiler version). If the field is left empty, or the requested compiler cannot be found, a list of compilers will be checked until one is found in the system. If none is found, the script will abort. 

2. The codebase relies on `libpthread` and `libomp`.
In general, `libomp` may not be provided by default in your system. A possible solution is to install the LLVM framework, which can be done via your preferred package manager: for macOS more information [here](https://formulae.brew.sh/formula/llvm), for linux systems more information [here](https://apt.llvm.org/). Even after proper installation, it might happen that the linker is not able to find the library. This can happen, for instance, if the library has been installed in a directory not included in your `PATH`, or in a format like `libomp.so.<version>`. In the latter case, if the library is in a directory already in the `PATH`, a symlink `libomp.so` must be created by the user. In the former case, the configuration script can be informed of the path to the library via 
   ```sh
   LIBOMP_PATH="<path-to-libomp.so>" ./configure.sh <opt>
   ```
   The configuration script will then try to compile a test program using the provided path to `libomp.so`. **In case of a failure, please check carefully that the provided path is correct.** If no error is encountered, `compile.sh` and a `Makefile` are generated, which can be used to compile the library and the tests by running  `./compile.sh`
   The compilation will happen on all the available logical cores of the machine.
   If the compilation succeeds, the static library `libtw3ev.a` will be generated in the `$lib` folder.