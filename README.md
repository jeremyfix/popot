POPulation based Optimization Toolbox

### Installation

Let's call ROOT_DIR the directory where the sources have been extracted so that you should already have the subdirectories ROO
T_DIR/doc , ROOT_DIR/src , ROOT_DIR/examples ...

Compilation/Installation using cmake is done with :

- cd ROOT_DIR
- mkdir Build
- cd Build
- cmake .. -DCMAKE_INSTALL_PREFIX=/usr
- make
- sudo make install

This will build :
- the library libpopot
- the python wrapper libPyPopot
- the pkf-config file popot.pc
- the binaries for the examples popot-example-xxxx
- the documentation

and install everything in $PREFIX_INSTALL/lib, $PREFIX_INSTALL/lib/pkgconfig, $PREFIX_INSTALL/bin, $PREFIX_INSTALL/share/popot as well as the headers in $PREFIX_INSTALL/include/popot. The python wrapper shall be installed in your $PREFIX_INSTALL/python-xxx/site-packages and the python usage example in $PREFIX_INSTALL/bin/popot-xxx.py