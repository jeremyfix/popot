POPulation based Optimization Toolbox

YOTO

### Installation

Let's call ROOT_DIR the directory where the sources have been extracted so that you should already have the subdirectories ROO
T_DIR/doc , ROOT_DIR/src , ROOT_DIR/examples ...

Compilation/Installation using cmake is done with :

- cd ROOT_DIR
- mkdir build
- cd build
- cmake .. -DCMAKE_INSTALL_PREFIX=/usr
- make
- sudo make install

This will build :
- the library libpopot
- the python wrapper libPyPopot
- the pkg-config file popot.pc
- the binaries for the examples popot-example-xxxx
- the documentation

and install :
- the library in prefix/lib, 
- a pkg-config file popot.pc in prefix/lib/pkgconfig, 
- the examples in prefix/bin,
- the documentation in prefix/share/popot 
- the headers in prefix/include/popot. 
- the python wrapper in prefix/python-xxx/site-packages 
- the python usage example in prefix/bin/popot-xxx.py

### Testing

The examples have been installed in prefix/bin and are named popot-example-xxx