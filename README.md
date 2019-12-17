# WARNING
This code is in active development. It is expected that nothing works/compile by now. First release might come in some months.


## Compiling
The following third parties should be already installed:
```
cmake 3.1+
Boost
```


Then, clone the repository:
```
git clone https://gitlab.inria.fr/lishisoa/EYTA
```
and build the package:
```
cd EYTA && mkdir build && cd build && cmake -DABUNDANCE_TYPE=u_int32_t -DSKIP_DISCRETIZATION=1 .. &&  make -j8 && make package
```
Then, you should get a package `EYTA-<version>-Linux.tar.gz`.

## Debugging
For debugging, please compile as:
```
    cd EYTA && mkdir build_debug && cd build_debug && cmake -DABUNDANCE_TYPE=u_int32_t -DSKIP_DISCRETIZATION=1 -DCMAKE_BUILD_TYPE=Debug .. &&  make -j8
```

