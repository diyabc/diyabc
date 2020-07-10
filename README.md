# DIYABC RF V1.0

## Command line help
```plain
$ diyabc -h
USAGE :
-p <directory of header.txt>
-r <number of required data sets in the reftable>
-i <name given to the analysis
-c <coverage percentage for poolseq poisson law
-g <minimum number of particles simulated in a single bunch (default=100)>
-m <multithreaded version of the program
-q to merge all reftable_$j.bin 
-Q to merge all reftableRF_$j.bin 
-s <seed for the random generator (deprecated!!!)>
-t <required number of threads>
-w <computer's number if in a cluster (0 by default) >
          (each computer in the cluster is numbered between 0 and the maximal number of computers in the cluster.)
-l for producing a csv reftable file (reftable.csv) and a scenario file (scenario.txt) used by program rf
-x for translating a reftable.bin in reftable.txt
-y for translating a reftable.bin in reftable.txt writing summary statistics before parameters

-z <path/RNGfilename.bin> write the number of streams of the RNG into path/RNGfilename_cores.txt

-n for INITIALIZATION OF THE PARALLEL RNG'S (with parameters as a string including the following options separated par a semi-colon)
           t:<maximal number of the threads (per computers if cluster, 16 by default)>
           c:<maximal number of computers in the cluster (1 by default)>
           s:<seed of the first RNG (1 by default)>
           f:<forcing creation of new RNG's and overriding the old ones>

-d for ABC PRIOR/SCENARIO CHECKING (idem)
           a:<p for PCA, l for locate observed, pl for both>

-k to SIMULATE DATA FILES

-o to simulate summary statistics from a text file containing scenario and parameter values

-R <activate all stats (if empty) or selected stats families
option -p is compulsory
```

## Prebuilt binaries for Windows/MacOS/Linux

There are avalaible at the [Releases](./releases) page.
Please note that those are only x64 binaries. 

## Build

Prerequisites : 

- a C++17 (gcc-7 and ulterior version for example) compiler
- [OpenMp](https://en.wikipedia.org/wiki/OpenMP) for vastly improved simulator performance on multicore machines
- [CMake](https://cmake.org/)
  
### Linux & MacOS
```sh
git clone --recurse-submodules https://github.com/diyabc/diyabc.git
cd diyabc
mkdir build
cd build
cmake ../
cmake --build . --config Release
```

### Windows
``` Windows
git clone --recurse-submodules https://github.com/diyabc/diyabc.git
cd diyabc
mkdir build
cd build
cmake ..\
cmake --build . --config Release
```
