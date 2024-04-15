# transient-sparam-gnucap
Device plugin for Gnucap which uses vector fitting and recursive convolution to allow pure transient simulations with S-parameter networks which are described by tabulated frequency domain data.

## Description
Initially created by Seán Higginbotham as an M.A.I project at Trinity College Dublin in the Department of Electronic and Electrical Engineering, 2023-24. Supervisor: Dr. Justin King.

## File Descriptions

## Dependancies
Gnucap documentation and installation instructions found here: <https://savannah.gnu.org/projects/gnucap/>.
This project used 'master 2021.01.07', though later versions should work just fine.

LAPACK version 3.12.0 was used for linear algebra computations, see <https://netlib.org/lapack/> for LAPACK homepage. The  compressed ```.tar.gz``` file for LAPACK is given in the ```/rc``` directory. It will need to be extracted and built from source.

## Compilation
Instructions for building Gnucap are included in the downloaded source from its repository; from the top directory simply type
```
./configure
make install
```
Building LAPACK is more involved, but instructions are also given in the downloaded source. It will need to be compiled into ```.so``` dynamic library files for the compilation of ```rc_model.cc``` to work.
For simplicity, a ```make.inc``` file is provided in this repository, which should contain the required changes to compile LAPACK properly. Place this file into the ```.../lapack-3.12.0/``` top directory, and then run the following commands from this same directory.
```
make
cd LAPACKE
make
cd ../CBLAS
make
```
After all the compiling is done, there should be five dynamic library files in the top directory: ```libcblas.so```, ```liblapack.so```, ```libtmglib.so```, ```liblapacke.so```, and ```librefblas.so```.

To compile the Gnucap plugin, run the following in the directory where  ```vf.h``` and  ```rc_model.cc``` are located. It's assumed that LAPACK is installed in  ```~/code/```. Change the linking and include paths as necessary; both the Gnucap and LAPACK headers and  ```.so``` library files are needed.
```
g++ -shared -fPIC rc_model.cc ~/code/lapack-3.12.0/liblapack.so ~/code/lapack-3.12.0/libtmglib.so ~/code/lapack-3.12.0/liblapacke.so ~/code/lapack-3.12.0/librefblas.so ~/code/lapack-3.12.0/libcblas.so -o rc.so -I~/code/lapack-3.12.0/LAPACKE/include -I~/code/lapack-3.12.0/CBLAS/include -I/usr/local/include/gnucap -L~/code/lapack-3.12.0 -L/usr/local/lib/ -l lapack -l lapacke -l gfortran -lcblas -l tmglib -l refblas
```
The final result is a file, ```rc.so```, which is linked into Gnucap by typing ```load rc.so``` inside a Gnucap runtime. Note that ```rc.so``` should ideally be placed somewhere that is searched by the ```LD_LIBRARY_PATH``` environment variable.

## Licenses & Acknowledgements

**Gnucap**

Gnucap is the creation of Albert Davis and is developed by him and others. It is provided under
the GNU GPLv3, which is also the license that this project code is provided under.

See https://www.gnu.org/licenses/gpl-3.0.html. For the GNU GPLv3 license.

Additionally, see the Gnucap repository here https://savannah.gnu.org/projects/gnucap/.

**LAPACK**

LAPACK is a co-creation of The University of Tennessee and The University of Tennessee
Research Foundation, The University of California Berkeley, and The University of Colorado
Denver.

See the user guide here https://netlib.org/lapack/.

The LAPACKE C bindings are the creation of Intel Corp.

The relevant licensing files are found within the source code provided.
