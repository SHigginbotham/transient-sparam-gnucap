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

...

Compiling ```rc_model.cc``` ...

## Licenses & Acknowledgements
I would like to thank my M.A.I supervisor Dr. Justin King, whose previous work was the basis for this project. He provided invaluable insights and guidance which made the project both possible and an enjoyable experience, instilling curiosity at each discussion.

Relevant academic references are included in the M.A.I dissertation (see Trinity College Dublin Library, <https://www.tcd.ie/library/>).

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
