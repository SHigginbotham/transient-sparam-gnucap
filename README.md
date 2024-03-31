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

Gnucap is provided under the GNU GPLv3, which is also the license that this project code is provided under. See ```LICENSE``` file.

**LAPACK**
-----------
The relevant licensing files are found within ```lapack-3.12.0.tar.gz``` as appropriate.

**THE FOLLOWING COPYRIGHTS APPLY TO LAPACK:**

Copyright (c) 1992-2013 The University of Tennessee and The University
                        of Tennessee Research Foundation.  All rights
                        reserved.
                        
Copyright (c) 2000-2013 The University of California Berkeley. All
                        rights reserved.
                        
Copyright (c) 2006-2013 The University of Colorado Denver.  All rights
                        reserved.

**THE FOLLOWING COPYRIGHTS APPLY TO THE LAPACKE C BINDINGS:**
                        
Copyright (c) 2012, Intel Corp. All rights reserved.
