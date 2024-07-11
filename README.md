<h1 align="center">H<sub>2</sub>O + H<sub>2</sub>O Collisional Rates Coefficients</h1>

## Objective:
This program computes rate coefficients for state-to-state rotational transitions in H<sub>2</sub>O + H<sub>2</sub>O system as a function of Rotational and Kinetic temperatures (T<sub>rot</sub> & T<sub>kin</sub>). The code is primarily developed to be used for astronomical modeling in atmospheres of comets, and other astrophysical environments where H<sub>2</sub>O is discovered.

## Installation:
This code is written in Fortran language. The only requirement to compile and use this code is to install a Fortran compiler, such as gfortran or ifort. Here are the steps to compile this code and use for your astronomical modeling.
1. First, you download this project by from the [GitHub website](https://github.com/bikramaditya-mandal/Water_Rate_Coefficients.git). Alternatively, you can clone this project using CLI and the following commands.

```sh
   git clone https://github.com/bikramaditya-mandal/Water_Rate_Coefficients.git
   cd Water_Rate_Coefficients```

2. Then, you edit the first line of the Makefile to incorporate the appropriate installed compiler by modifying "FC=gfortran". You need to replace "gfortran" with your own choice of compiler.
3. Then you clean the directory to remove old and unnecessary module and object files using "make clean"
4. Finally, you use "make" to compile and create the executable.

## Citing this work:
For more details and to cite this work, please refer to:
1. Bikramaditya Mandal et al, 2024, Astronomy & Astrophysics
2. Bikramaditya Mandal and Dmitri Babikov, 2023, Astronomy & Astrophysics, 671, A51.
3. Bikramaditya Mandal and Dmitri Babikov, 2023, Astronomy & Astrophysics, 678, A51.


User should change the values of input variables Temp_rot and Temp_kin to the desired ones only in the file "Generate_Rates.f90".
For example, the current values in the main driver file are set to 250 and 300 K, respectively.
The unit of input temterature is K and the unit of output rate is cm^3/s.

