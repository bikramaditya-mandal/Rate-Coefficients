<h1 align="center">H<sub>2</sub>O + H<sub>2</sub>O Collisional Rates Coefficients</h1>

## Objective:
This program computes rate coefficients for state-to-state rotational transitions in H<sub>2</sub>O + H<sub>2</sub>O system as a function of Rotational and Kinetic temperatures (T<sub>rot</sub> & T<sub>kin</sub>). The code is primarily developed to be used for astronomical modeling in atmospheres of comets, and other astrophysical environments where H<sub>2</sub>O is discovered. The **temperatures (both T<sub>rot</sub> & T<sub>kin</sub>)** are used as input in the unit of **kelvin** and the rotational state-to-state transitions **rate coefficients (*k*)** is computed as ***cm<sup>3</sup>s<sup>-1</sup>***.

## Installation:
This code is written in Fortran language. The only requirement to compile and use this code is to install a Fortran compiler, such as gfortran or ifort. Here are the steps to compile this code and use for your astronomical modeling.

1. First, you download this project by from the [GitHub website](https://github.com/bikramaditya-mandal/Water_Rate_Coefficients.git). Alternatively, you can clone this project using CLI and the following commands.

```sh
   git clone https://github.com/bikramaditya-mandal/Water_Rate_Coefficients.git
   cd Water_Rate_Coefficients
```

2. Then, you edit the first line of the Makefile to incorporate the appropriate installed compiler by modifying 

```sh
    FC=gfortran
```

and  replace **gfortran** with your own choice of compiler.

3. Then you clean the directory to remove old and unnecessary module and object files using 

```sh
    make clean
```

4. Finally, you use "make" to compile and create the executable.
```sh
    make
```

## Modify temperatures:
There is only one file which needs to be modified to get rate coefficients for the desired Temperatures (both T<sub>rot</sub> & T<sub>kin</sub>). The file name is [**Generate_Rates.f90**](Generate_Rates.f90). Here are the steps to insert desired temperatures.

1. First, you open the file [**Generate_Rates.f90**](Generate_Rates.f90) using an editor. 
2. Then, you scroll down a bit to locate the line

```sh
    Temp_rot = 250.d0
```

and you set the value of the desired **Rotational Temperature** to this variable.

**Note:** You can modify it further to run a loop over ranges of temperatures. For example:

```sh
    do i = 1, 10
        Temp_rot = i * 10.0d0
    end do
```

In this example, the code will compute 10 values of **Rotational Temperature** of 
## Citing this work:
For more details and to cite this work, please refer to:
1. Bikramaditya Mandal et al, 2024, Astronomy & Astrophysics
2. Bikramaditya Mandal and Dmitri Babikov, 2023, Astronomy & Astrophysics, 671, A51.
3. Bikramaditya Mandal and Dmitri Babikov, 2023, Astronomy & Astrophysics, 678, A51.


User should change the values of input variables Temp_rot and Temp_kin to the desired ones only in the file "Generate_Rates.f90".
For example, the current values in the main driver file are set to 250 and 300 K, respectively.
The unit of input temterature is K and the unit of output rate is cm^3/s.

