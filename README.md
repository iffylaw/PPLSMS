PPLSMS
=====
PPLSMS is a software package with multiple scripting language for Pre-processing and Post-Processing of Land Surface Models.

Introduction
----

In order to process different data sets and comprehensive evaluate Land Surface Models (LSMs), PPLSMS integrates different program language into shell interface with modern land surface modeling capabilities, and establishes a set of methods for easy processing of data, and offers model-data fusion technologies, evaluation and visualization for land surface models. The Pre-processing and Post-processing of Land Surface Models using different scripting languages, including R, NCL, NCO, CDO, SHELL, LSF, GraDS GrADS and so on.

![alt tag](https://raw.github.com/iffylaw/PPLSMS/master/Img/Figure1.png "The flow diagram of PPLSMS")

  - LSMs: ORCHIDEE, JULES, CoLM, CLM, Noah LSM and JSBACH... 
  - SHELL script to integrated and control CDO (Climate Data Operators)(Schulzweida and Kornblueh, 2006), NCO (NetCDF Operators)(Zender, 2008), GrADS (Grid Analysis and Display System)(Doty and Kinter III, 1995), R(Gentleman et al., 1997), Platform LSF, and NCL(NCAR Command Language)(UCAR/NCAR/CISL/VETS, 2013)
  - also SHELL script can configure the model and interconnect different programming functions of the other scripting languages. 
  - GrADS, CDO and NCO are mainly used to read and write data for pre-processing and post-processing of LSMs. 
  - R and NCL are mainly used for the visualization and statistical analysis, which NCL is for data analysis and visualization of the surface, and R is for the visualization of the single point, station and surface.

![alt tag](https://raw.github.com/iffylaw/PPLSMS/master/Img/Figure2.png "The architecture of PPLSMS based on scripting languages")

ScreenShots
----
> Some figures of simulations of LSMS
![alt tag](https://raw.github.com/iffylaw/PPLSMS/master/Img/china-4plot-2004-01.gif "China FSH Jan 2004")

Version
----

1.0

Installation
--------------

```sh
chmod 755 *.R||sh||ncl
./*.R Arg1 Arg2 Arg3
./*.sh Arg1 Arg2 Arg3
./*.ncl Arg1 Arg2 Arg3
```

License
----

GNU GENERAL PUBLIC LICENSE


****Cold and Arid Regions Environmental and Engineering
Research Institute, Chinese Academy of Sciences****

[CAREERI, CAS]:http://www.careeri.cas.cn
