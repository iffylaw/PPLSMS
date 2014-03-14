PPLSMS
=====
PPLSMS is a software package with multiple scripting language for Pre-processing and Post-Processing of Land Surface Models.

Introduction
----
Land surface models applications often require data processing and extensive analysis of their input/output in order to complete the simulation or produce a visual image or animation with statistical analysis. Often this analysis and visualization cannot be done in situ because it requires calculating time-series statistics from state sampled over the entire length of the run or analyzing the relationship between similar time series from previous simulations or observations. 

In order to process different data sets and comprehensive evaluate Land Surface Models (LSMs), PPLSMS integrates different program language into shell interface with modern land surface modeling capabilities, and establishes a set of methods for easy processing of data, and offers model-data fusion technologies, evaluation and visualization for land surface models. The Pre-processing and Post-processing of Land Surface Models using different scripting languages, including R, NCL, NCO, CDO, SHELL, LSF, GraDS GrADS and so on.

> Land surface models (LSMs) are the components of global climate models (GCMs) that simulate land surface processes, such as Bucket, BATS, SiB, BIOME，LSM, CoLM, CLM, ORCHIDEE, JULES, JSBACH, IAP94 and so on. 

![alt tag](https://raw.github.com/iffylaw/PPLSMS/master/Img/Figure1.png "The flow diagram of PPLSMS")

![alt tag](https://raw.github.com/iffylaw/PPLSMS/master/Output/LH_Animation.gif "The gif animation")

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

0.1

Installation
--------------

```sh
dos2unix *.R||sh||ncl
chmod 755 *.R||sh||ncl
./*.R Arg1 Arg2 Arg3
./*.sh Arg1 Arg2 Arg3
./*.ncl Arg1 Arg2 Arg3
```

Discussion
----
Pre-processing and post-processing of LSMs are crucial components of the scientific process in the Earth sciences. To build the software package, we have leveraged several well-engineered script languages and software libraries that provide key capabilities in the area of data processing, statistical analysis and visualization in different data sets. Using different scripts languages for this purpose not only simplifies construction of an integrated data analysis capability but also makes the pre-processing and post-processing operations more convenience with well-established visualization and analysis scripts. This software package supports the standard input and output of LSMs with the NetCDF or text file format, although each models are differences, but it is capable of nine core variables of different data for statistical analysis and visualization. The software package provides high-level functions that hide the details of the processing and analyzing from the end user, potentially allowing the familiar script-based analysis approach to custom and modify according to your requirements. Thus most tools and programming languages need the user having sophisticated analysis, display and task schedule skills in programming to generate products in batch, while the software package allow normal user do most products customization for generation using visualization application conveniently, and the implementation of the system changes the way of data products automation generation. 

License
----

GNU GENERAL PUBLIC LICENSE


> Cold and Arid Regions Environmental and Engineering Research Institute, 
> Chinese Academy of Sciences

[CAREERI, CAS]:http://www.careeri.cas.cn
