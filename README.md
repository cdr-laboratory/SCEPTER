# PyWeath 
1 D weathering model of pyrite and albite

*pysil.f90*  
Main code for simulation of 1D weathering.    

---- Prerequisite: BLAS library -----

>>> BLAS (OpenBLAS) install 
(Google yourself or)
1. download: http://www.openblas.net/ -> TAR
2. tar zxvf OpenBLAS-0.2.20.tar.gz
3. cd OpenBLAS-0.2.20
4. make BINARY=64 CC="gcc -m64" FC="gfortran -m64"
5. su
6. make PREFIX=/usr/local install

If you get the error 'libopenblas.so.0: cannot open shared object file: No such file or directory' then type
1. sudo apt-get install libopenblas-base
2. export LD_LIBRARY_PATH=/usr/lib/openblas-base/

Compile & run:  
```gfortran -cpp -Dporoevol -Dsurfevol1 pysil.f90 -lopenblas  ```   
```./a```

Input files:  
*frame.in*   
A file to define major parameters including the simulation name which becomes the directory name where simulation results are stored.     
*slds.in*   
A file to specify solid speces to explicitly simulate parameters including the simulation name which becomes the directory name where simulation results are stored.     
*solutes.in*   
A file to specify aqueous species.     
*gases.in*   
A file to specify gaseous species.     
*extrxns.in*   
A file to specify reactions to be included other than mineral reactions.     
*parentrock.in*   
A file to specify parent rock composition (or lower boundary conditions for solid phases).      
*rain.in*   
A file to specify rainwater composition (or upper boundary conditions for aqueous species).      
*atm.in*   
A file to specify atmospheric composition (or upper boundary conditions for gaseous species).      
*dust.in*   
A file to specify dust composition (this specifys the mineral compositions when solid phase is added to profile; the total amount of dust is specified in *frame.in*.  
*OM_rain.in*   
A file to specify OM rain composition (this specifys the fraction of different classes of OM; the total amount of OM rain is specified in *frame.in*.      
*switches.in*   
A file to enable/disable some aspects of simulations, including the type of bio-mixing in the profile.      

Output files:  
Solid (mol/m3), gaseous (atm) and aqueous (mol/L) concentrations are recorded in prof_sld_xxx, prof_aq_xxx and prof_gas_xxx text files.  
xxx denotes number of results recording which is from 000 to 019.  
In case of solid phase, wt% converted results are also stored in file of prof_sld(wt%)_xxx.txt.  
All primary and secondary aqueous phase is recorded in file of charge_balance_xxx.txt where you can check the charge balance of porewater.   
Fluxes of individual species are recorded in files of flx_sld/aq/gas-yyy.txt where yyy is the species of solid/aq/gas phases.   

Concentration profiles can be plotted by a python script basalt_exp.py (Python 2.7).    
When you run the script, you are asked to choose the simulation you want to plot.  
