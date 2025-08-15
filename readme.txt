NOTE: Outside joint collaborative agreements, there is absolutely no support with respect to model download, setup, development, or application.

scepter

scepter.f90 	- source code  

makefile_RE     - compile the source code to run soil hydrology benchmark
makefile_ex11   - compile the source code to run cation exchange benchmark
makefile_AMD    - compile the source code to run acid rock drainage (ARD) benchmark

spinup.py       		- run soil hydrology benchmark
test_phreeqc_ex11.py	- run cation exchange benchmark
spinup_AMD.py			- run ARD benchmark-1
spinup_AMD2.py			- run ARD benchmark-2
spinup_AMD3.py			- run ARD benchmark-3


$$ running benchmark simulations in GMD paper $$

*** correct the "outdir" parameter for all Python scripts

[A] soil hydrology 
*** one need openRE repository installed parallel to SCEPTER directory as it contains input 

(1) $ make --file=makefile_RE clean
(2) $ make --file=makefile_RE 
(3) $ python3 spinup.py


[B] cation exchange 

(1) $ make --file=makefile_ex11 clean
(2) $ make --file=makefile_ex11 
(3) make sure " runtype = 'init' " in test_phreeqc_ex11.py
(4) $ python3 test_phreeqc_ex11.py
(5) make sure " runtype = 'rstrt' " in test_phreeqc_ex11.py
(6) $ python3 test_phreeqc_ex11.py

[C] ARD 

(1) $ make --file=makefile_AMD clean
(2) $ make --file=makefile_AMD 
(3) make sure " runtype = 'init' " in spinup_AMD.py
(4) $ python3 spinup_AMD.py
(5) make sure " runtype = 'rstrt' " in spinup_AMD.py
(6) $ python3 spinup_AMD.py
(7) repeat (3)-(6) with spinup_AMD2.py and spinup_AMD3.py