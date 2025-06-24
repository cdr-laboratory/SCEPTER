NOTE: Outside joint collaborative agreements, there is absolutely no support with respect to model download, setup, development, or application.

scepter

scepter.f90             - source code  
makefile                - compile (data dir in L52 needs to be specified)


scripts:
get_int_prof.py         - contain functions to get output data 
get_int_prof_time.py    - same above but with specifing the model time  
get_int_prof_time_dep.py
                        - same above but with specifing the model depth 
get_soilpH_time.py      - contain functions to calculate soil pH 
get_soilpH_time_dep.py  - same above but with specifing the model depth
get_inputs.py           - contain functions to retrieve input data
make_inputs.py          - contain functions to make input data
tunespin_3_newton_inert_buff_v2_clean.py  
                        - conduct 3 variable field-run iterations 
                            (output dir needs to be specified)
basalt_buff_tunespin_bisec_v2.py    
                        - conduct field-run iteration to get to a target soil/pw pH 
                            (output dir needs to be specified)
spinup.py               - run a spin-up run
spinup_inert.py         - run series of run with bulk speces varying CECs
spinup_inrt2.py         - run series of run with bulk and OM speces varying CECs etc.
test_phreeqc_ex11.py	- run series of run with an exchanger bulk species 
							initially equilibrated with porewater of 1 mM Na and 0.2 mM K and 1.2 mM NO3
							and replaced by porewater with a different composition with 1.2 mM CaCl2
							through advection and dispersion characterized with the Peclet number 40 
							as in Appelo (1994)
test_phreeqc_ex11_init.py
						- run series of run with an exchanger bulk species 
							with 1 mM Na and 0.2 mM K and 1-15 mM NO3



$$ running example simulations in GMD paper $$

1. test for Sikora buffer used in in Section 3
1.1. in-silico field samples 
    a. modify outdir in L20 of spinup_inert.py 
    b. type: python3 spinup_inert.py
1.1. buffer pH for field samples in silico
    a. modify outdir in L528 and 
        undo comments-out in L849-852 of get_soilpH_time.py
    b. type: python3 get_soilpH_time.py
    
2. cation exchange experiment in Section 3 
2.1. cation exchange equilibrium simulation
    a. modify outdir in L176 of test_phreeqc_ex11_init.py
    b. type: python3 test_phreeqc_ex11_init.py
2.2. cation exchange + advection + dispersion simulation
	2.2.1. spin-up
		a. modify outdir in L176 of test_phreeqc_ex11.py
		b. type: python3 test_phreeqc_ex11.py
	2.1.2. dynamic simulation 
		a. undo comment-out in L16,17 of makefile
		b. type: make 
		c. undo comment-out in L277-280 of test_phreeqc_ex11.py
		d. type: python3 test_phreeqc_ex11.py
	
3. mesocosm experiment in Section 3 
3.1. field simulation
    a. modify outdir in L159 of spinup_inrt2.py 
    b. type: python3 spinup_inrt2.py

3.2. laboratory simulation
    a. modify outdir in L528, L888 of get_soilpH_time.py 
    b. type: python3 get_soilpH_time.py
    
4. alkalinity requiremnt for ERW in Section 4 
4.1. field + laboratory simulation
4.1.1. spin/tune-up 
    a. modify outdir in L148 of tunespin_3_newton_inert_buff_v2_clean.py 
    b. type: python3 tunespin_3_newton_inert_buff_v2_clean.py spinup_run_name 
        21.103289732688683 6.058006742238197 20.980309042371502 2.0516666666666667 
        8.222189843654622 0.282726679550165 0.35136107875550837 0.0010131311683626316 
        1.005952418781816
    [ note that the runtime inputs are to specify: run ID, CEC (cmol/kg), target soil pH, 
        target exchange acidity (%CEC), target soil OM (wt%), temperature (oC),
        moisture, runoff (m/yr), erosion rate (m/yr), nitrification rate (gN/m2/yr) ]
4.1.2. basalt application 
    a. modify outdir in L93 of basalt_buff_tunespin_bisec_v2.py and
        option of using soil/porewater pH (phnorm_pw=False/True, L25-26)
    b. type: python3 basalt_buff_tunespin_bisec_v2.py 6.2 1 21.103289732688683 
        basalt_run_name spinup_run_name
    [ note that the runtime inputs are to specify: 
        target pH, duration of run, CEC (cmol/kg), basalt run name, spin/tune-up run name ]
4.2. calculating soil pH prodiles
    a. modify outdir in L49, L299 of get_soilpH_time_dep.py
    b. type: get_soilpH_time_dep.py basalt_run_name
