scepter

scepter.f90                     - source code  
makefile                        - compile (data dir needs to be specified)


scripts:
get_int_prof.py                 - contain functions to get output data 
make_inputs.py                  - contain functions to make input data
tunespin_3_newton_inert_buff.py - run 3 variable iterations (output dir needs to be specified)
water_amb.py                    - sample porewater and leave it in the beaker in lab (output dir needs to be specified; no use for pipeline?)
basalt_buff_tunespin_bisec.py   - run basalt application for a target pH etc. (output dir needs to be specified)
sub_jobs.py                     - submit multiple jobs (for spinup and basalt exp)
run_a_shell.sbatch              - submit a single job (for GT cluster)