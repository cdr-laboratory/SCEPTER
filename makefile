# Start of the makefile
# Defining variables

FC            = gfortran
# FC            = ifort

CPFLAGS       = 
# CPFLAGS       += -Dno_intr_findloc # need to use in cluster
# CPFLAGS       += -Dshow_PSDiter # showing iteration process during PSD calculation
# CPFLAGS       += -Dparallel_ON # testing parallelization
# CPFLAGS       += -Dnpar_in=1 # number of threads for parallelization 
# CPFLAGS       += -Dnpar_in=1 # number of threads for parallelization 
# CPFLAGS       += -Dparpsd_chk # checking parallelization results
CPFLAGS       += -Dksld_chk # checking rate consts for sld species
# CPFLAGS       += -Ddebug_newton # printing out more during main newton iteration
# CPFLAGS       += -Ddebug_phcalc # printing out more to check pH calculation
# CPFLAGS       += -Ddebug_season # printing out more to check seasonaility forcing calculation
CPFLAGS       += -Dsat_AMD_banch # use saturation depth profile for benchmarking
CPFLAGS       += -Dredox_TST_ON # use saturation depth profile for benchmarking
CPFLAGS       += -DAMD_benchmark # use parameterization for AMD benchmarking
# CPFLAGS       += -DAMD_benchmark_grid # use parameterization for AMD benchmarking
# CPFLAGS       += -DAMD_boundary # use boundary conditions for AMD benchmarking
# CPFLAGS       += -DAMD_initial # use initial conditions for AMD benchmarking
CPFLAGS       += -DAMD_benchmark_grid # use grid for AMD benchmarking
CPFLAGS       += -Dredox_eq # assume equilibrium for redox reactions
# CPFLAGS       += -Ddisp_cnst=1e-2 # forcing constant and common dispersion for aqueous species
# CPFLAGS       += -Ddisp_cnst=0 # forcing constant and common dispersion for aqueous species
# CPFLAGS       += -DGoldberg_Sikora # using Goldberg et al. 2002 data for Sikora buffer dissociation constants
# CPFLAGS       += -DolddustPSD # using old PSD for dust (not user input but prescribed one)
# CPFLAGS       += -Derrmtx_printout # 
# CPFLAGS       += -Dmtx_printout # 
CPFLAGS       += -Dmod_basalt_cmp # using basalt composition defined in <basalt_define.h>
# CPFLAGS       += -Dlocate_sb # display subroutine names when into and out of the subroutines 
# CPFLAGS       += -Ddef_flx_save_alltime # flux reported each integration (costs lots of bites)
# CPFLAGS       += -Dfull_flux_report # output all cumulative flux
# CPFLAGS       += -Ddisp_lim # limiting the display of results
# CPFLAGS       += -Ddiss_only # not allowing precipitation of minerals
# CPFLAGS       += -Dlim_minsld # limiting mineral lowest conc. 
# CPFLAGS       += -Dporoiter # do iteration for porosity  
# CPFLAGS       += -Dcalcw_full # fully coupled w calcuation  
# CPFLAGS       += -Ddispiter # showing flux/concs in each iteration   
# CPFLAGS       += -DdispPSDiter # showing PSD flux in each iteration   
# CPFLAGS       += -Dcalcporo_full # fully coupled poro calcuation  
# CPFLAGS       += -Diwtypein=0 # uplift type 0--cnst w, 1-- cnst poro*w, 2-- cnst (1-poro)*w, 3--- w-flexible, if not defined 0 is taken
# CPFLAGS       += -Dimixtype_background_in=0 # background mixing type 0--no mix w, 1-- fickian, 2-- homogeneous, 3-- tilling, 4-- LABS, if not defined 0 is taken
CPFLAGS       += -Dnrec_prof_in=200

ifeq ($(FC),gfortran)
  # CFLAGS        = -fcheck=all -g -O3  
  # CFLAGS        = -Wall -O3 -g -fcheck=all -ffpe-trap=invalid,zero,overflow -fbacktrace
  # CFLAGS        = -Wall -O3 -g -fcheck=all -fbacktrace
  CFLAGS        = -fimplicit-none  -Wall  -Wline-truncation  -Wcharacter-truncation  -Wsurprising  \
	  -Waliasing  -Wimplicit-interface  -Wunused-parameter  -fwhole-file  -fcheck=all  -std=gnu  -pedantic  -fbacktrace -O3

endif

ifeq ($(FC),ifort)
  CFLAGS        = -O3 -heap-arrays -g -traceback -check bounds -fp-stack-check -gen-interfaces -warn interfaces -check arg_temp_created 
endif 

# LDFLAGS       = -L/usr/local/lib
LDFLAGS       = 

LIBS          = -lopenblas

ifneq (,$(findstring -Dmod_basalt_cmp,$(CPFLAGS)))
  # Found -Dmod_basalt_cmp
  INC          = -I/storage/coda1/p-creinhard3/0/ykanzaki3/PyWeath/data 
  # INC          = -I/home/ykanz/PyWeath/data 
else
  # Not found
  INC          = 
endif

ifneq (,$(findstring -Dparallel_ON,$(CPFLAGS)))
  # Found -Dparallel_ON
  CFLAGS        += -fopenmp
else
  # Not found
endif

codeid = _H_part
codeid = _H_part_MOD
codeid = _H_part_SAVE_ATTEMPT_IS_CLEAN
# codeid = _H
# codeid = 

# OBJS          = scepter.o
# SRC           = scepter.f90
# OBJS          = scepter_PREV_STRCT.o
# SRC           = scepter_PREV_STRCT.f90
# OBJS          = scepter_PREV_STRCT_DEV.o
# SRC           = scepter_PREV_STRCT_DEV.f90
# OBJS          = scepter_PREV_STRCT_DEV_H.o
# SRC           = scepter_PREV_STRCT_DEV_H.f90
# OBJS          = scepter_PREV_STRCT_DEV_H_part.o
# SRC           = scepter_PREV_STRCT_DEV_H_part.f90
OBJS          = scepter_PREV_STRCT_DEV$(codeid).o
SRC           = scepter_PREV_STRCT_DEV$(codeid).f90
# OBJS          = scepter_PREV_STRCT_DEV_MAQFT.o
# SRC           = scepter_PREV_STRCT_DEV_MAQFT.f90
# OBJS          = scepter_PREV_STRCT_DEV_GHOST.o
# SRC           = scepter_PREV_STRCT_DEV_GHOST.f90
# OBJS          = scepter_AMD.o
# SRC           = scepter_AMD.f90
# OBJS          = scepter_AMD_DEV.o
# SRC           = scepter_AMD_DEV.f90
# OBJS          = sb_ph_calc.o scepter_AMD_test.o 
# SRC           = sb_ph_calc.f90 scepter_AMD_test.f90  # added (po2)^ss_add2 to charge balance calculation
# OBJS          = sb_ph_calc.o scepter_AMD_test_dz_boundary.o 
# OBJS          = sb_ph_calc.o scepter_AMD_test_H.o 
# SRC           = sb_ph_calc.f90 scepter_AMD_test_H.f90  # added (po2)^ss_add2 to charge balance calculation
# OBJS          = sb_ph_calc.o scepter_AMD_test_H_simple.o 
# SRC           = sb_ph_calc.f90 scepter_AMD_test_H_simple.f90  # added (po2)^ss_add2 to charge balance calculation
# SRC           = sb_ph_calc.f90 scepter_AMD_test_dz_boundary.f90  # added (po2)^ss_add2 to charge balance calculation
# OBJS          = sb_ph_calc.o scepter_AMD_test_dz_boundary_DEV.o 
# SRC           = sb_ph_calc.f90 scepter_AMD_test_dz_boundary_DEV.f90  # added (po2)^ss_add2 to charge balance calculation
# OBJS          = sb_ph_calc_DEV.o scepter_AMD_test_dz_boundary_DEVDEV.o 
# SRC           = sb_ph_calc_DEV.f90 scepter_AMD_test_dz_boundary_DEVDEV.f90  # added (po2)^ss_add2 to charge balance calculation
# OBJS          = scepter_v1_SAVE-11-8-2024.o
# SRC           = scepter_v1_SAVE-11-8-2024.f90
# OBJS          = scepter_v1.o
# SRC           = scepter_v1.f90
# OBJS          = scepter_v1_modparvals_fail_SAVE-11-11-2024.o
# SRC           = scepter_v1_modparvals_fail_SAVE-11-11-2024.f90
# OBJS          = scepter_PyWwath_SAVE-11-7-2024.o
# SRC           = scepter_PyWwath_SAVE-11-7-2024.f90
# OBJS          = scepter_PyWwath_SAVE-11-10-2024.o
# SRC           = scepter_PyWwath_SAVE-11-10-2024.f90
# OBJS          = scepter_SAVE-11-6-2024.o
# SRC           = scepter_SAVE-11-6-2024.f90
# OBJS          = scepter_DEVSAVE-9-16-2024.o
# SRC           = scepter_DEVSAVE-9-16-2024.f90
                            
PROGRAM       = scepter_AMD$(codeid)
# PROGRAM       = scepter

all:            $(PROGRAM)

$(PROGRAM):     $(OBJS)
	$(FC) $(OBJS) -o $(PROGRAM) -cpp $(CPFLAGS) $(CFLAGS) $(LIBS) $(LDFLAGS) $(INC)

$(OBJS):        $(SRC) 
	$(FC) $(SRC) -c -cpp $(CPFLAGS) $(CFLAGS) $(LIBS) $(LDFLAGS) $(INC)

clean:;         rm -f *.o  *~ $(PROGRAM)
blank:;         truncate -s 0 *.out
cleanall:;         rm -f *.o *.out *~ $(PROGRAM)

