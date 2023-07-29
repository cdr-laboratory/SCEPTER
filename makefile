# Start of the makefile
# Defining variables

FC            = gfortran
# FC            = ifort

CPFLAGS       = 
CPFLAGS       += -Dno_intr_findloc # need to use in cluster
# CPFLAGS       += -Dshow_PSDiter # showing iteration process during PSD calculation
# CPFLAGS       += -Dparallel_ON # testing parallelization
# CPFLAGS       += -Dnpar_in=1 # number of threads for parallelization 
# CPFLAGS       += -Dnpar_in=1 # number of threads for parallelization 
# CPFLAGS       += -Dparpsd_chk # checking parallelization results
# CPFLAGS       += -Dksld_chk # checking rate consts for sld species
CPFLAGS       += -Ddebug_phcalc # printing out more to check pH calculation
CPFLAGS       += -Ddebug_season # printing out more to check seasonaility forcing calculation
# CPFLAGS       += -Ddisp_cnst=1e-2 # forcing constant and common dispersion for aqueous species
# CPFLAGS       += -Ddisp_cnst=0 # forcing constant and common dispersion for aqueous species
# CPFLAGS       += -DGoldberg_Sikora # using Goldberg et al. 2002 data for Sikora buffer dissociation constants
# CPFLAGS       += -DolddustPSD # using old PSD for dust (not user input but prescribed one)
# CPFLAGS       += -Derrmtx_printout # 
CPFLAGS       += -Dmod_basalt_cmp # using basalt composition defined in <basalt_define.h>
# CPFLAGS       += -Ddef_flx_save_alltime # flux reported each integration (costs lots of bites)
# CPFLAGS       += -Dfull_flux_report # output all cumulative flux
# CPFLAGS       += -Ddisp_lim # limiting the display of results
# CPFLAGS       += -Ddiss_only # not allowing precipitation of minerals
# CPFLAGS       += -Dlim_minsld # limiting mineral lowest conc. 
# CPFLAGS       += -Dporoiter # do iteration for porosity  
# CPFLAGS       += -Dcalcw_full # fully coupled w calcuation  
# CPFLAGS       += -Ddispiter # showing PSD flux in each iteration   
# CPFLAGS       += -DdispPSDiter # showing PSD flux in each iteration   
# CPFLAGS       += -Dcalcporo_full # fully coupled poro calcuation  
# CPFLAGS       += -Diwtypein=0 # uplift type 0--cnst w, 1-- cnst poro*w, 2-- cnst (1-poro)*w, 3--- w-flexible, if not defined 0 is taken
# CPFLAGS       += -Dimixtype_background_in=0 # background mixing type 0--no mix w, 1-- fickian, 2-- homogeneous, 3-- tilling, 4-- LABS, if not defined 0 is taken

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

OBJS          = scepter.o
SRC           = scepter.f90
                            
# PROGRAM       = scepter_DEV
PROGRAM       = scepter

all:            $(PROGRAM)

$(PROGRAM):     $(OBJS)
	$(FC) $(OBJS) -o $(PROGRAM) -cpp $(CPFLAGS) $(CFLAGS) $(LIBS) $(LDFLAGS) $(INC)

$(OBJS):        $(SRC) 
	$(FC) $(SRC) -c -cpp $(CPFLAGS) $(CFLAGS) $(LIBS) $(LDFLAGS) $(INC)

clean:;         rm -f *.o  *~ $(PROGRAM)
blank:;         truncate -s 0 *.out
cleanall:;         rm -f *.o *.out *~ $(PROGRAM)

