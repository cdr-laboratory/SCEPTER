# Start of the makefile
# Defining variables

FC            = gfortran

CPFLAGS       = 
CPFLAGS       += -Dno_intr_findloc # need to use in cluster
# CPFLAGS       += -Dshow_PSDiter # showing iteration process during PSD calculation
# CPFLAGS       += -Dksld_chk # checking rate consts for sld species
# CPFLAGS       += -Ddef_flx_save_alltime # flux reported each integration (costs lots of bites)
# CPFLAGS       += -Dfull_flux_report # output all cumulative flux
# CPFLAGS       += -Ddisp_lim # limiting the display of results
# CPFLAGS       += -Ddiss_only # not allowing precipitation of minerals
# CPFLAGS       += -Dlim_minsld # limiting mineral lowest conc. 
# CPFLAGS       += -Dporoiter # do iteration for porosity  
# CPFLAGS       += -Dcalcw_full # fully coupled w calcuation  
# CPFLAGS       += -DdispPSDiter # showing PSD flux in each iteration   
# CPFLAGS       += -Dcalcporo_full # fully coupled poro calcuation  
# CPFLAGS       += -Diwtypein=0 # uplift type 0--cnst w, 1-- cnst poro*w, 2-- cnst (1-poro)*w, 3--- w-flexible, if not defined 0 is taken

# CFLAGS        = -fcheck=all -g -O3  
CFLAGS        = -fimplicit-none  -Wall  -Wline-truncation  -Wcharacter-truncation  -Wsurprising  \
	-Waliasing  -Wimplicit-interface  -Wunused-parameter  -fwhole-file  -fcheck=all  -std=gnu  -pedantic  -fbacktrace -O3

# LDFLAGS       = -L/usr/local/lib
LDFLAGS       = 

LIBS          = -lopenblas
# OBJS          = pysil_dev.o 
# SRC           = pysil_dev.f90 
# OBJS          = pysil_sent.o 
# SRC           = pysil_sent.f90 
# OBJS          = pysil.o 
# SRC           = pysil.f90  
OBJS          = pysil_dev_dev.o 
SRC           = pysil_dev_dev.f90 
# OBJS          = pysil_dev_dev_SSV_NOT_WORKING.o 
# SRC           = pysil_dev_dev_SSV_NOT_WORKING.f90 
PROGRAM       = weathering

all:            $(PROGRAM)

$(PROGRAM):     $(SRC)
	$(FC) $(SRC) -o $(PROGRAM) -cpp $(CPFLAGS) $(CFLAGS) $(LIBS) $(LDFLAGS)


clean:;         rm -f *.o  *~ $(PROGRAM)
blank:;         truncate -s 0 *.out
cleanall:;         rm -f *.o *.out *~ $(PROGRAM)

