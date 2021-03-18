# Start of the makefile
# Defining variables

FC            = gfortran

CPFLAGS       = 
CPFLAGS       += -Dno_intr_findloc # need to use in cluster
# CPFLAGS       += -Dfull_flux_report # output all cumulative flux
CPFLAGS       += -Ddisp_lim # limiting the display of results
# CPFLAGS       += -Ddiss_only # not allowing precipitation of minerals

# CFLAGS        = -fcheck=all -g -O3  
CFLAGS        = -fimplicit-none  -Wall  -Wline-truncation  -Wcharacter-truncation  -Wsurprising  \
	-Waliasing  -Wimplicit-interface  -Wunused-parameter  -fwhole-file  -fcheck=all  -std=gnu  -pedantic  -fbacktrace -O3

# LDFLAGS       = -L/usr/local/lib
LDFLAGS       = 

LIBS          = -lopenblas
OBJS          = pysil_dev.o 
SRC           = pysil_dev.f90 
# OBJS          = pysil.o 
# SRC           = pysil.f90  
# OBJS          = pysil_save02242021.o 
# SRC           = pysil_save02242021.f90 
# OBJS          = pysil_dev_save_tmp.o 
# SRC           = pysil_dev_save_tmp.f90 
PROGRAM       = weathering

all:            $(PROGRAM)

$(PROGRAM):     $(SRC)
	$(FC) $(SRC) -o $(PROGRAM) -cpp $(CPFLAGS) $(CFLAGS) $(LIBS) $(LDFLAGS)


clean:;         rm -f *.o *~ $(PROGRAM)

