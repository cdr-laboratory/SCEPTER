#!/bin/bash
# gfortran -cpp -o weathering pysil_dev.f90  -lopenblas -g -fcheck=all -O3
# gfortran -cpp -o test1.exe pysil_dev.f90   -lopenblas -g -fcheck=all -O3
# gfortran -cpp -Dphtest -o weathering pysil_dev.f90   -lopenblas -fimplicit-none  -Wall  -Wline-truncation  -Wcharacter-truncation  -Wsurprising \
  # -Waliasing  -Wimplicit-interface  -Wunused-parameter  -fwhole-file  -fcheck=all  -std=gnu  -pedantic  -fbacktrace -O3
make

FILE="./slds.in"
echo "** choose and list from [fo, ab, an, cc, ka, gb, py, ct, fa, gt, cabd, dp, hb, kfs, om, omb, amsi, arg, dlm, hm, ill, anl, nph, qtz, gps, by, olg, and, tm, la, cpx, en, fer, opx, g1, g2, g3] " > $FILE
echo "fo" >> $FILE
echo "ab" >> $FILE
echo "an" >> $FILE
echo "cc" >> $FILE
echo "ka" >> $FILE
echo "fa" >> $FILE
echo "gt" >> $FILE
# echo "cabd" >> $FILE
echo "dp" >> $FILE
echo "hb" >> $FILE
echo "kfs" >> $FILE
# echo "g1" >> $FILE
echo "g2" >> $FILE
# echo "g3" >> $FILE


FILE="./solutes.in"
echo "** choose and list from [mg, si, na, ca, al, fe2, fe3, k, so4, no3]" > $FILE
echo "mg" >> $FILE
echo "si" >> $FILE
echo "na" >> $FILE
echo "ca" >> $FILE
echo "al" >> $FILE
echo "fe2" >> $FILE
echo "fe3" >> $FILE
echo "k" >> $FILE
echo "no3" >> $FILE


FILE="./gases.in"
echo "** choose and list from [pco2, po2, pnh3, pn2o]" > $FILE
echo "po2" >> $FILE
echo "pco2" >> $FILE
echo "pnh3" >> $FILE
echo "pn2o" >> $FILE


FILE="./extrxns.in"
echo "** choose and list from [resp, fe2o2, omomb, ombto, pyfe3, amo2o, g2n0, g2n21, g2n22]" > $FILE
echo "fe2o2" >> $FILE
echo "amo2o" >> $FILE
# echo "g2n0" >> $FILE
echo "g2n21" >> $FILE
echo "g2n22" >> $FILE


FILE="./parentrock.in"
echo "** choose and list from [fo, ab, an, cc, ka, gb, py, ct, fa, gt, cabd, dp, hb, kfs, om, omb, amsi, arg, dlm, hm, ill, anl, nph, qtz, gps, by, olg, and, tm, la, cpx, en, fer, opx, g1, g2, g3] " > $FILE
echo "ka        1e-3" >> $FILE
echo "gt        1e-3" >> $FILE
# echo "cabd      1e-4" >> $FILE


FILE="./rain.in"
echo "** solute concs. of rain in mol/L (if not specified assumed 1e-20)" > $FILE


FILE="./atm.in"
echo "** atmospheric composition in atm (if not specified assumed 1 PAL)" > $FILE
echo "pco2      3.16e-4" >> $FILE
echo "po2       0.21" >> $FILE
echo "pnh3      1.0e-9" >> $FILE
echo "pn2o      270e-9" >> $FILE


FILE="./switches.in"
echo "** switches on and off [true if on, false if off]" > $FILE
echo "false           no bioturbation" >> $FILE
echo "false           Fickian mixng" >> $FILE
echo "false           homogeneous mixng" >> $FILE
echo "false           LABS mixing" >> $FILE
echo "true           tilling" >> $FILE
echo "true            display results at runtime" >> $FILE
echo "false           restart from a previous run" >> $FILE
echo "true            include roughness in mineral surface area" >> $FILE
echo "false           include Al-inhibition in mineral reaction kinetics " >> $FILE
echo "false           time step fixed" >> $FILE
echo "false            pre-calculation before attempting to obtain flly coupled results" >> $FILE
echo "true            adopting a regular grid      " >> $FILE
echo "false            enforcing solid profile (only aq and gas are explicitly simulated)" >> $FILE
echo "true            enabling porosity evolution " >> $FILE
echo "true            enabling SA evolution 1 (SA decreases as porosity increases)" >> $FILE
echo "false            enabling SA evolution 2 (SA increases with porosity)" >> $FILE



FILE="./frame.in"

echo "** values to determine the boundary conditions" > $FILE
echo "2			total depth of weathering profile [m]" >> $FILE
echo "30			number of grids into which calculation domain is divided" >> $FILE
echo "100.0			total duration of simulation [yr]" >> $FILE
echo "15			temperature [oC]" >> $FILE
echo "40e2			amounts of dusts [g/m2/yr]" >> $FILE
echo "5e2      OM [g C/m2/yr]" >> $FILE
echo "0.6			scale depth of e-folding decrease for supplied dust distribution [m]" >> $FILE
echo "0.5			initial porosity" >> $FILE
echo "0.1			water saturation at the surface of profile" >> $FILE
echo "1.0			depth of water table [m]" >> $FILE
echo "1.0			depth of mixed layer [m]" >> $FILE
echo "5e-05			uplift rate [m/yr]" >> $FILE
echo "1.0e-2			net water flux [m/yr]" >> $FILE
echo "1e-06			radius of particles [m]" >> $FILE
echo "1000			interations needed to increase time step by a factor of 10" >> $FILE
echo "test_rain-0.40E+04_p80r-0.10E-05_q-0.10E+00_zsat-1" >> $FILE
echo "^ directory name (only used when restart from a previous run; switch is in switches.in)" >> $FILE
echo "test_v8_N" >> $FILE
echo "^ simulation name" >> $FILE

./weathering

