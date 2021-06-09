#!/bin/bash
workdir="../pyweath_output/"
# workdir="/storage/scratch1/0/ykanzaki3/pyweath_output/"
workname=$1 
mkdir -p "$workdir$workname"

make

cp weathering "$workdir$workname"

FILE="slds.in"
echo "** choose and list from [fo, ab, an, cc, ka, gb, py, ct, fa, gt, cabd, dp, hb, kfs, om, omb, amsi, arg, dlm, hm, ill, anl, nph, qtz, gps, by, olg, and, tm, la, cpx, en, fer, opx, g1, g2, g3, zrc] " > $FILE
# echo "fo" >> $FILE
echo "ab" >> $FILE
# echo "an" >> $FILE
# echo "cc" >> $FILE
# echo "agt" >> $FILE
echo "ka" >> $FILE
# echo "fa" >> $FILE
# echo "cabd" >> $FILE
# echo "dp" >> $FILE
# echo "hb" >> $FILE
# echo "kfs" >> $FILE
# echo "gps" >> $FILE
# echo "py" >> $FILE
# echo "gt" >> $FILE
# echo "g1" >> $FILE
# echo "g2" >> $FILE
# echo "g3" >> $FILE
# echo "zrc" >> $FILE
# echo "amsi" >> $FILE

cp "$FILE" "$workdir$workname"

FILE="solutes.in"
echo "** choose and list from [mg, si, na, ca, al, fe2, fe3, k, so4, no3]" > $FILE
# echo "mg" >> $FILE
echo "si" >> $FILE
echo "na" >> $FILE
# echo "ca" >> $FILE
echo "al" >> $FILE
# echo "fe2" >> $FILE
# echo "fe3" >> $FILE
# echo "k" >> $FILE
# echo "so4" >> $FILE
# echo "no3" >> $FILE

cp "$FILE" "$workdir$workname"

FILE="gases.in"
echo "** choose and list from [pco2, po2, pnh3, pn2o]" > $FILE
# echo "po2" >> $FILE
# echo "pco2" >> $FILE
# echo "pnh3" >> $FILE
# echo "pn2o" >> $FILE


cp "$FILE" "$workdir$workname"

FILE="extrxns.in"
echo "** choose and list from [resp, fe2o2, omomb, ombto, pyfe3, amo2o, g2n0, g2n21, g2n22]" > $FILE
# echo "fe2o2" >> $FILE
# echo "pyfe3" >> $FILE
# echo "amo2o" >> $FILE
# echo "g2n0" >> $FILE
# echo "g2n21" >> $FILE
# echo "g2n22" >> $FILE

cp "$FILE" "$workdir$workname"


FILE="kinspc.in"
echo "** specify rate const in [mol/m2/yr] except for OMs which should be presented as turnover year [yr] (e.g., g2   1.0)" > $FILE
# echo "g2     1" >> $FILE

cp "$FILE" "$workdir$workname"


FILE="parentrock.in"
echo "** parent rock wt fraction (e.g., 'ab      0.2' in one line and 'ka     0.001' in the next) (if not specified assumed 1e-20)" > $FILE
echo "ka        1e-5" >> $FILE
# echo "ab        1.0" >> $FILE
echo "ab        0.1" >> $FILE
# echo "ab        0.03" >> $FILE
# echo "ab        0.01" >> $FILE
# echo "ab        0.003" >> $FILE
# echo "fa        0.1" >> $FILE
# echo "amsi        0.5" >> $FILE
# echo "agt        0.2" >> $FILE
# echo "agt        0.02" >> $FILE
# echo "gps        1e-3" >> $FILE
# echo "gps        0.1140" >> $FILE
# echo "gps        0.1141" >> $FILE
# echo "gps        0.1135" >> $FILE
# echo "gps        0.005" >> $FILE
# echo "cc        0.03" >> $FILE
# echo "cc        0.004" >> $FILE
# echo "py        0.0056" >> $FILE
# echo "py        0.0280" >> $FILE
# echo "py        0.0560" >> $FILE
# echo "py        0.0112" >> $FILE
# echo "gt        1e-5" >> $FILE
# echo "cabd      1e-4" >> $FILE


cp "$FILE" "$workdir$workname"

FILE="rain.in"
echo "** solute concs. of rain in mol/L (if not specified assumed 1e-20)" > $FILE

cp "$FILE" "$workdir$workname"

FILE="atm.in"
echo "** atmospheric composition in atm (if not specified assumed 1 PAL)" > $FILE
echo "pco2      3.16e-4" >> $FILE
echo "po2       0.21" >> $FILE
# echo "pnh3      1.0e-9" >> $FILE
# echo "pn2o      270e-9" >> $FILE
echo "pnh3      1e-100" >> $FILE
echo "pn2o      1e-100" >> $FILE

cp "$FILE" "$workdir$workname"


FILE="OM_rain.in"
echo "** OM rain fraction wrt the value in frame.in (if not specified assumed 0)" > $FILE
echo "g1      0.0" >> $FILE
echo "g2      1." >> $FILE
echo "g3      0.0" >> $FILE

cp "$FILE" "$workdir$workname"


FILE="dust.in"
FILE2="dust_basalt.in"
cp "$FILE2" "$FILE"
cp "$FILE" "$workdir$workname"


FILE="T_temp.in"
cp "$FILE" "$workdir$workname"
FILE="q_temp.in"
cp "$FILE" "$workdir$workname"
FILE="Wet_temp.in"
cp "$FILE" "$workdir$workname"


FILE="switches.in"
echo "** switch number or on/off [true if on, false if off]" > $FILE
echo "0           erosion scheme: 0-- cnst w, 1-- cnst poro*w, 2-- cnst (1-poro)*w, 3--- w-flexible, if not defined 0 is taken" >> $FILE
echo "0           bio-mixing style: 0-- no mixing, 1-- fickian mixing, 2-- homogeneous mixng, 3--- tilling, 4--- LABS mixing, if not defined 0 is taken" >> $FILE
echo "true            display results at runtime" >> $FILE
echo "true            limited results display" >> $FILE
echo "false           restart from a previous run" >> $FILE
echo "false            include roughness in mineral surface area" >> $FILE
echo "false           include Al-inhibition in mineral reaction kinetics " >> $FILE
echo "false           time step fixed" >> $FILE
echo "false            pre-calculation before attempting to obtain flly coupled results" >> $FILE
echo "true            adopting a regular grid      " >> $FILE
echo "false            enforcing solid profile (only aq and gas are explicitly simulated)" >> $FILE
echo "true            enabling porosity evolution " >> $FILE
echo "true            enabling SA evolution 1 (SA decreases as porosity increases)" >> $FILE
echo "false            enabling SA evolution 2 (SA increases with porosity)" >> $FILE
echo "false            enabling PSD tracking" >> $FILE

cp "$FILE" "$workdir$workname"

FILE="frame.in"

echo "** values to determine the boundary conditions" > $FILE
echo "2			total depth of weathering profile [m]" >> $FILE
echo "30			number of grids into which calculation domain is divided" >> $FILE
echo "1e5			total duration of simulation [yr]" >> $FILE
echo "15			temperature [oC]" >> $FILE
echo "0e2			amounts of dusts [g/m2/yr]" >> $FILE
echo "0e2      OM [g C/m2/yr]" >> $FILE
echo "0.6			scale depth of e-folding decrease for supplied dust distribution [m]" >> $FILE
echo "0.1			initial porosity" >> $FILE
echo "0.1			water saturation at the surface of profile" >> $FILE
echo "1.0			depth of water table [m]" >> $FILE
echo "1.0			depth of mixed layer [m]" >> $FILE
echo "5e-05			uplift rate [m/yr]" >> $FILE
echo "0.01			net water flux [m/yr]" >> $FILE
echo "1e-05			radius of particles [m]" >> $FILE
echo "10			interations needed to increase time step by a factor of 10" >> $FILE
echo "gps_base_rain-0.00E+00_pevol_sevol1_p80r-0.10E-05_q-0.10E-01_zsat-1" >> $FILE
# echo "self" >> $FILE
echo "^ directory name (only used when restart from a previous run; switch is in switches.in)" >> $FILE
# echo "gps_base" >> $FILE
echo "$workname" >> $FILE
echo "^ simulation name" >> $FILE

cp "$FILE" "$workdir$workname"


if [ -z "$2" ]; then
# when directly running from this shell
$workdir$workname/weathering  
else
# else submit as a job in pbs file 
FILE="weather_test_pbs.pbs"

echo "#PBS -N $workname           # job name" > $FILE
echo "#PBS -A GT-creinhard3-FY20Phase2               # account to which job is charged" >> $FILE
echo "#PBS -l nodes=1:ppn=1           # number of nodes and cores per node required" >> $FILE
echo "#PBS -l pmem=2gb                # memory per core" >> $FILE
echo "#PBS -l walltime=05:20:00          # duration of the job (ex: 15 min)" >> $FILE
echo "#PBS -j oe                      # combine output and error messages into 1 file" >> $FILE
echo "#PBS -o $workname.out       # output file name" >> $FILE
echo "  " >> $FILE
echo "cd \$PBS_O_WORKDIR" >> $FILE
echo "module load gcc/8.3.0" >> $FILE
echo "module load openblas/0.3.7" >> $FILE
echo "$workdir$workname/weathering  " >> $FILE

cp "$FILE" "$workdir$workname"

echo "${workdir}${workname}/${FILE}"
qsub "${workdir}${workname}/${FILE}"
fi


