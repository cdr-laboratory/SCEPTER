#!/bin/bash

# --- OM rain 
declare -a array=(1)

# ---- particle radius (um)
declare -a array2=(100)

# ---- porosity 
declare -a array3=(25)

# ---- OM decay turn over year
declare -a array4=(16)

# ---- erosion/uplift (+ 5 cm/kyr)
declare -a array5=(0)


# ---- mineralogy 
# plagioclase
declare -a plglist=("ab")

# pyroxene
declare -a pxnlist=("dp")

# hornblende
declare -a hbldlist=("tm")

# 10-A clay
declare -a clay10list=("ill")

# ==== cation exchange
# cation_exchange=true
cation_exchange=false

# ==== choose how to spin-up
opensystem=true

exe_src="../scepter_gmd22"

arraylength=${#array[@]}
arraylength2=${#array2[@]}
arraylength3=${#array3[@]}
arraylength4=${#array4[@]}
arraylength5=${#array5[@]}

plgn=${#plglist[@]}
pxnn=${#pxnlist[@]}
hbldn=${#hbldlist[@]}
clay10n=${#clay10list[@]}

for (( i=1; i<${arraylength}+1; i++ )); do
for (( j=1; j<${arraylength2}+1; j++ )); do
for (( k=1; k<${arraylength3}+1; k++ )); do
for (( m=1; m<${arraylength4}+1; m++ )); do
for (( n=1; n<${arraylength5}+1; n++ )); do



for (( iplg=1; iplg<${plgn}+1; iplg++ )); do
for (( ipxn=1; ipxn<${pxnn}+1; ipxn++ )); do
for (( ihbld=1; ihbld<${hbldn}+1; ihbld++ )); do
for (( iclay10=1; iclay10<${clay10n}+1; iclay10++ )); do


echo ${array[$i-1]}
echo ${array2[$j-1]}
echo ${array3[$k-1]}
echo ${array4[$m-1]}
echo ${array5[$n-1]}

# Mineralogy variation 
echo ${plglist[$iplg-1]}
echo ${pxnlist[$ipxn-1]}
echo ${hbldlist[$ihbld-1]}
echo ${clay10list[$iclay10-1]}

workdir="/storage/project/r-creinhard3-0/ykanzaki3/scepter_output/GMD22/"

file="./pristine_soil_input.txt"

parent=$(basename "$(dirname "$PWD")")
echo "$parent"

branch=$(git rev-parse --abbrev-ref HEAD)
echo "You’re on branch: $branch"

version="_v_${1}"

declare -a f1 f2 f3 f4 f5 f6 f7 f8 f9 f10 f11 f12 f13 f14 f15 f16 f17 f18 f19 f20 f21 f22 f23 f24 f25 f26
while read -r f1 f2 f3 f4 f5 f6 f7 f8 f9 f10 f11 f12 f13 f14 f15 f16 f17 f18 f19 f20 f21 f22 f23 f24 f25 f26;do

if [ "$cation_exchange" = true ] ; then
workname="${parent}${branch}${version}_example_wce_5-1-A1_site${f1}"
else
workname="${parent}${branch}${version}_example_5-1-A1_site${f1}"
fi

if [ "$f1" -ne 1 ] ; then
continue
fi

mkdir -p "$workdir$workname"
# make
cp $exe_src $workdir$workname/scepter

echo 
echo $f1
echo $f2
echo $f3
echo $f4
echo $f5
echo $f6
echo $f7
echo $f8
echo $f9
echo $f10
echo $f11
echo $f12
echo $f13
echo $f14
echo $f15
echo $f16
echo $f17
echo $f18
echo $f19
echo $f20
echo $f21
echo $f22
echo $f23
echo $f24
echo $f25
echo $f26
echo




FILE4="$workdir$workname/solutes.in"
echo "** choose and list from [mg, si, na, ca, al, fe2, fe3, k, so4, no3]" > $FILE4
echo "mg" >> $FILE4
echo "si" >> $FILE4
echo "na" >> $FILE4
echo "ca" >> $FILE4
echo "al" >> $FILE4
echo "fe3" >> $FILE4
echo "k" >> $FILE4

# cp "$FILE4" "$workdir$workname"


FILE3="$workdir$workname/gases.in"
echo "** choose and list from [pco2, po2, pnh3, pn2o]" > $FILE3
echo "pco2" >> $FILE3
echo "po2" >> $FILE3

# cp "$FILE3" "$workdir$workname"

FILE5="$workdir$workname/extrxns.in"
echo "** choose and list from [resp, fe2o2, omomb, ombto, pyfe3, amo2o, g2n0, g2n21, g2n22]" > $FILE5



FILE3="$workdir$workname/kinspc.in"
echo "** specify rate const in [mol/m2/yr] except for OMs which should be presented as turnover year [yr] (e.g., g2   1.0)" > $FILE3
echo "g2     5" >> $FILE3
echo "g3     85" >> $FILE3

# cp "$FILE3" "$workdir$workname"



FILE3="$workdir$workname/a.in"
echo "** parent rock particle radii in meter (e.g., 'ab      1e-5' in one line and 'ka     1e-6' in the next) (if not specified value in frame.in is used for all sld sp.)" > $FILE3

# cp "$FILE3" "$workdir$workname"


FILE3="$workdir$workname/OM_rain.in"
echo "** OM rain fraction wrt the value in frame.in (if not specified assumed 0)" > $FILE3
echo "g1      0.0" >> $FILE3
echo "g2      1." >> $FILE3
echo "g3      0.1016" >> $FILE3

# cp "$FILE3" "$workdir$workname"


FILE3="$workdir$workname/dust.in"
echo "** dusts wt fraction (if not specified assumed 0)" > $FILE3

# cp "$FILE3" "$workdir$workname"

FILE3="$workdir$workname/dust_2nd.in"
echo "** dusts wt fraction (if not specified assumed 0)" > $FILE3

# cp "$FILE3" "$workdir$workname"


FILE3="$workdir$workname/switches.in"
echo "** switch number or on/off [true if on, false if off]" > $FILE3
echo "1           erosion scheme: 0-- cnst w, 1-- cnst poro*w, 2-- cnst (1-poro)*w, 3--- w-flexible(cnst porosity prof), if not defined 0 is taken" >> $FILE3
echo "1           bio-mixing style: 0-- no mixing, 1-- fickian mixing, 2-- homogeneous mixng, 3--- tilling, 4--- LABS mixing, if not defined 0 is taken" >> $FILE3
echo "false            porosity  iteration" >> $FILE3
echo "false            limiting mineral lowest conc." >> $FILE3
# echo "true            display results at runtime" >> $FILE3
# echo "true            limited results display" >> $FILE3
echo "1	display results at runtime: 0-- none, 1-- only reporting time, 2-- every time iteration, if not defined 1 is taken" >> $FILE3
echo "0	report files: 0-- basics, 1-- +saturation time series" >> $FILE3
echo "false           restart from a previous run" >> $FILE3
echo "false            include roughness in mineral surface area" >> $FILE3
echo "false           enabling activity coefficients" >> $FILE3
# echo "true           enabling activity coefficients" >> $FILE3
echo "false           time step fixed" >> $FILE3
echo "${cation_exchange}            enabling adsorption for cation exchange" >> $FILE3
echo "true            adopting a regular grid      " >> $FILE3
if [ "$opensystem" = true ] ; then
# open system spinup
echo '!! open system !!'
echo "false            enforcing solid profile (only aq and gas are explicitly simulated)" >> $FILE3
echo "true            enabling porosity evolution " >> $FILE3
echo "true            enabling SA evolution 1 (SA decreases as porosity increases)" >> $FILE3
echo "false            enabling SA evolution 2 (SA increases with porosity)" >> $FILE3                                          
echo "false            enabling PSD tracking" >> $FILE3
echo "false            enabling PSD tracking for individual solid species" >> $FILE3
echo "false            enabling full seasonality" >> $FILE3
else 
#close (porewater-only) spinup
echo '!! closed (porewater-only) system !!'
echo "true            enforcing solid profile (only aq and gas are explicitly simulated)" >> $FILE3
echo "false            enabling porosity evolution " >> $FILE3
echo "false            enabling SA evolution 1 (SA decreases as porosity increases)" >> $FILE3
echo "false            enabling SA evolution 2 (SA increases with porosity)" >> $FILE3
echo "false            enabling PSD tracking" >> $FILE3
fi

# cp "$FILE3" "$workdir$workname"

dstp="$(echo "1.0*${array[$i-1]}*$f22" | bc -l)"
echo $dstp

wmod_tmp="$(echo "5.0*${f25}+${array5[$n-1]}" | bc -l)"
wmod_min=0.5
if [ `echo "$wmod_tmp < $wmod_min" | bc` == 1 ]; then
wmod=$wmod_min
else
wmod=$wmod_tmp
fi

echo $wmod

varporo="$(echo "0.01*${array3[$k-1]}" | bc -l)"
echo $varporo

qmod="$(echo "1.0*${f19}" | bc -l)"

pmod="$(echo "1.0*${array2[$j-1]}" | bc -l)"

FILE2="$workdir$workname/frame.in"

echo "** values to determine the boundary conditions" > $FILE2
echo "2.0			total depth of weathering profile [m]" >> $FILE2
echo "60			number of grids into which calculation domain is divided" >> $FILE2
echo "100000.0			total duration of simulation [yr]" >> $FILE2
echo "${f21}			temperature [oC]" >> $FILE2
echo "0e2			amounts of dusts [g/m2/yr]" >> $FILE2
echo "0e2			amounts of 2nd dusts [g/m2/yr]" >> $FILE2
echo "0			duration of dust app [yr]" >> $FILE2
echo "${dstp}      OM [g C/m2/yr]" >> $FILE2
echo "2.0			depth of mixed layer for soil [m]" >> $FILE2
echo "${varporo}			initial porosity" >> $FILE2
echo "${f20}e-3			water saturation at the surface of profile" >> $FILE2
echo "1000			depth of water table [m]" >> $FILE2
echo "2.0			depth of mixed layer for dust [m]" >> $FILE2
echo "${wmod}e-05			uplift rate [m/yr]" >> $FILE2
echo "${qmod}e-3			net water flux [m/yr]" >> $FILE2
echo "${pmod}e-6			radius of particles [m]" >> $FILE2
echo "10			interations needed to increase time step by a factor of 10" >> $FILE2
echo "self" >> $FILE2
echo "^ directory name for restart (switch is in switches.in) (type 'self' when restart from the same directory)" >> $FILE2
echo "$workname" >> $FILE2
echo "^ simulation name" >> $FILE2

# cp "$FILE2" "$workdir$workname"


FILE2="$workdir$workname/slds.in"
echo "** choose and list from [fo, ab, an, cc, ka, gb, py, ct, fa, gt, cabd, dp, hb, kfs, om, omb, amsi, arg, dlm, hm, ill, anl, nph, qtz, gps, by, olg, and, tm, la, cpx, en, fer, opx, g1, g2, g3] " > $FILE2

FILE3="$workdir$workname/parentrock.in"
echo "** parent rock wt fraction (e.g., 'ab      0.2' in one line and 'ka     0.001' in the next) (if not specified assumed 1e-20)" > $FILE3


FILE6="$workdir$workname/rain.in"
echo "** solute concs. of rain in mol/L (if not specified assumed 1e-20)" > $FILE6
echo "mg	1e-20" >> $FILE6
echo "si	1e-20" >> $FILE6
echo "na	1e-20" >> $FILE6
echo "ca	1e-20" >> $FILE6
echo "al	1e-20" >> $FILE6
echo "fe3	1e-20" >> $FILE6
echo "k		1e-20" >> $FILE6

threshold=0.00001

if [ `echo "$f2 > $threshold" | bc` == 1 ]; then
echo "amsi" >> $FILE2
echo "amsi    ${f2}e-2" >> $FILE3
fi
if [ `echo "$f3 > $threshold" | bc` == 1 ]; then
echo "arg" >> $FILE2
echo "arg     ${f3}e-2" >> $FILE3
fi
if [ `echo "$f4 > $threshold" | bc` == 1 ]; then
echo "cc" >> $FILE2
echo "cc      ${f4}e-2" >> $FILE3
fi
if [ `echo "$f5 > $threshold" | bc` == 1 ]; then
echo "dlm" >> $FILE2
echo "dlm     ${f5}e-2" >> $FILE3
fi
if [ `echo "$f6 > $threshold" | bc` == 1 ]; then
echo "gb" >> $FILE2
echo "gb     ${f6}e-2" >> $FILE3
fi
if [ `echo "$f7 > $threshold" | bc` == 1 ]; then
echo "gt" >> $FILE2
echo "gt     ${f7}e-2" >> $FILE3
fi
if [ `echo "$f8 > $threshold" | bc` == 1 ]; then
echo "gps" >> $FILE2
echo "gps     ${f8}e-2" >> $FILE3
echo "so4" >> $FILE4
echo "so4	1e-20" >> $FILE6
fi
if [ `echo "$f9 > $threshold" | bc` == 1 ]; then
echo "hm" >> $FILE2
echo "hm      ${f9}e-2" >> $FILE3
fi
if [ `echo "$f10 > $threshold" | bc` == 1 ]; then

if [ "${hbldlist[$ihbld-1]}" == "tm" ]; then 
echo "tm" >> $FILE2
echo "tm      ${f10}e-2" >> $FILE3
elif [ "${hbldlist[$ihbld-1]}" == "antp" ]; then 
echo "antp" >> $FILE2
echo "antp      ${f10}e-2" >> $FILE3
fi

fi
if [ `echo "$f11 > $threshold" | bc` == 1 ]; then
echo "ka" >> $FILE2
echo "ka      ${f11}e-2" >> $FILE3
fi
if [ `echo "$f12 > $threshold" | bc` == 1 ]; then
echo "${pxnlist[$ipxn-1]}" >> $FILE2
echo "${pxnlist[$ipxn-1]}      ${f12}e-2" >> $FILE3
if [ "${pxnlist[$ipxn-1]}" == "dp" ]; then 
echo "*** not including Fe2+ ***"
else 
echo "*** Including Fe2+***"
echo "fe2" >> $FILE4
echo "fe2	1e-20" >> $FILE6
echo "fe2o2" >> $FILE5
fi 
fi
if [ `echo "$f13 > $threshold" | bc` == 1 ]; then
echo "qtz" >> $FILE2
echo "qtz     ${f13}e-2" >> $FILE3
fi
if [ `echo "$f14 > $threshold" | bc` == 1 ]; then
echo "${clay10list[$iclay10-1]}" >> $FILE2
echo "${clay10list[$iclay10-1]}     ${f14}e-2" >> $FILE3
fi
if [ `echo "$f15 > $threshold" | bc` == 1 ]; then
echo "cabd" >> $FILE2
echo "cabd    ${f15}e-2" >> $FILE3
fi
if [ `echo "$f16 > $threshold" | bc` == 1 ]; then
echo "kfs" >> $FILE2
echo "kfs     ${f16}e-2" >> $FILE3
fi
if [ `echo "$f17 > $threshold" | bc` == 1 ]; then
echo "${plglist[$iplg-1]}" >> $FILE2
echo "${plglist[$iplg-1]}      ${f17}e-2" >> $FILE3
fi
if [ `echo "$f18 > $threshold" | bc` == 1 ]; then
echo "anl" >> $FILE2
echo "anl     ${f18}e-2" >> $FILE3
fi
if [ `echo "$f23 > $threshold" | bc` == 1 ]; then
echo "g2" >> $FILE2
echo "g3" >> $FILE2
# enforce OM if closed system
if [ "$opensystem" = false ] ; then
echo "g2     ${f23}e-2" >> $FILE3
fi

fi


# cp "$FILE2" "$workdir$workname"
# cp "$FILE3" "$workdir$workname"
# cp "$FILE4" "$workdir$workname"
# cp "$FILE5" "$workdir$workname"

# cp "$FILE3" "$workdir$workname"


FILE3="$workdir$workname/atm.in"
echo "** atmospheric composition in atm (if not specified assumed 1 PAL)" > $FILE3
echo "pco2      3.16e-4" >> $FILE3
echo "po2       0.21" >> $FILE3
echo "pnh3      1e-50" >> $FILE3
echo "pn2o      1e-50" >> $FILE3

# cp "$FILE3" "$workdir$workname"


# specify secondary phase 
cp "../data/2ndslds_def.in" "$workdir$workname/2ndslds.in"



# when directly running from this shell
$workdir$workname/scepter  


rm -f *.in

done <$file

# BCs
done 
done 
done 
done
done 

# Mineralogy 
done 
done 
done 
done 

