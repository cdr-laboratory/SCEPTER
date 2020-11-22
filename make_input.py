# -*- coding: utf-8 -*-
# import os

ztot = 2            #total depth of weathering profile [m]
nz = 30             #number of grids into which calculation domain is divided
ttot = 1e2          #total duration of simulation [yr]
tc = 10.05826792             #temperature [oC]
dust = 40e2         #amounts of dusts [g/m2/yr]
OM_rain = 5e2       #OM [g C/m2/yr]
zdust = 0.6             #scale depth of e-folding decrease for supplied dust distribution [m]
poroi = 0.5            #initial porosity
sati = 446.4186198e-3             #water saturation at the surface of profile
zsat = 1e3              #depth of water table [m]
zml = 1.             #depth of mixed layer [m]
w = 5e-5            #uplift rate [m/yr]
q = 0.963             #net water flux [m/yr]
p80 = 1e-6           #radius of particles [m]
iter = 1000            #interations needed to increase time step by a factor of 10
sim_name_restart = 'site98_ex_v5_rain-0.00E+00_p80r-0.10E-05_q-0.10E+01_zsat-1000'
#^ directory name (only used when restart from a previous run; switch is in switches.in)
sim_name = 'site98_basalt_OM-'+str(int(OM_rain))
#^ simulation name



f = open('frame.in','w')
f.write('** values to determine the boundary conditions'+'\n')
f.write(str(ztot)+'\t\t\t'+ 'total depth of weathering profile [m]'+'\n')
f.write(str(int(nz))+'\t\t\t'+'number of grids into which calculation domain is divided'+'\n')
f.write(str(ttot)+'\t\t\t' + 'total duration of simulation [yr]'+'\n')
f.write(str(tc)+'\t\t\t' + 'temperature [oC]'+'\n')
f.write(str(dust)+'\t\t\t' + 'amounts of dusts [g/m2/yr]'+'\n')
f.write(str(OM_rain)+'\t\t\t' + 'OM [g C/m2/yr]'+'\n')
f.write(str(zdust)+'\t\t\t' + 'scale depth of e-folding decrease for supplied dust distribution [m]'+'\n')
f.write(str(poroi)+'\t\t\t' + 'initial porosity'+'\n')
f.write(str(sati)+'\t\t\t' + 'water saturation at the surface of profile'+'\n')
f.write(str(zsat)+'\t\t\t' + 'depth of water table [m]'+'\n')
f.write(str(zml)+'\t\t\t' + 'depth of mixed layer [m]'+'\n')
f.write(str(w)+'\t\t\t' + 'uplift rate [m/yr]'+'\n')
f.write(str(q)+'\t\t\t' + 'net water flux [m/yr]'+'\n')
f.write(str(p80)+'\t\t\t' + 'radius of particles [m]'+'\n')
f.write(str(int(iter))+'\t\t\t' + 'interations needed to increase time step by a factor of 10'+'\n')
f.write(sim_name_restart+'\n')
f.write('^ directory name (only used when restart from a previous run; switch is in switches.in)'+'\n')
f.write(sim_name+'\n')
f.write('^ simulation name')
f.close()

# os.system('./a.exe')
# os.system('C:/cygwin64/home/YK/PyWeath/a.exe')

