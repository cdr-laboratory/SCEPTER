import os
import sys
import time
import numpy as np
import shutil
import get_int_prof
import make_inputs
import random

include_N = False
include_N = True

include_Al = False
# include_Al = True

alphase = 'gb'
alphase = 'amal'

use_CaCl2 = False
# use_CaCl2 = True

include_DIC = False
include_DIC = True

phnorm_pw = True
phnorm_pw = False

use_local_storage = True
use_local_storage = False

cec=float(sys.argv[3])
# targetpH = 6.0
targetpH = float(sys.argv[1])
tau = float(sys.argv[2])

dep_sample = 0.18
dep_sample = 0.25
dep_sample = 0.15

added_sp = 'gbas'
# added_sp = 'cc'

# mxing style

imix = 2 # homogeneous
imix = 3 # tilling

    
ztot_lab = 0.5
ttot_lab = 1000

ztot_lab = 0.05
ttot_lab = 100

water_frac = 5.
water_frac = 1.
# water_frac = 2.5

catlist = ['ca','mg','k','na']

spinid = 'test_inert_fert_buff_v2'
spinid = 'test_inert_fert_buff_v2_omx4'
spinid = 'test_inert_fert_buff_v2_omx9'
spinid = 'test_inert_fert_buff_3v'
spinid = 'test_inert_fert_buff_3v_v5'
spinid = 'test_inert_buff_3v_v5_pw'
# spinid = 'test_inert_buff_3v_v5_sw'
spinid = sys.argv[5]

# spinup = 'test_iter_excl2nd'
# spinup = 'test_inert_spintuneup'
# spinup = 'test_inert_spintuneup_incl2nd'
# spinup_field    = 'test_inert_buff_spintuneup_field'
# spinup_lab      = 'test_inert_buff_spintuneup_lab'
spinup_field    = spinid+'_spintuneup_field'
spinup_lab      = spinid+'_spintuneup_lab'

expid = 'test_inert_fert_buff_excl2nd_v2'
expid = 'test_inert_fert_buff_excl2nd_v2_omx4'
expid = 'test_inert_fert_buff_excl2nd_v2_omx9'
expid = 'test_inert_fert_buff_excl2nd_3v'
expid = 'test_inert_fert_buff_excl2nd_3v_v5'
expid = 'test_inert_buff_excl2nd_3v_v5_pw'
expid = sys.argv[4]
# expid = 'test_inert_buff_excl2nd_3v_v5_sw'

# runname_field   = 'test_iter_buff_excl2nd_basalt_field_tpH'+sys.argv[1].replace('.','p')+'_tau'+sys.argv[2].replace('.','p')
# runname_lab     = 'test_iter_buff_excl2nd_basalt_lab_tpH'+sys.argv[1].replace('.','p')+'_tau'+sys.argv[2].replace('.','p')
runname_field   = expid+'_'+added_sp+'_field_tpH'+sys.argv[1].replace('.','p')+'_tau'+sys.argv[2].replace('.','p')
runname_lab     = expid+'_'+added_sp+'_lab_tpH'+sys.argv[1].replace('.','p')+'_tau'+sys.argv[2].replace('.','p')
# runname = 'chk_iter_incl2nd_basalt_tpH'+sys.argv[1].replace('.','p')+'_tau'+sys.argv[2].replace('.','p')

outdir = '../scepter_output/'
if use_local_storage:  
    outdir_src = '../scepter_output/'
    outdir = os.environ['TMPDIR'] + '/scepter_output/'
    for runname in [spinup_field,spinup_lab]:

        src = outdir_src + runname
        dst = outdir + runname

        if not os.path.exists(dst): 
            shutil.copytree(src, dst)
        else:
            shutil.rmtree(dst)
            shutil.copytree(src, dst)
datadir='./data/'

# dupricate directories from spinups

for (runname,spinup) in [(runname_field,spinup_field),(runname_lab,spinup_lab)]:

    src = outdir + spinup
    dst = outdir + runname

    if not os.path.exists(dst): 
        shutil.copytree(src, dst)
    else:
        shutil.rmtree(dst)
        shutil.copytree(src, dst)


error = 1e4
tol = 1e-4

fdust = 5000
taudust =0.005

# if added_sp == 'cc':
    # fdust = 10000

exename = 'scepter'
# exename_src = 'scepter_test'
exename_src = 'scepter'
to = ' '
where = '/'

# ============ common input file modification wrt spinup: for field run
filename = '/switches.in'
src = outdir + spinup_field + filename
dst = outdir + runname_field  + filename
with open(src, 'r') as file:
    data = file.readlines()
data[2] = '{:d}\tbio-mixing style: 0-- no mixing, 1-- fickian mixing, 2-- homogeneous mixng, 3--- tilling, 4--- LABS mixing, if not defined 0 is taken\n'.format(imix)
data[7] = 'true\trestart from a previous run\n'
data[-3] = 'true\tenabling PSD tracking\n'
data[-2] = 'true\tenabling PSD tracking for individual solid species\n'
with open(dst, 'w') as file:
    file.writelines(data)

if added_sp == 'gbas': dustsrc = './data/dust_gbasalt.in'
if added_sp == 'cc': dustsrc = './data/dust_lime.in'
dustdst = 'dust.in'

os.system('cp ' + dustsrc + to + outdir + runname_field + where + dustdst) 


filename = '/slds.in'
src = outdir + spinup_field + filename
dst = outdir + runname_field  + filename
with open(src, 'r') as file:
    data = file.readlines()

data.insert(1, added_sp+'\t\n')
with open(dst, 'w') as file:
    file.writelines(data)
    
# ============ adding Fe(II) as tracer and its oxidation =================
# filename = '/solutes.in'
# src = outdir + spinup + filename
# dst = outdir + runname  + filename
# with open(src, 'r') as file:
    # data = file.readlines()

# data.insert(1, 'fe2\t\n')
# with open(dst, 'w') as file:
    # file.writelines(data)
    
# filename = '/gases.in'
# src = outdir + spinup + filename
# dst = outdir + runname  + filename
# with open(src, 'r') as file:
    # data = file.readlines()

# data.insert(1, 'po2\t\n')
# with open(dst, 'w') as file:
    # file.writelines(data)
    
# filename = '/extrxns.in'
# src = outdir + spinup + filename
# dst = outdir + runname  + filename
# with open(src, 'r') as file:
    # data = file.readlines()

# data.insert(1, 'fe2o2\t\n')
# with open(dst, 'w') as file:
    # file.writelines(data)
    
# ============ common input file modification wrt spinup: for lab run ============

filename = '/slds.in'
src = outdir + spinup_lab + filename
dst = outdir + runname_lab  + filename
with open(src, 'r') as file:
    data = file.readlines()

data.insert(1, added_sp+'\t\n')
with open(dst, 'w') as file:
    file.writelines(data)
    
# filename = '/solutes.in'
# src = outdir + spinup_lab + filename
# dst = outdir + runname_lab  + filename
# with open(src, 'r') as file:
    # data = file.readlines()

# data.insert(1, 'k\t\nna\t\nmg\t\n')
# with open(dst, 'w') as file:
    # file.writelines(data)
    
filename = '/kinspc.in'
src = outdir + spinup_lab + filename
dst = outdir + runname_lab  + filename
with open(src, 'r') as file:
    data = file.readlines()
    
data.insert(1, added_sp+'\t0\n')
with open(dst, 'w') as file:
    file.writelines(data)
    
filename = '2ndslds.in'
srcfile = './data/2ndslds_def.in'
make_inputs.get_input_sld_properties(
    outdir=outdir
    ,runname=runname_field
    ,filename = filename
    ,srcfile = srcfile
    )
make_inputs.get_input_sld_properties(
    outdir=outdir
    ,runname=runname_lab
    ,filename = filename
    ,srcfile = srcfile
    )
    
# compile 
# os.system('make')
# os.system('make --file=makefile_test')
# for runname in [runname_field,runname_lab]:
    # if not os.path.exists( outdir + runname) : 
        # os.system('mkdir -p ' + outdir + runname)
        # os.system('cp ' + exename_src + to + outdir + runname + where + exename)

maxiter = 50

res_list = []
cnt = 0

# get ph of spin-up
phint = get_int_prof.get_ph_int_site(outdir,spinup_lab,dep_sample)
phint_field = get_int_prof.get_ph_int_site(outdir,spinup_field,dep_sample)
ymx = phint - targetpH
if phnorm_pw: ymx = phint_field - targetpH
res_list.append([cnt, phint_field,phint, targetpH, 0, abs( ymx/targetpH ) ])

cnt += 1

while (error > tol):
    
    
    ## /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    
    ## --- setup for field run --- ##
    
    filename = '/frame.in'
    src = outdir + spinup_field + filename
    dst = outdir + runname_field  + filename

    with open(src, 'r') as file:
        data = file.readlines()
    data[3]     = '{:.8f}\ttotal duration of simulation [yr]\n'.format(tau)
    data[5]     = '{:.8f}\tamounts of dusts [g/m2/yr]\n'.format(fdust)
    data[7]     = '{:.8f}\tduration of dust application [yr]\n'.format(taudust)
    data[18]    = '{}\n'.format(spinup_field)
    data[20]    = '{}\n'.format(runname_field)
    with open(dst, 'w') as file:
        file.writelines(data)
        
        
    ## --- run field run --- ##
    
    print(outdir+runname_field+'/scepter')
    os.system(outdir+runname_field+'/scepter')

    ## --- getting data from field run --- ##
    
    phint_field = get_int_prof.get_ph_int_site(outdir,runname_field,dep_sample)
    dense_lab = get_int_prof.get_rhobulk_int_site(outdir,runname_field,dep_sample)
    sps = ['g2','inrt',added_sp]
    if include_Al:  sps.append( alphase )
    sldwt_list = []
    for sp in sps:
        sldwt = get_int_prof.get_sldwt_int_site(outdir,runname_field,dep_sample,[sp])
        sldwt_list.append(sldwt)
    
    aqsps,btmconcs,dep = get_int_prof.get_totsave_site(outdir,runname_field,dep_sample)  # returning mol/ solid m3 depth averaged value 
    print(aqsps,btmconcs,dep)
    
    dic,dep = get_int_prof.get_ave_DIC_save(outdir,runname_field,dep_sample) # returning DIC in mol/ solid m3
    
    ## /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    
    ## --- setup for lab run --- ##
    
    poro_lab = water_frac/(1./dense_lab+water_frac)
    
    oxide_ctnm_list = ['ca','mg','na','k'] 
    oxide_oxnm_list = ['cao','mgo','na2o','k2o'] 
    oxide_stch_list = [1,1,2,2] 
    oxide_mass_list = [56.1 ,40.3, 62, 94.2]
    
    if include_N:
        oxide_ctnm_list.append( 'no3' )
        oxide_oxnm_list.append( 'amnt' )
        oxide_stch_list.append( 2 )
        oxide_mass_list.append( 80 )
    if include_Al:
        oxide_ctnm_list.append( 'al' )
        oxide_oxnm_list.append( 'al2o3' )
        oxide_stch_list.append( 2 )
        oxide_mass_list.append( 1.02E+02 )
    
    fdust_list = []
    fdust_nm_list = []

    for sp in aqsps:
        isp = aqsps.index(sp)
        iox = oxide_ctnm_list.index(sp)
        conc = btmconcs[isp]
        fdust_tmp = ztot_lab*(1-poro_lab)*conc* oxide_mass_list[iox]/oxide_stch_list[iox]
        fdust_list.append(fdust_tmp)
        fdust_nm_list.append(oxide_oxnm_list[iox])
        
        
    fdust_lab = fdust_list[aqsps.index('ca')]
    
    if use_CaCl2:
        cacl2_conc = 0.01
        cacl2_wt2  = 110.98
        fdust_cacl2 = ztot_lab*poro_lab*cacl2_conc*1e3*cacl2_wt2
    
        if fdust_cacl2 > 0:
            fdust_list.append(fdust_cacl2)
            fdust_nm_list.append('cacl2')
            
            
    if include_DIC:                                     # (added 3.23.2023)
        fdust_dic = ztot_lab*(1-poro_lab)*dic* 30.      # (added 3.22.2023)
        fdust_list.append(fdust_dic)                    # (added 3.22.2023)
        fdust_nm_list.append('g1')                      # (added 3.22.2023)

    fdust_list = [fdust_tmp/fdust_lab  for fdust_tmp in fdust_list  ]  
    
    filename = '/frame.in'
    src = outdir + spinup_lab + filename
    dst = outdir + runname_lab  + filename

    with open(src, 'r') as file:
        data = file.readlines()
    data[1]     = '{:.8f}\ttotal depth of weathering profile [m]\n'.format(ztot_lab)
    data[3]     = '{:.8f}\ttotal duration of simulation [yr]\n'.format(ttot_lab)
    data[5]     = '{:.8f}\tamounts of dusts [g/m2/yr]\n'.format(fdust_lab)
    data[10]     = '{:.8f}\tinitial porosity\n'.format(poro_lab)
    with open(dst, 'w') as file:
        file.writelines(data)
        
        
    filename = 'dust.in'
    sld_varlist = [ ( fdust_nm_list[i], fdust_list[i]) for i in range(len(fdust_nm_list)) ] 
    make_inputs.get_input_sld_properties(
        outdir=outdir
        ,runname=runname_lab
        ,filename = filename
        ,sld_varlist=sld_varlist
        )
        
    pr_list_lab = [(sp,sldwt_list[sps.index(sp)]/100.) for sp in sps]
    atm_list_lab = [('pco2',3.16e-4),('po2',0.21),('pnh3',1e-50),('pn2o',1e-50)]
    if include_DIC: atm_list_lab = [('pco2',1e-20),('po2',0.21),('pnh3',1e-50),('pn2o',1e-50)] # (added 3.22.2023)
    make_inputs.get_input_tracer_bounds(
        outdir=outdir
        ,runname=runname_lab
        ,pr_list = pr_list_lab
        # ,rain_list=rain_list
        ,atm_list=atm_list_lab
        )
    
    # break
    
    ## --- run lab run --- ##
    
    print(outdir+runname_lab+'/scepter')
    os.system(outdir+runname_lab+'/scepter')
    
    ## --- getting data from lab run --- ##
    
    phint = get_int_prof.get_ph_int_site(outdir,runname_lab,dep_sample)
    print(phint_field,phint)
    
    ## /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    
    ## --- Newton + bisection iteration --- ## 
    
    time.sleep(5)
    

    ymx = phint - targetpH
    if phnorm_pw: ymx = phint_field - targetpH
    ymx = -ymx
        

    emx = ymx/targetpH
    error = np.max(np.abs(emx))

    print(fdust)
    
    print(' ')
    print(' ')
    print(' ')
    print(' ')
    print("******* {:d}'s iteration, error = {:.8f}".format(cnt,error))
    print(' ')
    print(' ')
    print(' ')
    print(' ')
    
    res_list.append([cnt, phint_field, phint, targetpH, fdust, abs( ymx/targetpH ) ])
    
    cnt += 1

    if np.isnan(ymx):
        print('nan detected in solution')
        error = 1e-99
    else:
        data_tmp = np.array(res_list)
        iph = 1 
        idust = 3 
        if not phnorm_pw: iph = 2
        idust = 4 
        if np.max(data_tmp[:,iph]) < targetpH: 
            fdust = np.max(data_tmp[:,idust])*1.5
            print(np.max(data_tmp[:,iph]), targetpH, fdust)
        else:
            # if cnt >2:  fdust_old_old = fdust_old
            # fdust_old = fdust
            imin = np.argmin(1./(data_tmp[:,iph]-targetpH))
            imax = np.argmax(1./(data_tmp[:,iph]-targetpH))
            # fdust = 0.5 * (data_tmp[imin,idust] +  data_tmp[imax,idust])
            # define slope = DpHmax_min/Dfdustmax_min
            # and Dfmax_x = DpHmax_x/slope
            # and then fx = fmax - Dfmax_x
            slope = (data_tmp[imax,iph] - data_tmp[imin,iph])/(data_tmp[imax,idust] - data_tmp[imin,idust])
            fdust = data_tmp[imax,idust] - ( data_tmp[imax,iph] - targetpH )/slope
            
            # detecting if the updated fdust is already tested
            # if it is, it means non-linearlity
            iclose = np.argmin(np.abs(data_tmp[:,idust]-fdust))
            fdust_old =  data_tmp[iclose,idust]
            
            # just avoid repeated solution by randomly picking up from most likely range 
            while abs((fdust- fdust_old)/fdust) < 1e-6:  # the updated fdust is already tested  
                fdust = random.uniform(data_tmp[imin,idust],data_tmp[imax,idust])
                iclose = np.argmin(np.abs(data_tmp[:,idust]-fdust))
                fdust_old =  data_tmp[iclose,idust]
            
            if abs((fdust- fdust_old)/fdust) < 1e-6:  # the updated fdust is already tested  
                if (data_tmp[iclose,iph] - targetpH)* (data_tmp[imax,iph] - targetpH) >= 0.: 
                    # if the tested value is above target pH, new slope and approximation is calculated switching imax with iclose
                    # Newton-ish new solution
                    # slope = (data_tmp[iclose,iph] - data_tmp[imin,iph])/(data_tmp[iclose,idust] - data_tmp[imin,idust])
                    # fdust = data_tmp[iclose,idust] - ( data_tmp[iclose,iph] - targetpH )/slope
                    
                    # Bisection-ish new solution
                    fdust = 0.5*(data_tmp[iclose,idust] + data_tmp[imin,idust])
                elif (data_tmp[iclose,iph] - targetpH)* (data_tmp[imin,iph] - targetpH) >= 0.:
                    # if the tested value is below target pH, new slope and approximation is calculated switching imin with iclose
                    # Newton-ish new solution
                    # slope = (data_tmp[imax,iph] - data_tmp[iclose,iph])/(data_tmp[imax,idust] - data_tmp[iclose,idust])
                    # fdust = data_tmp[imax,idust] - ( data_tmp[imax,iph] - targetpH )/slope
                    
                    # Bisection-ish new solution
                    fdust = 0.5*(data_tmp[iclose,idust] + data_tmp[imax,idust])
                else:
                    # there must be something is wrong
                    print('there must be something wrong')
                    fdust = np.nan
                # fdust = fdust_old * 0.5
            
            # if cnt >2 and abs((fdust- fdust_old_old)/fdust) < 1e-6:  fdust = fdust_old_old * 0.5
            print(data_tmp[imin,iph],data_tmp[imax,iph], targetpH 
                ,data_tmp[imin,idust],data_tmp[imax,idust], fdust)
    
    time.sleep(5)
    
    if cnt > maxiter: break
    
    for runname in [runname_field,runname_lab]:
        np.savetxt(outdir + runname + where + 'iteration_tmp.res',np.array(res_list))
    
    

name_list = [
    'iter.'
    ,'porewater_pH[-]'
    ,'soil_pHw[-]'
    ,'target_pH[-]'
    ,'basalt[g/m2/yr]'
    ,'error'
    ]

for runname in [runname_field,runname_lab]:
    dst = outdir + runname + where + 'iteration.res'

    with open(dst, 'w') as file:
        for item in name_list:
            if name_list.index(item)==len(name_list)-1:
                file.write('{}\n'.format(item))
            else:
                file.write('{}\t'.format(item))
        for j in range(len(res_list)):
            item_list = res_list[j]
            for i in range(len(item_list)):
                item = item_list[i]
                if i==0:
                    file.write('{:d}\t'.format(item))
                elif i==len(item_list)-1:
                    file.write('{:.6e}\n'.format(item))
                else:
                    file.write('{:.6e}\t'.format(item))
            
    print(res_list)
    
    
if use_local_storage:
    for runname in [runname_field,runname_lab]:
        src = outdir + runname 
        dst = outdir_src + runname 
        
        if not os.path.exists(dst): 
            shutil.copytree(src, dst)
        else:
            shutil.rmtree(dst)
            shutil.copytree(src, dst)