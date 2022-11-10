import os
import sys
import time
import numpy as np
import shutil
import get_int_prof
import make_inputs

phnorm_pw = True
phnorm_pw = False

cec=float(sys.argv[3])
# targetpH = 6.0
targetpH = float(sys.argv[1])
tau = float(sys.argv[2])

dep_sample = 0.18
dep_sample = 0.25
# dep_sample = 0.15

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
runname_field   = expid+'_basalt_field_tpH'+sys.argv[1].replace('.','p')+'_tau'+sys.argv[2].replace('.','p')
runname_lab     = expid+'_basalt_lab_tpH'+sys.argv[1].replace('.','p')+'_tau'+sys.argv[2].replace('.','p')
# runname = 'chk_iter_incl2nd_basalt_tpH'+sys.argv[1].replace('.','p')+'_tau'+sys.argv[2].replace('.','p')

outdir='/storage/scratch1/0/ykanzaki3/scepter_output/'
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

exename = 'scepter'
to = ' '
where = '/'

# ============ common input file modification wrt spinup: for field run
filename = '/switches.in'
src = outdir + spinup_field + filename
dst = outdir + runname_field  + filename
with open(src, 'r') as file:
    data = file.readlines()
data[2] = '2\tbio-mixing style: 0-- no mixing, 1-- fickian mixing, 2-- homogeneous mixng, 3--- tilling, 4--- LABS mixing, if not defined 0 is taken\n'
data[7] = 'true\trestart from a previous run\n'
with open(dst, 'w') as file:
    file.writelines(data)

dustsrc = 'dust_gbasalt.in'
dustdst = 'dust.in'

os.system('cp ' + dustsrc + to + outdir + runname_field + where + dustdst) 


filename = '/slds.in'
src = outdir + spinup_field + filename
dst = outdir + runname_field  + filename
with open(src, 'r') as file:
    data = file.readlines()

data.insert(1, 'gbas\t\n')
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

# data.insert(1, 'k2o\t\nna2o\t\nmgo\t\n')
data.insert(1, 'gbas\t\nk2o\t\nna2o\t\nmgo\t\n')
with open(dst, 'w') as file:
    file.writelines(data)
    
filename = '/solutes.in'
src = outdir + spinup_lab + filename
dst = outdir + runname_lab  + filename
with open(src, 'r') as file:
    data = file.readlines()

data.insert(1, 'k\t\nna\t\nmg\t\n')
with open(dst, 'w') as file:
    file.writelines(data)
    
filename = '/kinspc.in'
src = outdir + spinup_lab + filename
dst = outdir + runname_lab  + filename
with open(src, 'r') as file:
    data = file.readlines()
    
data.insert(1, 'gbas\t0\n')
with open(dst, 'w') as file:
    file.writelines(data)
    
# compile 
os.system('make')
for runname in [runname_field,runname_lab]:
    if not os.path.exists( outdir + runname) : os.system('mkdir -p ' + outdir + runname)
    os.system('cp ' + exename + to + outdir + runname)

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
        
        
    print(outdir+runname_field+'/scepter')
    os.system(outdir+runname_field+'/scepter')


    phint_field = get_int_prof.get_ph_int_site(outdir,runname_field,dep_sample)
    dense_lab = get_int_prof.get_rhobulk_int_site(outdir,runname_field,dep_sample)
    sps = ['g2','inrt','gbas']
    sldwt_list = []
    for sp in sps:
        sldwt = get_int_prof.get_sldwt_int_site(outdir,runname_field,dep_sample,[sp])
        sldwt_list.append(sldwt)
    exchanger = sldwt_list[sps.index('g2')]  + sldwt_list[sps.index('inrt')]
    
    catsat_list = []
    for sp in catlist:
        catsat = get_int_prof.get_spex_int_site(outdir,runname_field,dep_sample,sp)
        catsat_list.append(catsat)
    
    caint = catsat_list[catlist.index('ca')]
    mgint = catsat_list[catlist.index('mg')]
    kint = catsat_list[catlist.index('k')]
    naint = catsat_list[catlist.index('na')]
    ztot = 0.5
    # cec=4
    poro_lab = 5./(1./dense_lab+5.)
    fdust_cao_lab   = ztot*(1-poro_lab)*dense_lab*1e3*exchanger/100.*cec*1e-2*caint/100./2. * 56.1 
    fdust_mgo_lab   = ztot*(1-poro_lab)*dense_lab*1e3*exchanger/100.*cec*1e-2*mgint/100./2. * 40.3
    fdust_na2o_lab  = ztot*(1-poro_lab)*dense_lab*1e3*exchanger/100.*cec*1e-2*naint/100./2. * 62
    fdust_k2o_lab   = ztot*(1-poro_lab)*dense_lab*1e3*exchanger/100.*cec*1e-2*kint /100./2. * 94.2
    
    filename = '/frame.in'
    src = outdir + spinup_lab + filename
    dst = outdir + runname_lab  + filename

    with open(src, 'r') as file:
        data = file.readlines()
    data[5]     = '{:.8f}\tamounts of dusts [g/m2/yr]\n'.format(fdust_cao_lab)
    with open(dst, 'w') as file:
        file.writelines(data)
        
        
    filename = 'dust.in'
    sld_varlist = [ ('cao',1) , ('mgo',fdust_mgo_lab/fdust_cao_lab) 
        , ('na2o',fdust_na2o_lab/fdust_cao_lab) , ('k2o',fdust_k2o_lab/fdust_cao_lab)  ] 
    make_inputs.get_input_sld_properties(
        outdir=outdir
        ,runname=runname_lab
        ,filename = filename
        ,sld_varlist=sld_varlist
        )
        
    pr_list_lab = [(sp,sldwt_list[sps.index(sp)]/100.) for sp in sps]
    atm_list_lab = [('pco2',3.16e-4),('po2',0.21),('pnh3',1e-50),('pn2o',1e-50)]
    make_inputs.get_input_tracer_bounds(
        outdir=outdir
        ,runname=runname_lab
        ,pr_list = pr_list_lab
        # ,rain_list=rain_list
        ,atm_list=atm_list_lab
        )
        
    print(outdir+runname_lab+'/scepter')
    os.system(outdir+runname_lab+'/scepter')
    
    phint = get_int_prof.get_ph_int_site(outdir,runname_lab,dep_sample)

    print(phint_field,phint)
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
            fdust_old = fdust
            imin = np.argmin(1./(data_tmp[:,iph]-targetpH))
            imax = np.argmax(1./(data_tmp[:,iph]-targetpH))
            # fdust = 0.5 * (data_tmp[imin,idust] +  data_tmp[imax,idust])
            # define slope = DpHmax_min/Dfdustmax_min
            # and Dfmax_x = DpHmax_x/slope
            # and then fx = fmax - Dfmax_x
            slope = (data_tmp[imax,iph] - data_tmp[imin,iph])/(data_tmp[imax,idust] - data_tmp[imin,idust])
            fdust = data_tmp[imax,idust] - ( data_tmp[imax,iph] - targetpH )/slope
            if abs((fdust- fdust_old)/fdust) < 1e-6:  fdust = fdust_old * 0.5
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