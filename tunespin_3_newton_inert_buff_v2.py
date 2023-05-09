import os
import numpy as np
import get_int_prof
import make_inputs
import time
import sys

print(sys.argv)

include_N = False
include_N = True

include_Al = False
# include_Al = True

alphase = 'gb'
alphase = 'amal'

use_CaCl2 = False
# use_CaCl2 = True

phnorm_pw = True
phnorm_pw = False

water_frac = 5.
water_frac = 1.

iter_max = 120

cec=4
cec=float(sys.argv[2])
targetpH = 5.8
targetpH = float(sys.argv[3])
basesat = 75
acidsat = 25
acidsat = float(sys.argv[4])
targetOM = 5
targetOM = float(sys.argv[5])
z = 2 
alpha = 3.4
kh = 10**3.5
logkh = np.log10(kh)
ca=500 # uM

dense_lab = 2.59296482412060000 # only inertphase

dep_sample = 0.15
# dep_sample = 0.18
# dep_sample = 0.25


# outdir = '/storage/scratch1/0/ykanzaki3/scepter_output/'
outdir = '/storage/coda1/p-creinhard3/0/ykanzaki3/scepter_output/'
# runname_field = 'test_inert_buff_spintuneup_field'
# runname_lab = 'test_inert_buff_spintuneup_lab'
simid = 'test_inert_fert_buff_3v'
simid = 'test_inert_fert_buff_3v_v2'
simid = 'test_inert_fert_buff_3v_v3'
simid = 'test_inert_fert_buff_3v_v4'
simid = 'test_inert_fert_buff_3v_v5'
simid = 'test_inert_fert_buff_3v_v5_chk'
simid = 'test_inert_fert_buff_3v_v5_chk2'
simid = 'test_inert_buff_3v_v5_pw'
simid = 'test_inert_buff_3v_v5_sw'
simid = sys.argv[1]
runname_field   = simid+'_spintuneup_field'
runname_lab     = simid+'_spintuneup_lab'

# compile 
exename = 'scepter'
to = ' '
where = '/'
os.system('make')
if not os.path.exists( outdir + runname_field) : os.system('mkdir -p ' + outdir + runname_field)
if not os.path.exists( outdir + runname_lab) : os.system('mkdir -p ' + outdir + runname_lab)
os.system('cp ' + exename + to + outdir + runname_field)
os.system('cp ' + exename + to + outdir + runname_lab)

ztot=0.5
nz=30
ttot_field=10000
ttot_lab=1000
temp_field=15
temp_field=18
temp_field=float(sys.argv[6])
temp_lab=25
fdust_field=0
fdust_lab=0
fdust2=0
taudust_field=0
taudust_lab=0.01
# omrain_field=100
# omrain_field=400
omrain_field=900
omrain_lab=0
zom=0.5
poro_field=0.5
poro_lab=0.928391508
moistsrf_field=0.5
moistsrf_field=float(sys.argv[7])
moistsrf_lab=1.0
zwater=100000
# zdust_field=0.18
zdust_field=0.25
zdust_lab=0.15
w_field=100e-5
w_field=float(sys.argv[9])
w_lab=0
q_field=1200e-3
q_field=float(sys.argv[8])
q_lab=0
p=10e-6
nstep=10
rstrt='self'
runid_field=runname_field
runid_lab=runname_lab

N_rain = 0  # gN/m2/yr
N_rain = 8.406375  # gN/m2/yr ( <---> 75 lbs/acre/year)
N_rain = 24.6587   # gN/m2/yr ( <---> 220 lbs/acre/year)
N_rain = float(sys.argv[10])   # gN/m2/yr ( <---> 220 lbs/acre/year)
N_rain = N_rain/14  # mol N/m2/yr
N_rain = N_rain*80  # g NH4NO3/m2/yr
N_rain = N_rain/2.  # only half is required as 1 mol NH4NO3 contains 2 moles of N
N_rain = N_rain/omrain_field  # normalize against OC rain 

make_inputs.get_input_frame(
    outdir=outdir
    ,runname=runname_field
    ,ztot=ztot
    ,nz=nz
    ,ttot=ttot_field
    ,temp=temp_field
    ,fdust=fdust_field
    ,fdust2=fdust2
    ,taudust=taudust_field
    ,omrain=omrain_field
    ,zom=zom
    ,poro=poro_field
    ,moistsrf=moistsrf_field
    ,zwater=zwater
    ,zdust=zdust_field
    ,w=w_field
    ,q=q_field
    ,p=p
    ,nstep=nstep
    ,rstrt=rstrt
    ,runid=runname_field
    )

w_scheme_field=1
mix_scheme_field=1
poro_iter_field='false' 
sldmin_lim ='true'
display='true'
disp_lim='true'
restart ='false'
rough_field      ='true'
al_inhib ='false'
dt_fix='false'
cec_on='false'
cec_on='true'
dz_fix='true'
close_aq_field='false'
poro_evol='true'
sa_evol_1 ='true'
sa_evol_2='false'
psd_bulk_field='true'
psd_full_field='true'
season='false'

w_scheme_lab=0
mix_scheme_lab=0 
poro_iter_lab='true' 
rough_lab      ='false'
close_aq_lab='true'
psd_bulk_lab='false'
psd_full_lab ='false'
    
make_inputs.get_input_switches(
    outdir=outdir
    ,runname=runname_field
    ,w_scheme=w_scheme_field
    ,mix_scheme=mix_scheme_field
    ,poro_iter=poro_iter_field
    ,sldmin_lim=sldmin_lim 
    ,display=display
    ,disp_lim=disp_lim
    ,restart=restart 
    ,rough=rough_field
    ,al_inhib=al_inhib 
    ,dt_fix=dt_fix
    ,cec_on=cec_on
    ,dz_fix=dz_fix
    ,close_aq=close_aq_field
    ,poro_evol=poro_evol
    ,sa_evol_1=sa_evol_1 
    ,sa_evol_2=sa_evol_2
    ,psd_bulk=psd_bulk_field
    ,psd_full=psd_full_field
    ,season=season
    )

sld_list_field=['inrt','g2']
aq_list_field = ['ca','k','mg','na']
if include_N: 
    sld_list_field.append('amnt')
    aq_list_field.append('no3')
if include_Al: 
    sld_list_field.append(alphase)
    aq_list_field.append('al')

gas_list_field = ['pco2']
exrxn_list_field = []
make_inputs.get_input_tracers(
    outdir=outdir
    ,runname=runname_field
    ,sld_list = sld_list_field
    ,aq_list = aq_list_field
    ,gas_list = gas_list_field
    ,exrxn_list = exrxn_list_field
    )
    
pr_list_field = [('inrt',1.0)]
rain_list_field = [('ca',5.0e-6)]
atm_list_field = [('pco2',3.16e-4),('po2',0.21),('pnh3',1e-50),('pn2o',1e-50)]
if include_Al:  pr_list_field.append( (alphase,1e-2) )
make_inputs.get_input_tracer_bounds(
    outdir=outdir
    ,runname=runname_field
    ,pr_list = pr_list_field
    ,rain_list=rain_list_field
    ,atm_list=atm_list_field
    )
    
filename = 'dust.in'
srcfile = '/storage/coda1/p-creinhard3/0/ykanzaki3/PyWeath/data/dust_gbasalt.in'
sld_varlist =[]
make_inputs.get_input_sld_properties(
    outdir=outdir
    ,runname=runname_field
    ,filename = filename
    # ,srcfile = srcfile
    ,sld_varlist=sld_varlist
    )
    
filename = 'cec.in'
sld_varlist = [('inrt',4,4.1,alpha) ,('g2',4,4.1,alpha) ] 
if include_Al:  sld_varlist.append( (alphase,4,4.1,alpha)  )
make_inputs.get_input_sld_properties(
    outdir=outdir
    ,runname=runname_field
    ,filename = filename
    ,sld_varlist=sld_varlist
    )
    
filename = 'OM_rain.in'
sld_varlist = [ ('g2',1) ] 
if include_N: sld_varlist.append( ('amnt',N_rain) ) 
make_inputs.get_input_sld_properties(
    outdir=outdir
    ,runname=runname_field
    ,filename = filename
    ,sld_varlist=sld_varlist
    )
    
filename = '2ndslds.in'
srcfile = '/storage/coda1/p-creinhard3/0/ykanzaki3/PyWeath/data/2ndslds_def.in'
make_inputs.get_input_sld_properties(
    outdir=outdir
    ,runname=runname_field
    ,filename = filename
    ,srcfile = srcfile
    )
    
    
# ============= buffer exp. setup ============= 
make_inputs.get_input_frame(
    outdir=outdir
    ,runname=runname_lab
    ,ztot=ztot
    ,nz=nz
    ,ttot=ttot_lab
    ,temp=temp_lab
    ,fdust=fdust_lab
    ,fdust2=fdust2
    ,taudust=taudust_lab
    ,omrain=omrain_lab
    ,zom=zom
    ,poro=poro_lab
    ,moistsrf=moistsrf_lab
    ,zwater=zwater
    ,zdust=zdust_lab
    ,w=w_lab
    ,q=q_lab
    ,p=p
    ,nstep=nstep
    ,rstrt=rstrt
    ,runid=runname_lab
    )
    
make_inputs.get_input_switches(
    outdir=outdir
    ,runname=runname_lab
    ,w_scheme=w_scheme_lab
    ,mix_scheme=mix_scheme_lab
    ,poro_iter=poro_iter_lab
    ,sldmin_lim=sldmin_lim 
    ,display=display
    ,disp_lim=disp_lim
    ,restart=restart 
    ,rough=rough_lab
    ,al_inhib=al_inhib 
    ,dt_fix=dt_fix
    ,cec_on=cec_on
    ,dz_fix=dz_fix
    ,close_aq=close_aq_lab
    ,poro_evol=poro_evol
    ,sa_evol_1=sa_evol_1 
    ,sa_evol_2=sa_evol_2
    ,psd_bulk=psd_bulk_lab
    ,psd_full=psd_full_lab
    ,season=season
    )

sld_list_lab = ['inrt','g2','cao','mgo','na2o','k2o']
aq_list_lab = ['ca','mg','na','k'] 
if include_N: 
    sld_list_lab.append('amnt')
    aq_list_lab.append('no3')
if include_Al: 
    sld_list_lab.append(alphase)
    aq_list_lab.append('no3')
if use_CaCl2:
    sld_list_lab.append('cacl2')
    aq_list_lab.append('cl')
make_inputs.get_input_tracers(
    outdir=outdir
    ,runname=runname_lab
    ,sld_list = sld_list_lab
    ,aq_list = aq_list_lab
    # ,gas_list = gas_list
    # ,exrxn_list = exrxn_list
    )
    
# pr_list = [('inrt',1.0)]
pr_list_lab = [('inrt',0.95),('g2',0.05)]
atm_list_lab = [('pco2',3.16e-4),('po2',0.21),('pnh3',1e-50),('pn2o',1e-50)]
make_inputs.get_input_tracer_bounds(
    outdir=outdir
    ,runname=runname_lab
    ,pr_list = pr_list_lab
    # ,rain_list=rain_list
    ,atm_list=atm_list_lab
    )
    
filename = 'dust.in'
srcfile = '/storage/coda1/p-creinhard3/0/ykanzaki3/PyWeath/data/dust_cao.in'
sld_varlist =[]
make_inputs.get_input_sld_properties(
    outdir=outdir
    ,runname=runname_lab
    ,filename = filename
    ,srcfile = srcfile
    # ,sld_varlist=sld_varlist
    )
    
filename = 'cec.in'
# sld_varlist = [('inrt',4,4.1)  ] 
sld_varlist = [('inrt',4,4.1,alpha) ,('g2',4,4.1,alpha) ] 
if include_Al: sld_varlist.append( (alphase,4,4.1,alpha) )
make_inputs.get_input_sld_properties(
    outdir=outdir
    ,runname=runname_lab
    ,filename = filename
    ,sld_varlist=sld_varlist
    )
    
filename = 'OM_rain.in'
sld_varlist = [ ('g2',1) ] 
make_inputs.get_input_sld_properties(
    outdir=outdir
    ,runname=runname_lab
    ,filename = filename
    ,sld_varlist=sld_varlist
    )
    
filename = 'kinspc.in'
sld_varlist = [ ('g2',0) ] 
if include_Al:  sld_varlist.append(  (alphase,0)  )
make_inputs.get_input_sld_properties(
    outdir=outdir
    ,runname=runname_lab
    ,filename = filename
    ,sld_varlist=sld_varlist
    )
    
filename = '2ndslds.in'
srcfile = '/storage/coda1/p-creinhard3/0/ykanzaki3/PyWeath/data/2ndslds_def.in'
make_inputs.get_input_sld_properties(
    outdir=outdir
    ,runname=runname_lab
    ,filename = filename
    ,srcfile = srcfile
    )



error = 1e4
tol = 1e-4

cnt = 1
res_list = []

while (error > tol):

    ymx = np.zeros(3)
    emx = np.zeros(3)
    amx = np.zeros((3,3),dtype=np.float64)
    
    facts = [1e-4]*3
    # facts = [1e-2]*3

    # command = ' {:.8f} {:.8f} {:.8f} {}'.format(cec,np.log10(kh),ca,runname)  
    
    # =========== field sim (1st) =========== 
    
    rain_list_field = [ ('ca',ca*1e-6 )]
    make_inputs.get_input_tracer_bounds(
        outdir=outdir
        ,runname=runname_field
        ,pr_list = pr_list_field
        ,rain_list=rain_list_field
        ,atm_list=atm_list_field
        )

    filename = 'cec.in'
    sld_varlist = [('inrt',cec, np.log10(kh),alpha) ,('g2',cec,np.log10(kh),alpha) ]
    if include_Al:  sld_varlist.append(  (alphase,cec,np.log10(kh),alpha)  )
    make_inputs.get_input_sld_properties(
        outdir=outdir
        ,runname=runname_field
        ,filename = filename
        ,sld_varlist=sld_varlist
        )
        
    make_inputs.get_input_frame(
        outdir=outdir
        ,runname=runname_field
        ,ztot=ztot
        ,nz=nz
        ,ttot=ttot_field
        ,temp=temp_field
        ,fdust=fdust_field
        ,fdust2=fdust2
        ,taudust=taudust_field
        ,omrain=omrain_field 
        ,zom=zom
        ,poro=poro_field
        ,moistsrf=moistsrf_field
        ,zwater=zwater
        ,zdust=zdust_field
        ,w=w_field
        ,q=q_field
        ,p=p
        ,nstep=nstep
        ,rstrt=rstrt
        ,runid=runname_field
        )

    # N_rain = 8.406375  # gN/m2/yr
    N_rain = float(sys.argv[10])  # gN/m2/yr
    N_rain = N_rain/14  # mol N/m2/yr
    N_rain = N_rain*80  # g NH4NO3/m2/yr
    N_rain = N_rain/2.  # only half is required as 1 mol NH4NO3 contains 2 moles of N
    N_rain = N_rain/omrain_field  # normalize against OC rain 
    
    filename = 'OM_rain.in'
    sld_varlist = [ ('g2',1) ] 
    if include_N: sld_varlist.append(  ('amnt',N_rain)  ) 
    make_inputs.get_input_sld_properties(
        outdir=outdir
        ,runname=runname_field
        ,filename = filename
        ,sld_varlist=sld_varlist
        )

    # run 
    os.system(outdir+runname_field+where+exename)


    phint_field = get_int_prof.get_ph_int_site(outdir,runname_field,dep_sample)
    # acint = get_int_prof.get_ac_int_site(outdir,runname,dep_sample)
    acint = get_int_prof.get_ac_int_site_v2(outdir,runname_field,dep_sample)
    
    dense_lab = get_int_prof.get_rhobulk_int_site(outdir,runname_field,dep_sample)
    sps = ['g2','inrt']
    if include_Al:   sps.append( alphase ) 
    sldwt_list = []
    for sp in sps:
        sldwt = get_int_prof.get_sldwt_int_site(outdir,runname_field,dep_sample,[sp])
        sldwt_list.append(sldwt)
    exchanger = sum(sldwt_list)
    
    omint = sldwt_list[sps.index('g2')]
    
    aqsps,btmconcs,dep = get_int_prof.get_totsave_site(outdir,runname_field,dep_sample)  # returning mol/ solid m3 depth averaged value 
    
    print(aqsps,btmconcs,dep)
    
    # =========== lab sim (1st) =========== 
    poro_lab = water_frac/(1./dense_lab+water_frac)
    # fdust_lab = ztot*(1-poro_lab)*dense_lab*1e3*exchanger/100*cec*1e-2*(1.-acint/100.)/2. * 56.1 
    
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
        fdust = ztot*(1-poro_lab)*conc* oxide_mass_list[iox]/oxide_stch_list[iox]
        fdust_list.append(fdust)
        fdust_nm_list.append(oxide_oxnm_list[iox])
        
        
    fdust_lab = fdust_list[aqsps.index('ca')]
    
    if use_CaCl2:
        cacl2_conc = 0.01
        cacl2_wt2  = 110.98
        fdust_cacl2 = ztot*poro_lab*cacl2_conc*1e3*cacl2_wt2
    
        if fdust_cacl2 > 0:
            fdust_list.append(fdust_cacl2)
            fdust_nm_list.append('cacl2')

    fdust_list = [fdust/fdust_lab  for fdust in fdust_list  ]  
    
    make_inputs.get_input_frame(
        outdir=outdir
        ,runname=runname_lab
        ,ztot=ztot
        ,nz=nz
        ,ttot=ttot_lab
        ,temp=temp_lab
        ,fdust=fdust_lab
        ,fdust2=fdust2
        ,taudust=taudust_lab
        ,omrain=omrain_lab
        ,zom=zom
        ,poro=poro_lab
        ,moistsrf=moistsrf_lab
        ,zwater=zwater
        ,zdust=zdust_lab
        ,w=w_lab
        ,q=q_lab
        ,p=p
        ,nstep=nstep
        ,rstrt=rstrt
        ,runid=runname_lab
        )

    pr_list_lab = [('inrt',sldwt_list[sps.index('inrt')]/100.),('g2',sldwt_list[sps.index('g2')]/100.)]
    if include_Al:  pr_list_lab.append(  (alphase,sldwt_list[sps.index(alphase)]/100.) )
    make_inputs.get_input_tracer_bounds(
        outdir=outdir
        ,runname=runname_lab
        ,pr_list = pr_list_lab
        # ,rain_list=rain_list
        ,atm_list=atm_list_lab
        )

    filename = 'cec.in'
    sld_varlist = [('inrt',cec, np.log10(kh),alpha),('g2',cec,np.log10(kh),alpha)  ]
    if include_Al:   sld_varlist.append( (alphase,cec,np.log10(kh),alpha) )
    make_inputs.get_input_sld_properties(
        outdir=outdir
        ,runname=runname_lab
        ,filename = filename
        ,sld_varlist=sld_varlist
        )
        
    filename = 'dust.in'
    sld_varlist =[ ( fdust_nm_list[i], fdust_list[i]) for i in range(len(fdust_nm_list)) ]
    make_inputs.get_input_sld_properties(
        outdir=outdir
        ,runname=runname_lab
        ,filename = filename
        ,sld_varlist=sld_varlist
        )
    
    os.system(outdir+runname_lab+where+exename)
    
    phint = get_int_prof.get_ph_int_site(outdir,runname_lab,dep_sample)
    
    # break

    print(phint_field,phint,acint)

    dca = facts[0]
    dca = ca*facts[0]

    # command = ' {:.8f} {:.8f} {:.8f} {}'.format(cec,np.log10(kh),ca+dca,runname)  
    # print('./'+shell + command)
    # os.system('./'+shell + command)
    
    # =========== field sim (2nd) =========== 
    
    rain_list_field = [('ca',(ca+dca)*1e-6)]
    make_inputs.get_input_tracer_bounds(
        outdir=outdir
        ,runname=runname_field
        ,pr_list = pr_list_field
        ,rain_list=rain_list_field
        ,atm_list=atm_list_field
        )
        
    # run 
    os.system(outdir+runname_field+where+exename)
    
    
    dphint_field = get_int_prof.get_ph_int_site(outdir,runname_field,dep_sample)
    # dacint = get_int_prof.get_ac_int_site(outdir,runname,dep_sample)
    dacint = get_int_prof.get_ac_int_site_v2(outdir,runname_field,dep_sample)
    # dense = get_int_prof.get_rhobulk_int_site(outdir,runname_field,dep_sample)
    dense_lab = get_int_prof.get_rhobulk_int_site(outdir,runname_field,dep_sample)
    sps = ['g2','inrt']
    if include_Al:   sps.append( alphase ) 
    sldwt_list = []
    for sp in sps:
        sldwt = get_int_prof.get_sldwt_int_site(outdir,runname_field,dep_sample,[sp])
        sldwt_list.append(sldwt)
    exchanger = sum(sldwt_list)
    
    domint = sldwt_list[sps.index('g2')]
    
    aqsps,btmconcs,dep = get_int_prof.get_totsave_site(outdir,runname_field,dep_sample)  # returning mol/ solid m3 depth averaged value 
    
    print(aqsps,btmconcs,dep)
    
    # =========== lab sim (2nd) =========== 
    poro_lab = water_frac/(1./dense_lab+water_frac)
    
    fdust_list = []
    fdust_nm_list = []

    for sp in aqsps:
        isp = aqsps.index(sp)
        iox = oxide_ctnm_list.index(sp)
        conc = btmconcs[isp]
        fdust = ztot*(1-poro_lab)*conc* oxide_mass_list[iox]/oxide_stch_list[iox]
        fdust_list.append(fdust)
        fdust_nm_list.append(oxide_oxnm_list[iox])
        
        
    fdust_lab = fdust_list[aqsps.index('ca')]
    
    if use_CaCl2:
        cacl2_conc = 0.01
        cacl2_wt2  = 110.98
        fdust_cacl2 = ztot*poro_lab*cacl2_conc*1e3*cacl2_wt2
    
        if fdust_cacl2 > 0:
            fdust_list.append(fdust_cacl2)
            fdust_nm_list.append('cacl2')

    fdust_list = [fdust/fdust_lab  for fdust in fdust_list  ]  
    
    make_inputs.get_input_frame(
        outdir=outdir
        ,runname=runname_lab
        ,ztot=ztot
        ,nz=nz
        ,ttot=ttot_lab
        ,temp=temp_lab
        ,fdust=fdust_lab
        ,fdust2=fdust2
        ,taudust=taudust_lab
        ,omrain=omrain_lab
        ,zom=zom
        ,poro=poro_lab
        ,moistsrf=moistsrf_lab
        ,zwater=zwater
        ,zdust=zdust_lab
        ,w=w_lab
        ,q=q_lab
        ,p=p
        ,nstep=nstep
        ,rstrt=rstrt
        ,runid=runname_lab
        )
        
    pr_list_lab = [('inrt',sldwt_list[sps.index('inrt')]/100.),('g2',sldwt_list[sps.index('g2')]/100.)]
    if include_Al:  pr_list_lab.append(  (alphase,sldwt_list[sps.index(alphase)]/100.) )
    make_inputs.get_input_tracer_bounds(
        outdir=outdir
        ,runname=runname_lab
        ,pr_list = pr_list_lab
        # ,rain_list=rain_list
        ,atm_list=atm_list_lab
        )
        
    filename = 'dust.in'
    sld_varlist =[ ( fdust_nm_list[i], fdust_list[i]) for i in range(len(fdust_nm_list)) ]
    make_inputs.get_input_sld_properties(
        outdir=outdir
        ,runname=runname_lab
        ,filename = filename
        ,sld_varlist=sld_varlist
        )
        
    # run 
    os.system(outdir+runname_lab+where+exename)

    dphint = get_int_prof.get_ph_int_site(outdir,runname_lab,dep_sample)
    
    
    # if phnorm_pw:       dphint_dca = (dphint_field-phint_field)/dca * ca
    # if not phnorm_pw:   dphint_dca = (dphint-phint)/dca * ca
    if phnorm_pw:       dphint_dca = (10**-dphint_field-10**-phint_field)/dca * ca
    if not phnorm_pw:   dphint_dca = (10**-dphint-10**-phint)/dca * ca
    dacint_dca = (dacint-acint)/dca * ca
    domint_dca = (domint-omint)/dca * ca

    print(dphint_dca,dacint_dca,domint_dca)


    
    # =========== field sim (3rd) =========== 
    
    dkh = facts[1]
    dkh = kh*facts[1]
    
    rain_list_field = [('ca',ca*1e-6)]
    make_inputs.get_input_tracer_bounds(
        outdir=outdir
        ,runname=runname_field
        ,pr_list = pr_list_field
        ,rain_list=rain_list_field
        ,atm_list=atm_list_field
        )

    filename = 'cec.in'
    sld_varlist = [('inrt',cec, np.log10(kh+dkh),alpha) ,('g2',cec,np.log10(kh+dkh),alpha) ] 
    if include_Al:   sld_varlist.append( (alphase,cec,np.log10(kh+dkh),alpha) )
    make_inputs.get_input_sld_properties(
        outdir=outdir
        ,runname=runname_field
        ,filename = filename
        ,sld_varlist=sld_varlist
        )
    
    # run 
    os.system(outdir+runname_field+where+exename)
    
    dphint_field = get_int_prof.get_ph_int_site(outdir,runname_field,dep_sample)
    # dacint = get_int_prof.get_ac_int_site(outdir,runname,dep_sample)
    dacint = get_int_prof.get_ac_int_site_v2(outdir,runname_field,dep_sample)
    # dense = get_int_prof.get_rhobulk_int_site(outdir,runname_field,dep_sample)
    dense_lab = get_int_prof.get_rhobulk_int_site(outdir,runname_field,dep_sample)
    sps = ['g2','inrt']
    if include_Al:   sps.append( alphase ) 
    sldwt_list = []
    for sp in sps:
        sldwt = get_int_prof.get_sldwt_int_site(outdir,runname_field,dep_sample,[sp])
        sldwt_list.append(sldwt)
    exchanger = sum(sldwt_list)
    
    domint = sldwt_list[sps.index('g2')]
    
    aqsps,btmconcs,dep = get_int_prof.get_totsave_site(outdir,runname_field,dep_sample)  # returning mol/ solid m3 depth averaged value 
    
    print(aqsps,btmconcs,dep)
    
    # =========== lab sim (3rd) =========== 
    poro_lab = water_frac/(1./dense_lab+water_frac)
    
    fdust_list = []
    fdust_nm_list = []

    for sp in aqsps:
        isp = aqsps.index(sp)
        iox = oxide_ctnm_list.index(sp)
        conc = btmconcs[isp]
        fdust = ztot*(1-poro_lab)*conc* oxide_mass_list[iox]/oxide_stch_list[iox]
        fdust_list.append(fdust)
        fdust_nm_list.append(oxide_oxnm_list[iox])
        
        
    fdust_lab = fdust_list[aqsps.index('ca')]
    
    if use_CaCl2:
        cacl2_conc = 0.01
        cacl2_wt2  = 110.98
        fdust_cacl2 = ztot*poro_lab*cacl2_conc*1e3*cacl2_wt2
    
        if fdust_cacl2 > 0:
            fdust_list.append(fdust_cacl2)
            fdust_nm_list.append('cacl2')

    fdust_list = [fdust/fdust_lab  for fdust in fdust_list  ]  
    
    make_inputs.get_input_frame(
        outdir=outdir
        ,runname=runname_lab
        ,ztot=ztot
        ,nz=nz
        ,ttot=ttot_lab
        ,temp=temp_lab
        ,fdust=fdust_lab
        ,fdust2=fdust2
        ,taudust=taudust_lab
        ,omrain=omrain_lab
        ,zom=zom
        ,poro=poro_lab
        ,moistsrf=moistsrf_lab
        ,zwater=zwater
        ,zdust=zdust_lab
        ,w=w_lab
        ,q=q_lab
        ,p=p
        ,nstep=nstep
        ,rstrt=rstrt
        ,runid=runname_lab
        )
        
    pr_list_lab = [('inrt',sldwt_list[sps.index('inrt')]/100.),('g2',sldwt_list[sps.index('g2')]/100.)]
    if include_Al:  pr_list_lab.append(  (alphase,sldwt_list[sps.index(alphase)]/100.) )
    make_inputs.get_input_tracer_bounds(
        outdir=outdir
        ,runname=runname_lab
        ,pr_list = pr_list_lab
        # ,rain_list=rain_list
        ,atm_list=atm_list_lab
        )

    filename = 'cec.in'
    sld_varlist = [('inrt',cec, np.log10(kh+dkh),alpha),('g2',cec,np.log10(kh+dkh),alpha)  ] 
    if include_Al:   sld_varlist.append( (alphase,cec,np.log10(kh+dkh),alpha) )
    make_inputs.get_input_sld_properties(
        outdir=outdir
        ,runname=runname_lab
        ,filename = filename
        ,sld_varlist=sld_varlist
        )
        
    filename = 'dust.in'
    sld_varlist =[ ( fdust_nm_list[i], fdust_list[i]) for i in range(len(fdust_nm_list)) ]
    make_inputs.get_input_sld_properties(
        outdir=outdir
        ,runname=runname_lab
        ,filename = filename
        ,sld_varlist=sld_varlist
        )
        
    # run 
    os.system(outdir+runname_lab+where+exename)

    dphint = get_int_prof.get_ph_int_site(outdir,runname_lab,dep_sample)
    
    

    # dphint_dlogkh = (dphint-phint)/dlogkh
    # dacint_dlogkh = (dacint-acint)/dlogkh
    # if phnorm_pw:       dphint_dlogkh = (dphint_field-phint_field)/dkh * kh
    # if not phnorm_pw:   dphint_dlogkh = (dphint-phint)/dkh * kh
    if phnorm_pw:       dphint_dlogkh = (10**-dphint_field-10**-phint_field)/dkh * kh
    if not phnorm_pw:   dphint_dlogkh = (10**-dphint-10**-phint)/dkh * kh
    dacint_dlogkh = (dacint-acint)/dkh * kh
    domint_dlogkh = (domint-omint)/dkh * kh

    print(dphint_dlogkh,dacint_dlogkh,domint_dlogkh)
    
    # =========== field sim (4th) =========== 
    
    domrain = facts[2]
    domrain = omrain_field*facts[2]
    
    rain_list_field = [('ca',ca*1e-6)]
    make_inputs.get_input_tracer_bounds(
        outdir=outdir
        ,runname=runname_field
        ,pr_list = pr_list_field
        ,rain_list=rain_list_field
        ,atm_list=atm_list_field
        )

    filename = 'cec.in'
    sld_varlist = [('inrt',cec, np.log10(kh),alpha) ,('g2',cec,np.log10(kh),alpha) ] 
    if include_Al:   sld_varlist.append( (alphase,cec,np.log10(kh),alpha) )
    make_inputs.get_input_sld_properties(
        outdir=outdir
        ,runname=runname_field
        ,filename = filename
        ,sld_varlist=sld_varlist
        )
        
    make_inputs.get_input_frame(
        outdir=outdir
        ,runname=runname_field
        ,ztot=ztot
        ,nz=nz
        ,ttot=ttot_field
        ,temp=temp_field
        ,fdust=fdust_field
        ,fdust2=fdust2
        ,taudust=taudust_field
        ,omrain=omrain_field + domrain
        ,zom=zom
        ,poro=poro_field
        ,moistsrf=moistsrf_field
        ,zwater=zwater
        ,zdust=zdust_field
        ,w=w_field
        ,q=q_field
        ,p=p
        ,nstep=nstep
        ,rstrt=rstrt
        ,runid=runname_field
        )

    # N_rain = 8.406375  # gN/m2/yr
    N_rain = float(sys.argv[10])  # gN/m2/yr
    N_rain = N_rain/14  # mol N/m2/yr
    N_rain = N_rain*80  # g NH4NO3/m2/yr
    N_rain = N_rain/2.  # only half is required as 1 mol NH4NO3 contains 2 moles of N
    N_rain = N_rain/(omrain_field + domrain)  # normalize against OC rain 
    
    filename = 'OM_rain.in'
    sld_varlist = [ ('g2',1) ] 
    if include_N: sld_varlist.append( ('amnt',N_rain) ) 
    make_inputs.get_input_sld_properties(
        outdir=outdir
        ,runname=runname_field
        ,filename = filename
        ,sld_varlist=sld_varlist
        )
    
    # run 
    os.system(outdir+runname_field+where+exename)
    
    dphint_field = get_int_prof.get_ph_int_site(outdir,runname_field,dep_sample)
    # dacint = get_int_prof.get_ac_int_site(outdir,runname,dep_sample)
    dacint = get_int_prof.get_ac_int_site_v2(outdir,runname_field,dep_sample)
    # dense = get_int_prof.get_rhobulk_int_site(outdir,runname_field,dep_sample)
    dense_lab = get_int_prof.get_rhobulk_int_site(outdir,runname_field,dep_sample)
    sps = ['g2','inrt']
    if include_Al:   sps.append( alphase ) 
    sldwt_list = []
    for sp in sps:
        sldwt = get_int_prof.get_sldwt_int_site(outdir,runname_field,dep_sample,[sp])
        sldwt_list.append(sldwt)
    exchanger = sum(sldwt_list)
    
    domint = sldwt_list[sps.index('g2')]
    
    aqsps,btmconcs,dep = get_int_prof.get_totsave_site(outdir,runname_field,dep_sample)  # returning mol/ solid m3 depth averaged value 
    
    print(aqsps,btmconcs,dep)
    
    # =========== lab sim (3rd) =========== 
    poro_lab = water_frac/(1./dense_lab+water_frac)
    
    fdust_list = []
    fdust_nm_list = []

    for sp in aqsps:
        isp = aqsps.index(sp)
        iox = oxide_ctnm_list.index(sp)
        conc = btmconcs[isp]
        fdust = ztot*(1-poro_lab)*conc* oxide_mass_list[iox]/oxide_stch_list[iox]
        fdust_list.append(fdust)
        fdust_nm_list.append(oxide_oxnm_list[iox])
        
        
    fdust_lab = fdust_list[aqsps.index('ca')]
    
    if use_CaCl2:
        cacl2_conc = 0.01
        cacl2_wt2  = 110.98
        fdust_cacl2 = ztot*poro_lab*cacl2_conc*1e3*cacl2_wt2
    
        if fdust_cacl2 > 0:
            fdust_list.append(fdust_cacl2)
            fdust_nm_list.append('cacl2')

    fdust_list = [fdust/fdust_lab  for fdust in fdust_list  ]
    
    make_inputs.get_input_frame(
        outdir=outdir
        ,runname=runname_lab
        ,ztot=ztot
        ,nz=nz
        ,ttot=ttot_lab
        ,temp=temp_lab
        ,fdust=fdust_lab
        ,fdust2=fdust2
        ,taudust=taudust_lab
        ,omrain=omrain_lab
        ,zom=zom
        ,poro=poro_lab
        ,moistsrf=moistsrf_lab
        ,zwater=zwater
        ,zdust=zdust_lab
        ,w=w_lab
        ,q=q_lab
        ,p=p
        ,nstep=nstep
        ,rstrt=rstrt
        ,runid=runname_lab
        )
        
    pr_list_lab = [('inrt',sldwt_list[sps.index('inrt')]/100.),('g2',sldwt_list[sps.index('g2')]/100.)]
    if include_Al:  pr_list_lab.append(  (alphase,sldwt_list[sps.index(alphase)]/100.) )
    make_inputs.get_input_tracer_bounds(
        outdir=outdir
        ,runname=runname_lab
        ,pr_list = pr_list_lab
        # ,rain_list=rain_list
        ,atm_list=atm_list_lab
        )

    filename = 'cec.in'
    sld_varlist = [('inrt',cec, np.log10(kh),alpha),('g2',cec,np.log10(kh),alpha)  ] 
    if include_Al:   sld_varlist.append( (alphase,cec,np.log10(kh),alpha) )
    make_inputs.get_input_sld_properties(
        outdir=outdir
        ,runname=runname_lab
        ,filename = filename
        ,sld_varlist=sld_varlist
        )
        
    filename = 'dust.in'
    sld_varlist =[ ( fdust_nm_list[i], fdust_list[i]) for i in range(len(fdust_nm_list)) ]
    make_inputs.get_input_sld_properties(
        outdir=outdir
        ,runname=runname_lab
        ,filename = filename
        ,sld_varlist=sld_varlist
        )
        
    # run 
    os.system(outdir+runname_lab+where+exename)

    dphint = get_int_prof.get_ph_int_site(outdir,runname_lab,dep_sample)
    
    

    # dphint_dlogkh = (dphint-phint)/dlogkh
    # dacint_dlogkh = (dacint-acint)/dlogkh
    # if phnorm_pw:       dphint_domrain = (dphint_field-phint_field)/domrain * omrain_field
    # if not phnorm_pw:   dphint_domrain = (dphint-phint)/domrain * omrain_field
    if phnorm_pw:       dphint_domrain = (10**-dphint_field-10**-phint_field)/domrain * omrain_field
    if not phnorm_pw:   dphint_domrain = (10**-dphint-10**-phint)/domrain * omrain_field
    dacint_domrain = (dacint-acint)/domrain * omrain_field
    domint_domrain = (domint-omint)/domrain * omrain_field

    print(dphint_domrain,dacint_domrain,domint_domrain)
    
    

    # updating logkh and ca by newton method 
    # f1 = phint - targetpH = 0
    # f2 = acint - acidsat = 0

    # if phnorm_pw:       ymx[0] = phint_field - targetpH
    # if not phnorm_pw:   ymx[0] = phint - targetpH
    if phnorm_pw:       ymx[0] = 10**-phint_field - 10**-targetpH
    if not phnorm_pw:   ymx[0] = 10**-phint - 10**-targetpH
    ymx[1] = 1e5*(acint - acidsat) 
    ymx[2] = omint - targetOM 

    amx[0,0] = dphint_dca
    amx[0,1] = dphint_dlogkh
    amx[0,2] = dphint_domrain
    amx[1,0] = 1e5*dacint_dca
    amx[1,1] = 1e5*dacint_dlogkh
    amx[1,2] = 1e5*dacint_domrain
    amx[2,0] = domint_dca
    amx[2,1] = domint_dlogkh
    amx[2,2] = domint_domrain

    ymx = -ymx

    dx = np.linalg.solve(amx, ymx)

    print(ca,logkh,omrain_field)
    # print(ca*np.exp(dx[0]),logkh+dx[1])
    print(ca*np.exp(dx[0]),np.log10(kh*np.exp(dx[1])),omrain_field*np.exp(dx[2]) )

    if np.isnan(ymx).any():
        print('nan detected in solution')
        error = 1e-99
    else:
        thrs = [5]*3
        # thrs = [3]*3
        # thrs = [7]*3
        if dx[0]>thrs[0]:
            ca = 1.5*ca
        elif dx[0]<-1*thrs[0]:
            ca = ca/1.5
        else:
            ca = ca*np.exp(dx[0])
        
        if dx[1] > thrs[1]:
            kh = 1.5*kh
        elif dx[1] < -1*thrs[1]:
            kh = kh/1.5
        else:
            kh = kh*np.exp(dx[1])
        
        if dx[2] > thrs[2]:
            omrain_field = 1.5*omrain_field
        elif dx[2] < -1*thrs[2]:
            omrain_field = omrain_field/1.5
        else:
            omrain_field = omrain_field*np.exp(dx[2])
        # logkh = logkh+dx[1]    

        # if phnorm_pw:       emx[0] = ymx[0]/targetpH
        # if not phnorm_pw:   emx[0] = ymx[0]/targetpH
        if phnorm_pw:       emx[0] = ymx[0]/10**-targetpH
        if not phnorm_pw:   emx[0] = ymx[0]/10**-targetpH
        
        if acidsat!=0: 
            emx[1] = ymx[1]/acidsat 
        else:
            if acint!=0: 
                emx[1] = ymx[1]/acint
            else:
                emx[1] = ymx[1]
                # emx[1] = ymx[1]*1e2
                
                
        if targetOM!=0:
            emx[2] = ymx[2]/targetOM 
        else:
            if omrain_field!=0:
                emx[2] = ymx[2]/omrain_field
            else:
                emx[2] = ymx[2]
                # emx[2] = ymx[2]*1e2

        error = np.max(np.abs(emx))
        
        
        # when log10 kh is below 1, calculation can become difficult
        # if kh < 10**1:
            # kh = 10**1
            # error = 100
        
        # if omrain_field > 2000:
            # omrain_field = 2000
            
    
    print(' ')
    print(' ')
    print(' ')
    print(' ')
    print('******* error = ',error)
    print(' ')
    print(' ')
    print(' ')
    print(' ')
    
    time.sleep(5)
    
    res_list.append([cnt,phint_field,phint,omint,acint,error,targetpH,targetOM,acidsat,ca,omrain_field,np.log10(kh)])
    cnt += 1
    
    if cnt > iter_max: break
    
    name_list = [
        'iter.'
        ,'porewater_pH[-]'
        ,'soil_pHw[-]'
        ,'OM[wt%]'
        ,'exchangeable_acidity[%CEC]'
        ,'error'
        ,'target_pH[-]'
        ,'target_OM[wt%]'
        ,'target_exchangeable_acidity[%CEC]'
        ,'Ca[uM]'
        ,'OM_rain[gC/m2/yr]'
        ,'log10(KH/Na)'
        ]
    for runname in [runname_field,runname_lab]:
        # np.savetxt(outdir + runname + where + 'iteration_tmp.res',np.array(res_list))
        dst = outdir + runname + where + 'iteration_tmp.res'

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


# name_list = [
    # 'iter.'
    # ,'porewater_pH[-]'
    # ,'soil_pHw[-]'
    # ,'error'
    # ]
name_list = [
    'iter.'
    ,'porewater_pH[-]'
    ,'soil_pHw[-]'
    ,'OM[wt%]'
    ,'exchangeable_acidity[%CEC]'
    ,'error'
    ,'target_pH[-]'
    ,'target_OM[wt%]'
    ,'target_exchangeable_acidity[%CEC]'
    ,'Ca[uM]'
    ,'OM_rain[gC/m2/yr]'
    ,'log10(KH/Na)'
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