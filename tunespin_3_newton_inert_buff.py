import os
import numpy as np
import get_int_prof
import make_inputs
import time


shell='inert_tune_spinup_newton.sh'

cec=4
targetpH = 5.8
basesat = 75
acidsat = 25
targetOM = 5
z = 2 
alpha = 3.4
kh = 10**3.5
logkh = np.log10(kh)
ca=35 # uM

dense_lab = 2.59296482412060000 # only inertphase

dep_sample = 0.18


outdir = '/storage/scratch1/0/ykanzaki3/pyweath_output/'
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
runname_field   = simid+'_spintuneup_field'
runname_lab     = simid+'_spintuneup_lab'

# compile 
exename = 'weathering'
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
moistsrf_lab=1.0
zwater=100000
zdust_field=0.18
zdust_lab=0.15
w_field=100e-5
w_lab=0
q_field=1200e-3
q_lab=0
p=10e-6
nstep=10
rstrt='self'
runid_field=runname_field
runid_lab=runname_lab

N_rain = 0  # gN/m2/yr
N_rain = 8.406375  # gN/m2/yr
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
precalc='false'
precalc='true'
dz_fix='true'
sld_fix_field='false'
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
sld_fix_lab='true'
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
    ,precalc=precalc
    ,dz_fix=dz_fix
    ,sld_fix=sld_fix_field
    ,poro_evol=poro_evol
    ,sa_evol_1=sa_evol_1 
    ,sa_evol_2=sa_evol_2
    ,psd_bulk=psd_bulk_field
    ,psd_full=psd_full_field
    ,season=season
    )

sld_list_field=['inrt','g2']
aq_list_field = ['ca','k','mg','na']
# sld_list_field =['inrt','g2','amnt']
# aq_list_field = ['ca','k','mg','na','no3']
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
sld_varlist = [('inrt',4,4.1) ,('g2',4,4.1) ] 
make_inputs.get_input_sld_properties(
    outdir=outdir
    ,runname=runname_field
    ,filename = filename
    ,sld_varlist=sld_varlist
    )
    
filename = 'OM_rain.in'
sld_varlist = [ ('g2',1) ] 
# sld_varlist = [ ('g2',1), ('amnt',N_rain) ] 
make_inputs.get_input_sld_properties(
    outdir=outdir
    ,runname=runname_field
    ,filename = filename
    ,sld_varlist=sld_varlist
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
    ,precalc=precalc
    ,dz_fix=dz_fix
    ,sld_fix=sld_fix_lab
    ,poro_evol=poro_evol
    ,sa_evol_1=sa_evol_1 
    ,sa_evol_2=sa_evol_2
    ,psd_bulk=psd_bulk_lab
    ,psd_full=psd_full_lab
    ,season=season
    )

make_inputs.get_input_tracers(
    outdir=outdir
    ,runname=runname_lab
    # ,sld_list = ['inrt','cao']
    ,sld_list = ['inrt','cao','g2']
    ,aq_list = ['ca']
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
sld_varlist = [('inrt',4,4.1) ,('g2',4,4.1) ] 
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
make_inputs.get_input_sld_properties(
    outdir=outdir
    ,runname=runname_lab
    ,filename = filename
    ,sld_varlist=sld_varlist
    )



error = 1e4
tol = 1e-4

cnt = 1
res_list = []

iter_max = 100

while (error > tol):

    ymx = np.zeros(3)
    emx = np.zeros(3)
    amx = np.zeros((3,3),dtype=np.float64)
    
    facts = [1e-4]*3

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
    sld_varlist = [('inrt',cec, np.log10(kh)) ,('g2',cec,np.log10(kh)) ] 
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

    N_rain = 8.406375  # gN/m2/yr
    N_rain = N_rain/14  # mol N/m2/yr
    N_rain = N_rain*80  # g NH4NO3/m2/yr
    N_rain = N_rain/2.  # only half is required as 1 mol NH4NO3 contains 2 moles of N
    N_rain = N_rain/omrain_field  # normalize against OC rain 
    
    filename = 'OM_rain.in'
    sld_varlist = [ ('g2',1) ] 
    # sld_varlist = [ ('g2',1), ('amnt',N_rain) ] 
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
    sldwt_list = []
    for sp in sps:
        sldwt = get_int_prof.get_sldwt_int_site(outdir,runname_field,dep_sample,[sp])
        sldwt_list.append(sldwt)
    exchanger = sum(sldwt_list)
    
    omint = sldwt_list[sps.index('g2')]
    
    # =========== lab sim (1st) =========== 
    poro_lab = 5./(1./dense_lab+5.)
    fdust_lab = ztot*(1-poro_lab)*dense_lab*1e3*exchanger/100*cec*1e-2*(1.-acint/100.)/2. * 56.1 
    
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
    make_inputs.get_input_tracer_bounds(
        outdir=outdir
        ,runname=runname_lab
        ,pr_list = pr_list_lab
        # ,rain_list=rain_list
        ,atm_list=atm_list_lab
        )

    filename = 'cec.in'
    sld_varlist = [('inrt',cec, np.log10(kh)),('g2',cec,np.log10(kh))  ] 
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
    sldwt_list = []
    for sp in sps:
        sldwt = get_int_prof.get_sldwt_int_site(outdir,runname_field,dep_sample,[sp])
        sldwt_list.append(sldwt)
    exchanger = sum(sldwt_list)
    
    domint = sldwt_list[sps.index('g2')]
    
    # =========== lab sim (2nd) =========== 
    poro_lab = 5./(1./dense_lab+5.)
    # fdust_lab = ztot*(1-poro_lab)*dense_lab*1e3*cec*1e-2*(1.-dacint/100.)/2. * 56.1 
    fdust_lab = ztot*(1-poro_lab)*dense_lab*1e3*exchanger/100*cec*1e-2*(1.-dacint/100.)/2. * 56.1 
    
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
    make_inputs.get_input_tracer_bounds(
        outdir=outdir
        ,runname=runname_lab
        ,pr_list = pr_list_lab
        # ,rain_list=rain_list
        ,atm_list=atm_list_lab
        )
    # run 
    os.system(outdir+runname_lab+where+exename)

    dphint = get_int_prof.get_ph_int_site(outdir,runname_lab,dep_sample)
    
    
    # dphint_dca = (dphint_field-phint_field)/dca * ca
    # dphint_dca = (10**-dphint_field-10**-phint_field)/dca * ca
    # dphint_dca = (dphint-phint)/dca * ca
    dphint_dca = (10**-dphint-10**-phint)/dca * ca
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
    sld_varlist = [('inrt',cec, np.log10(kh+dkh)) ,('g2',cec,np.log10(kh+dkh)) ] 
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
    sldwt_list = []
    for sp in sps:
        sldwt = get_int_prof.get_sldwt_int_site(outdir,runname_field,dep_sample,[sp])
        sldwt_list.append(sldwt)
    exchanger = sum(sldwt_list)
    
    domint = sldwt_list[sps.index('g2')]
    
    # =========== lab sim (3rd) =========== 
    poro_lab = 5./(1./dense_lab+5.)

    # fdust_lab = ztot*(1-poro_lab)*dense_lab*1e3*cec*1e-2*(1.-dacint/100.)/2. * 56.1 
    fdust_lab = ztot*(1-poro_lab)*dense_lab*1e3*exchanger/100*cec*1e-2*(1.-dacint/100.)/2. * 56.1 
    
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
    make_inputs.get_input_tracer_bounds(
        outdir=outdir
        ,runname=runname_lab
        ,pr_list = pr_list_lab
        # ,rain_list=rain_list
        ,atm_list=atm_list_lab
        )

    filename = 'cec.in'
    sld_varlist = [('inrt',cec, np.log10(kh+dkh)),('g2',cec,np.log10(kh+dkh))  ] 
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
    # dphint_dlogkh = (dphint_field-phint_field)/dkh * kh
    # dphint_dlogkh = (10**-dphint_field-10**-phint_field)/dkh * kh
    # dphint_dlogkh = (dphint-phint)/dkh * kh
    dphint_dlogkh = (10**-dphint-10**-phint)/dkh * kh
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
    sld_varlist = [('inrt',cec, np.log10(kh)) ,('g2',cec,np.log10(kh)) ] 
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

    N_rain = 8.406375  # gN/m2/yr
    N_rain = N_rain/14  # mol N/m2/yr
    N_rain = N_rain*80  # g NH4NO3/m2/yr
    N_rain = N_rain/2.  # only half is required as 1 mol NH4NO3 contains 2 moles of N
    N_rain = N_rain/(omrain_field + domrain)  # normalize against OC rain 
    
    filename = 'OM_rain.in'
    sld_varlist = [ ('g2',1) ] 
    # sld_varlist = [ ('g2',1), ('amnt',N_rain) ] 
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
    sldwt_list = []
    for sp in sps:
        sldwt = get_int_prof.get_sldwt_int_site(outdir,runname_field,dep_sample,[sp])
        sldwt_list.append(sldwt)
    exchanger = sum(sldwt_list)
    
    domint = sldwt_list[sps.index('g2')]
    
    # =========== lab sim (3rd) =========== 
    poro_lab = 5./(1./dense_lab+5.)

    # fdust_lab = ztot*(1-poro_lab)*dense_lab*1e3*cec*1e-2*(1.-dacint/100.)/2. * 56.1 
    fdust_lab = ztot*(1-poro_lab)*dense_lab*1e3*exchanger/100*cec*1e-2*(1.-dacint/100.)/2. * 56.1 
    
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
    make_inputs.get_input_tracer_bounds(
        outdir=outdir
        ,runname=runname_lab
        ,pr_list = pr_list_lab
        # ,rain_list=rain_list
        ,atm_list=atm_list_lab
        )

    filename = 'cec.in'
    sld_varlist = [('inrt',cec, np.log10(kh+dkh)),('g2',cec,np.log10(kh+dkh))  ] 
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
    # dphint_domrain = (dphint_field-phint_field)/domrain * omrain_field
    # dphint_domrain = (10**-dphint_field-10**-phint_field)/domrain * omrain_field
    # dphint_domrain = (dphint-phint)/domrain * omrain_field
    dphint_domrain = (10**-dphint-10**-phint)/domrain * omrain_field
    dacint_domrain = (dacint-acint)/domrain * omrain_field
    domint_domrain = (domint-omint)/domrain * omrain_field

    print(dphint_domrain,dacint_domrain,domint_domrain)
    
    

    # updating logkh and ca by newton method 
    # f1 = phint - targetpH = 0
    # f2 = acint - acidsat = 0

    # ymx[0] = phint_field - targetpH
    # ymx[0] = 10**-phint_field - 10**-targetpH
    # ymx[0] = phint - targetpH
    ymx[0] = 10**-phint - 10**-targetpH
    ymx[1] = acint - acidsat 
    ymx[2] = omint - targetOM 

    amx[0,0] = dphint_dca
    amx[0,1] = dphint_dlogkh
    amx[0,2] = dphint_domrain
    amx[1,0] = dacint_dca
    amx[1,1] = dacint_dlogkh
    amx[1,2] = dacint_domrain
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

        # emx[0] = ymx[0]/targetpH
        emx[0] = ymx[0]/10**-targetpH
        emx[1] = ymx[1]/acidsat 
        emx[2] = ymx[2]/targetOM 

        error = np.max(np.abs(emx))
    
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
    
    res_list.append([cnt,phint_field,phint,error])
    cnt += 1
    
    if cnt > iter_max: break
    
    # break
    



name_list = [
    'iter.'
    ,'porewater_pH[-]'
    ,'soil_pHw[-]'
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