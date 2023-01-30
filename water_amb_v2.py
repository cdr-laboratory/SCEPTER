import os
import numpy as np
import get_int_prof
import make_inputs
import time
import shutil
import copy
import sys

def labsim_aq_close(runname_org):

    dep_sample = 0.5
    # dep_sample = 0.0
            
    exename = 'scepter'
    to = ' '
    where = '/'
            
    srcdir = '/storage/scratch1/0/ykanzaki3/scepter_archive/'
    dstdir = '/storage/coda1/p-creinhard3/0/ykanzaki3/scepter_output/'

    # runname_org   = 'glob_ag_site'+str(i)+'_lith_fick_noiter_test_psdfull_RCP85'
    runname_lab   = runname_org+'_water_'+str(dep_sample).replace('.','p')

    if 'archive' in srcdir:
        if os.path.exists(dstdir + runname_lab):  
            shutil.rmtree(dstdir + runname_lab)
            os.system('tar xzf ' + srcdir + runname_org + '.tar.gz -C ' + dstdir )
            os.system('mv ' + dstdir + runname_org + to + dstdir + runname_lab)
        else:
            os.system('tar xzf ' + srcdir + runname_org + '.tar.gz -C ' + dstdir )
            os.system('mv ' + dstdir + runname_org + to + dstdir + runname_lab)

    else:
        src = srcdir + runname_org
        dst = dstdir + runname_lab

        if not os.path.exists(dst): 
            shutil.copytree(src, dst)
        else:
            shutil.rmtree(dst)
            shutil.copytree(src, dst)

    os.system('cp ' + exename + to + dstdir + runname_lab)

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
    poro_lab=0.9999
    # poro_lab=0.9999999
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
    runid_field=runname_org
    runid_lab=runname_lab

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
    # cec_on='true'
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


    # sps,btmconcs,dep = get_int_prof.get_btmwater_site(outdir,runname_org)
    sps,btmconcs,dep = get_int_prof.get_water_site_OLD(dstdir,runname_lab,dep_sample)

    print(sps,btmconcs,dep)

    incl_SO4 = False

    if 'so4' in sps:
        print('so4 in the original run')
        incl_SO4 = True

    DIC_btm_OLD,dep,i_save = get_int_prof.get_DIC_save_get_site_OLD(dstdir,runname_lab,dep_sample)

    time.sleep(5)

    # oxide_ctnm_list = ['ca','mg','na','k','no3'] 
    # oxide_oxnm_list = ['cao','mgo','na2o','k2o','amnt'] 
    # oxide_stch_list = [1,1,2,2,2] 
    # oxide_mass_list = [56.1 ,40.3, 62, 94.2, 80]

    oxide_ctnm_list = ['ca','mg','na','k','al','si'] 
    oxide_oxnm_list = ['cao','mgo','na2o','k2o','al2o3','sio2'] 
    oxide_stch_list = [1,1,2,2,2,1] 
    oxide_mass_list = [56.079 ,40.304, 61.979, 94.195, 101.962, 60.085]


    if incl_SO4:
        oxide_ctnm_list.append('so4')
        oxide_oxnm_list.append('caso4')
        oxide_stch_list.append(1)
        oxide_mass_list.append(136.138)
        
    fdust_list = []
    fdust_nm_list = []

    for sp in sps:
        if sp  not in oxide_ctnm_list: continue
        isp = sps.index(sp)
        iox = oxide_ctnm_list.index(sp)
        conc = btmconcs[isp]
        fdust = ztot*poro_lab*1e3*conc* oxide_mass_list[iox]/oxide_stch_list[iox]
        fdust_list.append(fdust)
        fdust_nm_list.append(oxide_oxnm_list[iox])
        
    if incl_SO4:
        ica = oxide_ctnm_list.index('ca')
        iso4 = oxide_ctnm_list.index('so4')
        if fdust_list[ica]/oxide_mass_list[ica] > fdust_list[iso4]/oxide_mass_list[iso4]: 
            # when mole flux of cao is larger than those of caso4, cao mole flux must be reduced to account for Ca mole flux from caso4
            fdust_list[ica] = ( fdust_list[ica]/oxide_mass_list[ica] - fdust_list[iso4]/oxide_mass_list[iso4] ) * oxide_mass_list[ica]
        else:
            # when mole flux of caso4 is larger than those of cao, cao mole flux must be reduced to account for Ca mole flux from caso4
            # ca is going to be added too much 
            # and cannot be subtracted ....
            print ('Too much SO4')
    
            os.system('tar -czf '+dstdir+runname_lab+'.tar.gz ' + dstdir+runname_lab)
            os.system('mv '+dstdir+runname_lab+'.tar.gz ' + srcdir)
            shutil.rmtree(dstdir+runname_lab)
            return np.nan

    fdust_lab = fdust_list[oxide_ctnm_list.index('ca')]

    fdust_list = [fdust/fdust_lab  for fdust in fdust_list  ]  


    # ============= buffer exp. setup ============= 
    make_inputs.get_input_frame(
        outdir=dstdir
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
        outdir=dstdir
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

    sld_list = copy.copy(oxide_oxnm_list)
    sld_list.append('inrt')
    print(sld_list)
    make_inputs.get_input_tracers(
        outdir=dstdir
        ,runname=runname_lab
        # ,sld_list = ['inrt','cao']
        ,sld_list = sld_list
        ,aq_list = oxide_ctnm_list
        # ,gas_list = gas_list
        # ,exrxn_list = exrxn_list
        )
        
    pr_list_lab = [('inrt',1.0)]
    # pr_list_lab = [('inrt',0.95),('g2',0.05)]
    atm_list_lab = [('pco2',3.16e-4),('po2',0.21),('pnh3',1e-50),('pn2o',1e-50)]
    make_inputs.get_input_tracer_bounds(
        outdir=dstdir
        ,runname=runname_lab
        ,pr_list = pr_list_lab
        # ,rain_list=rain_list
        ,atm_list=atm_list_lab
        )
        
    filename = 'dust.in'
    srcfile = '/storage/coda1/p-creinhard3/0/ykanzaki3/PyWeath/data/dust_cao.in'
    sld_varlist =[ ( fdust_nm_list[i], fdust_list[i]) for i in range(len(fdust_nm_list)) ]
    make_inputs.get_input_sld_properties(
        outdir=dstdir
        ,runname=runname_lab
        ,filename = filename
        # ,srcfile = srcfile
        ,sld_varlist=sld_varlist
        )
        
    filename = 'cec.in'
    # sld_varlist = [('inrt',4,4.1)  ] 
    # sld_varlist = [('inrt',0,0)  ] 
    sld_varlist = [('inrt',0,0,3.4)  ] 
    # sld_varlist = [('inrt',4,4.1) ,('g2',4,4.1) ] 
    make_inputs.get_input_sld_properties(
        outdir=dstdir
        ,runname=runname_lab
        ,filename = filename
        ,sld_varlist=sld_varlist
        )
        
    filename = 'OM_rain.in'
    sld_varlist = [ ('g2',1) ] 
    make_inputs.get_input_sld_properties(
        outdir=dstdir
        ,runname=runname_lab
        ,filename = filename
        ,sld_varlist=sld_varlist
        )
        
    filename = 'kinspc.in'
    sld_varlist = [ ('g2',0) ] 
    make_inputs.get_input_sld_properties(
        outdir=dstdir
        ,runname=runname_lab
        ,filename = filename
        ,sld_varlist=sld_varlist
        )
        
        
    os.system(dstdir+runname_lab+where+exename)


    DIC_btm_NEW,dep = get_int_prof.get_DIC_save_get_site_NEW(dstdir,runname_lab,dep_sample,i_save)

    print(DIC_btm_NEW/DIC_btm_OLD)
    
    os.system('tar -czf '+dstdir+runname_lab+'.tar.gz ' + dstdir+runname_lab)
    os.system('mv '+dstdir+runname_lab+'.tar.gz ' + srcdir)
    shutil.rmtree(dstdir+runname_lab)
    
    return DIC_btm_NEW/DIC_btm_OLD


def main():
    # runname_org   = 'glob_ag_site1_lith_fick_noiter_test_psdfull_RCP85'
    runname_org   = sys.argv[1]
    labsim_aq_close(runname_org)
    
if __name__ == '__main__':
    main()