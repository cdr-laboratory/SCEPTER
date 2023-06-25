import os
import numpy as np
import get_int_prof_time_dep
import get_int_prof_time
import get_int_prof
import make_inputs
import get_inputs
import time
import sys
import shutil

def calc_soilpH(inputlist):
    
    idep = inputlist[0]
    itime = inputlist[1]

    include_N = False
    include_N = True

    include_Al = False
    # include_Al = True

    use_CaCl2 = False
    # use_CaCl2 = True
    
    include_pw = False
    include_pw = True

        
    water_frac = 5.
    water_frac = 1.

    dep_sample = 0.15
    # dep_sample = 0.18
    # dep_sample = 0.25

    ttot_lab=1000
    temp_lab=25


    outdir = '../scepter_output/'

    exename = 'scepter'
    to = ' '
    where = '/'

    labstr = '_lab'
    if use_CaCl2: labstr = '_cacl2'
    if include_pw: labstr += '_wtot'
    else: labstr += '_w'
    labstr += '_'+str(int(water_frac))
    labstr += '_'+str(int(idep))
    labstr += '_'+str(int(itime))

    runname_field   = sys.argv[1] #'US_cropland_311_sph_N_spintuneup_field'
    runname_lab     = runname_field+labstr

        
    # if not os.path.exists( outdir + runname_lab) : os.system('mkdir -p ' + outdir + runname_lab)

    src = outdir + runname_field
    dst = outdir + runname_lab

    if not os.path.exists(dst): 
        shutil.copytree(src, dst)
    else:
        shutil.rmtree(dst)
        shutil.copytree(src, dst)

    # """
    runname = runname_field
    ztot,nz,ttot,temp,fdust,fdust2,taudust,omrain,zom,poro,moistsrf,zwater,zdust,w,q,p,nstep,rstrt,runid \
        = get_inputs.get_input_frame(outdir,runname)
    print(
        ztot,nz,ttot,temp,fdust,fdust2,taudust,omrain,zom,poro,moistsrf,zwater,zdust,w,q,p,nstep,rstrt,runid
        )
        
    w_scheme,mix_scheme,poro_iter,sldmin_lim,display,disp_lim,restart,rough,al_inhib,dt_fix \
        ,cec_on,dz_fix,close_aq,poro_evol,sa_evol_1,sa_evol_2,psd_bulk,psd_full,season \
        = get_inputs.get_input_switches(outdir,runname)
    print(
        w_scheme,mix_scheme,poro_iter,sldmin_lim,display,disp_lim,restart,rough,al_inhib,dt_fix 
        ,cec_on,dz_fix,close_aq,poro_evol,sa_evol_1,sa_evol_2,psd_bulk,psd_full,season 
        )
        
    sld_list,aq_list,gas_list,exrxn_list \
        = get_inputs.get_input_tracers(outdir,runname)
    print(
        sld_list,aq_list,gas_list,exrxn_list
        )
        
        
    # """
    # get data from field run 

    if include_pw:
        try:
            aqsps,btmconcs,dep,time = get_int_prof_time_dep.get_totsave_site(outdir,runname_field,idep,itime)  # returning mol/ solid m3 depth averaged value 
        except:
            print('data not available: ', runname_field)
            return
    else:
        try:
            aqsps,btmconcs,dep,time = get_int_prof_time_dep.get_adssave_site(outdir,runname_field,idep,itime)  # returning mol/ solid m3 depth averaged value 
        except:
            print('data not available: ', runname_field)
            return

    print(aqsps,btmconcs,dep)

    phint_field = get_int_prof_time_dep.get_ph_int_site(outdir,runname_field,idep,itime)
    # acint = get_int_prof.get_ac_int_site(outdir,runname,dep_sample)
    acint = get_int_prof_time_dep.get_ac_int_site_v2(outdir,runname_field,idep,itime)

    dense_lab = get_int_prof_time_dep.get_rhobulk_int_site(outdir,runname_field,idep,itime)
    sps = list(sld_list)
    sldwt_list = []
    for sp in sps:
        sldwt = get_int_prof_time_dep.get_sldwt_int_site(outdir,runname_field,idep,[sp],itime)
        sldwt_list.append(sldwt)



    # """

    # setting up conditions for lab experiments 
        
    poro_lab = water_frac/(1./dense_lab+water_frac)
    oxide_ctnm_list = ['ca','mg','na','k'] 
    oxide_oxnm_list = ['cao','mgo','na2o','k2o'] 
    oxide_stch_list = [1,1,2,2] 
    oxide_mass_list = [56.1 ,40.3, 62, 94.2]
    
    if include_pw:
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
        if sp  not in oxide_ctnm_list: continue
        isp = aqsps.index(sp)
        iox = oxide_ctnm_list.index(sp)
        conc = btmconcs[isp]
        fdust = ztot*(1-poro_lab)*conc* oxide_mass_list[iox]/oxide_stch_list[iox]
        fdust_list.append(fdust)
        fdust_nm_list.append(oxide_oxnm_list[iox])
        
    # print (fdust_nm_list)
    # print (fdust_list)
    # fdust_lab = fdust_list[aqsps.index('ca')]
    fdust_lab = fdust_list[fdust_nm_list.index('cao')]

    if use_CaCl2:
        cacl2_conc = 0.01
        cacl2_wt2  = 110.98
        fdust_cacl2 = ztot*poro_lab*cacl2_conc*1e3*cacl2_wt2

        if fdust_cacl2 > 0:
            fdust_list.append(fdust_cacl2)
            fdust_nm_list.append('cacl2')
            oxide_ctnm_list.append('cl')

    fdust_list = [fdust/fdust_lab  for fdust in fdust_list  ]  


    sld_list_lab = [sld for sld in sld_list if sld not in fdust_nm_list ]
    sld_list_lab.extend(fdust_nm_list)
    aq_list_lab = [aq for aq in aq_list if aq not in oxide_ctnm_list ]
    aq_list_lab.extend(oxide_ctnm_list)

    # """
    make_inputs.get_input_frame(
        outdir=outdir
        ,runname=runname_lab
        ,ztot=ztot
        ,nz=nz
        ,ttot=ttot_lab
        ,temp=temp_lab
        ,fdust=fdust_lab
        ,fdust2=fdust2
        ,taudust=0.01
        ,omrain=0
        ,zom=zom
        ,poro=poro_lab
        ,moistsrf=1.0
        ,zwater=zwater
        ,zdust=0.15
        ,w=0
        ,q=0
        ,p=p
        ,nstep=10
        ,rstrt=rstrt
        ,runid=runname_lab
        )
        

    w_scheme_lab=0
    mix_scheme_lab=0 
    poro_iter_lab='true' 
    rough_lab      ='false'
    sld_fix_lab='true'
    psd_bulk_lab='false'
    psd_full_lab ='false'

    make_inputs.get_input_switches(
        outdir=outdir
        ,runname=runname_lab
        ,w_scheme=w_scheme_lab
        ,mix_scheme=mix_scheme_lab
        ,poro_iter=poro_iter_lab
        ,sldmin_lim=sldmin_lim 
        ,display=display
        ,disp_lim=disp_lim
        ,restart='false' 
        ,rough=rough_lab
        ,al_inhib=al_inhib 
        ,dt_fix=dt_fix
        ,cec_on=cec_on
        ,dz_fix=dz_fix
        ,close_aq='true'
        ,poro_evol=poro_evol
        ,sa_evol_1=sa_evol_1 
        ,sa_evol_2=sa_evol_2
        ,psd_bulk='false'
        ,psd_full='false'
        ,season=season
        )

    make_inputs.get_input_tracers(
        outdir=outdir
        ,runname=runname_lab
        ,sld_list = sld_list_lab
        ,aq_list = aq_list_lab
        )
    # print(sps)
    # print(sldwt_list)
    pr_list_lab = [(sp,sldwt_list[sps.index(sp)]/100.) for sp in sps if sp not in oxide_oxnm_list]
    make_inputs.get_input_tracer_bounds(
        outdir=outdir
        ,runname=runname_lab
        ,pr_list = pr_list_lab
        )
        
    filename = 'dust.in'
    sld_varlist =[ ( fdust_nm_list[i], fdust_list[i]) for i in range(len(fdust_nm_list)) ]
    make_inputs.get_input_sld_properties(
        outdir=outdir
        ,runname=runname_lab
        ,filename = filename
        ,sld_varlist=sld_varlist
        )


    filename = 'kinspc.in'
    sld_varlist = [ (sld,0) for sld in sld_list if sld not in oxide_oxnm_list ] 
    make_inputs.get_input_sld_properties(
        outdir=outdir
        ,runname=runname_lab
        ,filename = filename
        ,sld_varlist=sld_varlist
        )

    os.system(outdir+runname_lab+where+exename)
    
    phint_lab = get_int_prof.get_ph_int_site(outdir,runname_lab,ztot)
    
    shutil.rmtree(dst)
    
    return dep,phint_lab,time


def main():
    
    n_dep = 30
    n_time = 20
    
    inputlists = [[idep,itime+1] for itime in range(n_time) for idep in range(n_dep)  ]
    
    dep,phint_lab,time = calc_soilpH([0,1])
    res_list = []
    # n_sample = n_time
    # for i in range(n_sample):
    for inputlist in inputlists:
        dep,phint_lab,time = calc_soilpH(inputlist)
        res_list.append([dep,phint_lab,time])
    
    field_dir = sys.argv[1]
    
    outdir = '../scepter_output/'
    
    for itime in range(n_time):
        file = outdir + field_dir + '/prof/soil_ph_{:03d}.txt'.format(itime+1)
        with open(file,'w') as f:
            f.write('z\tsoil_ph\ttime\n')
            for idep in range(n_dep):
                f.write('{:.6e}\t{:.6e}\t{:.6e}\n'.format(res_list[idep+n_dep*itime][0],res_list[idep+n_dep*itime][1],res_list[idep+n_dep*itime][2]))
    
    print(res_list)
        
        
        
if __name__ == '__main__':
    main()