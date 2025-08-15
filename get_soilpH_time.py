import os
import numpy as np
import get_int_prof
import make_inputs
import get_inputs
import print_control
import time
import sys
import shutil
import inspect
import subprocess
import platform
import re
from sys import exit
import json

disply_glbl = print_control.get_global_display()

def my_rmdir(dir_to_remove):
    
    for _ in range(3):  # Try up to 3 times
        try:
            shutil.rmtree(dir_to_remove)
            break
        except Exception as e:
            print(f"Retrying deletion: {e}")
            time.sleep(1)


def calc_soilpH(outdir, runname_field, dep_sample, itime, **kwargs):
    
    include_N           = kwargs.get('include_N',           True)
    include_Cl          = kwargs.get('include_Cl',          True)
    include_Al          = kwargs.get('include_Al',          False)
    include_DIC         = kwargs.get('include_DIC',         False)
    use_CaCl2           = kwargs.get('use_CaCl2',           False)
    use_local_storage   = kwargs.get('use_local_storage',   True)
    save_lab_dir        = kwargs.get('save_lab_dir',        False)
    use_Sikora_buffer   = kwargs.get('use_Sikora_buffer',   False)
    set_ExTimeLimit     = kwargs.get('set_ExTimeLimit',     True)
    add_salt            = kwargs.get('add_salt',            False)
    add_NO3             = kwargs.get('add_NO3',             False)
    no_depth_integral   = kwargs.get('no_depth_integral',   False)
    re_compile          = kwargs.get('re_compile',          False)
    water_frac          = kwargs.get('water_frac',          1)
    cacl2_conc          = kwargs.get('cacl2_conc',          0.01)
    salt_conc           = kwargs.get('salt_conc',           0.01)
    NO3_conc            = kwargs.get('NO3_conc',            0.01)
    # cl_salt             = kwargs.get('cl_salt',             'hcl')
    cl_salt             = kwargs.get('cl_salt',             'nacl')
    salt_sp             = kwargs.get('salt_sp',             'nacl')
    Ex_TimeLim          = kwargs.get('Ex_TimeLim',          20)
    labstr              = kwargs.get('labstr',              'lab')
    make_runlogfile     = kwargs.get('make_runlogfile',     False)
    sub_as_a_job        = kwargs.get('sub_as_a_job',        False)
    exename_src         = kwargs.get('exename_src',         'scepter_test')
        
    # water_frac = 5.
    # water_frac = 1.
    # water_frac = 2.5
    
    include_pw = False
    include_pw = True
    
    if use_Sikora_buffer:
        if use_CaCl2: use_CaCl2 = False
        water_frac = 2.
        # water_frac = 1.

    # dep_sample = 0.15
    # dep_sample = 0.18
    # dep_sample = 0.25
    
    flag_error = False
    
    if salt_sp not in ['nacl','kcl','cacl2']: 
        exit('salt_sp need be nacl, kcl or cacl2: stop')
    
    if cl_salt not in ['nacl','kcl','cacl2','mgcl2','hcl']: 
        exit('cl_salt need be nacl, kcl, cacl2, mgcl2, hcl: stop')
    else:
        if cl_salt=='nacl':
            cl_salt_wt = 58.443
            cl_salt_st = 1
        elif cl_salt=='kcl':
            cl_salt_wt = 74.551
            cl_salt_st = 1
        elif cl_salt=='cacl2':
            cl_salt_wt = 110.986
            cl_salt_st = 2
        elif cl_salt=='mgcl2':
            cl_salt_wt = 95.211
            cl_salt_st = 2
        elif cl_salt=='hcl':
            cl_salt_wt = 36.46
            cl_salt_st = 1
    
    
    ztot_lab=0.5
    ttot_lab=1000
    
    ztot_lab=0.05
    ttot_lab=100
    
    temp_lab=25


    # outdir = '../scepter_output/'

    exename = 'scepter'
    # exename_src = 'scepter_test'
    to = ' '
    where = '/'

    # labstr = '_labtmp'
    # labstr = '_labtmppy'
    
    outdir_tmp = outdir
    
    if use_local_storage: outdir_tmp = os.environ['TMPDIR'] + '/scepter_output/'

    # runname_field   = 'US_cropland_311_sph_N_spintuneup_field'
    if 'field' in runname_field:
        runname_lab = runname_field.replace('field',labstr)
    else:
        runname_lab= runname_field+'_'+labstr

        
    # if not os.path.exists( outdir + runname_lab) : os.system('mkdir -p ' + outdir + runname_lab)

    src = outdir + runname_field
    if not os.path.exists(src): 
        if use_local_storage: 
            outdir = outdir_tmp
            src = outdir + runname_field
            if not os.path.exists(src):
                sys.exit(f"error source directory ({src}) could not be found")
        else:
            sys.exit(f"error source directory ({src}) could not be found")
    
    
    dst = outdir_tmp + runname_lab

    if not os.path.exists(dst): 
        shutil.copytree(src, dst, dirs_exist_ok=True)
    else:
        print('field = {}'.format(src))
        print('lab = {}'.format(dst))
        # shutil.rmtree(dst)
        my_rmdir(dst)
        shutil.copytree(src, dst, dirs_exist_ok=True)
    
    # remove results directories (this is safe as long as experimrnts are not restart from self) | 10-27-2023 added 
    shutil.rmtree(dst+'/flx')
    shutil.rmtree(dst+'/prof')

    # ------------------------------------------------
    # get input data for field run 
    # ------------------------------------------------
    runname = runname_field
    # (1) frame
    ztot,nz,ttot,temp,fdust,fdust2,taudust,omrain,zom,poro,moistsrf,zwater,zdust,w,q,p,nstep,rstrt,runid \
        = get_inputs.get_input_frame(outdir,runname)
    
    if disply_glbl: 
        print(
            ztot,nz,ttot,temp,fdust,fdust2,taudust,omrain,zom,poro,moistsrf,zwater,zdust,w,q,p,nstep,rstrt,runid
            )
    # (2) switches
    w_scheme,mix_scheme,poro_iter,sldmin_lim,display,report,restart,rough,act_ON,dt_fix \
        ,cec_on,dz_fix,close_aq,poro_evol,sa_evol_1,sa_evol_2,psd_bulk,psd_full,season \
        = get_inputs.get_input_switches(outdir,runname)
    if disply_glbl: 
        print(
            w_scheme,mix_scheme,poro_iter,sldmin_lim,display,report,restart,rough,act_ON,dt_fix 
            ,cec_on,dz_fix,close_aq,poro_evol,sa_evol_1,sa_evol_2,psd_bulk,psd_full,season 
            )
    # (3) tracers
    sld_list,aq_list,gas_list,exrxn_list \
        = get_inputs.get_input_tracers(outdir,runname)
    if disply_glbl: 
        print(
            sld_list,aq_list,gas_list,exrxn_list
            )
    # (4) cec
    filename = 'cec.in'
    sld_data_list = get_inputs.get_input_sld_properties(outdir,runname,filename)
    sld_list_cecsp_field = [listtmp[0] for listtmp in sld_data_list]
        
        
    # ------------------------------------------------
    # get data from field run 
    # ------------------------------------------------
    
    # (1) cation and anion concs. 
    if disply_glbl: print('get chem comp')
    if include_pw:
        try:
            aqsps,btmconcs,dep,time_dum = get_int_prof.get_totsave_site(outdir,runname_field,dep_sample,itime,no_depth_integral=no_depth_integral)  # returning mol/ solid m3 depth averaged value 
        except:
            print('data not available: ', runname_field)
            return
    else:
        try:
            aqsps,btmconcs,dep,time_dum = get_int_prof.get_adssave_site(outdir,runname_field,dep_sample,itime,no_depth_integral=no_depth_integral)  # returning mol/ solid m3 depth averaged value 
        except:
            print('data not available: ', runname_field)
            return

    if disply_glbl: print(aqsps,btmconcs,dep)
    
    # exit()

    # (2) solid phase density 
    dense_lab = get_int_prof.get_rhobulk_int_site(outdir,runname_field,dep_sample,itime,no_depth_integral=no_depth_integral)
    
    # (3) solid phase wt%
    # if sld_list.shape ==(): sld_list = [str(sld_list)]  # passing solid phase list for field runs
    sps = list(sld_list)
    sldwt_list = []
    for sp in sps:
        sldwt = get_int_prof.get_sldwt_int_site(outdir,runname_field,dep_sample,[sp],itime,no_depth_integral=no_depth_integral)
        sldwt_list.append(sldwt)

    # (4) get DIC conc. (added 3.22.2023)
    dic,dep = get_int_prof.get_ave_DIC_save(outdir,runname_field,dep_sample,itime,no_depth_integral=no_depth_integral)



    # ------------------------------------------------
    # setting up conditions for lab experiments 
    # ------------------------------------------------
    
    # (1) dust composition
    poro_lab = water_frac/(1./dense_lab+water_frac)
    # poro_lab = 0.9999
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
        if include_Cl:
            oxide_ctnm_list.append( 'cl' )
            oxide_oxnm_list.append( cl_salt )
            oxide_stch_list.append( cl_salt_st )
            oxide_mass_list.append( cl_salt_wt )
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
        fdust = ztot_lab*(1-poro_lab)*conc* oxide_mass_list[iox]/oxide_stch_list[iox]
        fdust_list.append(fdust)
        fdust_nm_list.append(oxide_oxnm_list[iox])
       
    if include_pw:  
        if include_Cl:
            if disply_glbl: print(use_Sikora_buffer,use_CaCl2)
            
            if cl_salt =='nacl':
                inacl = fdust_nm_list.index('nacl')
                ina2o = fdust_nm_list.index('na2o')
                na_flx_nacl = fdust_list[inacl]/oxide_mass_list[inacl]*oxide_stch_list[inacl]
                na_flx_na2o = fdust_list[ina2o]/oxide_mass_list[ina2o]*oxide_stch_list[ina2o]
                if  na_flx_nacl <= na_flx_na2o: # na > cl   
                    fdust_list[ina2o] = (na_flx_na2o - na_flx_nacl) * oxide_mass_list[ina2o]/oxide_stch_list[ina2o]
                else:
                    print('*** ERROR: situation Na < Cl which necessitates change in Cl dust species from NaCl')
                    print('*** ERROR: for now, just stop simulation')
                    exit()
            elif cl_salt =='kcl':
                ikcl = fdust_nm_list.index('kcl')
                ik2o = fdust_nm_list.index('k2o')
                k_flx_kcl = fdust_list[ikcl]/oxide_mass_list[ikcl]*oxide_stch_list[ikcl]
                k_flx_k2o = fdust_list[ik2o]/oxide_mass_list[ik2o]*oxide_stch_list[ik2o]
                if  k_flx_kcl <= k_flx_k2o: # k > cl   
                    fdust_list[ik2o] = (k_flx_k2o - k_flx_kcl) * oxide_mass_list[ik2o]/oxide_stch_list[ik2o]
                else:
                    print('*** ERROR: situation K < Cl which necessitates change in Cl dust species from KCl')
                    print('*** ERROR: for now, just stop simulation')
                    exit()
            elif cl_salt =='cacl2':
                icacl2 = fdust_nm_list.index('cacl2')
                icao = fdust_nm_list.index('cao')
                ca_flx_cacl2 = fdust_list[icacl2]/oxide_mass_list[icacl2]*oxide_stch_list[icacl2]
                ca_flx_cao = fdust_list[icao]/oxide_mass_list[icao]*oxide_stch_list[icao]
                if  ca_flx_cacl2 <= ca_flx_cao: # cao > cacl2   
                    fdust_list[icao] = (ca_flx_cao - ca_flx_cacl2) * oxide_mass_list[icao]/oxide_stch_list[icao]
                else:
                    print('*** ERROR: situation Ca < 2Cl which necessitates change in Cl dust species from CaCl2')
                    print('*** ERROR: for now, just stop simulation')
                    exit()
            elif cl_salt =='mgcl2':
                imgcl2 = fdust_nm_list.index('mgcl2')
                imgo = fdust_nm_list.index('mgo')
                mg_flx_mgcl2 = fdust_list[imgcl2]/oxide_mass_list[imgcl2]*oxide_stch_list[imgcl2]
                mg_flx_mgo = fdust_list[imgo]/oxide_mass_list[imgo]*oxide_stch_list[imgo]
                if  mg_flx_mgcl2 <= mg_flx_mgo: # mgo > mgcl2   
                    fdust_list[imgo] = (mg_flx_mgo - mg_flx_mgcl2) * oxide_mass_list[imgo]/oxide_stch_list[imgo]
                else:
                    print('*** ERROR: situation Mg < 2Cl which necessitates change in Cl dust species from MgCl2')
                    print('*** ERROR: for now, just stop simulation')
                    exit()
            
            # for cl_salt, cat_oxide in zip(['nacl','kcl','cacl2'],['na2o','k2o','cao']):
                # isalt = fdust_nm_list.index(cl_salt)
                # ioxide = fdust_nm_list.index('cao')
                # cat_flx_salt = fdust_list[isalt]/oxide_mass_list[isalt]*oxide_stch_list[isalt]
                # cat_flx_oxide = fdust_list[ioxide]/oxide_mass_list[ioxide]*oxide_stch_list[ioxide]
                # if  cat_flx_salt <= cat_flx_oxide: # cao > cacl2   
                    # fdust_list[ioxide] = (cat_flx_oxide - cat_flx_salt) * oxide_mass_list[ioxide]/oxide_stch_list[ioxide]
                    # break
                # else:
                    # print('*** ERROR: situation K < Cl which necessitates change in Cl dust species from CaCl2')
                    # print('*** ERROR: for now, just stop simulation')
                    # exit()
            

    if use_CaCl2:
        # cacl2_conc = 0.01
        cacl2_wt2  = 110.986
        fdust_cacl2 = ztot_lab*poro_lab*cacl2_conc*1e3*cacl2_wt2

        if fdust_cacl2 > 0:
            if 'cacl2' not in fdust_nm_list:
                fdust_list.append(fdust_cacl2)
                fdust_nm_list.append('cacl2')
            else:
                fdust_list[fdust_nm_list.index('cacl2')] += fdust_cacl2
            
            if 'cl' not in oxide_ctnm_list: oxide_ctnm_list.append('cl')
            
    if add_salt:
        # nacl_conc = 0.01
        if salt_sp == 'nacl':
            salt_wt2  = 58.443
        elif salt_sp == 'kcl':
            salt_wt2  = 74.551
        elif salt_sp == 'cacl2':
            salt_wt2  = 110.986
        fdust_salt = ztot_lab*poro_lab*salt_conc*1e3*salt_wt2

        if fdust_salt > 0:
            if salt_sp not in fdust_nm_list:
                fdust_list.append(fdust_salt)
                fdust_nm_list.append(salt_sp)
            else:
                fdust_list[fdust_nm_list.index(salt_sp)] += fdust_salt
            
            if 'cl' in salt_sp:
                if 'cl' not in oxide_ctnm_list: oxide_ctnm_list.append('cl')

    if add_NO3:
        wtamnt = 80 # g/mol
        fdust_salt = ztot_lab*poro_lab*NO3_conc*1e3*wtamnt/2.

        if fdust_salt > 0:
            if 'amnt' not in fdust_nm_list:
                fdust_list.append(fdust_salt)
                fdust_nm_list.append('amnt')
            else:
                fdust_list[fdust_nm_list.index('amnt')] += fdust_salt
            
            if 'no3' not in oxide_ctnm_list: oxide_ctnm_list.append('no3')

    if include_DIC:                                     # (added 3.23.2023)
        fdust_dic = ztot_lab*(1-poro_lab)*dic* 30.      # (added 3.22.2023)
        fdust_list.append(fdust_dic)                    # (added 3.22.2023)
        fdust_nm_list.append('g1')                      # (added 3.22.2023)
    
    if use_Sikora_buffer:
        naoh_conc = 0.05  # as prescribed by original paper Sikora 2006
        naoh_conc = 0.0605430 # pH 7.51752
        if act_ON=='true': naoh_conc = 0.058 # pH 7.55725 # database = Sikora 2006
        # naoh_conc = 0.072 # pH 7.55725 # database = Goldberg et al. 2002
        # if act_ON=='true': naoh_conc = 0.068 # pH 7.55725 # database = Goldberg et al. 2002
        # naoh_conc = 0.066 # pH 7.661763
        # naoh_conc = 0.0674 # pH 7.69695
        sikora_ingredients_names = ['teas','ims','mesmh','gac','kcl','naoh']
        sikora_ingredients_mwts = [149.190,68.077,213.25,60.052,74.551,39.9971]
        sikora_ingredients_concs = [0.0696,0.0137,0.0314,0.0893,2.,naoh_conc]
        sikora_aqsps = ['tea','im','mes','ac','k','na','cl']
        
        for o in range(len(sikora_ingredients_names)):
            mwt_tmp = sikora_ingredients_mwts[o]
            conc_tmp = sikora_ingredients_concs[o]/water_frac
            nm_tmp = sikora_ingredients_names[o]
            fdust_tmp = ztot_lab*poro_lab*conc_tmp*1e3*mwt_tmp
            
            if nm_tmp in fdust_nm_list:
                j = fdust_nm_list.index(nm_tmp)
                fdust_list[j] = fdust_list[j]+fdust_tmp
            else:            
                fdust_list.append(fdust_tmp)
                fdust_nm_list.append(nm_tmp)
        
        for o in range(len(sikora_aqsps)):
            nm_tmp = sikora_aqsps[o]
            if nm_tmp not in oxide_ctnm_list:
                oxide_ctnm_list.append(nm_tmp)

    # fdust_lab = fdust_list[fdust_nm_list.index('cao')]
    fdust_lab = max(fdust_list)
    
    fdust_list = [fdust/fdust_lab  for fdust in fdust_list  ]  

    # ------------------------------------------------
    # making inputs for lab experiments 
    # ------------------------------------------------
    
    # (1) tracers 
    sld_list_lab = [sld for sld in sld_list if sld not in fdust_nm_list ]
    sld_list_lab.extend(fdust_nm_list)
    aq_list_lab = [aq for aq in aq_list if aq not in oxide_ctnm_list ]
    aq_list_lab.extend(oxide_ctnm_list)
    gas_list_lab = [] # (added 3.23.2023)
    if include_DIC: gas_list_lab.append('pco2')
    
    if disply_glbl: print(sld_list_lab)
    if disply_glbl: print(aq_list_lab)
    if disply_glbl: print(gas_list_lab)
    
    make_inputs.get_input_tracers(
        outdir=outdir_tmp
        ,runname=runname_lab
        ,sld_list = sld_list_lab
        ,aq_list = aq_list_lab
        ,gas_list = gas_list_lab
        )

    # (2) framework
    make_inputs.get_input_frame(
        outdir=outdir_tmp
        ,runname=runname_lab
        ,ztot=ztot_lab
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
        # ,nstep=200
        ,rstrt=rstrt
        ,runid=runname_lab
        )
        

    # (3) switches
    w_scheme_lab=0
    mix_scheme_lab=0 
    report_lab=0
    poro_iter_lab='true' 
    rough_lab      ='false'
    # dt_fix_lab='true'
    dt_fix_lab='false'
    close_aq_lab='true'
    psd_bulk_lab='false'
    psd_full_lab ='false'

    make_inputs.get_input_switches(
        outdir=outdir_tmp
        ,runname=runname_lab
        ,w_scheme=w_scheme_lab
        ,mix_scheme=mix_scheme_lab
        ,poro_iter=poro_iter_lab
        ,sldmin_lim=sldmin_lim 
        ,display=display
        ,report=report_lab
        ,restart='false' 
        ,rough=rough_lab
        ,act_ON=act_ON
        ,dt_fix=dt_fix_lab
        ,cec_on=cec_on
        ,dz_fix=dz_fix
        ,close_aq=close_aq_lab
        ,poro_evol=poro_evol
        ,sa_evol_1=sa_evol_1 
        ,sa_evol_2=sa_evol_2
        ,psd_bulk='false'
        ,psd_full='false'
        ,season='false'
        )
    
    # (4) boundary values
    # pr_list_lab = [(sp,sldwt_list[sps.index(sp)]/100.) for sp in sps if sp not in oxide_oxnm_list]
    # below; modified 8.28.2024 oxide existing in the field profile should be dissolved away 
    pr_list_lab = [(sp,sldwt_list[sps.index(sp)]/100.) for sp in sps]  
    atm_list_lab = [('pco2',3.16e-4),('po2',0.21),('pnh3',1e-50),('pn2o',1e-50)]
    if include_DIC: atm_list_lab = [('pco2',1e-20),('po2',0.21),('pnh3',1e-50),('pn2o',1e-50)] # (added 3.22.2023)
    rain_list_lab = [(sp,1e-20) for sp in aq_list_lab] # (added 8.10.2025)
    make_inputs.get_input_tracer_bounds(
        outdir=outdir_tmp,
        runname=runname_lab,
        pr_list = pr_list_lab,
        atm_list = atm_list_lab,
        rain_list = rain_list_lab,
        )
    
    # (5) dust composition
    filename = 'dust.in'
    sld_varlist =[ ( fdust_nm_list[i], fdust_list[i]) for i in range(len(fdust_nm_list)) ]
    make_inputs.get_input_sld_properties(
        outdir=outdir_tmp
        ,runname=runname_lab
        ,filename = filename
        ,sld_varlist=sld_varlist
        )


    # (6) suppressing reactions for solid phase except for those in dust
    # (6-i) kinetics
    filename = 'kinspc.in'
    sld_varlist = [ (sld,0) for sld in sld_list if sld not in oxide_oxnm_list ] 
    # sld_varlist = [ (sld,0) for sld in sld_list if sld not in oxide_oxnm_list + ['cc'] ]  # modifed so that cc can dissolve in lab | 10-11-2023
    make_inputs.get_input_sld_properties(
        outdir=outdir_tmp
        ,runname=runname_lab
        ,filename = filename
        ,sld_varlist=sld_varlist
        )
        
    # (6-ii) cec    
    filename = 'cec.in'
    sld_varlist = [(listtmp[0],listtmp[1],listtmp[2],listtmp[3],listtmp[4],listtmp[5],listtmp[6],listtmp[7])  for listtmp in sld_data_list] 
    sld_varlist += [ (sld,0,0,0,0,0,0,0) for sld in fdust_nm_list if sld not in sld_list_cecsp_field ] 
    make_inputs.get_input_sld_properties(
        outdir=outdir_tmp
        ,runname=runname_lab
        ,filename = filename
        ,sld_varlist=sld_varlist
        )
    
    # ------------------------------------------------
    # run lab experiments 
    # ------------------------------------------------
    
    if re_compile:
        # os.system('make')
        # os.system('make --file=makefile_test')
        os.system('cp ' + exename_src + to + outdir_tmp + runname_lab + where + exename)
    
    run_success = False
    
    print('chmod u+x '+outdir_tmp+runname_lab+where+exename)
    os.system('chmod u+x '+outdir_tmp+runname_lab+where+exename)
    
    if not sub_as_a_job:
        
        if not set_ExTimeLimit:
            print()
            
            if os.path.exists(outdir_tmp+runname_lab+'/run_complete.txt'): os.remove(outdir_tmp+runname_lab+'/run_complete.txt')
            
            if make_runlogfile:
                os.system(outdir_tmp+runname_lab+where+exename+' > '+outdir_tmp+runname_lab+'/logfile.txt')
            else:
                os.system(outdir_tmp+runname_lab+where+exename)
            
            run_success = True
            
            if not os.path.exists(outdir_tmp+runname_lab+'/run_complete.txt'): 
                run_success = False
                flag_error = True
                print( '*** run is not completed ' )
                return time_dum, np.nan, np.nan, np.nan, flag_error, runname_lab
            
        else:
            
            if os.path.exists(outdir_tmp+runname_lab+'/run_complete.txt'): os.remove(outdir_tmp+runname_lab+'/run_complete.txt')
            
            if make_runlogfile:
                logf = open(outdir_tmp+runname_lab+'/logfile.txt', 'w')
                proc = subprocess.Popen([outdir_tmp+runname_lab+where+exename], stdout=logf)
            else:
                proc = subprocess.Popen([outdir_tmp+runname_lab+where+exename])
            
            my_timeout =60*Ex_TimeLim
            
            try:
                proc.wait(my_timeout)
                print('run finished within {:f} min'.format(int(my_timeout/60.)))
            
                run_success = True
                
                if not os.path.exists(outdir_tmp+runname_lab+'/run_complete.txt'): 
                    run_success = False
                    flag_error = True
                    print( '*** run is not completed ' )
                    return time_dum, np.nan, np.nan, np.nan, flag_error, runname_lab

            except subprocess.TimeoutExpired:
                proc.kill()
                flag_error = True
                print('run UNfinished within {:f} min'.format(int(my_timeout/60.)))
                return time_dum, np.nan, np.nan, np.nan, flag_error, runname_lab
        
        if make_runlogfile and run_success:
            os.remove(outdir_tmp+runname_lab+'/logfile.txt')
    
    else:
        if os.path.exists(outdir_tmp+runname_lab+'/run_complete.txt'): os.remove(outdir_tmp+runname_lab+'/run_complete.txt')
        slurm_cmd = 'sbatch  --time=0-24:00 --account=gts-creinhard3 --nodes=1 --ntasks=1  -qinferno  --mem-per-cpu=4G run_a_shell.sbatch '
        print(slurm_cmd+outdir_tmp+runname_lab+where+exename + ' > ' + outdir_tmp+runname_lab+'/joblog.txt' + ' 2> ' + outdir_tmp+runname_lab+'/joberr.txt') # submit a job
        os.system(slurm_cmd+outdir_tmp+runname_lab+where+exename + ' > ' + outdir_tmp+runname_lab+'/joblog.txt' + ' 2> ' + outdir_tmp+runname_lab+'/joberr.txt') # submit a job
        
        # get job id 
        while(True):
            try:
                f = open( outdir_tmp+runname_lab+'/joblog.txt', mode = 'r' )
                jobidstr = re.sub( "[^0-9]", "",  f.read() )
                f.close()
                break
            except:
                pass
        
        st = time.time() # time before loop
        cnt = 0
        while(cnt <= Ex_TimeLim):
            time.sleep(60)
            if os.path.exists(outdir_tmp+runname_lab+'/run_complete.txt'):
                run_success = True
                break
            cnt += 1
        et = time.time() 

        if os.path.exists(outdir_tmp+runname_lab+'/run_complete.txt'):
            run_success = True
        
        if run_success:
            print('run finished in {:f} min'.format(( et - st )/60.))
            print('slurm-{}.out'.format(jobidstr))
            os.remove('slurm-{}.out'.format(jobidstr))
            # kill job just in case 
            print('scancel ' + jobidstr )
            os.system('scancel ' + jobidstr )
        else:
            print('run UNfinished within {:f} min'.format(int( Ex_TimeLim )))
        
            # kill job
            print('scancel ' + jobidstr )
            os.system('scancel ' + jobidstr )

            flag_error = True
            return time_dum, np.nan, np.nan, np.nan, flag_error, runname_lab
    
    # ------------------------------------------------
    # retrieve data from lab experiments  
    # ------------------------------------------------
    
    # (1) ph
    phint_lab = get_int_prof.get_ph_int_site(outdir_tmp,runname_lab,dep_sample,20,no_depth_integral=False)
    
    # (2) IS 
    IS_ave,dep = get_int_prof.get_ave_IS(outdir_tmp,runname_lab,dep_sample,20,no_depth_integral=False)
    
    # (3) exchangeable acidity (%CEC)
    acint = get_int_prof.get_ac_int_site_v2(outdir_tmp,runname_lab,dep_sample,20,no_depth_integral=False)
    
    
    # ------------------------------------------------
    # remove or copy lab experiments  
    # ------------------------------------------------
    
    # :: (A) case to remove 
    if not save_lab_dir:
        if not use_local_storage: 
            print('>>> save_lab_dir = {} and use_local_storage = {}; so going to remove lab dir [{}]'.format(save_lab_dir,use_local_storage,dst))
            shutil.rmtree( dst ) 
    
    # :: (B) case to save 
    if save_lab_dir:
        if use_local_storage and outdir_tmp!=outdir:
            src = outdir_tmp + runname_lab 
            dst = outdir + runname_lab 
            
            print('>>> save_lab_dir = {} and use_local_storage = {}; so going to mv lab tmp dir [{}] to working dir [{}]'.format(save_lab_dir,use_local_storage,src,dst))
            
            if not os.path.exists(dst): 
                shutil.copytree(src, dst, dirs_exist_ok=True)
            else:
                shutil.rmtree(dst)
                shutil.copytree(src, dst, dirs_exist_ok=True)
    
    return time_dum,IS_ave,acint,phint_lab, flag_error, runname_lab


def timeseries_soilpH():
    
    outdir = '/storage/coda1/p-creinhard3/0/ykanzaki3/scepter_output/'
    runname_field   = 'US_cropland_297_sph_N_DIC_v2_spintuneup_field'
    runname_field   = 'US_cropland_311_sph_N_DIC_v2_spintuneup_field'
    runname_field   = 'US_cropland_311_sph_N_DIC_cacl2_2p5_all_spintuneup_field'
    runname_field   = 'US_cropland/sph_N_cacl2_2p5_pco2_ps_hmix_b/alpha2p0/spintuneup_deptest_v3/128_field'
    # runname_field   = 'US_cropland/sph_N_cacl2_2p5_pco2_ps_hmix_b/alpha2p0/spintuneup_deptest/311_field'
    # runname_field   = 'US_cropland/sph_N_cacl2_2p5_pco2_ps_hmix_b/alpha2p0/spintuneup_depcnst/128_field'
    # runname_field   = 'test_inert_spinup'
    # runname_field   = sys.argv[1]
    
    dep_sample = 0.15
    # dep_sample = 0.18
    # dep_sample = 0.25
    
    time_dum,IS_lab,acint_lab,phint_lab,flag_error,runname_lab = calc_soilpH(
        outdir,runname_field,dep_sample, 20, 
        include_N           = True,
        include_Cl          = False,
        include_Al          = False,
        include_DIC         = True,
        use_CaCl2           = True,
        # use_local_storage   = True,
        use_local_storage   = False,
        save_lab_dir        = True,
        use_Sikora_buffer   = False,
        water_frac          = 2.5,
        cacl2_conc          = 0.01,
        )
        
    # res_list = []
    # n_sample = 20
    # for i in range(n_sample):
        # time,IS_lab,phint_lab = calc_soilpH(outdir,runname_field,dep_sample, i+1)
        # res_list.append([time,IS_lab,phint_lab])
    
    if disply_glbl: print(time_dum,IS_lab,acint_lab,phint_lab,flag_error)
    
    # np.savetxt(outdir+runname_field+'/timeseries_soil_ph_'+str(int(dep_sample*100)) +'cm.res',np.array(res_list))


def depseries_soilpH(runname_field,cl_salt):
    
    outdir = '/storage/coda1/p-creinhard3/0/ykanzaki3/scepter_output/'
    if 'CYGWIN' in platform.system(): outdir = '/home/ykanz/scepter_output/'    
    
    # runname_field   = 'US_cropland_297_sph_N_DIC_v2_spintuneup_field'
    # runname_field   = 'US_cropland_311_sph_N_DIC_v2_spintuneup_field'
    # runname_field   = 'US_cropland_311_sph_N_DIC_cacl2_2p5_all_spintuneup_field'
    # runname_field   = 'US_cropland/sph_N_cacl2_2p5_pco2_ps_hmix_b/alpha2p0/spintuneup_deptest_v2/128_field'
    # runname_field   = 'US_cropland/sph_N_cacl2_2p5_pco2_ps_hmix_b/alpha2p0/spintuneup_depcnst_v2/128_field'
    # runname_field   = 'US_cropland/sph_N_cacl2_2p5_pco2_ps_hmix_b/alpha2p0/spintuneup_deptest_v3/128_field'
    # runname_field   = 'US_cropland/sph_N_cacl2_2p5_pco2_ps_hmix_b/alpha2p0/spintuneup_depcnst/128_field'
    # runname_field   = 'test_inert_spinup'
    # runname_field   = 'JCU/test'
    # runname_field   = sys.argv[1]
    
    infile  = outdir+runname_field+'/prof/prof_aq-{:03d}.txt'.format(20)
    data    = np.loadtxt(infile,skiprows=1)
    depths = data[:,0]
    
    dep_sample = 0.15
    # dep_sample = 0.18
    # dep_sample = 0.25
    
    res_list = []
    for dep_sample in depths:
        time_dum,IS_lab,acint_lab,phint_lab,flag_error,runname_lab = calc_soilpH(
            outdir,runname_field,dep_sample, 20, 
            include_N           = True,
            # include_Cl          = False,
            include_Cl          = True,
            # cl_salt             = 'kcl',
            cl_salt             = cl_salt,
            include_Al          = False,
            include_DIC         = True,
            # use_CaCl2           = True,
            use_CaCl2           = False,
            no_depth_integral   = True,
            # use_local_storage   = True,
            use_local_storage   = False,
            # save_lab_dir        = True,
            save_lab_dir        = False,
            set_ExTimeLimit     = False,
            use_Sikora_buffer   = False,
            water_frac          = 5,
            cacl2_conc          = 0.01,
            add_salt            = True,
            salt_sp             = 'nacl',
            salt_conc           = 0.04,
            )
        if disply_glbl: print(time_dum,IS_lab,acint_lab,phint_lab,flag_error)
        res_list.append([dep_sample,IS_lab,acint_lab,phint_lab,time_dum])
    
    filename = outdir+runname_field+'/prof/prof_soilph-{:03d}.txt'.format(20)
    np.savetxt(filename,np.array(res_list))

    with open(filename, "r") as f:
        contents = f.readlines()

    contents.insert(0, "{}\t{}\t{}\t{}\t{}\n".format('depth','IS','Hexch','soilpH','time'))

    with open(filename, "w") as f:
        contents = "".join(contents)
        f.write(contents)

def CECseries_bufferpH(act_ON,database_name,alpha_str):
    
    
    cec_list = [0,1,2,4,8,16,32,64]
    cec_list = list(range(21))
    # cec_list = list(range(21)) + [0.001,0.003,0.01,0.03,0.1,0.3] + [50,100,200]
    cec_list = list(range(41)) + [0.001,0.003,0.01,0.03,0.1,0.3] + [50,100,200]
    cec_list.sort()
    # cec_list = list(range(13))
    
    outdir = '/storage/coda1/p-creinhard3/0/ykanzaki3/scepter_output/'
    
    res_list = []
    
    for cec in cec_list:
        if not act_ON: 
            if cec<1 and cec>0:
                runname_field   = 'Sikora/alpha{}/test_inert_spinup_noact_cec_{}'.format( alpha_str, '{:.3f}'.format(cec).replace('.','p') ) 
            else:
                runname_field   = 'Sikora/alpha{}/test_inert_spinup_noact_cec_{:d}'.format(alpha_str, cec) 
        else: 
            runname_field   = 'Sikora/alpha{}/test_inert_spinup_act_cec_{:d}'.format(alpha_str, cec) 
        
        runname_field += '_'+database_name
        
        dep_sample = 0.15
        # dep_sample = 0.18
        # dep_sample = 0.25
        
        time_dum,IS_lab,acint_lab,phint_lab,flag_error,runname_lab = calc_soilpH(
            outdir,runname_field,dep_sample, 20, 
            include_N           = True,
            include_Cl          = False,
            include_Al          = False,
            include_DIC         = True,
            use_CaCl2           = True,
            # use_local_storage   = True,
            use_local_storage   = False,
            use_Sikora_buffer   = True,
            water_frac          = 2.5,
            cacl2_conc          = 0.01,
            save_lab_dir        = True,
            )
        
        res_list.append([cec,IS_lab,acint_lab,phint_lab])
    
    actstr = 'no' if not act_ON else ''
    datadir = f'/storage/coda1/p-creinhard3/0/ykanzaki3/scripts/pyweath/data/Sikora/alpha{alpha_str}/'
    if not os.path.exists(datadir): os.makedirs(datadir)
    outfile_name = f'{actstr}act_cecseries_bufferpH_{database_name}.res'
    
    np.savetxt(datadir+outfile_name,np.array(res_list))

def calc_bufferpH(runname_field,dep_sample):
    
    
    outdir = '/storage/coda1/p-creinhard3/0/ykanzaki3/scepter_output/'
    
        
    # dep_sample = 0.15
    # dep_sample = 0.18
    # dep_sample = 0.25
    
    time_dum,IS_lab,acint_lab,phint_lab,flag_error,runname_lab = calc_soilpH(
        outdir,runname_field,dep_sample, 20, 
        include_N           = True,
        include_Cl          = False,
        include_Al          = False,
        include_DIC         = True,
        use_CaCl2           = True,
        # use_local_storage   = True,
        use_local_storage   = False,
        use_Sikora_buffer   = True,
        water_frac          = 2.5,
        cacl2_conc          = 0.01,
        # save_lab_dir        = True,
        save_lab_dir        = False,
        )
        

    excludes = ['excludes',]

    output_dict = {k: v for k, v in locals().items() if k not in excludes}

    print(output_dict)
    
    # ---- save output
    dst = outdir + runname_field + '/Sikora_buffer_pH.res'
    with open(dst, 'w') as f:
        json.dump(output_dict, f, sort_keys=True, indent=4)

    return output_dict
    
    
    
    


def CationSeries_soilpH():
    
    # outdir = '/storage/coda1/p-creinhard3/0/ykanzaki3/scepter_output/'
    outdir = '/storage/coda1/p-creinhard3/0/ykanzaki3/scepter_output/tests/'
    
    cec_list = [5,15,30]
    cats = ['ca','mg','na','k']
    catrains = [1e-7,1e-6,1e-5,1e-4,1e-3,1e-2] 
    
    # cec_list = [5]
    # cats = ['mg']
    # catrains = [1e-5] 
    
    # cec_list = [30]
    # cats = ['ca']
    # catrains = [1e-6] 
    
    # cec_list = [5]
    # cats = ['na']
    # catrains = [1e-6] 
    
    for cat in cats:
        res_list = []
        for cec in cec_list:
            for catrain in catrains:
                
                # runname_field = 'test_inert_spinup_cec_{:d}_{:}_p{:d}'.format(cec,cat,int(-np.log10(catrain)))
                runname_field = 'inert_spinup_cec_{:d}_{:}_p{:d}'.format(cec,cat,int(-np.log10(catrain)))

                dep_sample = 0.15
                # dep_sample = 0.18
                # dep_sample = 0.25

                # (1) get porewater pH 
                phint = get_int_prof.get_ph_int_site(outdir,runname_field,dep_sample,20)
                
                time_dum,IS_lab,acint_lab,phint_lab,flag_error,runname_lab = calc_soilpH(
                    outdir,runname_field,dep_sample, 20, 
                    include_N           = False,
                    include_Cl          = False,
                    include_Al          = False,
                    include_DIC         = True,
                    use_CaCl2           = True,
                    # use_local_storage   = True,
                    use_local_storage   = False,
                    use_Sikora_buffer   = False,
                    water_frac          = 2.5,
                    cacl2_conc          = 0.01,
                    # save_lab_dir        = False,
                    save_lab_dir        = True,
                    )
                
                res_list.append([cec,catrain,phint,phint_lab,acint_lab])
    
    
        outfile_name = '{}_conc_cecseries_soilpH_py.res'.format(cat.title())
        
        np.savetxt(outfile_name,np.array(res_list))
        

def S2006_bufferpH(
    outdir,
    runname_field,
    dep_sample,
    include_N,
    include_Cl,
    include_Al,
    include_DIC,
    add_salt,
    add_NO3,
    use_local_storage,
    exename_src,
    re_compile,
    ):
    
    # outdir = '../scepter_output/'
    # runname_field   = 'US_cropland_297_sph_N_DIC_v2_spintuneup_field'
    # runname_field   = 'US_cropland_311_sph_N_DIC_v2_spintuneup_field'
    # runname_field   = 'US_cropland_297_sph_N_DIC_act_all_spintuneup_field'
    # runname_field   = 'test_inert_spinup'
    # runname_field   = sys.argv[1]
    
    # dep_sample = 0.15
    # dep_sample = 0.18
    # dep_sample = 0.25
    
    res_list = []
    
    water_frac_list = [5,2,1,0.5]
    cacl2_conc      = 0.0
    
    save_lab_dir = True
    # save_lab_dir = False
    
    salt_sp         = 'kcl'
    
    water_frac = 2
    
    salt_conc       = 0.02995515189123628  / water_frac
    NO3_conc        = 0.00876723282429123 / water_frac # + 0.0003118665211289568*2
    
    time_dum,IS_lab,acint_lab,phint_lab,flag_error,runname_lab = calc_soilpH(
        outdir,runname_field,dep_sample, 20, 
        include_N           = include_N,
        include_Cl          = include_Cl,
        include_Al          = include_Al,
        include_DIC         = include_DIC,
        use_CaCl2           = False,
        add_salt            = add_salt,
        add_NO3             = add_NO3,
        use_local_storage   = use_local_storage,
        save_lab_dir        = save_lab_dir,
        use_Sikora_buffer   = True,
        set_ExTimeLimit     = True,
        water_frac          = water_frac,
        cacl2_conc          = cacl2_conc,
        salt_conc           = salt_conc,
        NO3_conc            = NO3_conc,
        salt_sp             = salt_sp,
        Ex_TimeLim          = 20,
        exename_src         = exename_src,
        re_compile          = re_compile,
        )
        
    res_list.append([cacl2_conc,water_frac,IS_lab,acint_lab,phint_lab])
    
    frame = inspect.currentframe()
    if include_DIC: np.savetxt(outdir + runname_field + '/' + frame.f_code.co_name+'_DIC_H2O.res',np.array(res_list))
    else: np.savetxt(outdir + runname_field + '/' + frame.f_code.co_name+'_H2O.res',np.array(res_list))
    


def MK2010_series_soilpH(
    outdir,
    runname_field,
    dep_sample,
    include_N,
    include_Cl,
    include_Al,
    include_DIC,
    add_salt,
    add_NO3,
    use_local_storage,
    exename_src,
    re_compile,
    ):
    
    # outdir = '../scepter_output/'
    # runname_field   = 'US_cropland_297_sph_N_DIC_v2_spintuneup_field'
    # runname_field   = 'US_cropland_311_sph_N_DIC_v2_spintuneup_field'
    # runname_field   = 'US_cropland_297_sph_N_DIC_act_all_spintuneup_field'
    # runname_field   = 'test_inert_spinup'
    # runname_field   = sys.argv[1]
    
    # dep_sample = 0.15
    # dep_sample = 0.18
    # dep_sample = 0.25
    
    res_list = []
    
    water_frac_list = [5,2,1,0.5]
    cacl2_conc      = 0.0
    
    # save_lab_dir = True
    save_lab_dir = False
    
    salt_sp         = 'kcl'
    
    for water_frac in water_frac_list:
    
        salt_conc       = 0.02995515189123628  / water_frac
        NO3_conc        = 0.00876723282429123 / water_frac # + 0.0003118665211289568*2
    
        time_dum,IS_lab,acint_lab,phint_lab,flag_error,runname_lab = calc_soilpH(
            outdir,runname_field,dep_sample, 20, 
            include_N           = include_N,
            include_Cl          = include_Cl,
            include_Al          = include_Al,
            include_DIC         = include_DIC,
            use_CaCl2           = False,
            add_salt            = add_salt,
            add_NO3             = add_salt,
            use_local_storage   = use_local_storage,
            save_lab_dir        = save_lab_dir,
            use_Sikora_buffer   = False,
            set_ExTimeLimit     = True,
            water_frac          = water_frac,
            cacl2_conc          = cacl2_conc,
            salt_conc           = salt_conc,
            NO3_conc            = NO3_conc,
            salt_sp             = salt_sp,
            Ex_TimeLim          = 20,
            exename_src         = exename_src,
            re_compile          = re_compile,
            )
        
        res_list.append([cacl2_conc,water_frac,IS_lab,acint_lab,phint_lab])
    
    frame = inspect.currentframe()
    if include_DIC: np.savetxt(outdir + runname_field + '/' + frame.f_code.co_name+'_DIC_H2O.res',np.array(res_list))
    else: np.savetxt(outdir + runname_field + '/' + frame.f_code.co_name+'_H2O.res',np.array(res_list))
    
    
    
    res_list = []
    
    cacl2_conc_list = [0.0025,0.005,0.01]
    # cacl2_conc_list = [0.0025,0.005,0.01, 0.05, 0.1, 0.5]
    water_frac_list = [1., 2.]
    water_frac_list = [1.]
    
    for cacl2_conc in cacl2_conc_list:
        for water_frac in water_frac_list:
    
            salt_conc       = 0.02995515189123628  / water_frac
            NO3_conc        = 0.00876723282429123 / water_frac # + 0.0003118665211289568*2
        
            time_dum,IS_lab,acint_lab,phint_lab,flag_error,runname_lab = calc_soilpH(
                outdir,runname_field,dep_sample, 20, 
                include_N           = include_N,
                include_Cl          = include_Cl,
                include_Al          = include_Al,
                include_DIC         = include_DIC,
                use_CaCl2           = True,
                add_salt            = add_salt,
                add_NO3             = add_salt,
                use_local_storage   = use_local_storage,
                save_lab_dir        = save_lab_dir,
                use_Sikora_buffer   = False,
                set_ExTimeLimit     = True,
                water_frac          = water_frac,
                cacl2_conc          = cacl2_conc,
                salt_conc           = salt_conc,
                NO3_conc            = NO3_conc,
                salt_sp             = salt_sp,
                Ex_TimeLim          = 20,
                exename_src         = exename_src,
                re_compile          = re_compile,
                )
            
            res_list.append([cacl2_conc,water_frac,IS_lab,acint_lab,phint_lab])
    
    frame = inspect.currentframe()
    if include_DIC: np.savetxt(outdir + runname_field + '/' + frame.f_code.co_name+'_DIC_CaCl2.res',np.array(res_list))
    else: np.savetxt(outdir + runname_field + '/' + frame.f_code.co_name+'_CaCl2.res',np.array(res_list))
        

def MK2010_series_soilpH_2(
    outdir,
    runname_field,
    dep_sample,
    include_N,
    include_Cl,
    include_Al,
    include_DIC,
    add_salt,
    add_NO3,
    use_local_storage,
    ):
    
    # outdir = '../scepter_output/'
    # runname_field   = 'US_cropland_297_sph_N_DIC_v2_spintuneup_field'
    # runname_field   = 'US_cropland_311_sph_N_DIC_v2_spintuneup_field'
    # runname_field   = 'US_cropland_297_sph_N_DIC_act_all_spintuneup_field'
    # runname_field   = 'test_inert_spinup'
    # runname_field   = sys.argv[1]
    
    # dep_sample = 0.15
    # dep_sample = 0.18
    # dep_sample = 0.25
    
    save_lab_dir = True
    save_lab_dir = False
    
    res_list = []
    
    nacl_conc_list  = [1e-2,2e-2,3e-2]
    # nacl_conc_list  = [0.5e-2,1e-2,1.5e-2]
    cacl2_conc      = 0.01
    water_frac      = 1
    salt_sp         = 'kcl'
    NO3_conc        = 0.00876723282429123 # + 0.0003118665211289568*2
    
    # 0.02995515189123628 0.00876723282429123 0.0003118665211289568
    
    for nacl_conc in nacl_conc_list:
    
        time_dum,IS_lab,acint_lab,phint_lab,flag_error = calc_soilpH(
            outdir,runname_field,dep_sample, 20, 
            include_N           = include_N,
            include_Cl          = include_Cl,
            include_Al          = include_Al,
            include_DIC         = include_DIC,
            use_CaCl2           = False,
            add_salt            = add_salt,
            add_NO3             = add_NO3,
            use_local_storage   = use_local_storage,
            save_lab_dir        = save_lab_dir,
            use_Sikora_buffer   = False,
            set_ExTimeLimit     = True,
            water_frac          = water_frac,
            cacl2_conc          = cacl2_conc,
            salt_conc           = nacl_conc,
            NO3_conc            = NO3_conc,
            salt_sp             = salt_sp,
            Ex_TimeLim          = 20,
            )
        
        res_list.append([cacl2_conc,nacl_conc,IS_lab,acint_lab,phint_lab])
    
    frame = inspect.currentframe()
    if include_DIC: np.savetxt(outdir + runname_field + '/' + frame.f_code.co_name+'_DIC_H2O_NaCl_ws1.res',np.array(res_list))
    else: np.savetxt(outdir + runname_field + '/' + frame.f_code.co_name+'_H2O_NaCl_ws1.res',np.array(res_list))
    
    
    
    res_list = []
    
    for nacl_conc in nacl_conc_list:
    
        time_dum,IS_lab,acint_lab,phint_lab,flag_error = calc_soilpH(
            outdir,runname_field,dep_sample, 20, 
            include_N           = include_N,
            include_Cl          = include_Cl,
            include_Al          = include_Al,
            include_DIC         = include_DIC,
            use_CaCl2           = True,
            add_salt            = True,
            add_NO3             = True,
            use_local_storage   = use_local_storage,
            save_lab_dir        = save_lab_dir,
            use_Sikora_buffer   = False,
            set_ExTimeLimit     = True,
            water_frac          = water_frac,
            cacl2_conc          = cacl2_conc,
            salt_conc           = nacl_conc,
            NO3_conc            = NO3_conc,
            salt_sp             = salt_sp,
            Ex_TimeLim          = 20,
            )
            
        res_list.append([cacl2_conc,nacl_conc,IS_lab,acint_lab,phint_lab])
    
    frame = inspect.currentframe()
    if include_DIC: np.savetxt(outdir + runname_field + '/' + frame.f_code.co_name+'_DIC_0p01MCaCl2_NaCl_ws1.res',np.array(res_list))
    else: np.savetxt(outdir + runname_field + '/' + frame.f_code.co_name+'_0p01MCaCl2_NaCl_ws1.res',np.array(res_list))
    
    
def main():   
    
    
    # dep_sample=0.15
    # runname_field=sys.argv[1]
    # calc_bufferpH(runname_field,dep_sample)
    # exit("end run")
    
    # CationSeries_soilpH()
    # exit("end run")
    
    # timeseries_soilpH()
    # exit("end run")
    
    # runname_field = 'JCU/test'
    # depseries_soilpH(runname_field)
    # exit("end run")
    
    
    # act_ON = False
    # database_name = 'Sikora'
    # alpha = 0
    # alpha = 2
    # alpha = 5
    # alpha = 10
    # alpha_str = '{:.1f}'.format(alpha).replace('.','p')
    # CECseries_bufferpH(act_ON,database_name,alpha_str)
    # exit("end run")
    
    
    
    
    
    # i_parallel = int(sys.argv[1])
    # n_parallel = int(sys.argv[2])
    
    # name_base  = sys.argv[3]
    # name_base  = 'test_Pot7_25C_v2_alpha'
    name_base  = 'test_Pot7_25C_v2_poro0p4_alpha'
    name_base  = 'test_Pot7_25C_v2_densqrtz_alpha'
    name_base  = 'test_Pot7_25C_v2_CEC10p0'
    name_base  = 'test_Pot7_25C_v2_amal'
    name_base  = 'test_Pot7_25C_v2_noact_alpha'
    name_base  = 'test_Pot7_25C_v2_act_alpha'
    # name_base  = 'test_Pot7_25C_v3_noact_alpha'
    name_base  = 'test_Pot7_25C_v2_act_2alpha'
    
    ph_pw = 6.3
    ph_pw = 6.0
    ph_pw = 6.1
    # ph_pw = 5.8
    
    om_tgt = 4.9
    cec = 8.9
    cec_2 = 120
    cec_1 = (100 *cec - om_tgt*cec_2)/(100 - om_tgt)
    
    name_base  = 'test_Pot7_25C_v2_act_pH{:.2f}_alpha'.format(ph_pw).replace('.','p')
    
    name_base  = 'test_Pot7_25C_v2_act_2CEC{:.1f}-{:.1f}_pH{:.1f}_alpha'.format(cec_1,cec_2,ph_pw).replace('.','p')
    name_base  = 'chk_Pot7_25C_v2_act_2CEC{:.1f}-{:.1f}_pH{:.1f}_alpha'.format(cec_1,cec_2,ph_pw).replace('.','p')
    name_base  = 'chk2_Pot7_25C_v2_act_2CEC{:.1f}-{:.1f}_pH{:.1f}_alpha'.format(cec_1,cec_2,ph_pw).replace('.','p')
    name_base  = 'GMD24/bnch_Pot7_25C_v2_act_2CEC{:.1f}-{:.1f}_pH{:.1f}_alpha'.format(cec_1,cec_2,ph_pw).replace('.','p')
    
    outdir  = '../scepter_output/'
    
    alpha_list = list(np.linspace(0.1,5.5,19))
    alpha_list = [3.4]
    alpha_list = [0.1]
    alpha_list = [0.0]
    alpha_list = [0.5]
    alpha_list = [1.0]
    alpha_list = [1.1]
    alpha_list = [1.2]
    alpha_list = [1.3]
    # alpha_list = [1.4]
    # alpha_list = [1.7]
    # alpha_list = [5.5]
    # alpha_list = [(0, 6)]
    
    exename_src = 'scepter_gmd24'
    
    use_local_storage  = False
    
    re_compile  = True
    
    for i in range(len(alpha_list)):
                   
        # if i%n_parallel!=i_parallel: continue
        
        alpha = alpha_list[i]
        
        runname_field = name_base+'{:.1f}'.format(alpha).replace('.','p')
        
        
        # alpha_1 = alpha_list[i][0]
        # alpha_2 = alpha_list[i][1]
        # runname_field = name_base+'{:.1f}-{:.1f}'.format(alpha_1,alpha_2).replace('.','p')
        
        # runname_field = 'test_Pot7_v2'
        # runname_field = 'test_Pot7_25C_alpha0p0'
        # runname_field = 'test_Pot7_25C_alpha3p4'
        # runname_field = 'test_Pot7_25C_alpha3p7'
        # runname_field = 'test_Pot7_25C_alpha4p0'
        # runname_field = 'test_Pot7_25C_alpha4p3'
        # runname_field = 'test_Pot7_25C_alpha4p9'
        # runname_field = 'test_Pot7_25C_alpha5p2'
        # runname_field = 'test_Pot7_25C_alpha5p5'
        
        dep_sample = 0.15
        include_N = True
        include_Cl  = True
        # include_Cl  = False
        # include_Al  = True
        include_Al  = False
        # use_local_storage  = True
        include_DIC = True
                    
        # MK2010_series_soilpH_2(
            # outdir,
            # runname_field,
            # dep_sample,
            # include_N,
            # include_Cl,
            # include_Al,
            # include_DIC,
            # add_salt,
            # add_NO3,
            # use_local_storage,
            # )

        dep_sample = 0.15
        include_N = True
        include_Cl  = True
        # include_Cl  = False
        # include_Al  = True
        include_Al  = False
        # use_local_storage  = True
        include_DIC = True
        
        add_salt = True
        add_NO3 = True
        
        MK2010_series_soilpH(
            outdir,
            runname_field,
            dep_sample,
            include_N,
            include_Cl,
            include_Al,
            include_DIC,
            add_salt,
            add_NO3,
            use_local_storage,
            exename_src,
            re_compile,
            )
            
        S2006_bufferpH(
            outdir,
            runname_field,
            dep_sample,
            include_N,
            include_Cl,
            include_Al,
            include_DIC,
            add_salt,
            add_NO3,
            use_local_storage,
            exename_src,
            re_compile,
            )
        
    exit("end run")
    
    
    
    
                
    
if __name__ == '__main__':
    main()