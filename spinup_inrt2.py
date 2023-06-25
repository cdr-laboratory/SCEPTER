import os,shutil,sys
import spinup,get_int_prof
import numpy as np
   
   
def spin_inert(outdir_src,runname,
    sld_list,
    aq_list,
    gas_list,
    pr_list,
    rain_list,
    omrain,
    sld_varlist_omrain,
    sld_varlist_cec,
    act_ON,
    ):

    # runname = 'test'
    
    # cec     = 10.0
    # logkh   = 5.9
    # alpha   = 3.4
    
    # na      = 1e-5
    
    # ---- frame.in ----
    ztot                = 0.5
    nz                  = 30
    ttot                = 1e5
    temp                = 25
    fdust               = 0
    fdust2              = 0
    taudust             = 0
    # omrain              = 300
    zom                 = 0.25
    poro                = 0.5
    moistsrf            = 0.22
    zwater              = 1000
    zdust               = 0.25
    w                   = 1e-3
    q                   = 0.55
    p                   = 1e-5
    nstep               = 10
    rstrt               = 'self'
    runid               = runname
    # ---- switches.in ----
    w_scheme            = 1 
    mix_scheme          = 1 # 1 --Fickian, 2 --Homogeneous
    poro_iter           = 'false' 
    sldmin_lim          = 'true'
    display             = 'true'
    disp_lim            = 'true'
    restart             = 'false'
    rough               = 'true'
    # act_ON              = 'true'
    # act_ON              = 'false'
    dt_fix              = 'false'
    cec_on              = 'true'
    dz_fix              = 'true'
    close_aq            = 'false'
    poro_evol           = 'true'
    sa_evol_1           = 'true'
    sa_evol_2           = 'false'
    psd_bulk            = 'true'
    psd_full            = 'true'
    season              = 'false'
    # ---- tracers ----
    # sld_list            = ['inrt','g2','nacl']
    # aq_list             = ['ca','k','mg','na','no3']
    # gas_list            = ['pco2']
    exrxn_list          = []
    # ---- boundary values ----    
    # pr_list             = [('inrt',1.0)]
    # rain_list           = [('ca',ca)]
    atm_list            = [('pco2',3.16e-4),('po2',0.21),('pnh3',1e-50),('pn2o',1e-50)]
    # ---- sopecify solid phase properties ----    
    sld_varlist_dust    = []
    # sld_varlist_cec     = [('inrt', cec, logkh, alpha), ('g2', cec, logkh, alpha) ]
    # sld_varlist_omrain  = [('g2',1.0),('nacl',nacl/omrain)]
    sld_varlist_kinspc  = []
    sld_varlist_2ndslds = []
    srcfile_dust        = None
    srcfile_omrain      = None
    srcfile_cec         = None
    srcfile_kinspc      = None
    srcfile_2ndslds     = './data/2ndslds_def.in'
    # ---- python stuff ----
    use_local_storage   = True
    
    spinup.run_a_scepter_run(
        runname,outdir_src,
        # ---- frame.in ----
        ztot                = ztot,
        nz                  = nz,
        ttot                = ttot,
        temp                = temp,
        fdust               = fdust,
        fdust2              = fdust2,
        taudust             = taudust,
        omrain              = omrain,
        zom                 = zom,
        poro                = poro,
        moistsrf            = moistsrf,
        zwater              = zwater,
        zdust               = zdust,
        w                   = w,
        q                   = q,
        p                   = p,
        nstep               = nstep,
        rstrt               = rstrt,
        runid               = runid,
        # ---- switches.in ----
        w_scheme            = w_scheme,
        mix_scheme          = mix_scheme,
        poro_iter           = poro_iter,
        sldmin_lim          = sldmin_lim,
        display             = display,
        disp_lim            = disp_lim,
        restart             = restart,
        rough               = rough,
        act_ON              = act_ON,
        dt_fix              = dt_fix,
        cec_on              = cec_on,
        dz_fix              = dz_fix,
        close_aq            = close_aq,
        poro_evol           = poro_evol,
        sa_evol_1           = sa_evol_1,
        sa_evol_2           = sa_evol_2,
        psd_bulk            = psd_bulk,
        psd_full            = psd_full,
        season              = season,
        # ---- tracers ----
        sld_list            = sld_list,
        aq_list             = aq_list,
        gas_list            = gas_list,
        exrxn_list          = exrxn_list,
        # ---- boundary values ----    
        pr_list             = pr_list,
        rain_list           = rain_list,
        atm_list            = atm_list,
        # ---- sopecify solid phase properties ----    
        sld_varlist_dust    = sld_varlist_dust,
        sld_varlist_cec     = sld_varlist_cec,
        sld_varlist_omrain  = sld_varlist_omrain,
        sld_varlist_kinspc  = sld_varlist_kinspc,
        sld_varlist_2ndslds = sld_varlist_2ndslds,
        srcfile_dust        = srcfile_dust,
        srcfile_omrain      = srcfile_omrain,
        srcfile_cec         = srcfile_cec,
        srcfile_kinspc      = srcfile_kinspc,
        srcfile_2ndslds     = srcfile_2ndslds,
        # ---- python stuff ----
        use_local_storage   = use_local_storage,
        )
        
        
        
def main():
    outdir_src = '../scepter_output/'
    
    runname = 'test_Pot7_25C_alpha0p0'
    runname = 'test_Pot7_25C_alpha3p1'
    # runname = 'test_Pot7_25C_alpha3p4'
    # runname = 'test_Pot7_25C_alpha3p7'
    # runname = 'test_Pot7_25C_alpha4p0'
    # runname = 'test_Pot7_25C_alpha4p3'
    # runname = 'test_Pot7_25C_alpha4p9'
    # runname = 'test_Pot7_25C_alpha5p2'
    # runname = 'test_Pot7_25C_alpha5p5'
    
    cl_tmp = 2.9e-4  # poro = 0.5
    # cl_tmp = 2.98e-4  # poro = 0.5 no act
    # cl_tmp = 4.285e-4  # poro = 0.5 no act, target pw pH = 5.41
    # cl_tmp = 2.38e-4  # poro = 0.4
    # cl_tmp = 2.8e-4  # CEC=10
    cl_tmp = 1.6e-4  # double CEC
    cl_tmp = 5.1e-5  # double CEC
    cl_tmp = 3.6e-4  # double CEC
    cl_tmp = 3.65e-4  # double CEC OM - 150, pwpH=5.9, alpha = 0 
    cl_tmp = 3.7e-4  # double CEC OM - 160, pwpH=5.8, alpha = 0 
    cl_tmp = 4.e-4  # double CEC OM - 120, pwpH=6.0, alpha = 0 
    cl_tmp = 3.7e-4  # double CEC OM - 120, pwpH=6.3, alpha = 0 
    cl_tmp = 3.1e-4  # double CEC OM - 120, pwpH=6.3, alpha = 0.5 
    cl_tmp = 2.6e-4  # double CEC OM - 120, pwpH=6.3, alpha = 1.0 
    cl_tmp = 3.05e-4  # double CEC OM - 120, pwpH=6.0, alpha = 1.0 
    cl_tmp = 3.e-4  # double CEC OM - 120, pwpH=6.1, alpha = 1.0 
    cl_tmp = 2.9e-4  # double CEC OM - 120, pwpH=6.1, alpha = 1.1 
    cl_tmp = 2.8e-4  # double CEC OM - 120, pwpH=6.1, alpha = 1.2 
    cl_tmp = 2.68e-4  # double CEC OM - 120, pwpH=6.1, alpha = 1.3 
    # cl_tmp = 2.55e-4  # double CEC OM - 120, pwpH=6.1, alpha = 1.4 
    
    ph_pw = 6.68
    ph_pw = 6.3
    ph_pw = 6.0
    ph_pw = 6.1
    # ph_pw = 6.2
    # ph_pw = 5.9
    # ph_pw = 5.8
        
    na_frac = 0.006
    k_frac = 0.016
    ca_frac = 0.56
    mg_frac = 0.168
    
    na_pw = 9.59485E-05
    k_pw = 0.000715796
    ca_pw = 0.001362459
    mg_pw = 0.000192035
    
    h_pw = 10.**-ph_pw
    
    na_val = 1
    k_val = 1
    ca_val = 2
    mg_val = 2
    al_val = 3
    
    sld_list = ['inrt','g2','amnt']
    aq_list = ['na','k','ca','mg','no3','cl']
    gas_list = ['pco2']
    pr_list = [('inrt',1.0)]
    rain_list = [('na',na_pw),('k',k_pw),('ca',ca_pw),('mg',mg_pw),('cl',cl_tmp)]
    
    omrain = 480 # 15C
    omrain = 1338 # 25C
    no3rain =  216 # lbs N/acre 
    no3rain =  no3rain/8.92179 # converting to g/m2 using 1 g/m2 = 8.92179 lbs/acre  
    no3rain =  no3rain/14  # mol N/m2/yr
    no3rain =  no3rain*80  # g NH4NO3/m2/yr
    no3rain =  no3rain/2.  # only half is required as 1 mol NH4NO3 contains 2 moles of N
    no3rain =  no3rain/omrain  # normalize against OC rain 
    # no3rain =  no3rain * 1.112  # randomly multiplying to match porewater pH 
    # no3rain =  no3rain * 1.117  # randomly multiplying to match porewater pH 
    # no3rain =  no3rain * 1.11  # randomly multiplying to match porewater pH (15C) 
    # no3rain =  no3rain * 1.08  # randomly multiplying to match porewater pH (25C)
    # no3rain =  no3rain * 1.085  # randomly multiplying to match porewater pH (25C)
    # no3rain =  no3rain * .97  # randomly multiplying to match porewater pH (25C)
    
    sld_varlist_omrain = [('g2',1.0),('amnt',no3rain)]
    
    
    alpha = 0
    alpha = 1.7
    # alpha = 3.1
    # alpha = 3.4
    alpha = 3.7
    # alpha = 4
    # alpha = 4.3
    # alpha = 4.9
    # alpha = 5.2
    # alpha = 5.5
        
    om_tgt = 4.9
        
    cec = 8.9
    
    cec_2 = 330
    cec_2 = 150
    cec_2 = 180
    cec_2 = 120
    # cec_2 = 100
    # cec_2 = 127.4
    # cec_2 = 50
    cec_1 = (100 *cec - om_tgt*cec_2)/(100 - om_tgt)
    # cec_1 = 1.
    # cec = 10.0
    
    if cec_1<0: exit("cec_1<0")
    
    
    # i_parallel = int(sys.argv[1])
    # n_parallel = int(sys.argv[2])
    
    # name_base  = sys.argv[3]
    # name_base  = 'test_Pot7_25C_v2_CEC10p0'
    name_base  = 'test_Pot7_25C_v2_noact_alpha'
    name_base  = 'test_Pot7_25C_v2_act_alpha'
    name_base  = 'test_Pot7_25C_v2_act_2alpha'
    name_base  = 'test_Pot7_25C_v2_act_2CEC{:.1f}-{:.1f}_pH{:.1f}_alpha'.format(cec_1,cec_2,ph_pw).replace('.','p')
    name_base  = 'chk_Pot7_25C_v2_act_2CEC{:.1f}-{:.1f}_pH{:.1f}_alpha'.format(cec_1,cec_2,ph_pw).replace('.','p')
    name_base  = 'chk2_Pot7_25C_v2_act_2CEC{:.1f}-{:.1f}_pH{:.1f}_alpha'.format(cec_1,cec_2,ph_pw).replace('.','p')
    # name_base  = 'test_Pot7_25C_v3_noact_alpha'
    
    act_ON = 'true'
    
    alpha_list = list(np.linspace(0.1,5.5,19))
    alpha_list = [3.4]
    alpha_list = [5.5]
    alpha_list = [(0.5, 0.5)]
    alpha_list = [(1.0, 1.0)]
    alpha_list = [(1.1, 1.1)]
    alpha_list = [(1.2, 1.2)]
    alpha_list = [(1.3, 1.3)]
    # alpha_list = [(1.4, 1.4)]
    
    for i in range(len(alpha_list)):
                   
        # if i%n_parallel!=i_parallel: continue
        
        alpha_1 = alpha_list[i][0]
        alpha_2 = alpha_list[i][1]
        
        fH_base = 0.25
        
        
        runname = name_base+'{:.1f}'.format(alpha_1).replace('.','p')
        
        logkhna_1 = ( - np.log10( na_frac**(1./na_val)*h_pw/fH_base/na_pw**(1./na_val) ) + alpha_1 * ( fH_base) ) * na_val
        logkhk_1  = ( - np.log10(  k_frac**(1./ k_val)*h_pw/fH_base/ k_pw**(1./ k_val) ) + alpha_1 * ( fH_base) ) *  k_val
        logkhca_1 = ( - np.log10( ca_frac**(1./ca_val)*h_pw/fH_base/ca_pw**(1./ca_val) ) + alpha_1 * ( fH_base) ) * ca_val
        logkhmg_1 = ( - np.log10( mg_frac**(1./mg_val)*h_pw/fH_base/mg_pw**(1./mg_val) ) + alpha_1 * ( fH_base) ) * mg_val
        logkhal_1 = 16.47       + alpha_1 * ( fH_base) * 3
        
        logkhna_2 = ( - np.log10( na_frac**(1./na_val)*h_pw/fH_base/na_pw**(1./na_val) ) + alpha_2 * ( fH_base) ) * na_val
        logkhk_2  = ( - np.log10(  k_frac**(1./ k_val)*h_pw/fH_base/ k_pw**(1./ k_val) ) + alpha_2 * ( fH_base) ) *  k_val
        logkhca_2 = ( - np.log10( ca_frac**(1./ca_val)*h_pw/fH_base/ca_pw**(1./ca_val) ) + alpha_2 * ( fH_base) ) * ca_val
        logkhmg_2 = ( - np.log10( mg_frac**(1./mg_val)*h_pw/fH_base/mg_pw**(1./mg_val) ) + alpha_2 * ( fH_base) ) * mg_val
        logkhal_2 = 16.47       + alpha_2 * ( fH_base) * 3
        
        # logkhna_1 = 4.281827003 + alpha_1 * ( fH_base)
        # logkhk_1  = 4.728609281 + alpha_1 * ( fH_base)
        # logkhca_1 = 9.542015311 + alpha_1 * ( fH_base) * 2
        # logkhmg_1 = 9.213951903 + alpha_1 * ( fH_base) * 2
        # logkhal_1 = 16.47       + alpha_1 * ( fH_base) * 3
        
        # logkhna_2 = 4.281827003 + alpha_2 * ( fH_base)
        # logkhk_2  = 4.728609281 + alpha_2 * ( fH_base)
        # logkhca_2 = 9.542015311 + alpha_2 * ( fH_base) * 2
        # logkhmg_2 = 9.213951903 + alpha_2 * ( fH_base) * 2
        # logkhal_2 = 16.47       + alpha_2 * ( fH_base) * 3
        
        # logkhna = 3.011827003 + alpha * ( fH_base)
        # logkhk  = 3.458609281 + alpha * ( fH_base)
        # logkhca = 7.002015311 + alpha * ( fH_base) * 2
        # logkhmg = 6.673951903 + alpha * ( fH_base) * 2
        # logkhal = 16.47       + alpha * ( fH_base) * 3
        
        sld_varlist_cec = [
            ('inrt',cec_1,logkhna_1,logkhk_1,logkhca_1,logkhmg_1,logkhal_1,alpha_1),
            ('g2',cec_2,logkhna_2,logkhk_2,logkhca_2,logkhmg_2,logkhal_2,alpha_2),
            ]
        
        spin_inert(outdir_src,runname,
            sld_list,
            aq_list,
            gas_list,
            pr_list,
            rain_list,
            omrain,
            sld_varlist_omrain,
            sld_varlist_cec,
            act_ON,
        )


        dep_sample = 0.15
        # (1) get porewater pH 
        phint = get_int_prof.get_ph_int_site(outdir_src,runname,dep_sample)
        
        # (2) get acidity in %
        acint = get_int_prof.get_ac_int_site_v2(outdir_src,runname,dep_sample)
        
        # (3) get bulk density
        dense = get_int_prof.get_rhobulk_int_site(outdir_src,runname,dep_sample)
        
        # (4) get field solid wt%
        sps = ['g2','inrt']
        sldwt_list = []
        for sp in sps:
            sldwt = get_int_prof.get_sldwt_int_site(outdir_src,runname,dep_sample,[sp])
            sldwt_list.append(sldwt)
        
        # (5) get SOM wt%
        omint = sldwt_list[sps.index('g2')]
        
        print(phint,acint,omint)
        
        cl_ppm = 1062
        cl_frc = 1062e-6 # g/g
        cl_wtmol = 35.453 # g/mol
        
        n_ppm = 120 
        n_frc = 122.8e-6 # g/g 
        n_wtmol = 14.0067 # g/mol
        
        s_ppm = 10 
        s_frc = 10e-6 # g/g 
        s_wtmol = 32.065 # g/mol
        
        
        for water_frac in [5,2,1,0.5]:
            # water_frac = 1.
            ztot_lab = 0.05
            poro_lab = water_frac/(1./dense+water_frac)
            cl_add = ztot_lab*(1.0 - poro_lab)*dense*1e6*cl_frc  # g/m2
            # fdust = ztot_lab*poro_lab*cl_conc*1e3*cl_wtmol  # g/m2
            cl_conc = ztot_lab*(1.0 - poro_lab)*dense*1e6*cl_frc  / (ztot_lab*poro_lab*1e3*cl_wtmol)  # mol/L
            no3_conc = ztot_lab*(1.0 - poro_lab)*dense*1e6*n_frc  / (ztot_lab*poro_lab*1e3*n_wtmol)  # mol/L
            s_conc = ztot_lab*(1.0 - poro_lab)*dense*1e6*s_frc  / (ztot_lab*poro_lab*1e3*s_wtmol)  # mol/L
            
            print(water_frac,cl_conc,no3_conc,s_conc)

   
if __name__ == '__main__':
    main()
    
    