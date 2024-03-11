import os,shutil,sys
import spinup,get_int_prof,get_soilpH_time
import numpy as np
   
   
def spin_inert(outdir_src,runname,rstrt,ttot,
    sld_list,
    aq_list,
    gas_list,
    pr_list,
    rain_list,
    omrain,
    sld_varlist_omrain,
    sld_varlist_cec,
    act_ON,
    use_local_storage,
    ):

    # outdir_src = '/storage/coda1/p-creinhard3/0/ykanzaki3/scepter_output/'
    # runname = 'test'
    
    # cec     = 10.0
    # logkh   = 5.9
    # alpha   = 3.4
    
    # na      = 1e-5
    
    # replacement of porewater 
    # volume of porewater in 1 m2 of soil = 1 (m2)* ztot(m) * poro(m3/m3)*sat(m3/m3)
    # flux of water to 1 m2 of soil = 1 (m2) * q (m/m2/yr)
    # 1 cycle of water requires (ztot * poro * sat)/q (yr)
    # N cycle = N* (ztot * poro * sat)/q (yr)
    
    # ---- frame.in ----
    ztot                = 0.5
    # nz                  = 30
    nz                  = 100
    # nz                  = 300
    ttot                = ttot
    temp                = 25
    fdust               = 0
    fdust2              = 0
    taudust             = 0
    # omrain              = 300
    zom                 = 0.5
    poro                = 0.5
    # moistsrf            = 0.28
    moistsrf            = 1.0
    zwater              = 1000
    zdust               = 0.25
    # w                   = 1e-3
    w                   = 0
    # q                   = 0.01
    # q                   = 0.1
    q                   = 1
    # q                   = 10
    p                   = 1e-5
    nstep               = 10
    rstrt               = rstrt
    runid               = runname
    # ---- switches.in ----
    w_scheme            = 0 
    # mix_scheme          = 2 # 1 --Fickian, 2 --Homogeneous
    mix_scheme          = 0 # 1 --Fickian, 2 --Homogeneous
    poro_iter           = 'true' 
    sldmin_lim          = 'true'
    display             = 'true'
    disp_lim            = 'true'
    restart             = 'false' if rstrt == 'self' else 'true'
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
    # atm_list            = [('pco2',3.16e-4),('po2',0.21),('pnh3',1e-50),('pn2o',1e-50)]
    atm_list            = [('pco2',3.16e-40),('po2',0.21),('pnh3',1e-50),('pn2o',1e-50)]
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
    # use_local_storage   = True
    
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
        
        
        
def ex_11_series():
    outdir_src = '/storage/coda1/p-creinhard3/0/ykanzaki3/scepter_output/'
    
    
    use_local_storage   = True
    use_local_storage   = False
    
    ## filling solution 
    ca_pw   = 0.6e-3
    cl_pw   = 1.2e-3
    
    # initial solution 
    na_pw   = 1.0e-3
    k_pw    = 0.2e-3
    no3_pw  = 1.2e-3
    
    # ---- from Tipping_Hurley.dat ---- 
    # Na+ + X- = NaX
    # log_k           0.0
    #
    # K+ + X- = KX
    # log_k           0.7
    #
    # H+ + X- = HX
    # log_k           1.0
    #
    # Ca+2 + 2X- = CaX2
    # log_k           0.8
    #
    # Mg+2 + 2X- = MgX2
    # log_k           0.6
    #
    # Al+3 + 3X- = AlX3
    # log_k           0.67
    # ---------------------------------
    
    # ---- how SCEPTER parameterizes ---- 
    # X-Na + H+ = Na+(aq) + X-H (logkhna)
    # where 
    #       logkhna   = logkhna_0 * gamma
    #       gamma = 10^(  alpha * f[X-H] ) 
    # 
    # assume database is defined with alpha = 0
    # logkhna   = 1.0 - 0.0 = 1.0
    # logkhk    = 1.0 - 0.7 = 0.3
    # 
    # multi valent cations 
    # X2-Ca + 2H+ = Ca++(aq) + 2X-H 
    # 
    # logkhca   = (1.0 - 0.8)*2 = 0.4
    # logkhmg   = (1.0 - 0.6)*2 = 0.8
    # logkhal   = (1.0 - 0.67)*3 = 0.99
    # ---------------------------------- 
    
    alpha = 0.0
    logkh = 1.0
    
    logkhna = logkh + alpha
    logkhk  = logkh - 0.7 + alpha
    logkhca = (logkh - 0.8 + alpha)*2.
    logkhmg = (logkh - 0.6 + alpha)*2.
    logkhal = (logkh - 0.67 + alpha)*3.
        
    
    # alpha = 50
    
    # 0.0011 mol/1kgw of exchanger 
    # assume poro m3/m3 porosity, sat as m3/m3 water saturation, 1 - poro m3/m3 of solid phase
    # poro * sat * 1,000 kg/m3 of solution assuming 1,000 kg/m3 = 1 g/cm3 solution density
    # (1-poro) * rho kg/m3 of solid given the density rho in kg/solid kg
    # assuming cec in units of cmol/kg, total exchange sites in soil is (1-poro) * rho * cec cmol/m3
    # exchange sites for a given solution mass is then given as 
    #   (1-poro) * rho * cec * 0.01 / ( poro * sat * 1,000 )  mol/kgw 
    # and this must be equal to 0.0011 mol/1kgw
    #   (1-poro) * rho * cec * 0.01 / ( poro * sat * 1,000 ) = 0.0011
    
    rho = 258.162/99.52 * 1000 # kg/m3
    # cec = 0.0011 * (0.5 * 0.28 * 1000 )/ ( (1.-0.5) * rho * 0.01 )
    cec = 0.0011 * (0.5 * 1 * 1000 )/ ( (1.-0.5) * rho * 0.01 )
    
    sld_list = ['inrt']
    aq_list = ['na','k','ca','mg','no3','cl']
    # gas_list = ['pco2']
    gas_list = []
    pr_list = [('inrt',1.0)]
    
    omrain = 0 
    
    no3rain =  0.0 
    
    sld_varlist_omrain = [('g2',1.0),('amnt',no3rain)]
    
    
    act_ON = 'true'
    act_ON = 'false'
    
    rstrt = 'self'
    runname = 'ex11_init_nofh_disp_n100_alpha0'
    ttot = 1.
    rain_list = [('na',na_pw), ('k',k_pw), ('no3',no3_pw)]
    
    
    # rstrt = 'ex11_init_nofh_disp_n100_alpha0'
    # runname = 'ex11_replace_nofh_disp_n100_alpha0_noact_q1p0'
    # ttot = 0.63
    # rain_list = [('ca',ca_pw), ('cl',cl_pw)]
    
    sld_varlist_cec = [
        ('inrt',cec,logkhna,logkhk,logkhca,logkhmg,logkhal,alpha),
        ]
    
    spin_inert(outdir_src,runname,rstrt,ttot,
        sld_list,
        aq_list,
        gas_list,
        pr_list,
        rain_list,
        omrain,
        sld_varlist_omrain,
        sld_varlist_cec,
        act_ON,
        use_local_storage,
        )
   
def main():
    ex_11_series()
   
if __name__ == '__main__':
    main()
    
    