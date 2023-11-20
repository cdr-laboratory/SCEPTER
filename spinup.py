import os,shutil,sys,subprocess
import make_inputs
import numpy as np
    
def run_a_scepter_run(
    runname,outdir_src,
    **kwargs
    ):
    # ---- frame.in ----
    ztot                = kwargs.get('ztot',                0.5)
    nz                  = kwargs.get('nz',                  30)
    ttot                = kwargs.get('ttot',                1e5)
    temp                = kwargs.get('temp',                15)
    fdust               = kwargs.get('fdust',               0)
    fdust2              = kwargs.get('fdust2',              0)
    taudust             = kwargs.get('taudust',             0)
    omrain              = kwargs.get('omrain',              900)
    zom                 = kwargs.get('zom',                 0.25)
    poro                = kwargs.get('poro',                0.5)
    moistsrf            = kwargs.get('moistsrf',            0.5)
    zwater              = kwargs.get('zwater',              1000)
    zdust               = kwargs.get('zdust',               0.25)
    w                   = kwargs.get('w',                   1e-3)
    q                   = kwargs.get('q',                   0.3)
    p                   = kwargs.get('p',                   1e-5)
    nstep               = kwargs.get('nstep',               10)
    rstrt               = kwargs.get('rstrt',               'self')
    runid               = kwargs.get('runid',               runname)
    # ---- switches.in ----
    w_scheme            = kwargs.get('w_scheme',            0)
    mix_scheme          = kwargs.get('mix_scheme',          0) 
    poro_iter           = kwargs.get('poro_iter',           'true') 
    sldmin_lim          = kwargs.get('sldmin_lim',          'true')
    display             = kwargs.get('display',             1)
    report              = kwargs.get('report',              0)
    restart             = kwargs.get('restart',             'false')
    rough               = kwargs.get('rough',               'true')
    act_ON              = kwargs.get('act_ON',              'true')
    dt_fix              = kwargs.get('dt_fix',              'false')
    cec_on              = kwargs.get('cec_on',              'true')
    dz_fix              = kwargs.get('dz_fix',              'true')
    close_aq            = kwargs.get('close_aq',            'false')
    poro_evol           = kwargs.get('poro_evol',           'true')
    sa_evol_1           = kwargs.get('sa_evol_1',           'true')
    sa_evol_2           = kwargs.get('sa_evol_2',           'false')
    psd_bulk            = kwargs.get('psd_bulk',            'true')
    psd_full            = kwargs.get('psd_full',            'true')
    season              = kwargs.get('season',              'false')
    # ---- tracers ----
    sld_list            = kwargs.get('sld_list',            ['inrt','g2'])
    aq_list             = kwargs.get('aq_list',             ['ca','k','mg','na'])
    gas_list            = kwargs.get('gas_list',            ['pco2'])
    exrxn_list          = kwargs.get('exrxn_list',          [])
    # ---- boundary values ----    
    pr_list             = kwargs.get('pr_list',             [('inrt',1.0)])
    rain_list           = kwargs.get('rain_list',           [('ca',5.0e-6)])
    atm_list            = kwargs.get('atm_list',            [('pco2',3.16e-4),('po2',0.21),('pnh3',1e-50),('pn2o',1e-50)])
    # ---- sopecify solid phase properties ----    
    sld_varlist_psdpr   = kwargs.get('sld_varlist_psdpr',   [])
    sld_varlist_dust    = kwargs.get('sld_varlist_dust',    [])
    sld_varlist_cec     = kwargs.get('sld_varlist_cec',     [('inrt', 10, 5.9, 4.8, 10.47, 10.786, 16.47, 3.4) ])
    sld_varlist_omrain  = kwargs.get('sld_varlist_omrain',  [('g2',1.0)])
    sld_varlist_kinspc  = kwargs.get('sld_varlist_kinspc',  [])
    sld_varlist_keqspc  = kwargs.get('sld_varlist_keqspc',  [])
    sld_varlist_sa      = kwargs.get('sld_varlist_sa',      [])
    sld_varlist_nopsd   = kwargs.get('sld_varlist_nopsd',   [])
    sld_varlist_2ndslds = kwargs.get('sld_varlist_2ndslds', [])
    srcfile_psdpr       = kwargs.get('srcfile_psdpr',       None)
    srcfile_dust        = kwargs.get('srcfile_dust',        None)
    srcfile_omrain      = kwargs.get('srcfile_omrain',      None)
    srcfile_cec         = kwargs.get('srcfile_cec',         None)
    srcfile_kinspc      = kwargs.get('srcfile_kinspc',      None)
    srcfile_keqspc      = kwargs.get('srcfile_keqspc',      None)
    srcfile_sa          = kwargs.get('srcfile_sa',          None)
    srcfile_nopsd       = kwargs.get('srcfile_nopsd',       None)
    srcfile_2ndslds     = kwargs.get('srcfile_2ndslds',     './data/2ndslds_def.in')
    # ---- seasonality properties ----
    T_temp              = kwargs.get('T_temp',              [ list(1./np.linspace(12,1,12)),[15]*12 ] )
    moist_temp          = kwargs.get('moist_temp',          [ list(1./np.linspace(12,1,12)),[0.3]*12 ] )
    q_temp              = kwargs.get('q_temp',              [ list(1./np.linspace(12,1,12)),[0.5]*12 ])
    dust_temp           = kwargs.get('dust_temp',           [ list(1./np.linspace(12,1,12)),[0]*12 ])
    
    
    # ---- python stuff ----
    use_local_storage   = kwargs.get('use_local_storage',   True)
    lim_calc_time       = kwargs.get('lim_calc_time',       False)
    max_calc_time       = kwargs.get('max_calc_time',       20)
    
    outdir = outdir_src
    if use_local_storage:  outdir = os.environ['TMPDIR'] + '/scepter_output/'
    
    # compile 
    exename = 'scepter'
    exename_src = 'scepter_DEV'
    # exename_src = 'scepter_test'
    to = ' '
    where = '/'
    os.system('make')
    # os.system('make --file=makefile_test')
    if not os.path.exists( outdir + runname) : os.system('mkdir -p ' + outdir + runname)
    os.system('cp ' + exename_src + to + outdir + runname + where + exename)

    make_inputs.get_input_frame(
        outdir=outdir
        ,runname=runname
        ,ztot=ztot
        ,nz=nz
        ,ttot=ttot
        ,temp=temp
        ,fdust=fdust
        ,fdust2=fdust2
        ,taudust=taudust
        ,omrain=omrain
        ,zom=zom
        ,poro=poro
        ,moistsrf=moistsrf
        ,zwater=zwater
        ,zdust=zdust
        ,w=w
        ,q=q
        ,p=p
        ,nstep=nstep
        ,rstrt=rstrt
        ,runid=runid
        )
        
    make_inputs.get_input_switches(
        outdir=outdir
        ,runname=runname
        ,w_scheme=w_scheme
        ,mix_scheme=mix_scheme  
        ,poro_iter=poro_iter 
        ,sldmin_lim=sldmin_lim 
        ,display=display
        ,report=report
        ,restart=restart 
        ,rough=rough      
        ,act_ON=act_ON 
        ,dt_fix=dt_fix
        ,cec_on=cec_on
        ,dz_fix=dz_fix
        ,close_aq=close_aq
        ,poro_evol=poro_evol
        ,sa_evol_1=sa_evol_1 
        ,sa_evol_2=sa_evol_2
        ,psd_bulk=psd_bulk
        ,psd_full=psd_full
        ,season=season
        )
        
    make_inputs.get_input_tracers(
        outdir=outdir
        ,runname=runname
        ,sld_list = sld_list
        ,aq_list = aq_list
        ,gas_list = gas_list
        ,exrxn_list = exrxn_list
        )
        
    make_inputs.get_input_tracer_bounds(
        outdir=outdir
        ,runname=runname
        ,pr_list = pr_list
        ,rain_list=rain_list
        ,atm_list=atm_list
        )
        
    sldvar_list_all = [
        ('psdpr.in',        srcfile_psdpr,      sld_varlist_psdpr),
        ('dust.in',         srcfile_dust,       sld_varlist_dust),
        ('cec.in',          srcfile_cec,        sld_varlist_cec),
        ('OM_rain.in',      srcfile_omrain,     sld_varlist_omrain),
        ('kinspc.in',       srcfile_kinspc,     sld_varlist_kinspc),
        ('keqspc.in',       srcfile_keqspc,     sld_varlist_keqspc),
        ('sa.in',           srcfile_sa,         sld_varlist_sa),
        ('nopsd.in',        srcfile_nopsd,      sld_varlist_nopsd),
        ('2ndslds.in',      srcfile_2ndslds,    sld_varlist_2ndslds),
        ]
    for i in range(len(sldvar_list_all)):
        filename    = sldvar_list_all[i][0]
        srcfile     = sldvar_list_all[i][1]
        sld_varlist = sldvar_list_all[i][2]
        if srcfile!=None:
            make_inputs.get_input_sld_properties(
                outdir=outdir
                ,runname=runname
                ,filename = filename
                ,srcfile = srcfile
                )
        else:
            make_inputs.get_input_sld_properties(
                outdir=outdir
                ,runname=runname
                ,filename = filename
                ,sld_varlist=sld_varlist
                )
    
    if season=='true':
        make_inputs.get_input_climate_temp(
            outdir=outdir,
            runname=runname,
            T_temp = T_temp,
            moist_temp = moist_temp,
            q_temp = q_temp,
            dust_temp = dust_temp,
            )
        
    # >>>> run 
    
    run_success = False
    
    if not lim_calc_time:
        os.system(outdir+runname+where+exename + ' > ' + outdir+runname+'/logfile.txt' + ' 2> ' + outdir+runname+'/err.txt')
        run_success = True
    else:
        logf = open(outdir+runname+'/logfile.txt', 'w')
        logerrf = open(outdir+runname+'/err.txt', 'w')
        proc = subprocess.Popen([outdir+runname+where+exename], stdout=logf, stderr=logerrf)

        my_timeout =60*max_calc_time
        
        try:
            proc.wait(my_timeout)
            print('run finished within {:f} min'.format(int(my_timeout/60.)))
            run_success = True
        except subprocess.TimeoutExpired:
            proc.kill()
            print('run UNfinished within {:f} min'.format(int(my_timeout/60.)))
    
        
    if use_local_storage:
        src = outdir + runname 
        dst = outdir_src + runname 
        
        if not os.path.exists(dst): 
            shutil.copytree(src, dst)
        else:
            shutil.rmtree(dst)
            shutil.copytree(src, dst)
   
    return run_success
        
def main():

    outdir_src = '../scepter_output/tests/'
    runname = 'test_richards'
    runname = 'test_richards_2'
    runname = 'test_richards_3'
    runname = 'test_richards_4'
    runname = 'test_richards_4_v2'
    runname = 'test_richards_4_v4'
    
    #  >>>> input variables of interests 
    cec     = 10.0
    logkhna = 5.9
    logkhk  = 4.8
    logkhca = 10.47
    logkhmg = 10.786
    logkhal = 16.47
    alpha   = 3.4
    
    ca      = 1e-5
    
    # >>>> define input variables written in input files 
    # ---- frame.in ----
    ztot                = 1.5
    nz                  = 30
    ttot                = 1e1
    temp                = 15
    fdust               = 0
    fdust2              = 0
    taudust             = 0
    omrain              = 300
    zom                 = 0.25
    # poro                = 0.5
    poro                = 0.396
    moistsrf            = 0.5
    zwater              = 1000
    zdust               = 0.25
    w                   = 1e-3
    q                   = 0.3
    p                   = 1e-5
    nstep               = 10
    rstrt               = 'self'
    runid               = runname
    # ---- switches.in ----
    w_scheme            = 1 
    mix_scheme          = 1 # 1 --Fickian
    poro_iter           = 'false' 
    sldmin_lim          = 'true'
    display             = 1
    report              = 0
    restart             = 'false'
    rough               = 'true'
    act_ON              = 'true'
    dt_fix              = 'false'
    cec_on              = 'true'
    dz_fix              = 'true'
    close_aq            = 'false'
    poro_evol           = 'true'
    sa_evol_1           = 'true'
    sa_evol_2           = 'false'
    psd_bulk            = 'true'
    psd_full            = 'true'
    # season              = 'false'
    season              = 'true'
    # ---- tracers ----
    sld_list            = ['inrt','g2']
    aq_list             = ['ca','k','mg','na']
    gas_list            = ['pco2']
    exrxn_list          = []
    # ---- boundary values ----    
    pr_list             = [('inrt',1.0)]
    rain_list           = [('ca',ca)]
    atm_list            = [('pco2',3.16e-4),('po2',0.21),('pnh3',1e-50),('pn2o',1e-50)]
    # ---- sopecify solid phase properties ----    
    sld_varlist_dust    = []
    sld_varlist_cec     = [('inrt', cec, logkhna, logkhk, logkhca, logkhmg, logkhal, alpha), ('g2', cec, logkhna, logkhk, logkhca, logkhmg, logkhal, alpha)]
    sld_varlist_omrain  = [('g2',1.0)]
    sld_varlist_kinspc  = []
    sld_varlist_2ndslds = []
    srcfile_dust        = None
    srcfile_omrain      = None
    srcfile_cec         = None
    srcfile_kinspc      = None
    srcfile_2ndslds     = './data/2ndslds_def.in'
    # ---- seasonality properties ----
    q_temp              = list(np.loadtxt('../openRE/infiltrationproblem/input/infiltration.dat',skiprows=1,delimiter=',',usecols=1)/1000.*365.25)
    time_list           = list((np.arange(len(q_temp)))/365.25)
    T_temp              = [time_list,[temp]*len(time_list)]
    moist_temp          = [time_list,[moistsrf]*len(time_list)]
    dust_temp           = [time_list,[0]*len(time_list)]
    q_temp              = [time_list,q_temp]
    # ---- python stuff ----
    use_local_storage   = False
    
    
    
    #  >>>> run the code 
    
    run_a_scepter_run(
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
        report              = report,
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
        # ---- seasonality properties ----
        T_temp              = T_temp,
        moist_temp          = moist_temp,
        q_temp              = q_temp,
        dust_temp           = dust_temp,
        # ---- python stuff ----
        use_local_storage   = use_local_storage,
        )


   
if __name__ == '__main__':
    main()
    
    