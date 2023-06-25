import os
import make_inputs
import shutil
    
def run_closed_inert_cec(cec,act_ON_IN,data):

    use_local_storage = True
    use_local_storage = False
    
    # cec = 0
    
    logkhna = 5.9
    logkhk  = 4.8
    logkhca = 10.47
    logkhmg = 10.786
    logkhal = 16.47
    
    alpha = 3.4
    
    outdir = '../scepter_output/'
    if use_local_storage:  outdir = os.environ['TMPDIR'] + '/scepter_output/'
    # runname = 'test_inert_spinup_cec_{:d}'.format(cec) 
    if act_ON_IN: runname = 'test_inert_spinup_act_cec_{:d}'.format(cec) + '_'+ data
    else: runname = 'test_inert_spinup_noact_cec_{:d}'.format(cec) + '_'+ data

    # compile 
    exename = 'scepter'
    exename_src = 'scepter'
    # exename_src = 'scepter_test'
    to = ' '
    where = '/'
    os.system('make')
    # os.system('make --file=makefile_test')
    if not os.path.exists( outdir + runname) : os.system('mkdir -p ' + outdir + runname)
    os.system('cp ' + exename_src + to + outdir + runname + where + exename)
    
    ztot=0.5
    nz=30
    ttot=1e5
    temp=15
    fdust=0
    fdust2=0
    taudust=0
    omrain=0
    zom=ztot
    # poro=0.5
    poro=0.999
    # moistsrf=0.1
    moistsrf=1
    zwater=1000
    zdust=ztot/2.
    # w=1e-5
    # q=0.1
    w=0
    q=0
    p=1e-5
    nstep=10
    rstrt='self'
    runid=runname

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
    
    w_scheme=0
    mix_scheme=0 
    poro_iter='true' 
    sldmin_lim ='true'
    display='true'
    disp_lim='true'
    restart ='false'
    rough      ='true'
    if act_ON_IN: act_ON ='true'
    else: act_ON ='false'
    dt_fix='false'
    cec_on='true'
    dz_fix='true'
    # close_aq='false'
    close_aq='true'
    poro_evol='true'
    sa_evol_1 ='true'
    sa_evol_2='false'
    psd_bulk='true'
    psd_full='true'
    season='false'
        
    make_inputs.get_input_switches(
        outdir=outdir
        ,runname=runname
        ,w_scheme=w_scheme
        ,mix_scheme=mix_scheme  
        ,poro_iter=poro_iter 
        ,sldmin_lim=sldmin_lim 
        ,display=display
        ,disp_lim=disp_lim
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
    
    sld_list=['inrt']
    aq_list = ['ca','k','mg','na']
    gas_list = []
    exrxn_list = []
    make_inputs.get_input_tracers(
        outdir=outdir
        ,runname=runname
        ,sld_list = sld_list
        ,aq_list = aq_list
        ,gas_list = gas_list
        ,exrxn_list = exrxn_list
        )
        
    pr_list = [('inrt',1.0)]
    # rain_list = [('ca',5.0e-6)]
    rain_list = []
    # atm_list = [('pco2',3.16e-4),('po2',0.21),('pnh3',1e-50),('pn2o',1e-50)]
    atm_list = [('pco2',1e-20),('po2',0.21),('pnh3',1e-50),('pn2o',1e-50)]
    make_inputs.get_input_tracer_bounds(
        outdir=outdir
        ,runname=runname
        ,pr_list = pr_list
        ,rain_list=rain_list
        ,atm_list=atm_list
        )
        
    filename = 'dust.in'
    srcfile = './data/dust_gbasalt.in'
    sld_varlist =[]
    make_inputs.get_input_sld_properties(
        outdir=outdir
        ,runname=runname
        ,filename = filename
        # ,srcfile = srcfile
        ,sld_varlist=sld_varlist
        )
        
    filename = 'cec.in'
    sld_varlist = [('inrt',cec, logkhna, logkhk, logkhca, logkhmg, logkhal, alpha) ] 
    make_inputs.get_input_sld_properties(
        outdir=outdir
        ,runname=runname
        ,filename = filename
        ,sld_varlist=sld_varlist
        )
        
    filename = 'OM_rain.in'
    sld_varlist = [ ('g2',1) ] 
    make_inputs.get_input_sld_properties(
        outdir=outdir
        ,runname=runname
        ,filename = filename
        ,sld_varlist=sld_varlist
        )
        
    # run 
    os.system(outdir+runname+where+exename)
    
        
    if use_local_storage:
        src = outdir + runname 
        dst = '../scepter_output/' + runname 
        
        if not os.path.exists(dst): 
            shutil.copytree(src, dst)
        else:
            shutil.rmtree(dst)
            shutil.copytree(src, dst)
   
def main():
    
    cec_list = list(range(21))
    act_ON_IN = False
    act_ON_IN = True
    
    data = 'Sikora'
    
    for cec in cec_list:
        run_closed_inert_cec(cec,act_ON_IN,data)
        
        
        
   
if __name__ == '__main__':
    main()
    
    