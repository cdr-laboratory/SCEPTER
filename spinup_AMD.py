import os
import numpy as np
import get_int_prof
import make_inputs
import get_inputs
import time
import sys
from datetime import datetime
import shutil
import subprocess
import json
import re
import copy
import random
from multiprocessing import Pool
# from multiprocessing.dummy import Pool
import multiprocessing, logging
mpl = multiprocessing.log_to_stderr()
# mpl.setLevel(logging.INFO)
mpl.setLevel(multiprocessing.SUBDEBUG)

def get_input(outdir,runtype,nz,scheme_act,ver): 
    
    true = 'true'
    false = 'false'
    
    # default values for cation exchange parameters
    cec = 0.
    logkh = 5.9 
    alpha = 3.4
    
    logkhna = logkh
    logkhk  = logkh - 1.1
    logkhca = (logkh - 0.665)*2.
    logkhmg = (logkh - 0.507)*2.
    logkhal = (logkh - 0.41)*3.
    
    # frame.in 
    ztot=5
    ttot=10
    temp=25
    fdust=0
    fdust2=0
    taudust=0
    omrain=600
    zom=0.25
    poro=0.5
    moistsrf=0.5
    zwater=2.5
    zdust=0.25
    w=0
    q=0.3 # restart 
    # q=0 # initial/boundary condition test
    p=1e-5
    nstep=10
    runname_restart='self'
    
    # switches.in 
    w_scheme=0
    mix_scheme=0
    # mix_scheme=1
    poro_iter=false
    sldmin_lim=true 
    display=1
    report=1
    restart=false 
    # restart=true # restart exp
    rough=true
    act_ON=true 
    # dt_fix=false
    dt_fix=true
    # cec_on=true
    cec_on=false
    dz_fix=true
    close_aq=false
    # close_aq=true # initial/boundary condition test
    # poro_evol=true
    poro_evol=false
    sa_evol_1=true
    # poro_evol=false
    # sa_evol_1=false
    sa_evol_2=false
    psd_bulk=false
    psd_full=false
    season=false
    
    # tracers 
    sld_list= ['inrt','py',]
    aq_list= ['fe2','so4','ca','k','al','cl','si',]
    gas_list= ['po2','pco2']
    exrxn_list= []
    
    # boundary values for tracers
    
    mv_py = 0.239400E+02 # cm3/mol
    mv_inrt = 0.995200E+02 # cm3/mol 
    
    py_mol_m3_0 = 0.002/(mv_py*1e-6) # m3_min/m3_bulk /(m3_min/mol_min) = mol/m3_bulk
    inrt_mol_m3_0 = (1.-0.002-0.5)/(mv_inrt*1e-6) # m3_min/m3_bulk /(m3_min/mol_min) = mol/m3_bulk
    
    mwt_py = 0.119967E+03 # g/mol
    mwt_inrt = 0.258162E+03 # g/mol
    
    py_wt_frc = py_mol_m3_0* mwt_py/(py_mol_m3_0* mwt_py + inrt_mol_m3_0*mwt_inrt)
    inrt_wt_frc = inrt_mol_m3_0*mwt_inrt/(py_mol_m3_0* mwt_py + inrt_mol_m3_0*mwt_inrt)
    
    pr_list=[('inrt',inrt_wt_frc),('py',py_wt_frc)]
    
    
    rain_list=[('ca',1.390e-03),('k',8.525e-03),('al',2.078e-09),('cl',1.14e-4),('si',1.99e-4),('so4',5.436592e-03),('fe2',8.453e-14),] # boundary condition obtained by phreeqc
    
    rain_list_bnd_bdot_v2_names=[
         'fe2',   'so4',   'ca',    'k',     'al',    'cl',    'si',
     ]
    rain_list_bnd_bdot_v4_concs=[
        8.4777990879914060E-014,   5.4675092162629665E-003,   1.3605088849070975E-003,   8.5101614181632614E-003,   1.9639800852907899E-009,   1.1400000000000001E-004,   1.9899663665549360E-004,
     ]
    rain_list_bnd_bdot_v4 = [ (nm,conc) for nm, conc in zip(rain_list_bnd_bdot_v2_names,rain_list_bnd_bdot_v4_concs) ]
    
    rain_list_init_bdot_v4_concs=[
        9.9203528766392195E-005,   1.1992520311594751E-002,   9.2252065676718051E-003,   8.6805043211892231E-003,   4.0251554217574911E-013,   1.1400000000000000E-003,   1.9265158561408754E-003, 
    ]
    rain_list_init_bdot_v4 = [ (nm,conc) for nm, conc in zip(rain_list_bnd_bdot_v2_names,rain_list_init_bdot_v4_concs) ]
    atm_list=[('pco2',3.17e-4),('po2',0.21),('pnh3',0),('pn2o',0)] # bounary conditions 
    
    atm_list_init_bdot_v4 = [('pco2',1.1192392171660575E-002),('po2',7.9974341617857271E-066),('pnh3',0),('pn2o',0)]
    
    # cec
    sld_cec_list=[
        ('inrt',cec, logkhna, logkhk, logkhca, logkhmg, logkhal, alpha) ,
        ]
    
    # OM rain 
    sld_omrain_list=[('g2',1)]
    
    # 2ndary phases
    # srcfile_2ndslds='./data/2ndslds_def.in'
    srcfile_2ndslds='./data/2ndslds_incl_py.in'
    
    # psd rain
    srcfile_psdrain='./data/psdrain_320um.in'
    
    # kinspc
    sld_kinspc_list=[('cc',1e-5)]
    # sld_kinspc_list=[('py',0)] # initial condition 
    
    if runtype == 'init':
        q=0 # initial/boundary condition test
        dt_fix=false
        # dt_fix=true
        sld_kinspc_list=[('py',0)] # initial condition 
        close_aq=true # initial/boundary condition test
        restart=false 
        rain_list=[('ca',9.626e-03),('k',8.660e-03),('al',3.798e-13),('cl',1.140e-03),('si',1.926e-03),('so4',1.260395e-02),('fe2',9.115e-05),] # initial condition obtained by phreeqc
        
        if scheme_act == 'davies':
            rain_list=[
                ('ca',9.626e-03),('k',8.660e-03),('al',3.798e-13),('cl',1.140e-03),
                ('si',1.926e-03),('so4',1.263405e-02),('fe2',9.115e-05),
                ] # initial condition obtained by phreeqc
        elif scheme_act == 'bdot':
            rain_list=rain_list_init_bdot_v4
            atm_list=atm_list_init_bdot_v4
            
        
    elif runtype == 'bnd':
        q=0 # initial/boundary condition test
        dt_fix=false
        # dt_fix=true
        sld_kinspc_list=[('py',0)] # initial condition 
        close_aq=true # initial/boundary condition test
        restart=false 
        rain_list=[('ca',1.390e-03),('k',8.525e-03),('al',2.078e-09),('cl',1.14e-4),('si',1.99e-4),('so4',5.436592e-03),('fe2',8.453e-14),] # boundary condition obtained by phreeqc
        
        if scheme_act == 'davies':
            rain_list=[
                ('ca',1.390e-03),('k',8.525e-03),('al',2.078e-09),('cl',1.14e-4),
                ('si',1.99e-4),('so4',5.5029e-03),('fe2',8.453e-14),
                ] 
        elif scheme_act == 'bdot':
            rain_list=rain_list_bnd_bdot_v4
                
        atm_list=[('pco2',3.17e-4),('po2',0.21),('pnh3',0),('pn2o',0)] # bounary conditions 
    if runtype == 'init_flow':
        q=1e70 # initial/boundary condition test
        ttot = 1e5
        dt_fix=false
        # dt_fix=true
        sld_kinspc_list=[('py',0)] # initial condition 
        # close_aq=true # initial/boundary condition test
        restart=false 
        rain_list=[('ca',9.626e-03),('k',8.660e-03),('al',3.798e-13),('cl',1.140e-03),('si',1.926e-03),('so4',1.260395e-02),('fe2',9.115e-05),] # initial condition obtained by phreeqc
        atm_list=[('pco2',1.1241006093342163E-002),('po2',7.9969554948625034E-066),('pnh3',0),('pn2o',0)] # initial conditions 
        
        if scheme_act == 'davies':
            rain_list=[
                ('ca',9.626e-03),('k',8.660e-03),('al',3.798e-13),('cl',1.140e-03),
                ('si',1.926e-03),('so4',1.263405e-02),('fe2',9.115e-05),
                ] # initial condition obtained by phreeqc
        elif scheme_act == 'bdot':
            rain_list=[
                ('ca',9.2826540170040262E-003),('k',8.6810294623506715E-003),('al',4.0291153578108660E-013),('cl',1.1400000000000000E-003),
                ('si',1.9265157003896770E-003),('so4',1.1978E-002),('fe2',9.9269505129143321E-005),
                ] 
            atm_list=[('pco2',1.1241006093342163E-002),('po2',7.585E-066),('pnh3',0),('pn2o',0)] # initial conditions
                
    elif runtype == 'bnd_flow':
        q=1e70 # initial/boundary condition test
        ttot = 1e5
        dt_fix=false
        # dt_fix=true
        sld_kinspc_list=[('py',0)] # initial condition 
        # close_aq=true # initial/boundary condition test
        restart=false 
        rain_list=[('ca',1.390e-03),('k',8.525e-03),('al',2.078e-09),('cl',1.14e-4),('si',1.99e-4),('so4',5.436592e-03),('fe2',8.453e-14),] # boundary condition obtained by phreeqc
        
        if scheme_act == 'davies':
            rain_list=[
                ('ca',1.390e-03),('k',8.525e-03),('al',2.078e-09),('cl',1.14e-4),
                ('si',1.99e-4),('so4',5.5029e-03),('fe2',8.453e-14),
                ] 
        elif scheme_act == 'bdot':
            rain_list=[
                ('ca',1.3605117072790210E-003),('k',8.5101618593205978E-003),('al',1.9631309897300844E-009),('cl',1.1400000000000001E-004),
                ('si',1.9899662924434806E-004),('so4',5.46751E-003),('fe2',8.4743613380706646E-014),
                ] 
                
        atm_list=[('pco2',3.17e-4),('po2',0.21),('pnh3',0),('pn2o',0)] # bounary conditions 
    elif runtype == 'rstrt_norxn':
        # q=0.3 
        # q=0.3 * 0.5
        # dt_fix=true
        dt_fix=false
        close_aq=false # initial/boundary condition test
        sld_kinspc_list=[('py',0)] # initial condition 
        rain_list=[('ca',1.390e-03),('k',8.525e-03),('al',2.078e-09),('cl',1.14e-4),('si',1.99e-4),('so4',5.436592e-03),('fe2',8.453e-14),] # boundary condition obtained by phreeqc
        
        if scheme_act == 'davies':
            rain_list=[
                ('ca',1.390e-03),('k',8.525e-03),('al',2.078e-09),('cl',1.14e-4),
                ('si',1.99e-4),('so4',5.5029e-03),('fe2',8.453e-14),
                ] 
        elif scheme_act == 'bdot':
            rain_list=rain_list_bnd_bdot_v4
        
        atm_list=[('pco2',3.17e-4),('po2',0.21),('pnh3',0),('pn2o',0)] # bounary conditions 
        restart=true 
        runname_restart=f'{outdir}AMD/test_init_{nz:d}{ver}'
    elif runtype == 'nonrstrt':
        # q=0.3 * 0.5
        dt_fix=true
        # dt_fix=false
        # nstep=1000
        close_aq=false # initial/boundary condition test
        rain_list=[('ca',1.390e-03),('k',8.525e-03),('al',2.078e-09),('cl',1.14e-4),('si',1.99e-4),('so4',5.436592e-03),('fe2',8.453e-14),] # boundary condition obtained by phreeqc
        
            
        if scheme_act == 'davies':
            rain_list=[
                ('ca',1.390e-03),('k',8.525e-03),('al',2.078e-09),('cl',1.14e-4),
                ('si',1.99e-4),('so4',5.5029e-03),('fe2',8.453e-14),
                ] 
        elif scheme_act == 'bdot':
            rain_list=rain_list_bnd_bdot_v4
                    
        atm_list=[('pco2',3.17e-4),('po2',0.21),('pnh3',0),('pn2o',0)] # bounary conditions 
        restart=false 
    elif runtype == 'rstrt':
        # q=0.3 
        # q=0.3 * 0.5
        dt_fix=true
        # dt_fix=false
        # nstep=1000
        close_aq=false # initial/boundary condition test
        rain_list=[('ca',1.390e-03),('k',8.525e-03),('al',2.078e-09),('cl',1.14e-4),('si',1.99e-4),('so4',5.436592e-03),('fe2',8.453e-14),] # boundary condition obtained by phreeqc
        
            
        if scheme_act == 'davies':
            rain_list=[
                ('ca',1.390e-03),('k',8.525e-03),('al',2.078e-09),('cl',1.14e-4),
                ('si',1.99e-4),('so4',5.5029e-03),('fe2',8.453e-14),
                ] 
        elif scheme_act == 'bdot':
            rain_list=rain_list_bnd_bdot_v4
                    
        atm_list=[('pco2',3.17e-4),('po2',0.21),('pnh3',0),('pn2o',0)] # bounary conditions 
        restart=true 
        runname_restart=f'{outdir}AMD/test_init_{nz:d}{ver}'
        
        
    # no psd phases 
    sld_nopsd_list=[('cc',true)]
    
    excludes = ['true','false','logkh','cec','logkhna', 'logkhk', 'logkhca', 'logkhmg', 'logkhal', 'alpha','excludes']
    
    input_dict = {k: v for k, v in locals().items() if k not in excludes}
    
    print(input_dict)
    
    return input_dict


def setup_run(input_dict):
    
    # input_dict should be disctionary 
    # for i in input_dict.keys():
        # if isinstance(input_dict[i],list): 
            # exec("%s = []" % (i), globals() )
            # for j in input_dict[i]:
                # if isinstance(j,str): exec("%s.append('%s')" % (i,j), globals() )
                # elif isinstance(j,int): exec("%s.append(%d)" % (i,j), globals() )
                # elif isinstance(j,float): exec("%s.append(%f)" % (i,j) , globals() )
        # else:
            # if isinstance(input_dict[i],str): exec("%s = '%s'" % (i,input_dict[i]), globals() )
            # elif isinstance(input_dict[i],float): exec("%s = %f" % (i,input_dict[i]), globals() )
            # elif isinstance(input_dict[i],int): exec("%s = %d" % (i,input_dict[i]) , globals() )
    
    # print(input_dict['ztot'])
    
    make_inputs.get_input_frame(
        outdir=input_dict['outdir'],
        runname=input_dict['runname'],
        ztot=input_dict['ztot'],
        nz=input_dict['nz'],
        ttot=input_dict['ttot'],
        temp=input_dict['temp'],
        fdust=input_dict['fdust'],
        fdust2=input_dict['fdust2'],
        taudust=input_dict['taudust'],
        omrain=input_dict['omrain'],
        zom=input_dict['zom'],
        poro=input_dict['poro'],
        moistsrf=input_dict['moistsrf'],
        zwater=input_dict['zwater'],
        zdust=input_dict['zdust'],
        w=input_dict['w'],
        q=input_dict['q'],
        p=input_dict['p'],
        nstep=input_dict['nstep'],
        rstrt=input_dict['runname_restart'],
        runid=input_dict['runname'],
        )
        
    make_inputs.get_input_switches(
        outdir=input_dict['outdir'],
        runname=input_dict['runname'],
        w_scheme=input_dict['w_scheme'],
        mix_scheme=input_dict['mix_scheme'],
        poro_iter=input_dict['poro_iter'],
        sldmin_lim=input_dict['sldmin_lim'], 
        display=input_dict['display'],
        report=input_dict['report'],
        restart=input_dict['restart'], 
        rough=input_dict['rough'],
        act_ON=input_dict['act_ON'], 
        dt_fix=input_dict['dt_fix'],
        cec_on=input_dict['cec_on'],
        dz_fix=input_dict['dz_fix'],
        close_aq=input_dict['close_aq'],
        poro_evol=input_dict['poro_evol'],
        sa_evol_1=input_dict['sa_evol_1'],
        sa_evol_2=input_dict['sa_evol_2'],
        psd_bulk=input_dict['psd_bulk'],
        psd_full=input_dict['psd_full'],
        season=input_dict['season'],
        )
        
    make_inputs.get_input_tracers(
        outdir=input_dict['outdir'],
        runname=input_dict['runname'],
        sld_list=input_dict['sld_list'],
        aq_list=input_dict['aq_list'],
        gas_list=input_dict['gas_list'],
        exrxn_list=input_dict['exrxn_list'],
        )
    
    make_inputs.get_input_tracer_bounds(
        outdir=input_dict['outdir'],
        runname=input_dict['runname'],
        pr_list=input_dict['pr_list'],
        rain_list=input_dict['rain_list'],
        atm_list=input_dict['atm_list'],
        )
        
    filename = 'cec.in'
    sld_varlist = input_dict['sld_cec_list'] 
    make_inputs.get_input_sld_properties(
        outdir=input_dict['outdir'],
        runname=input_dict['runname'],
        filename = filename,
        sld_varlist=sld_varlist,
        )
        
    filename = 'OM_rain.in'
    sld_varlist = input_dict['sld_omrain_list']
    make_inputs.get_input_sld_properties(
        outdir=input_dict['outdir'],
        runname=input_dict['runname'],
        filename = filename,
        sld_varlist=sld_varlist,
        )
        
    filename = '2ndslds.in'
    srcfile = input_dict['srcfile_2ndslds']
    make_inputs.get_input_sld_properties(
        outdir=input_dict['outdir'],
        runname=input_dict['runname'],
        filename = filename,
        srcfile = srcfile,
        )
        
    filename = 'psdrain.in'
    srcfile = input_dict['srcfile_psdrain']
    make_inputs.get_input_sld_properties(
        outdir=input_dict['outdir'],
        runname=input_dict['runname'],
        filename = filename,
        srcfile = srcfile,
        )
        
    filename = 'kinspc.in'
    sld_varlist = input_dict['sld_kinspc_list']
    make_inputs.get_input_sld_properties(
        outdir=input_dict['outdir'],
        runname=input_dict['runname'],
        filename = filename,
        sld_varlist=sld_varlist,
        )
    # changing thermodyniamcs
    # filename = 'keqspc.in'
    # sld_varlist = input_dict['sld_keqspc_list']
    # make_inputs.get_input_sld_properties(
        # outdir=input_dict['outdir'],
        # runname=input_dict['runname'],
        # filename = filename,
        # sld_varlist=sld_varlist,
        # )
        
    filename = 'nopsd.in'
    sld_varlist = input_dict['sld_nopsd_list']
    make_inputs.get_input_sld_properties(
        outdir=input_dict['outdir'],
        runname=input_dict['runname'],
        filename = filename,
        sld_varlist=sld_varlist,
        )


def run(input_list):
    
    outdir,runname,src_code =  input_list
    # >>>> run 
    
    run_success = False
    
    # outdir = input_dict['outdir']
    # runname = input_dict['runname']
    where = '/'
    exename = 'scepter'
    
    # update_code = input_dict['update_code']
    # src_code = input_dict['src_code']
    
    # sub_as_a_job = input_dict['sub_as_a_job']
    # lim_calc_time = input_dict['lim_calc_time']
    # max_calc_time = input_dict['max_calc_time']
    
    update_code = True
    # src_code = 'scepter_DEV'
    # src_code = 'scepter_AMD'
    
    # sub_as_a_job = True
    sub_as_a_job = False
    # lim_calc_time = True
    lim_calc_time = False
    max_calc_time = 40
    
    run_iterative = True
    # run_iterative = False
    
    if update_code:
        print( 'cp {} {}'.format(src_code, outdir+runname+where+exename) )
        os.system( 'cp {} {}'.format(src_code, outdir+runname+where+exename) )
    
    print('chmod u+x '+outdir+runname+where+exename)
    os.system('chmod u+x '+outdir+runname+where+exename)
    
    if not sub_as_a_job: 
        if not lim_calc_time:
            if not run_iterative:
                os.system(outdir+runname+where+exename + ' > ' + outdir+runname+'/logfile.txt' + ' 2> ' + outdir+runname+'/err.txt')
                run_success = True
            else:
                os.system(outdir+runname+where+exename )
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
        
        if run_success:
            if not run_iterative:
                os.remove(outdir+runname+'/logfile.txt')
                os.remove(outdir+runname+'/err.txt')
    else:
        if os.path.exists(outdir+runname+'/run_complete.txt'): os.remove(outdir+runname+'/run_complete.txt')
        slurm_cmd = 'sbatch  --time=0-24:00 --account=gts-creinhard3 --nodes=1 --ntasks=1  -qinferno  --mem-per-cpu=4G run_a_shell.sbatch '
        print(slurm_cmd+outdir+runname+where+exename + ' > ' + outdir+runname+'/joblog.txt' + ' 2> ' + outdir+runname+'/joberr.txt') # submit a job
        os.system(slurm_cmd+outdir+runname+where+exename + ' > ' + outdir+runname+'/joblog.txt' + ' 2> ' + outdir+runname+'/joberr.txt') # submit a job
        
        # get job id 
        while(True):
            try:
                f = open( outdir+runname+'/joblog.txt', mode = 'r' )
                jobidstr = re.sub( "[^0-9]", "",  f.read() )
                f.close()
                break
            except:
                pass
        
        st = time.time() # time before loop
        if lim_calc_time:
            cnt = 0        
            while(cnt <= max_calc_time):
                time.sleep(60)
                if os.path.exists(outdir+runname+'/run_complete.txt'):
                    run_success = True
                    break
                cnt += 1
        else:
            while(True):
                time.sleep(60)
                if os.path.exists(outdir+runname+'/run_complete.txt'):
                    run_success = True
                    break
        et = time.time() 

        if os.path.exists(outdir+runname+'/run_complete.txt'):
            run_success = True
        
        if run_success:
            print('run finished in {:f} min'.format(( et - st )/60.))
            print('slurm-{}.out'.format(jobidstr))
            os.remove('slurm-{}.out'.format(jobidstr))
            # kill job just in case 
            print('scancel ' + jobidstr )
            os.system('scancel ' + jobidstr )
        else:
            print('run UNfinished within {:f} min'.format(int( max_calc_time )))
        
            # kill job
            print('scancel ' + jobidstr )
            os.system('scancel ' + jobidstr )
    
    return run_success
    
    

def save_a_list(dst,res_all):

    with open(dst, 'w') as file:
        for j in range(len(res_all)):
            if j ==0:
                for item in res_all[j]:
                    if res_all[j].index(item)==0:
                        file.write('{:5}\t'.format(item))
                    elif res_all[j].index(item)==len(res_all[j])-1:
                        file.write('{:17}\n'.format(item))
                    else:
                        file.write('{:17}\t'.format(item))
            else:
                item_list = res_all[j]
                for i in range(len(item_list)):
                    item = item_list[i]
                    if i==0:
                        file.write('{:5d}\t'.format(item))
                    elif i==len(item_list)-1:
                        file.write('{:17.9e}\n'.format(item))
                    else:
                        file.write('{:17.9e}\t'.format(item))

def save_a_list_2(dst,res_all):

    with open(dst, 'w') as file:
        for j in range(len(res_all)):
            if j ==0:
                for item in res_all[j]:
                    if res_all[j].index(item)==len(res_all[j])-1:
                        file.write('{:17}\n'.format(item))
                    else:
                        file.write('{:17}\t'.format(item))
            else:
                item_list = res_all[j]
                for i in range(len(item_list)):
                    item = item_list[i]
                    if i==len(item_list)-1:
                        file.write('{:17.9e}\n'.format(item))
                    else:
                        file.write('{:17.9e}\t'.format(item))



def main():
    
    runtype = 'init'
    runtype = 'rstrt'
    
    nz = 101
    
    act_ON = True
    
    scheme_act = 'bdot'
    
    ver = '_v_1'
    
    outdir = '/storage/project/r-creinhard3-0/ykanzaki3/scepter_output/'
    
    runname = f'AMD/test_{runtype}_{nz:d}{ver}'
    
    
    # ---- getting default inputs 
    input_dict = get_input(outdir,runtype,nz,scheme_act,ver)
    
    
    # ---- modify input before starting run 
    input_dict['outdir'] = outdir
    input_dict['runname'] = runname
     
    # ---- making input files 
    setup_run(input_dict)
    
    input_dict['update_code'] = True
    
    input_dict['src_code'] = 'scepter_AMD'
    
    src_code = 'scepter_AMD'
    
    # input_dict['sub_as_a_job'] = True
    input_dict['sub_as_a_job'] = False
    input_dict['lim_calc_time'] = True
    # input_dict['lim_calc_time'] = False
    input_dict['max_calc_time'] = 40
    
    # ---- run the code for a field 
    # run_success = run(input_dict)
    run_success = run([outdir,runname,src_code])
    
    print('!!! run SUCCESS !!!' if run_success else '!!! run FAILED !!!' )
    
    
if __name__ == '__main__':
    main()
