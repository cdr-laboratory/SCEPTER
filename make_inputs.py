import os
import numpy as np
from shutil import copyfile
import time

def get_input_frame(**kwargs):
    outdir      = kwargs.get('outdir',  '../scepter_output/')
    runname     = kwargs.get('runname', 'test_input')
    ztot        = kwargs.get('ztot',    0.5)
    nz          = kwargs.get('nz',      30) 
    ttot        = kwargs.get('ttot',    1) 
    temp        = kwargs.get('temp',    15) 
    fdust       = kwargs.get('fdust',   0)
    fdust2      = kwargs.get('fdust2',  0)
    taudust     = kwargs.get('taudust', 0) 
    omrain      = kwargs.get('omrain',  0)      
    zom         = kwargs.get('zom',     ztot) 
    poro        = kwargs.get('poro',    0.5)
    moistsrf    = kwargs.get('moistsrf',0.1)
    zwater      = kwargs.get('zwater',  1000)
    zdust       = kwargs.get('zdust',   0.25)
    w           = kwargs.get('w',       1e-5)
    q           = kwargs.get('q',       0.01) 
    p           = kwargs.get('p',       10e-6)
    nstep       = kwargs.get('nstep',   10)
    rstrt       = kwargs.get('rstrt',   'self')
    runid       = kwargs.get('runid',   runname) 
    
    
    notes = [
        '** values to determine the boundary conditions'
        ,'total depth of weathering profile [m]'
        ,'number of grids into which calculation domain is divided'
        ,'total duration of simulation [yr]'
        ,'temperature [oC]'
        ,'amounts of dusts [g/m2/yr]'
        ,'amounts of 2nd dusts [g/m2/yr]'
        ,'duration of dust application [yr]'
        ,'OM [g C/m2/yr]'
        ,'depth of mixed layer for OM  [m]'
        ,'initial porosity'
        ,'water saturation at the surface of profile'
        ,'depth of water table [m]'
        ,'depth of mixed layer for dust [m]'
        ,'uplift rate [m/yr]'
        ,'net water flux [m/yr]'
        ,'radius of particles [m]'
        ,'interations needed to increase time step by a factor of 10'
        ,''
        ,'^ directory name for restart (switch is in switches.in) (type "self" when restart from the same directory)'
        ,''
        ,'^ simulation name'
        ]
    
    values = [
        ''
        ,ztot
        ,nz
        ,ttot
        ,temp
        ,fdust
        ,fdust2
        ,taudust
        ,omrain
        ,zom
        ,poro
        ,moistsrf
        ,zwater
        ,zdust
        ,w
        ,q
        ,p
        ,nstep
        ,rstrt
        ,''
        ,runid
        ,''
        ]
        
    if len(values) != len(notes): 
        print('error')
    
    n = len(values)
        
    input_text = ''
    for i in range(n):
        if values[i]=='':
            input_text += notes[i] + '\n'
        else:
            input_text += str(values[i]) + '\t' + notes[i] + '\n'
    
    if not os.path.exists(outdir + runname): os.makedirs(outdir + runname)
    
    input_file = outdir + runname + '/frame.in'
    
    with open(input_file, 'w') as file:
        file.write(input_text)
    
    print(input_text)

def get_input_switches(**kwargs):
    outdir      = kwargs.get('outdir',      '../scepter_output/')
    runname     = kwargs.get('runname',     'test_input')
    w_scheme    = kwargs.get('w_scheme',    1)
    mix_scheme  = kwargs.get('mix_scheme',  1) 
    poro_iter   = kwargs.get('poro_iter',   'false') 
    sldmin_lim  = kwargs.get('sldmin_lim',  'true') 
    display     = kwargs.get('display',     'true')
    disp_lim    = kwargs.get('disp_lim',    'true')
    restart     = kwargs.get('restart',     'false') 
    rough       = kwargs.get('rough',       'true')      
    act_ON      = kwargs.get('act_ON',      'false') 
    dt_fix      = kwargs.get('dt_fix',      'false')
    cec_on      = kwargs.get('cec_on',      'false')
    dz_fix      = kwargs.get('dz_fix',      'true')
    close_aq    = kwargs.get('close_aq',    'false')
    poro_evol   = kwargs.get('poro_evol',   'true')
    sa_evol_1   = kwargs.get('sa_evol_1',   'true') 
    sa_evol_2   = kwargs.get('sa_evol_2',   'false')
    psd_bulk    = kwargs.get('psd_bulk',    'true')
    psd_full    = kwargs.get('psd_full',    'true')
    season      = kwargs.get('season',      'false') 
    
    notes = [
        '** switch number or on/off [true if on, false if off]'
        ,'erosion scheme: 0-- cnst w, 1-- cnst poro*w, 2-- cnst (1-poro)*w, 3--- w-flexible(cnst porosity prof), if not defined 0 is taken'
        ,'bio-mixing style: 0-- no mixing, 1-- fickian mixing, 2-- homogeneous mixng, 3--- tilling, 4--- LABS mixing, if not defined 0 is taken'
        ,'porosity  iteration'
        ,'limiting mineral lowest conc.'
        ,'display results at runtime'
        ,'limited results display'
        ,'restart from a previous run'
        ,'include roughness in mineral surface area'
        ,'enabling activity coefficients'
        ,'time step fixed'
        ,'enabling adsorption for cation exchange'
        ,'adopting a regular grid'
        ,'closing system for aq phases'
        ,'enabling porosity evolution'
        ,'enabling SA evolution 1 (SA decreases as porosity increases)'
        ,'enabling SA evolution 2 (SA increases with porosity)'
        ,'enabling PSD tracking'
        ,'enabling PSD tracking for individual solid species'
        ,'enabling full seasonality'
        ]
    
    values = [
        ''
        ,w_scheme
        ,mix_scheme 
        ,poro_iter 
        ,sldmin_lim 
        ,display
        ,disp_lim
        ,restart 
        ,rough      
        ,act_ON 
        ,dt_fix
        ,cec_on
        ,dz_fix
        ,close_aq
        ,poro_evol
        ,sa_evol_1 
        ,sa_evol_2
        ,psd_bulk
        ,psd_full
        ,season
        ]
        
    if len(values) != len(notes): 
        print('error')
    
    n = len(values)
        
    input_text = ''
    for i in range(n):
        if values[i]=='':
            input_text += notes[i] + '\n'
        else:
            input_text += str(values[i]) + '\t' + notes[i] + '\n'
    
    if not os.path.exists(outdir + runname): os.makedirs(outdir + runname)
    
    input_file = outdir + runname + '/switches.in'
    
    with open(input_file, 'w') as file:
        file.write(input_text)
    
    print(input_text)

def get_input_tracers(**kwargs):
    outdir      = kwargs.get('outdir',      '../scepter_output/')
    runname     = kwargs.get('runname',     'test_input')
    sld_list    = kwargs.get('sld_list',    [])
    aq_list     = kwargs.get('aq_list',     []) 
    gas_list    = kwargs.get('gas_list',    []) 
    exrxn_list  = kwargs.get('exrxn_list',  []) 
    
    notes = [
        '** choose and list from [fo, ab, an, cc, ka, gb, py, ct, fa, gt, cabd, dp, hb, kfs ... ]'
        ,'** choose and list from [mg, si, na, ca, al, fe2, fe3, k, so4, no3 ... ]'
        ,'** choose and list from [pco2, po2, pnh3, pn2o]'
        ,'** choose and list from [resp, fe2o2, omomb, ombto, pyfe3, amo2o, g2n0, g2n21, g2n22]'
        ]
    
    sld_list.insert(0,'')
    aq_list.insert(0,'')
    gas_list.insert(0,'')
    exrxn_list.insert(0,'')
        
    trc_lists = [
        sld_list
        ,aq_list
        ,gas_list
        ,exrxn_list
        ]
        
    filenames = [
        'slds.in'
        ,'solutes.in'
        ,'gases.in'
        ,'extrxns.in'
        ]
    
    nn = 4
    
    if not os.path.exists(outdir + runname): os.makedirs(outdir + runname)
    
    for k in range(nn):
        values = trc_lists[k]
        filename = filenames[k]
        n = len(values)
        input_text = ''
        for i in range(n):
            if i==0:
                input_text += notes[k] + '\n'
            else:
                input_text += values[i] + '\n'
        
        input_file = outdir + runname + '/' + filename
        
        with open(input_file, 'w') as file:
            file.write(input_text)
        
        print(input_text)
        
    del sld_list[0]
    del aq_list[0]
    del gas_list[0]
    del exrxn_list[0]

def get_input_tracer_bounds(**kwargs):
    outdir      = kwargs.get('outdir',      '../scepter_output/')
    runname     = kwargs.get('runname',     'test_input')
    pr_list     = kwargs.get('pr_list',     [])
    rain_list   = kwargs.get('rain_list',   []) 
    atm_list    = kwargs.get('atm_list',    [('pco2',3.16e-4),('po2',0.21),('pnh3',1e-50),('pn2o',1e-50)]) 
    
    notes = [
        '** parent rock wt fraction (e.g., "ab      0.2" in one line and "ka     0.001" in the next) (if not specified assumed 1e-20)'
        ,'** solute concs. of rain in mol/L (if not specified assumed 1e-20)'
        ,'** atmospheric composition in atm (if not specified assumed 1 PAL)'
        ]
    
    pr_list.insert(0,'')
    rain_list.insert(0,'')
    atm_list.insert(0,'')
        
    values_lists = [
        pr_list
        ,rain_list
        ,atm_list
        ]
        
    filenames = [
        'parentrock.in'
        ,'rain.in'
        ,'atm.in'
        ]
    
    nn = 3
    
    if not os.path.exists(outdir + runname): os.makedirs(outdir + runname)
    
    for k in range(nn):
        values = values_lists[k]
        filename = filenames[k]
        n = len(values)
        input_text = ''
        for i in range(n):
            if i==0:
                input_text += notes[k] + '\n'
            else:
                try:
                    input_text += values[i][0] + '\t' + str(values[i][1]) + '\n'
                except:
                    print(n,i,values)
                    time.sleep(100)
        
        input_file = outdir + runname + '/' + filename
        
        with open(input_file, 'w') as file:
            file.write(input_text)
        
        print(input_text)
    
    
    del pr_list[0]
    del rain_list[0]
    del atm_list[0]

def get_input_sld_properties(**kwargs):
    outdir      = kwargs.get('outdir',      '../scepter_output/')
    runname     = kwargs.get('runname',     'test_input')
    sld_varlist = kwargs.get('sld_varlist', [])
    filename    = kwargs.get('filename',    'kinspc.in')
    srcfile     = kwargs.get('srcfile',     None)
    
    if filename == 'kinspc.in':
        note = '** specify rate const in [mol/m2/yr] except for OMs which should be presented as turnover year [yr] (e.g., g2   1.0)'
    elif filename == 'sa.in':
        note = '** parent rock particle radii in meter (e.g., "ab      1e-5") (if not specified value in frame.in is used for all sld sp.)'
    elif filename == 'OM_rain.in':
        note = '** OM rain fraction wrt the value in frame.in (if not specified assumed 0)'
    elif filename == 'dust.in' or filename == 'dust_2nd.in':
        note = '** dusts wt fraction (if not specified assumed 0)'
    elif filename == 'cec.in':
        note = '** cec [cmol/kg], log10(KH-X) [-] (X=Na,K,Ca,Mg,Al) and beta specified by users (e.g., "g2   90   5.9   4.8   10.47   10.786   16.47   3.4") (if not specified assumed code default values)'
    elif filename == 'nopsd.in':
        note = '** list of minerals whose PSDs are not tracked for some reason'
    elif filename == '2ndslds.in':
        note = '** list of minerals whose precipitation is allowed'
    elif filename == 'psdrain.in':
        note = '** mean radius [m], standard deviation in log10 [-], weight [-], gaussian parameters to define dust psd (e.g., 1e-5    0.2    1)'
    else:
        print('{} is not supposed to be input file'.format(filename))
    
    sld_varlist.insert(0,'')
        
    
    if not os.path.exists(outdir + runname): os.makedirs(outdir + runname)
    
    if srcfile != None:
        copyfile(srcfile, outdir + runname + '/' + filename)
    else:
        
        n = len(sld_varlist)
        input_text = ''
        for i in range(n):
            if i==0:
                input_text += note + '\n'
            else:
                if filename == 'nopsd.in' or filename == '2ndslds.in':
                    input_text += sld_varlist[i][0] + '\n'
                elif filename == 'cec.in':
                    for j in range(len(sld_varlist[i])):
                        if j==0:
                            input_text += sld_varlist[i][j] + '\t' 
                        elif j==len(sld_varlist[i])-1:
                            input_text += str(sld_varlist[i][j]) + '\n' 
                        else:
                            input_text += str(sld_varlist[i][j]) + '\t'
                elif filename == 'psdrain.in':
                    for j in range(len(sld_varlist[i])):
                        if j==len(sld_varlist[i])-1:
                            input_text += str(sld_varlist[i][j]) + '\n' 
                        else:
                            input_text += str(sld_varlist[i][j]) + '\t'
                else:
                    input_text += sld_varlist[i][0] + '\t' + str(sld_varlist[i][1]) + '\n'
        
        input_file = outdir + runname + '/' + filename
        
        with open(input_file, 'w') as file:
            file.write(input_text)
        
        print(input_text)
    
    
    del sld_varlist[0]
    
    
    
def get_input_climate_temp(**kwargs):
    outdir      = kwargs.get('outdir',      '/storage/scratch1/0/ykanzaki3/scepter_output/')
    runname     = kwargs.get('runname',     'test_input')
    T_ave       = kwargs.get('T_ave',       15)
    T_amp       = kwargs.get('T_amp',       0.3)
    moist_ave   = kwargs.get('moist_ave',   0.5)
    moist_amp   = kwargs.get('moist_amp',   0.3)
    q_ave       = kwargs.get('q_ave',       0.5)
    q_amp       = kwargs.get('q_amp',       0.3)
    tau         = kwargs.get('tau',         1.)
    timeline    = kwargs.get('timeline',    np.linspace(0,1,12,endpoint=False))

    N = timeline.shape[0]
    T = T_ave + T_ave * T_amp * np.sin(timeline/tau*2.*np.pi)
    q = q_ave + q_ave * q_amp * np.sin(timeline/tau*2.*np.pi)
    moist = moist_ave + moist_ave * moist_amp * np.sin(timeline/tau*2.*np.pi)
    

    notes = [
        '# time(yr) / runoff(mm/month)',
        '# time(yr) / T(oC)',
        '# time(yr) / moisture(mm/m)',
        ]
        
    filenames = [
        'q_temp.in',
        'T_temp.in',
        'Wet_temp.in',
        ]
        
    values_lists = [
        q,
        T,
        moist,
        ]
        
    factors = [
        1e3/12.,  # converting m/y to mm/month
        1,
        1e3,  # converting m/m to mm/m
        ]
    
    nn = 3
    
    if not os.path.exists(outdir + runname): os.makedirs(outdir + runname)
    
    for k in range(nn):
        values = values_lists[k]
        filename = filenames[k]
        factor = factors[k]
        input_text = ''
        input_text += '{}\n'.format(notes[k])
        for i in range(N):
            input_text += '{:f}\t{:f}\n'.format(timeline[i], values[i]*factor) 
        
        input_file = outdir + runname + '/' + filename
        
        with open(input_file, 'w') as file:
            file.write(input_text)
        
        print(input_text)
        
        
        
        
def main():

    outdir = '../scepter_output/'
    runname = 'test_input'
    
    ztot=0.5
    nz=30
    ttot=1e5
    temp=15
    fdust=0
    fdust2=0
    taudust=0
    omrain=100
    zom=ztot
    poro=0.5
    moistsrf=0.1
    zwater=1000
    zdust=ztot/2.
    w=1e-5
    q=0.01
    p=1e-5
    nstep=10
    rstrt='self'
    runid=runname

    get_input_frame(
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
    
    w_scheme=1
    mix_scheme=1 
    poro_iter='false' 
    sldmin_lim ='true'
    display='true'
    disp_lim='true'
    restart ='false'
    rough      ='true'
    act_ON ='false'
    dt_fix='false'
    precalc='false'
    dz_fix='true'
    sld_fix='false'
    poro_evol='true'
    sa_evol_1 ='true'
    sa_evol_2='false'
    psd_bulk='true'
    psd_full='true'
    season='false'
        
    get_input_switches(
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
        ,precalc=precalc
        ,dz_fix=dz_fix
        ,sld_fix=sld_fix
        ,poro_evol=poro_evol
        ,sa_evol_1=sa_evol_1 
        ,sa_evol_2=sa_evol_2
        ,psd_bulk=psd_bulk
        ,psd_full=psd_full
        ,season=season
        )
    
    sld_list=['inrt','g2']
    aq_list = ['ca','k','mg','na']
    gas_list = ['pco2']
    exrxn_list = []
    get_input_tracers(
        outdir=outdir
        ,runname=runname
        ,sld_list = sld_list
        ,aq_list = aq_list
        ,gas_list = gas_list
        ,exrxn_list = exrxn_list
        )
    pr_list = [('inrt',1.0)]
    rain_list = [('ca',5.0e-6)]
    atm_list = [('pco2',3.16e-4),('po2',0.21),('pnh3',1e-50),('pn2o',1e-50)]
    get_input_tracer_bounds(
        outdir=outdir
        ,runname=runname
        ,pr_list = pr_list
        ,rain_list=rain_list
        ,atm_list=atm_list
        )
    filename = 'dust.in'
    srcfile = './data/dust_gbasalt.in'
    sld_varlist =[]
    get_input_sld_properties(
        outdir=outdir
        ,runname=runname
        ,filename = filename
        # ,srcfile = srcfile
        ,sld_varlist=sld_varlist
        )
    filename = 'cec.in'
    sld_varlist = [('inrt',4,  5.9, 4.8, 10.47, 10.786, 16.47,  3.4) ,('g2',4,  5.9, 4.8, 10.47, 10.786, 16.47,  3.4) ] 
    get_input_sld_properties(
        outdir=outdir
        ,runname=runname
        ,filename = filename
        ,sld_varlist=sld_varlist
        )
    filename = 'OM_rain.in'
    sld_varlist = [ ('g2',1) ] 
    get_input_sld_properties(
        outdir=outdir
        ,runname=runname
        ,filename = filename
        ,sld_varlist=sld_varlist
        )
    filename = '2ndslds.in'
    srcfile = './data/2ndslds_def.in'
    srcfile = './data/2ndslds_rm_al2o3.in'
    sld_varlist =[]
    get_input_sld_properties(
        outdir=outdir
        ,runname=runname
        ,filename = filename
        ,srcfile = srcfile
        # ,sld_varlist=sld_varlist
        )
    filename = 'psdrain.in'
    srcfile = './data/psdrain_100um.in'
    sld_varlist = [ (5e-6,0.2,1), (20e-6,0.2,1), (50e-6,0.2,1), (70e-6,0.2,1) ] 
    get_input_sld_properties(
        outdir=outdir
        ,runname=runname
        ,filename = filename
        ,sld_varlist=sld_varlist
        # ,srcfile = srcfile
        )
    
    get_input_climate_temp(
        outdir=outdir,
        runname=runname,
        tau = 1,
        T_amp = 0,
        moist_amp = 0,
        )
   
if __name__ == '__main__':
    main()
    
    