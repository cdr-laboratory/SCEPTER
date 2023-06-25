import os
import numpy as np
from shutil import copyfile
import time

def get_input_frame(outdir,runname):
    infile = outdir + runname + '/frame.in'
    
    with open(infile) as f:
        lines = [line for line in f]
    
    ztot = float(lines[1].split()[0])
    nz = int(lines[2].split()[0])
    ttot = float(lines[3].split()[0])
    temp = float(lines[4].split()[0])
    fdust = float(lines[5].split()[0])
    fdust2 = float(lines[6].split()[0])
    taudust = float(lines[7].split()[0])
    omrain = float(lines[8].split()[0])
    zom = float(lines[9].split()[0])
    poro = float(lines[10].split()[0])
    moistsrf = float(lines[11].split()[0])
    zwater = float(lines[12].split()[0])
    zdust = float(lines[13].split()[0])
    w = float(lines[14].split()[0])
    q = float(lines[15].split()[0])
    p = float(lines[16].split()[0])
    nstep = int(lines[17].split()[0])
    rstrt = (lines[18].split()[0])
    runid = (lines[20].split()[0])
    
    return ztot,nz,ttot,temp,fdust,fdust2,taudust,omrain,zom,poro,moistsrf,zwater,zdust,w,q,p,nstep,rstrt,runid

def get_input_switches(outdir,runname):    
    infile = outdir + runname + '/switches.in'
    
    with open(infile) as f:
        lines = [line for line in f]
    
    w_scheme = int(lines[1].split()[0])
    mix_scheme = int(lines[2].split()[0])
    poro_iter = (lines[3].split()[0])
    sldmin_lim = (lines[4].split()[0])
    display = (lines[5].split()[0])
    disp_lim = (lines[6].split()[0])
    restart = (lines[7].split()[0])
    rough = (lines[8].split()[0])
    act_ON = (lines[9].split()[0])
    dt_fix = (lines[10].split()[0])
    cec_on = (lines[11].split()[0])
    dz_fix = (lines[12].split()[0])
    close_aq = (lines[13].split()[0])
    poro_evol = (lines[14].split()[0])
    sa_evol_1 = (lines[15].split()[0])
    sa_evol_2 = (lines[16].split()[0])
    psd_bulk = (lines[17].split()[0])
    psd_full = (lines[18].split()[0])
    season = (lines[19].split()[0])
    
    return w_scheme,mix_scheme,poro_iter,sldmin_lim,display,disp_lim,restart,rough,act_ON,dt_fix \
        ,cec_on,dz_fix,close_aq,poro_evol,sa_evol_1,sa_evol_2,psd_bulk,psd_full,season

def get_input_tracers(outdir,runname):
    
    sld_list = np.genfromtxt(outdir + runname + '/slds.in',dtype='str',skip_header=1)
    aq_list = np.genfromtxt(outdir + runname + '/solutes.in',dtype='str',skip_header=1)
    gas_list = np.genfromtxt(outdir + runname + '/gases.in',dtype='str',skip_header=1)
    exrxn_list = np.genfromtxt(outdir + runname + '/extrxns.in',dtype='str',skip_header=1)
    
    try:
        sld_list = list(sld_list)
    except:
        sld_list = [str(sld_list)]
    try:
        aq_list = list(aq_list)
    except:
        aq_list = [aq_list[0]]
    try:
        gas_list = list(gas_list)
    except:
        print(gas_list.shape)
        gas_list = [str(gas_list)]
    try:
        exrxn_list = list(exrxn_list)
    except:
        exrxn_list = [exrxn_list[0]]
    
    return sld_list,aq_list,gas_list,exrxn_list

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

def get_input_sld_properties(outdir,runname,filename):

    
    infile = outdir + runname + '/' + filename
    
    sld_data_list = []
    
    
    with open(infile) as f:
        cnt2 = 0
        for line in f:
            if cnt2 ==0: 
                cnt2 += 1
                continue 
            else:
                cnt2 += 1
            a = line.split()
            tmplist = []
            cnt = 0
            for item in a:
                if cnt ==0:  tmplist.append(item)
                else: tmplist.append(float(item))
                cnt += 1
            
            sld_data_list.append(tmplist)
    
    return sld_data_list
    
def main():

    outdir = '../scepter_output/'
    # runname = 'US_cropland_311_sph_N_spintuneup_field'
    runname = 'test_Pot7'
    
    ztot,nz,ttot,temp,fdust,fdust2,taudust,omrain,zom,poro,moistsrf,zwater,zdust,w,q,p,nstep,rstrt,runid = get_input_frame(outdir,runname)
    print(
        ztot,nz,ttot,temp,fdust,fdust2,taudust,omrain,zom,poro,moistsrf,zwater,zdust,w,q,p,nstep,rstrt,runid
        )
        
    w_scheme,mix_scheme,poro_iter,sldmin_lim,display,disp_lim,restart,rough,act_ON,dt_fix \
        ,cec_on,dz_fix,close_aq,poro_evol,sa_evol_1,sa_evol_2,psd_bulk,psd_full,season = get_input_switches(outdir,runname)
    print(
        w_scheme,mix_scheme,poro_iter,sldmin_lim,display,disp_lim,restart,rough,act_ON,dt_fix 
        ,cec_on,dz_fix,close_aq,poro_evol,sa_evol_1,sa_evol_2,psd_bulk,psd_full,season 
        )
        
    sld_list,aq_list,gas_list,exrxn_list = get_input_tracers(outdir,runname)
    print(
        sld_list,aq_list,gas_list,exrxn_list
        )
    
    filename = 'cec.in'
    sld_data_list = get_input_sld_properties(outdir,runname,filename)
    print(
        sld_data_list,
        sld_data_list[0][1]
        )
    
    sld_data_list = get_input_sld_properties('./','data','dust_gbasalt.in')
    print(
        sld_data_list
        )
   
if __name__ == '__main__':
    main()
    
    