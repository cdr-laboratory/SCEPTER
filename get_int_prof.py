import numpy as np
from scipy import interpolate
import scipy.integrate as integrate
import copy

def phint(dep,phdep,ztot):
    nx = 100
    if ztot>dep[-1]: ztot=dep[-1]
    func = interpolate.interp1d(dep,10.0**-phdep)
    z = np.linspace(dep,ztot,nx)
    result = integrate.quad(lambda x: func(x), dep[0], ztot)
    a = -np.log10( result[0]/(ztot-dep[0]))
    
    return a

def linave(dep,phdep,ztot):
    nx = 100
    if ztot>dep[-1]: ztot=dep[-1]
    func = interpolate.interp1d(dep,phdep)
    z = np.linspace(dep,ztot,nx)
    result = integrate.quad(lambda x: func(x), dep[0], ztot)
    a = ( result[0]/(ztot-dep[0]))
    
    return a

def get_ph_int_site(outdir,runname,dep_sample):
    # outdir = '../pyweath_output/'
    # runname = ''
    # runname = 'Sheldon_A_fick_noiter_test_psdfullpbe_w1_fit_r2000_sig0p2_loop_q0p2_chkall'

    # dep_sample = float(sys.argv[1])

    try:
        infile = outdir+runname+'/prof/prof_aq-020.txt'
        data = np.loadtxt(infile,skiprows=1)
        print('using aq-020')
    except:
        infile = outdir+runname+'/prof/prof_aq-save.txt'
        data = np.loadtxt(infile,skiprows=1)
        print('using aq-save')
        
    pH_dep = data[:,-2]
    dep = data[:,0]

    # print(dep)
    # print(pH_dep)
            
    phintval =  phint(dep,pH_dep, dep_sample)

    print(phintval)
    
    return phintval

def get_ac_int_site(outdir,runname,dep_sample):
    # outdir = '../pyweath_output/'
    # runname = ''
    # runname = 'Sheldon_A_fick_noiter_test_psdfullpbe_w1_fit_r2000_sig0p2_loop_q0p2_chkall'

    # dep_sample = float(sys.argv[1])

    infile = outdir+runname+'/prof/prof_aq(ads%cec)-020.txt'
    data = np.loadtxt(infile,skiprows=1)
    print('using prof_aq(ads%cec)-020')
        
    pH_dep = data[:,-2]
    dep = data[:,0]

    # print(dep)
    # print(pH_dep)
            
    phintval =  linave(dep,pH_dep, dep_sample)

    print(phintval)
    
    return phintval

def get_ac_int_site_v2(outdir,runname,dep_sample):

    infile  = outdir+runname+'/prof/prof_aq(ads%cec)-020.txt'
    data    = np.loadtxt(infile,skiprows=1)
    print('using prof_aq(ads%cec)-020')
    
    with open(infile) as f:
        first_line  = f.readline() 
        sp_list     = first_line.split()
        
    sps     = ['h','al']
    
    pH_dep = 0
    for sp in sps:
        try:
            isp     = sp_list.index(sp) 
        except:
            continue
        
        pH_dep  += data[:,isp]
    
    dep = data[:,0]
            
    phintval    =  linave(dep,pH_dep, dep_sample)

    print(phintval)
    
    return phintval

def get_bs_int_site(outdir,runname,dep_sample):

    infile  = outdir+runname+'/prof/prof_aq(ads%cec)-020.txt'
    data    = np.loadtxt(infile,skiprows=1)
    print('using prof_aq(ads%cec)-020')
    
    with open(infile) as f:
        first_line  = f.readline() 
        sp_list     = first_line.split()
        
    sps     = ['ca','mg','k','na']
    
    pH_dep = 0
    for sp in sps:
        try:
            isp     = sp_list.index(sp) 
        except:
            continue
        
        pH_dep  += data[:,isp]
    
    dep = data[:,0]
            
    phintval =  linave(dep,pH_dep, dep_sample)

    print(phintval)
    
    return phintval

def get_spex_int_site(outdir,runname,dep_sample,sp):

    infile  = outdir+runname+'/prof/prof_aq(ads%cec)-020.txt'
    data    = np.loadtxt(infile,skiprows=1)
    print('using prof_aq(ads%cec)-020')
    
    dep = data[:,0]
    
    with open(infile) as f:
        first_line  = f.readline() 
        sp_list     = first_line.split()
        
    # sps     = ['ca','mg','k','na']
    
    pH_dep = 0
    try:
        isp     = sp_list.index(sp) 
        pH_dep  += data[:,isp]
        phintval =  linave(dep,pH_dep, dep_sample)
    except:
        print('not exist',sp)
        phintval = 0
            

    print(phintval)
    
    return phintval

def get_rhobulk_int_site(outdir,runname,dep_sample):

    infile  = outdir+runname+'/prof/bsd-020.txt'
    data    = np.loadtxt(infile,skiprows=1)
    print('using bsd-020')
    
    dep = data[:,0]
    
    with open(infile) as f:
        first_line  = f.readline() 
        sp_list     = first_line.split()
        
    # sps     = ['ca','mg','k','na']
    
    pH_dep = 0
    try:
        isp     = sp_list.index('dens[g/cm3]') 
        pH_dep  += data[:,isp]
        phintval =  linave(dep,pH_dep, dep_sample)
    except:
        print('not exist dens[g/cm3]')
        phintval = 0
            

    print(phintval)
    
    return phintval

def get_sldwt_int_site(outdir,runname,dep_sample,sps):

    infile  = outdir+runname+'/prof/prof_sld(wt%)-020.txt'
    data    = np.loadtxt(infile,skiprows=1)
    print('using prof_sld(wt%)-020')
    
    dep = data[:,0]
    
    with open(infile) as f:
        first_line  = f.readline() 
        sp_list     = first_line.split()
        
    pH_dep = 0
    
    for sp in sps:
        try:
            isp     = sp_list.index(sp) 
        except:
            continue
        
        pH_dep  += data[:,isp]
            
    phintval =  linave(dep,pH_dep, dep_sample)

    print(phintval)
    
    return phintval

def get_btmwater_site(outdir,runname):

    infile  = outdir+runname+'/prof/prof_aq(tot)-020.txt'
    data    = np.loadtxt(infile,skiprows=1)
    print('using prof_aq(tot)-020')
    
    dep = data[-1,0]
    
    with open(infile) as f:
        first_line  = f.readline() 
        sp_list     = first_line.split()
        
    sps = copy.copy(sp_list)
    
    del sps[0]
    del sps[-1]
    del sps[-1]
    
    btmconcs = []
    
    for sp in sps:
        isp     = sp_list.index(sp)
        btmconcs.append(data[-1,isp])
    
    return sps,btmconcs,dep

def get_water_site(outdir,runname,dep_sample):

    infile  = outdir+runname+'/prof/prof_aq(tot)-020.txt'
    data    = np.loadtxt(infile,skiprows=1)
    print('using prof_aq(tot)-020')
    
    deps = data[:,0]
    deps_list = [dep for dep in deps]
    idep = 0
    for dep in deps:
        if dep_sample>=dep:
            idep = deps_list.index(dep)
    dep = deps[idep]
    
    with open(infile) as f:
        first_line  = f.readline() 
        sp_list     = first_line.split()
        
    sps = copy.copy(sp_list)
    
    del sps[0]
    del sps[-1]
    del sps[-1]
    
    btmconcs = []
    
    for sp in sps:
        isp     = sp_list.index(sp)
        btmconcs.append(data[idep,isp])
    
    return sps,btmconcs,dep