import numpy as np
from scipy import interpolate
import scipy.integrate as integrate
import copy
import os
import shutil

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

def get_ph_int_site(outdir,runname,idep,i):

    infile = outdir+runname+'/prof/prof_aq-{:03d}.txt'.format(i)
    data = np.loadtxt(infile,skiprows=1)
    print('using aq-{:03d}'.format(i))
        
    pH_dep = data[idep,-2]
    dep = data[idep,0]
    
    return pH_dep

def get_ac_int_site(outdir,runname,dep_sample):

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

def get_ac_int_site_v2(outdir,runname,idep,i):

    infile  = outdir+runname+'/prof/prof_aq(ads%cec)-{:03d}.txt'.format(i)
    data    = np.loadtxt(infile,skiprows=1)
    print('using prof_aq(ads%cec)-{:03d}'.format(i))
    
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
        
        pH_dep  += data[idep,isp]
    
    dep = data[idep,0]
    
    return pH_dep

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

def get_rhobulk_int_site(outdir,runname,idep,i):

    infile  = outdir+runname+'/prof/bsd-{:03d}.txt'.format(i)
    data    = np.loadtxt(infile,skiprows=1)
    print('using bsd-{:03d}'.format(i))
    
    dep = data[:,0]
    
    with open(infile) as f:
        first_line  = f.readline() 
        sp_list     = first_line.split()
        
    # sps     = ['ca','mg','k','na']
    
    pH_dep = 0
    try:
        isp     = sp_list.index('dens[g/cm3]') 
        pH_dep  += data[pH_dep,isp]
    except:
        print('not exist dens[g/cm3]')
            

    
    return pH_dep

def get_sldwt_int_site(outdir,runname,idep,sps,i):

    infile  = outdir+runname+'/prof/prof_sld(wt%)-{:03d}.txt'.format(i)
    data    = np.loadtxt(infile,skiprows=1)
    print('using prof_sld(wt%)-{:03d}'.format(i))
    
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
        
        pH_dep  += data[idep,isp]
    
    return pH_dep

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

def get_water_site_OLD(outdir,runname,dep_sample):

    for i in range(20,0,-1):
        infile  = outdir+runname+'/prof/prof_aq-{:03d}.txt'.format(i)
        if not os.path.exists(infile): continue
        else: break 
    infile  = outdir+runname+'/prof/prof_aq-{:03d}.txt'.format(i)
    data    = np.loadtxt(infile,skiprows=1)
    print('using prof_aq-{:03d}'.format(i))
    
    shutil.copyfile(infile, infile.replace('.txt','_org.txt'))
    # saving gas data
    infile2  = outdir+runname+'/prof/prof_gas-{:03d}.txt'.format(i)
    shutil.copyfile(infile2, infile2.replace('.txt','_org.txt'))
    
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

def get_DIC_save_get_site_OLD(outdir,runname,dep_sample):

    for i in range(20,0,-1):
        infile  = outdir+runname+'/prof/prof_aq-{:03d}.txt'.format(i)
        if not os.path.exists(infile): continue
        else: break 
    infile  = outdir+runname+'/prof/chrge_balance-{:03d}.txt'.format(i)
    data    = np.loadtxt(infile,skiprows=1)
    print('using chrge_balance-{:03d}'.format(i))
    
    shutil.copyfile(infile, infile.replace('.txt','_org.txt'))
    
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
        
    # sps = copy.copy(sp_list)
    
    # del sps[0]
    # del sps[-1]
    # del sps[-1]
    
    btmconcs = []
    
    co2sp_list = [
        'hco3','co3','mg(co3)','mg(hco3)','na(co3)','na(hco3)','ca(co3)','ca(hco3)','fe2(co3)','fe2(hco3)'
        ]
    
    hco3 = data[:,sp_list.index('hco3')]
    co3 = data[:,sp_list.index('co3')]
    h = data[:,sp_list.index('h')]
    k2 = h*co3/hco3
    Rg = 8.3e-3
    tempk_0 = 273
    k2_ref = 10**(-10.43)
    H_k2 = 17.00089
    temp_ref = 15
    # k2 = k2_ref*exp(-H_k2/Rg*(1/(temp + tempk_0)-1/(temp_ref + tempk_0)))
    # k2/k2_ref = exp(-H_k2/Rg*(1/(temp + tempk_0)-1/(temp_ref + tempk_0)))
    # log(k2/k2_ref) = -H_k2/Rg*(1/(temp + tempk_0)-1/(temp_ref + tempk_0))
    # log(k2/k2_ref)*(-Rg/H_k2) = 1/(temp + tempk_0)-1/(temp_ref + tempk_0)
    # log(k2/k2_ref)*(-Rg/H_k2)+1/(temp_ref + tempk_0) = 1/(temp + tempk_0)
    # 1/[log(k2/k2_ref)*(-Rg/H_k2)+1/(temp_ref + tempk_0)] = temp + tempk_0
    # 1/[log(k2/k2_ref)*(-Rg/H_k2)+1/(temp_ref + tempk_0)] - tempk_0 = temp
    temp =  1./(np.log(k2/k2_ref)*(-Rg/H_k2)+1./(temp_ref + tempk_0)) - tempk_0 
    
    k1_ref = 10**(-6.42)
    H_k1 = 11.94453
    k1 = k1_ref*np.exp(-H_k1/Rg*(1./(temp + tempk_0)-1./(temp_ref + tempk_0)))
    
    co2 = hco3*h/k1
    
    DIC = co2
    
    for co2sp in co2sp_list:
        isp     = sp_list.index(co2sp)
        DIC += data[:,isp]
        
    DIC_btm = DIC[idep]
    
    np.savetxt(outdir+runname+'/prof/prof_aq-DIC-{:03d}_OLD.txt'.format(i),np.transpose(np.array([list(deps),list(DIC)])))
    
    return DIC_btm,dep,i

def get_DIC_save_get_site_NEW(outdir,runname,dep_sample,i_save):

    for i in range(20,0,-1):
        infile  = outdir+runname+'/prof/prof_aq-{:03d}.txt'.format(i)
        if not os.path.exists(infile): continue
        else: break 
    infile  = outdir+runname+'/prof/charge_balance-{:03d}.txt'.format(i)
    data    = np.loadtxt(infile,skiprows=1)
    print('using charge_balance-{:03d}'.format(i))
    
    shutil.copyfile(infile, infile.replace('.txt','_org.txt'))
    
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
        
    # sps = copy.copy(sp_list)
    
    # del sps[0]
    # del sps[-1]
    # del sps[-1]
    
    btmconcs = []
    
    co2sp_list = [
        'hco3','co3','mg(co3)','mg(hco3)','na(co3)','na(hco3)','ca(co3)','ca(hco3)','fe2(co3)','fe2(hco3)'
        ]
    
    hco3 = data[:,sp_list.index('hco3')]
    co3 = data[:,sp_list.index('co3')]
    h = data[:,sp_list.index('h')]
    k2 = h*co3/hco3
    Rg = 8.3e-3
    tempk_0 = 273
    k2_ref = 10**(-10.43)
    H_k2 = 17.00089
    temp_ref = 15
    # k2 = k2_ref*exp(-H_k2/Rg*(1/(temp + tempk_0)-1/(temp_ref + tempk_0)))
    # k2/k2_ref = exp(-H_k2/Rg*(1/(temp + tempk_0)-1/(temp_ref + tempk_0)))
    # log(k2/k2_ref) = -H_k2/Rg*(1/(temp + tempk_0)-1/(temp_ref + tempk_0))
    # log(k2/k2_ref)*(-Rg/H_k2) = 1/(temp + tempk_0)-1/(temp_ref + tempk_0)
    # log(k2/k2_ref)*(-Rg/H_k2)+1/(temp_ref + tempk_0) = 1/(temp + tempk_0)
    # 1/[log(k2/k2_ref)*(-Rg/H_k2)+1/(temp_ref + tempk_0)] = temp + tempk_0
    # 1/[log(k2/k2_ref)*(-Rg/H_k2)+1/(temp_ref + tempk_0)] - tempk_0 = temp
    temp =  1./(np.log(k2/k2_ref)*(-Rg/H_k2)+1./(temp_ref + tempk_0)) - tempk_0 
    
    k1_ref = 10**(-6.42)
    H_k1 = 11.94453
    k1 = k1_ref*np.exp(-H_k1/Rg*(1./(temp + tempk_0)-1./(temp_ref + tempk_0)))
    
    co2 = hco3*h/k1
    
    DIC = co2
    
    for co2sp in co2sp_list:
        isp     = sp_list.index(co2sp)
        DIC += data[:,isp]
        
    DIC_btm = DIC[idep]
    
    np.savetxt(outdir+runname+'/prof/prof_aq-DIC-{:03d}_NEW.txt'.format(i_save),np.transpose(np.array([list(deps),list(DIC)])))
    
    return DIC_btm,dep

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

def get_waterave_site(outdir,runname,dep_sample):

    infile  = outdir+runname+'/prof/prof_aq(tot)-020.txt'
    data    = np.loadtxt(infile,skiprows=1)
    print('using prof_aq(tot)-020')
    
    infile2  = outdir+runname+'/prof/bsd-020.txt'
    data2    = np.loadtxt(infile2,skiprows=1)
    print('using bsd-020')
    
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
        
    with open(infile2) as f:
        first_line  = f.readline() 
        prop_list     = first_line.split()
        
    sat = data2[:,prop_list.index('sat')]
    poro = data2[:,prop_list.index('poro')]
    
    # print(data2.shape)
    # print(prop_list,prop_list.index('sat'),prop_list.index('poro'),sat,poro)
        
    sps = copy.copy(sp_list)
    
    del sps[0]
    del sps[-1]
    del sps[-1]
    
    btmconcs = []
    
    for sp in sps:
        isp     = sp_list.index(sp)
        btmconcs.append( 
            np.average( data[:idep+1,isp]*poro[:idep+1]*sat[:idep+1]*1e3 )  # mol/ soil m3
            /np.average( poro[:idep+1]*sat[:idep+1] )  # mol/ aq m3
            )
    
    return sps,btmconcs,dep

def get_totsave_site(outdir,runname,idep,itime):
    
    infile  = outdir+runname+'/prof/prof_ex(tot)-{:03d}.txt'.format(itime)
    if not os.path.exists(infile): 
        print('data does not exist: return nothing')
        return 
    
    data    = np.loadtxt(infile,skiprows=1)
    print('using prof_ex(tot)-{:03d}'.format(itime))
    
    infile2  = outdir+runname+'/prof/bsd-{:03d}.txt'.format(itime)
    data2    = np.loadtxt(infile2,skiprows=1)
    print('using bsd-{:03d}'.format(itime))
    
    time = data[-1,-1]
    deps = data[:,0]
    dep = deps[idep]
    
    with open(infile) as f:
        first_line  = f.readline() 
        sp_list     = first_line.split()
        
    with open(infile2) as f:
        first_line  = f.readline() 
        prop_list     = first_line.split()
        
    poro = data2[:,prop_list.index('poro')]
    
    # print(data2.shape)
    # print(prop_list,prop_list.index('sat'),prop_list.index('poro'),sat,poro)
        
    sps = copy.copy(sp_list)
    
    del sps[0]
    del sps[-1]
    del sps[-1]
    
    btmconcs = []
    
    for sp in sps:
        isp     = sp_list.index(sp)
        btmconcs.append( 
            ( data[idep,isp] ) # mol/ soil m3
            /( (1 - poro[idep]) )   # mol/ solid m3
            )  
    
    return sps,btmconcs,dep,time

def get_adssave_site(outdir,runname,idep,i):
    
    infile  = outdir+runname+'/prof/prof_aq(ads)-{:03d}.txt'.format(i)
    data    = np.loadtxt(infile,skiprows=1)
    print('using prof_aq(ads)-{:03d}'.format(i))
    
    infile2  = outdir+runname+'/prof/bsd-{:03d}.txt'.format(i)
    data2    = np.loadtxt(infile2,skiprows=1)
    print('using bsd-{:03d}'.format(i))
    
    time = data[-1,-1]
    
    deps = data[:,0]
    dep = deps[idep]
    
    with open(infile) as f:
        first_line  = f.readline() 
        sp_list     = first_line.split()
        
    with open(infile2) as f:
        first_line  = f.readline() 
        prop_list     = first_line.split()
        
    poro = data2[:,prop_list.index('poro')]
    dense = data2[:,prop_list.index('dens[g/cm3]')]
    
    # print(data2.shape)
    # print(prop_list,prop_list.index('sat'),prop_list.index('poro'),sat,poro)
        
    sps = copy.copy(sp_list)
    
    del sps[0]
    del sps[-1]
    del sps[-1]
    
    btmconcs = []
    
    for sp in sps:
        isp     = sp_list.index(sp)
        btmconcs.append( 
            ( 
                data[idep,isp]  # cmol/kg(solid)
                * 1e-2/1e3  # now mok/g(solid)
                * dense[idep]    # mol/ solid cm3
                * 1./1e-6
                )  
            )  
    
    return sps,btmconcs,dep,time

def get_dis_frac(outdir,runname,sp):
    
    infile  = outdir+runname+'/flx/int_flx_sld-'+sp+'.txt'
    data    = np.loadtxt(infile,skiprows=1)
    
    with open(infile) as f:
        first_line  = f.readline() 
        flx_list     = first_line.split()
        
    rain = data[-1,flx_list.index('rain')]
    time = data[-1,flx_list.index('time')]
    dis = data[-1,flx_list.index(sp)]
    
    frac = abs(dis/rain)*100.
    
    return frac,time