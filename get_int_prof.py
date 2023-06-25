import numpy as np
from scipy import interpolate
import scipy.integrate as integrate
import copy
import os
import shutil


debug_printout = False
debug_printout = True
 
def phint(dep,phdep,ztot):
    if ztot>dep[-1]: ztot=dep[-1]
    func = interpolate.interp1d(dep,10.0**-phdep)
    result = integrate.quad(lambda x: func(x), dep[0], ztot)
    a = -np.log10( result[0]/(ztot-dep[0]))
    
    return a

def linave(dep,var,ztrgt):
    if ztrgt>dep[-1]: ztrgt=dep[-1]
    func = interpolate.interp1d(dep,var)
    result = integrate.quad(lambda x: func(x), dep[0], ztrgt)
    a = ( result[0]/(ztrgt-dep[0]))
    
    return a

def get_ph_int_site(outdir,runname,dep_sample):
    # outdir = '../pyweath_output/'
    # runname = ''
    # runname = 'Sheldon_A_fick_noiter_test_psdfullpbe_w1_fit_r2000_sig0p2_loop_q0p2_chkall'

    # dep_sample = float(sys.argv[1])

    for i in range(20,0,-1):
        infile  = outdir+runname+'/prof/prof_aq(tot)-{:03d}.txt'.format(i)
        if not os.path.exists(infile): continue
        else: break 
    infile = outdir+runname+'/prof/prof_aq(tot)-{:03d}.txt'.format(i)
    data = np.loadtxt(infile,skiprows=1)
    print('using aq(tot)-{:03d}'.format(i))
        
    pH_dep = data[:,-2]
    dep = data[:,0]

    # print(dep)
    # print(pH_dep)
            
    phintval =  phint(dep,pH_dep, dep_sample)

    print(phintval)
    
    return phintval

def get_intph_int_site(outdir,runname,dep_sample):

    infile = outdir+runname+'/flx/int_ph.txt'
    data = np.loadtxt(infile,skiprows=1)
    
    with open(infile) as f:
        first_line  = f.readline() 
        dep_list     = first_line.split()
        
    dep = [  float(dep_list[i]) for i in range(1,len(dep_list))   ]
    
    for i in range(len(dep)):
        if dep[i] <=dep_sample:
            idep = i+1
        else:
            continue
    
    print(dep_list[idep])
            
    phintval =  data[-1,idep]

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

    for i in range(20,0,-1):
        infile  = outdir+runname+'/prof/prof_aq(ads%cec)-{:03d}.txt'.format(i)
        if not os.path.exists(infile): continue
        else: break 
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
        
        pH_dep  += data[:,isp]
    
    dep = data[:,0]
            
    phintval    =  linave(dep,pH_dep, dep_sample)

    print(phintval)
    
    return phintval

def get_ac_int_site_v3(outdir,runname,dep_sample):

    for i in range(20,0,-1):
        infile  = outdir+runname+'/prof/prof_aq(ads%cec)-{:03d}.txt'.format(i)
        if not os.path.exists(infile): continue
        else: break 
    infile  = outdir+runname+'/prof/prof_aq(ads%cec)-{:03d}.txt'.format(i)
    data    = np.loadtxt(infile,skiprows=1)
    print('using prof_aq(ads%cec)-{:03d}'.format(i))
    
    infile2  = outdir+runname+'/prof/bsd-{:03d}.txt'.format(i)
    data2    = np.loadtxt(infile2,skiprows=1)
    print('using bsd-{:03d} in {}'.format(i,outdir+runname))
    
    with open(infile) as f:
        first_line  = f.readline() 
        sp_list     = first_line.split()
    
    with open(infile2) as f:
        first_line  = f.readline() 
        prop_list     = first_line.split()
        
    dense = data2[:,prop_list.index('dens[g/cm3]')]
    poro = data2[:,prop_list.index('poro')]
    cec = data2[:,prop_list.index('cec[cmol/kg]')]
        
    sps     = ['h','al']
    
    cecfrac = 0
    for sp in sps:
        try:
            isp     = sp_list.index(sp) 
        except:
            continue
        
        cecfrac  += data[:,isp]
    
    dep = data[:,0]
            
    exchac_int    =  linave(dep,(1-poro)*dense*cec*cecfrac, dep_sample)/linave(dep,(1-poro)*dense*cec, dep_sample)

    print(exchac_int)
    
    return exchac_int

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

    for i in range(20,0,-1):
        infile  = outdir+runname+'/prof/bsd-{:03d}.txt'.format(i)
        if not os.path.exists(infile): continue
        else: break 
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
        pH_dep  += data[:,isp]
        phintval =  linave(dep,pH_dep, dep_sample)
    except:
        print('not exist dens[g/cm3]')
        phintval = 0
            

    print(phintval)
    
    return phintval

def get_sldwt_int_site(outdir,runname,dep_sample,sps):

    for i in range(20,0,-1):
        infile  = outdir+runname+'/prof/prof_sld(wt%)-{:03d}.txt'.format(i)
        if not os.path.exists(infile): continue
        else: break 
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
        
        pH_dep  += data[:,isp]
            
    phintval =  linave(dep,pH_dep, dep_sample)

    print(phintval)
    
    return phintval

def get_totsldwt_site(outdir,runname,dep_sample):

    for i in range(20,0,-1):
        infile  = outdir+runname+'/prof/prof_sld(wt%)-{:03d}.txt'.format(i)
        if not os.path.exists(infile): continue
        else: break 
    infile  = outdir+runname+'/prof/prof_sld(wt%)-{:03d}.txt'.format(i)
    data    = np.loadtxt(infile,skiprows=1)
    print('using prof_sld(wt%)-{:03d}'.format(i))
    
    infile2  = outdir+runname+'/prof/bsd-{:03d}.txt'.format(i)
    data2    = np.loadtxt(infile2,skiprows=1)
    print('using bsd-{:03d} in {}'.format(i,outdir+runname))
    
    dep = data[:,0]
    
    with open(infile) as f:
        first_line  = f.readline() 
        sp_list     = first_line.split()
    
    with open(infile2) as f:
        first_line  = f.readline() 
        prop_list     = first_line.split()
        
    dense = data2[:,prop_list.index('dens[g/cm3]')]
    poro = data2[:,prop_list.index('poro')]
        
    sps = copy.copy(sp_list)
    
    del sps[0]
    del sps[-1]
    
    concs = []
    
    for sp in sps:
        isp     = sp_list.index(sp)
        conc_int = linave(dep,(1-poro)*dense*data[:,isp], dep_sample)/ linave(dep,(1-poro)*dense, dep_sample)
        concs.append( conc_int )  
    
    return sps,concs

def get_gas_int_site(outdir,runname,dep_sample,sp):

    for i in range(20,0,-1):
        infile  = outdir+runname+'/prof/prof_gas-{:03d}.txt'.format(i)
        if not os.path.exists(infile): continue
        else: break 
    infile  = outdir+runname+'/prof/prof_gas-{:03d}.txt'.format(i)
    data    = np.loadtxt(infile,skiprows=1)
    
    if debug_printout: 
        print('using prof_gas-{:03d}'.format(i))
    
    dep = data[:,0]
    
    with open(infile) as f:
        first_line  = f.readline() 
        sp_list     = first_line.split()
        
    isp     = sp_list.index(sp) 
        
    var_dep  = data[:,isp]
            
    var_dep_int =  linave(dep,var_dep, dep_sample)

    if debug_printout: 
        print(var_dep_int)
    
    return var_dep_int

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

def get_ave_DIC_save(outdir,runname,dep_sample):

    for i in range(20,0,-1):
        infile  = outdir+runname+'/prof/prof_aq-{:03d}.txt'.format(i)
        if not os.path.exists(infile): continue
        else: break 
    infile  = outdir+runname+'/prof/charge_balance-{:03d}.txt'.format(i)
    data    = np.loadtxt(infile,skiprows=1)
    print('using charge_balance-{:03d}'.format(i))
    
    
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
        
    DIC_btm = (
        np.average( DIC[:idep+1]*poro[:idep+1]*sat[:idep+1]*1e3 )  # mol/ soil m3
        /np.average( (1 - poro[:idep+1]) )  # mol/ solid m3
        )
    
    np.savetxt(outdir+runname+'/prof/prof_aq-DIC-{:03d}.txt'.format(i),np.transpose(np.array([list(deps),list(DIC)])))
    
    return DIC_btm,dep

def get_ave_DIC_save_v2(outdir,runname,dep_sample):

    for i in range(20,0,-1):
        infile  = outdir+runname+'/prof/prof_aq-{:03d}.txt'.format(i)
        if not os.path.exists(infile): continue
        else: break 
    infile  = outdir+runname+'/prof/charge_balance-{:03d}.txt'.format(i)
    data    = np.loadtxt(infile,skiprows=1)
    print('using charge_balance-{:03d}'.format(i))
    
    
    infile2  = outdir+runname+'/prof/bsd-020.txt'
    data2    = np.loadtxt(infile2,skiprows=1)
    print('using bsd-020')
    
    dep = data[:,0]
    
    with open(infile) as f:
        first_line  = f.readline() 
        sp_list     = first_line.split()
        
    with open(infile2) as f:
        first_line  = f.readline() 
        prop_list     = first_line.split()
        
    sat = data2[:,prop_list.index('sat')]
    poro = data2[:,prop_list.index('poro')]
    
        
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
        
    DIC_ave = (
        linave(dep, DIC*poro*sat*1e3, dep_sample )  # mol/ soil m3
        /linave(dep, (1 - poro), dep_sample )  # mol/ solid m3
        )
    
    np.savetxt(outdir+runname+'/prof/prof_aq-DIC-{:03d}.txt'.format(i),np.transpose(np.array([list(dep),list(DIC)])))
    
    return DIC_ave

def get_ave_IS(outdir,runname,dep_sample):

    for i in range(20,0,-1):
        infile  = outdir+runname+'/prof/prof_aq-{:03d}.txt'.format(i)
        if not os.path.exists(infile): continue
        else: break 
    infile  = outdir+runname+'/prof/charge_balance-{:03d}.txt'.format(i)
    data    = np.loadtxt(infile,skiprows=1)
    print('using charge_balance-{:03d}'.format(i))
    
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
    
    
    IS = data[:,-2]
        
    IS_ave = (
        np.average( IS[:idep+1] )  # mol/ L
        )
    
    return IS_ave,dep

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

def get_totsave_site(outdir,runname,dep_sample):
    
    for i in range(20,0,-1):
        infile  = outdir+runname+'/prof/prof_ex(tot)-{:03d}.txt'.format(i)
        if not os.path.exists(infile): continue
        else: break 
    
    infile  = outdir+runname+'/prof/prof_ex(tot)-{:03d}.txt'.format(i)
    data    = np.loadtxt(infile,skiprows=1)
    print('using prof_ex(tot)-{:03d} in {}'.format(i,outdir+runname))
    
    infile2  = outdir+runname+'/prof/bsd-{:03d}.txt'.format(i)
    data2    = np.loadtxt(infile2,skiprows=1)
    print('using bsd-{:03d} in {}'.format(i,outdir+runname))
    
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
            np.average( data[:idep+1,isp] ) # mol/ soil m3
            /np.average( (1 - poro[:idep+1]) )   # mol/ solid m3
            )  
    
    return sps,btmconcs,dep

def get_totsave_site_v2(outdir,runname,dep_sample):
    
    for i in range(20,0,-1):
        infile  = outdir+runname+'/prof/prof_ex(tot)-{:03d}.txt'.format(i)
        if not os.path.exists(infile): continue
        else: break 
    
    infile  = outdir+runname+'/prof/prof_ex(tot)-{:03d}.txt'.format(i)
    data    = np.loadtxt(infile,skiprows=1)
    print('using prof_ex(tot)-{:03d} in {}'.format(i,outdir+runname))
    
    infile2  = outdir+runname+'/prof/bsd-{:03d}.txt'.format(i)
    data2    = np.loadtxt(infile2,skiprows=1)
    print('using bsd-{:03d} in {}'.format(i,outdir+runname))
    
    infile3  = outdir+runname+'/prof/prof_aq(tot)-{:03d}.txt'.format(i)
    data3    = np.loadtxt(infile3,skiprows=1)
    print('using prof_aq(tot)-{:03d} in {}'.format(i,outdir+runname))
    
    dep = data[:,0]
    
    with open(infile) as f:
        first_line  = f.readline() 
        sp_list     = first_line.split()
        
    with open(infile2) as f:
        first_line  = f.readline() 
        prop_list     = first_line.split()
        
    poro = data2[:,prop_list.index('poro')]
    sat = data2[:,prop_list.index('sat')]
        
    sps = copy.copy(sp_list)
    
    del sps[0]
    del sps[-1]
    del sps[-1]
    
    exchconcs = []
    aqconcs = []
    
    for sp in sps:
        isp     = sp_list.index(sp)
        exchconcs.append( 
            linave(dep, data[:,isp] ,dep_sample) # mol/ soil m3
            /linave(dep, (1 - poro), dep_sample )   # mol/ solid m3
            )  
        aqconcs.append( 
            linave(dep, poro*sat*1e3*data3[:,isp] ,dep_sample) # mol/ soil m3
            /linave(dep, (1 - poro), dep_sample )   # mol/ solid m3
            )  
    
    return sps,exchconcs,aqconcs

def get_adssave_site(outdir,runname,dep_sample):
    
    for i in range(20,0,-1):
        infile  = outdir+runname+'/prof/prof_aq(ads)-{:03d}.txt'.format(i)
        if not os.path.exists(infile): continue
        else: break 
    
    infile  = outdir+runname+'/prof/prof_aq(ads)-{:03d}.txt'.format(i)
    data    = np.loadtxt(infile,skiprows=1)
    print('using prof_aq(ads)-{:03d}'.format(i))
    
    infile2  = outdir+runname+'/prof/bsd-{:03d}.txt'.format(i)
    data2    = np.loadtxt(infile2,skiprows=1)
    print('using bsd-{:03d}'.format(i))
    
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
            np.average( 
                data[:idep+1,isp]  # cmol/kg(solid)
                * 1e-2/1e3  # now mok/g(solid)
                * dense[:idep+1]    # mol/ solid cm3
                * 1./1e-6
                )  
            )  
    
    return sps,btmconcs,dep

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


def main():
    outdir = '../scepter_output/'
    runname = 'US_cropland_251_pph_N_cacl2_2p5_homo2_mgo_field_tpH7p0_tau10'
    dep_sample = 0.15
    get_intph_int_site(outdir,runname,dep_sample)
    
    
    
if __name__=='__main__':
    main()