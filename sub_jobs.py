import os
import numpy as np
import time 

input_dir = './data/'
# input_dir = './'
input_file = 'friendship.txt'
input_file = 'GAdata.txt'
input_file = 'SL004810.txt'
skip_file = 'GA_fail.txt'
skip_file = 'SL004810_sph_fail.txt'
skip_file = 'SL004810_FAIL_sph_N_spintuneup_field.txt'

data = np.loadtxt(input_dir+input_file,skiprows=1)
skipdata = np.loadtxt(input_dir+skip_file)
skipdata = [int(i) for i in skipdata]

n_sample = data.shape[0]

runs_chosen = [
    # 2,6,7
    # 7
    ]

runtype = 'spinup'
runtype = 'basalt'

for i in range(n_sample):


    runtime = '12:20:00'
    soft    = 'python3'
    
    if 'friendship' in input_file:

        runid   = data[i,0]
        soilpH  = data[i,1]
        om      = data[i,2]
        cec     = data[i,3]
        acid    = data[i,4]
        buffpH  = data[i,5]
    
        if acid ==0:acid = 0.1 
        if soilpH  ==5.6: soilpH  = 5.7 
    
    elif 'GAdata' in input_file:
        runid   = data[i,0]
        soilpH  = data[i,8]
        om      = 5  # no data given
        cec     = data[i,1]
        acid    = data[i,6]
        buffpH  = 7  # no data given
    
    elif 'SL004810' in input_file:

        runid   = data[i,0]
        soilpH  = data[i,1]
        om      = data[i,2]
        cec     = data[i,3]
        acid    = data[i,4]
        buffpH  = 7  # no data given
    
    
    if len(runs_chosen)!=0 and runid not in runs_chosen: continue
    
    if runtype == 'spinup':
        
        code    = 'tunespin_3_newton_inert_buff.py'
        code    = 'tunespin_3_newton_inert_buff_v2.py'
        # code    = 'tunespin_3b_newton_inert_buff.py'
        
        if code == 'tunespin_3b_newton_inert_buff.py' and np.isnan(buffpH): continue
    
        # runname = 'potato_'+str(int(runid))
        # runname = 'potato_pw_'+str(int(runid))
        # runname = 'potato_fert_'+str(int(runid))
        # runname = 'potato_fert_pw_'+str(int(runid))
        runname = 'potato_buff_'+str(int(runid))
        runname = 'GA_OM5_sph_'+str(int(runid))
        runname = input_file.replace('.txt','')+'_'+str(int(runid))
        
        # runname += '_sph'
        runname += '_sph_N'

        jobname = runname

        inputs = [
            runname
            ,str(cec)
            ,str(soilpH)
            ,str(acid)
            ,str(om)
            ]
        
        if code == 'tunespin_3b_newton_inert_buff.py':
            inputs = [
                runname
                ,str(cec)
                ,str(soilpH)
                ,str(buffpH)
                ,str(om)
                ]
        
        input_runtime = ' '.join(inputs)

        command = 'sbatch  --output='+jobname+'.out --job-name='+jobname+' --time='+runtime + ' run_a_shell.sbatch ' + soft + ' ' + code  \
            + ' ' + input_runtime
            
        print(command)
        os.system(command)
        
        time.sleep(5)
    
    elif runtype == 'basalt':
    
        if acid <= 0.1 : continue
        
        # if 'GAdata' in input_file and int(runid) in skipdata: continue
        if 'SL004810' in input_file and int(runid) in skipdata: continue
    
        code    = 'basalt_buff_tunespin_bisec.py'
        code    = 'basalt_buff_tunespin_bisec_v2.py'
    
        # runname = 'potato_bas_'+str(int(runid))
        # runname = 'potato_bas_pw_'+str(int(runid))
        # runname = 'potato_bas_fert_50u_'+str(int(runid))
        # runname = 'potato_bas_50u_'+str(int(runid))
        # runname = 'potato_bas_50u_'+str(int(runid))
        # runname = 'GA_OM5_sph_bas_10u_'+str(int(runid))
        runname = input_file.replace('.txt','')+'_'+str(int(runid))
        
        # spinname = 'potato_'+str(int(runid))
        # spinname = 'potato_pw_'+str(int(runid))
        # spinname = 'potato_fert_'+str(int(runid))
        spinname = 'GA_OM5_sph_'+str(int(runid))
        spinname = input_file.replace('.txt','')+'_'+str(int(runid))
        
        # spinname += '_N'
        # runname  += '_N_10u_bas'
        # spinname += '_h'
        # runname  += '_10u_bas'
        # spinname += '_sph'
        # runname  += '_sph_10u_bas'
        spinname += '_sph_N'
        runname  += '_sph_N_10u_bas'
        
        targetpHs = [6.5, 6.8]
        targetpHs = [5.4, 5.6, 5.8]
        targettaus = [1]
        
        for targetpH in targetpHs:
            
            if soilpH >=targetpH: continue
            
            for targettau in targettaus:
        
                jobname = runname+'_'+'_'.join([str(targetpH),str(targettau)]).replace('.','p')
                
                
                inputs = [
                    str(targetpH)
                    ,str(targettau)
                    ,str(cec)
                    ,runname
                    ,spinname
                    ]

                input_runtime = ' '.join(inputs)

                command = 'sbatch  --output='+jobname+'.out --job-name='+jobname+' --time='+runtime + ' run_a_shell.sbatch ' + soft + ' ' + code  \
                    + ' ' + input_runtime
                    
                print(command)
                os.system(command)
                
                time.sleep(5)
    # break