import os
import numpy as np
import time 

input_dir = './'
input_file = 'friendship.txt'

data = np.loadtxt(input_file)

n_sample = data.shape[0]

runs_chosen = [
    # 2,6,7
    # 7
    ]

runtype = 'spinup'
# runtype = 'basalt'

for i in range(n_sample):


    runtime = '12:20:00'
    soft    = 'python3'

    runid   = data[i,0]
    soilpH  = data[i,1]
    om      = data[i,2]
    cec     = data[i,3]
    acid    = data[i,4]
    buffpH  = data[i,5]
    
    
    if len(runs_chosen)!=0 and runid not in runs_chosen: continue
    
    if acid ==0:acid = 0.1 
    if soilpH  ==5.6: soilpH  = 5.7 
    
    if runtype == 'spinup':
        
        code    = 'tunespin_3_newton_inert_buff.py'
        # code    = 'tunespin_3b_newton_inert_buff.py'
        
        if code == 'tunespin_3b_newton_inert_buff.py' and np.isnan(buffpH): continue
    
        # runname = 'potato_'+str(int(runid))
        # runname = 'potato_pw_'+str(int(runid))
        # runname = 'potato_fert_'+str(int(runid))
        # runname = 'potato_fert_pw_'+str(int(runid))
        runname = 'potato_buff_'+str(int(runid))

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
    
        code    = 'basalt_buff_tunespin_bisec.py'
    
        # runname = 'potato_bas_'+str(int(runid))
        # runname = 'potato_bas_pw_'+str(int(runid))
        # runname = 'potato_bas_fert_50u_'+str(int(runid))
        runname = 'potato_bas_50u_'+str(int(runid))
        
        spinname = 'potato_'+str(int(runid))
        # spinname = 'potato_pw_'+str(int(runid))
        # spinname = 'potato_fert_'+str(int(runid))
        
        targetpHs = [6.5, 6.8]
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