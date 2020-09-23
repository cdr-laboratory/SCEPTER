# -*- coding: utf-8 -*-
import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
import glob,os
import subprocess

outdir = '../pyweath_output/'

simlist = os.listdir(outdir)
intlist = [int(i+1) for i in range(len(simlist))]
simlist_show = [str(intlist[i])+'--'+simlist[i] for i in range(len(simlist))]
simlist_show_2='\n'.join(simlist_show)
str_tmp = input('What is your simulation?: choose from the list below\n'+simlist_show_2+'\n')
for i in intlist:
    if eval(str_tmp)== i: 
        runname = simlist[i-1]
        break

simlist = glob.glob(outdir+runname+'/o2profile-res-0*.txt')
datalist = []
for i in simlist:
    datalist.append(np.loadtxt(i))

simlist = glob.glob(outdir+runname+'/o2profile-res(rate)-0*.txt')
ratelist = []
for i in simlist:
    ratelist.append(np.loadtxt(i))

simlist = glob.glob(outdir+runname+'/o2profile-bsd-0*.txt')
baselist = []
for i in simlist:
    baselist.append(np.loadtxt(i))

agelist = []
for i in range(len(datalist)):
    agelist.append(datalist[i][0,datalist[i].shape[1]-1])

c = np.array(agelist)
mymap = mpl.colors.LinearSegmentedColormap.from_list('mycolors',['0.3','g'])
mymap = mpl.cm.turbo
norm = mpl.colors.Normalize(vmin=c.min(), vmax=c.max())
cmap = mpl.cm.ScalarMappable(norm=norm, cmap=mymap)
cmap.set_array([])

figsize = (7,5)
fig = plt.figure(figsize=figsize)

nx =3 
ny =2 

axes = [[plt.subplot2grid((ny,nx), (j,i)) for i in range(nx)] for j in range(ny)]

po2char = r'$\it{p}\rm{O}_2^{\rm{soil}}$'
pco2char = r'$\it{p}\rm{CO}_2^{\rm{soil}}$'
omegaabchar = r'$\mathregular{\Omega_{\rm{ab}}}$'
omegaanchar = r'$\mathregular{\Omega_{\rm{an}}}$'
omegafochar = r'$\mathregular{\Omega_{\rm{fo}}}$'
omegaccchar = r'$\mathregular{\Omega_{\rm{cc}}}$'
fe2char = 'Fe(II)'
fe3char = 'Fe(III)'
so4char = 'Sulfate'
nachar = 'Na'
mgchar = 'Mg'
cachar = 'Ca'
sichar = 'Si'
porochar = r'$\mathregular{\phi}$'
satchar = r'$\mathregular{\sigma}$'
advchar = r'$\mathregular{\it{v}}$'
diffchar = r'$\mathregular{\it{D}}$'

numplt = [
    [(12,':','Forsterite'),(10,'-','Albite'),(11,'-.','Anorthite'),(13,'--','Calcite')]
    ,[(14,'-',omegaabchar),(15,':',omegafochar),(16,'-.',omegaanchar),(17,'--',omegaccchar)]
    ,[(8,'--',mgchar),(6,':',nachar),(7,'-',cachar),(9,'-.',sichar)]
    ,[(19,'-','pH')]
    ,[(1,':',po2char),(18,'-',pco2char)]
    ,[(1,':',porochar),(5,'-','SA')]
    ]
    
xlabels = [
    'mol m'+r'${^{-3}}$'
    ,'dimensionless'
    ,'mol L'+r'${^{-1}}$'
    ,'pH'
    # ,'ppmv (CO'+r'$_2$'+') or '+u'$‰$'+ ' (O'+r'$_2$'+')'
    ,'PAL'
    ,'10'+r'${^6}$'+' m'+r'${^2}$'+' m'+r'${^{-3}}$' + ' or m' +r'${^3}$'+' m'+r'${^{-3}}$'
    ]

for base in baselist:
    base[:,5] = base[:,5]/1e6
    
for data in datalist:
    data[:,1] = data[:,1]/0.21
    data[:,18] = data[:,18]/(10**-3.5)

for k in range(nx*ny):
    i = int(k%nx)
    j = int((k-i)/nx)
    for o in range(len(datalist)):
        for p in numplt[k]:
            if o==0:
                if k==4:
                    axes[j][i].plot(datalist[o][:,p[0]],datalist[o][:,0],linestyle = p[1],c = cmap.to_rgba(c[o]),label = p[2])
                elif k==5:
                    axes[j][i].plot(baselist[o][:,p[0]],baselist[o][:,0],linestyle = p[1],c = cmap.to_rgba(c[o]),label = p[2])
                else:
                    axes[j][i].plot(datalist[o][:,p[0]],datalist[o][:,0],linestyle = p[1],c = cmap.to_rgba(c[o]),label = p[2])
            else:
                if k==4:
                    axes[j][i].plot(datalist[o][:,p[0]],datalist[o][:,0],linestyle = p[1],c = cmap.to_rgba(c[o]))
                elif k==5:
                    axes[j][i].plot(baselist[o][:,p[0]],baselist[o][:,0],linestyle = p[1],c = cmap.to_rgba(c[o]))
                else:
                    axes[j][i].plot(datalist[o][:,p[0]],datalist[o][:,0],linestyle = p[1],c = cmap.to_rgba(c[o]))
    axes[j][i].invert_yaxis()
    axes[j][i].set_xlabel(xlabels[k])
    # if k==5:axes[j][i].set_xscale('log')
    if i==0:axes[j][i].set_ylabel('Depth (m)')
    axes[j][i].legend(loc = 'lower center')


cbaxes = fig.add_axes([0.88, 0.175, 0.02, 0.75]) 
# cticks = [500,400,300,200,100,0]
cbar = fig.colorbar(cmap, cax = cbaxes, orientation = 'vertical')
cbar.set_label('Time (yr)', rotation=90)

fig.subplots_adjust(left=0.1,right= 0.85 ,bottom=0.12,top = 0.95, wspace = 0.2, hspace = 0.3)

plt.savefig(outdir+runname+'\profiles.svg', transparent=True)
plt.savefig(outdir+runname+'\profiles.pdf', transparent=True)
svgname = outdir+runname+'\profiles.svg'
subprocess.call('"C:\Program Files\Inkscape\inkscape.exe" -z -f ' 
    + svgname + ' --export-emf '+svgname[:-4] +
    '.emf',shell=True)

plt.show()
plt.clf()
plt.close()