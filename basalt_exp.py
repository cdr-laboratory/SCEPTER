# -*- coding: utf-8 -*-
import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
import glob
import subprocess

outdir = '../pyweath_output/'
runname = 'sil+ph_wet_iter---q1E-1_zsat5_basalt_test_cpl_high-rain'
simlist = glob.glob(outdir+runname+'/o2profile-res-01*.txt')
datalist = []
for i in simlist:
    datalist.append(np.loadtxt(i))

simlist = glob.glob(outdir+runname+'/o2profile-res(rate)-01*.txt')
ratelist = []
for i in simlist:
    ratelist.append(np.loadtxt(i))

simlist = glob.glob(outdir+runname+'/o2profile-bsd-01*.txt')
baselist = []
for i in simlist:
    baselist.append(np.loadtxt(i))

agelist = []
for i in range(len(datalist)):
    agelist.append(datalist[i][0,datalist[i].shape[1]-1])

c = np.array(agelist)
mymap = mpl.colors.LinearSegmentedColormap.from_list('mycolors',['black','red'])
mymap = mpl.cm.jet
norm = mpl.colors.Normalize(vmin=c.min(), vmax=c.max())
cmap = mpl.cm.ScalarMappable(norm=norm, cmap=mymap)
cmap.set_array([])

figsize = (7,5)
fig = plt.figure(figsize=figsize)

nx =3 
ny =2 

axes = [[plt.subplot2grid((ny,nx), (j,i)) for i in range(nx)] for j in range(ny)]

po2char = 'pO'+r'$\mathregular{_2}$'
omegaabchar = r'$\mathregular{\Omega_{\rm{ab}}}$'
omegafochar = r'$\mathregular{\Omega_{\rm{fo}}}$'
fe2char = 'Fe(II)'
fe3char = 'Fe(III)'
so4char = 'Sulfate'
nachar = 'Na'
mgchar = 'Mg'
porochar = r'$\mathregular{\phi}$'
satchar = r'$\mathregular{\sigma}$'
advchar = r'$\mathregular{\it{v}}$'
diffchar = r'$\mathregular{\it{D}}$'

numplt = [
    [(10,':','Forsterite'),(9,'-','Albite')]
    ,[(13,':',omegafochar),(12,'-',omegaabchar)]
    ,[(7,':',mgchar),(6,'-',nachar)]
    ,[(11,'-','pH')]
    ,[(6,':','Forsterite'),(5,'-','Albite')]
    ,[(1,':',porochar),(2,'dashdot',satchar),(3,'--',advchar),(4,'-',diffchar)]
    ]
    
xlabels = [
    'mol m'+r'${^{-3}}$'
    ,'dimensionless'
    ,'mol L'+r'${^{-1}}$'
    ,'pH'
    ,'mol m'+r'${^{-3}}$'+' yr'+r'${^{-1}}$'
    ,'m'+' yr'+r'${^{-1}}$'+ ', m'+r'${^2}$'+' yr'+r'${^{-1}}$' + ' or m' +r'${^3}$'+' m'+r'${^{-3}}$'
    ]

for k in range(nx*ny):
    i = k%nx
    j = (k-i)/nx
    for o in range(len(datalist)):
        for p in numplt[k]:
            if o==0:
                if k==4:
                    axes[j][i].plot(ratelist[o][:,p[0]],ratelist[o][:,0],linestyle = p[1],c = cmap.to_rgba(c[o]),label = p[2])
                elif k==5:
                    axes[j][i].plot(baselist[o][:,p[0]],baselist[o][:,0],linestyle = p[1],c = cmap.to_rgba(c[o]),label = p[2])
                else:
                    axes[j][i].plot(datalist[o][:,p[0]],datalist[o][:,0],linestyle = p[1],c = cmap.to_rgba(c[o]),label = p[2])
            else:
                if k==4:
                    axes[j][i].plot(ratelist[o][:,p[0]],ratelist[o][:,0],linestyle = p[1],c = cmap.to_rgba(c[o]))
                elif k==5:
                    axes[j][i].plot(baselist[o][:,p[0]],baselist[o][:,0],linestyle = p[1],c = cmap.to_rgba(c[o]))
                else:
                    axes[j][i].plot(datalist[o][:,p[0]],datalist[o][:,0],linestyle = p[1],c = cmap.to_rgba(c[o]))
    axes[j][i].invert_yaxis()
    axes[j][i].set_xlabel(xlabels[k])
    if k==5:axes[j][i].set_xscale('log')
    if i==0:axes[j][i].set_ylabel('Depth (m)')
    axes[j][i].legend()


cbaxes = fig.add_axes([0.88, 0.175, 0.01, 0.75]) 
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