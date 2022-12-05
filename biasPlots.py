""" Generate representative bias evolution plots 

    Cameron F. Abrams, cfa22@drexel.edu
    
"""
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import numpy as np

def v(z,k,zstar):
    return 0.5*k*(z-zstar)**2

fig,ax=plt.subplots(2,1,figsize=(4,4),sharex=True)

D=np.linspace(0,1,200)

nk=20
#k0,k1=4,4000
k=4*np.logspace(0,3,nk,base=10)
za,zb=0.2,0.8
nz=20
ax[0].set_xticks([zb,za])
ax[1].set_xticks([zb,za])
ax[1].set_xticklabels([r'$z_B$',r'$z_A$'],fontsize=14)
z=np.linspace(za,zb,nz)

for i in [0,1]:
    ax[i].set_ylim((0,10))
    ax[i].set_xlim((0,1))
    ax[i].tick_params(
        axis='y',          # changes apply to the x-axis
        which='both',      # both major and minor ticks are affected
        left=False,      # ticks along the bottom edge are off
        right=False,         # ticks along the top edge are off
        labelleft=False) # labels along the bottom edge are off
#ax.set_xlabel(r'$z$')
    ax[i].set_ylabel(r'$V^{\rm (b)}(z)$',fontsize=14)

ff=0.2
cmap=cm.get_cmap('Reds')
for i in range(nz):
    ax[0].plot(D,v(D,k[-1],z[i]),color=cmap((i+ff*nz)/((1+ff)*nz)),alpha=0.75)

cmap=cm.get_cmap('Blues')
for i in range(nk):
    ax[1].plot(D,v(D,k[i],z[-1]),color=cmap((i+ff*nk)/((1+ff)*nk)),alpha=0.5)

plt.savefig('BiasEvolutions.png',bbox_inches='tight')
