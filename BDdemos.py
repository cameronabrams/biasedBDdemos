""" Biased Brownian Dynamics Demonstrations 

    Cameron F Abrams, cfa22@drexel.edu

    

"""
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import numpy as np
import argparse as ap
from alive_progress import alive_bar
import pandas as pd

def vb(z,k,zstar):
    """ The bias potential """
    return 0.5*k*(z-zstar)**2,-k*(z-zstar)

def pes(x,*p):
    """ The potential energy (a quartic two-well) """
    a,b,c,d=p
    x2=x*x
    x3=x*x2
    x4=x*x3
    bb=b-2*a
    e=a*x4+bb*x2+c*x+d
    f=-(4*a*x3+2*bb*x+c)
    return e,f

def pes_roots(*p):
    """ returns the roots of the quartic potential """
    a,b,c,d=p
    coeff=[a,0.0,b-2*a,c,d]
    return np.roots(coeff)

def BD(well_positions,args,simtype='unbiased'):
    """ Runs Brownian dynamics; returns array of time points and variable values """
    wm,wp=well_positions
    h1=args.dt/args.gam
    h2=np.sqrt(6*args.T*h1)
    x=wm
    if simtype=='unbiased':
        k=0
        xstar=wp
    elif simtype=='steered':
        k=args.kmax
        xstar=wm
        dxstar=(wp-wm)/args.n
    elif simtype=='coupled':
        k=args.kmin
        dk=(args.kmax-args.kmin)/args.n
        xstar=wp
    elif simtype=='TAMD':
        k=args.kmax
        xstar=wm
        TAMDh1=args.dt/args.gamprime
        TAMDh2=np.sqrt(6*args.Tprime*TAMDh1)
    
    xsav=[]
    tsav=[]
    print(f'running {simtype}...')
    with alive_bar(args.n) as bar:
        for i in range(args.n):
            e_pes,f_pes=pes(x,*(args.pes))
            e_bias,f_bias=vb(x,k,xstar)
            force=f_pes+f_bias
            dx=h1*force+h2*2*(0.5-np.random.random())
            x+=dx
            # bias update
            if simtype=='steered':
                xstar+=dxstar
            elif simtype=='coupled':
                k+=dk
            elif simtype=='TAMD':
                TAMDforce=-f_bias
                dxstar=TAMDh1*TAMDforce+TAMDh2*2*(0.5-np.random.random())
                xstar+=dxstar
            if i==0 or i%args.logevery==0:
                tsav.append(i*args.dt)
                xsav.append(x)
            bar()
    print(f'done running {simtype}')
    return tsav,xsav

if __name__=='__main__':
    parser=ap.ArgumentParser()
    parser.add_argument('-n',type=int,default=20000,help='number of BD time-steps')
    parser.add_argument('-dt',type=float,default=1.e-4,help='BD time step')
    parser.add_argument('-T',type=float,default=1.0,help='reduced temperature')
    parser.add_argument('-gam',type=float,default=5.0,help='BD friction')
    parser.add_argument('-kmin',type=float,default=0.0,help='k min')
    parser.add_argument('-kmax',type=float,default=2.0,help='k max')
    parser.add_argument('-Tprime',type=float,default=20.0,help='TAMD temperature')
    parser.add_argument('-gamprime',type=float,default=500.0,help='TAMD friction')
    parser.add_argument('-pes',type=float,default=[0.02,-1.0,0.0,0.0],help='quartic PES parameters a b c',nargs='+')
    parser.add_argument('-logevery',type=int,default=100,help='log x every this many BD time steps')
    parser.add_argument('-data',type=str,default='bd-runs.csv',help='name of data file to write')
    parser.add_argument('--readonly',default=False,action=ap.BooleanOptionalAction)
    args=parser.parse_args()

    # compute the positions of the wells by finding roots of the derivative of the potential
    nc=np.array([4,3,2,1,0])*np.array([args.pes[0],0,*(args.pes[1:])])
    roots=np.roots(nc)
    wm=roots[1]
    wp=roots[0]
    mid=roots[2]
    fuzz=0.4*(wp-wm)
    # compute and report barrier height
    barrier=pes(mid,*(args.pes))[0]-pes(wm,*(args.pes))[0]
    print(f'Barrier height is {barrier:.2f}')

    data=pd.DataFrame()
    if not args.readonly:
        # Run the four simulations and save resulting data
        for simtype in ['unbiased','steered','coupled','TAMD']:
            tsav,xsav=BD([wm,wp],args,simtype)
            if data.empty:
                # time points are identical among runs so only save first one
                data['time']=tsav
            data[simtype]=xsav
        # write the data to a CSV file
        data.to_csv(args.data,header=True,index=False,float_format='%.5f')
    else:
        # Read existing data in from CSV file
        print(f'Reading {args.data}')
        data=pd.read_csv(args.data,header=0,index_col=None)
    
    # handles for each variable trajectory
    tsav=data['time']
    xsav=data['unbiased']
    xssav=data['steered']
    xcsav=data['coupled']
    xTsav=data['TAMD']

    # generate the plot
    fig,ax=plt.subplots(1,1,figsize=(6,4))
    ax.set_xlabel(r'$t$',fontsize=14)
    ax.tick_params(axis='x',which='both',bottom=False,top=False,labelbottom=False)
    ax.set_yticks([wm,wp])
    ax.set_yticklabels([r'$z_A$',r'$z_B$'],fontsize=14)
    ax.set_ylim([wm-fuzz,wp+2*fuzz])
    ax.plot(tsav,xsav,'k-',label=r'unbiased $T$=1 $\gamma$=5',alpha=0.6)
    ax.plot(tsav,xssav,'r-',label=r'steered $k$='+f'{args.kmax:.1f}',alpha=0.6)
    ax.plot(tsav,xcsav,'b-',label=r'coupled $k_{\rm max}$='+f'{args.kmax:.1f}',alpha=0.6)
    ax.plot(tsav,xTsav,'g-',label=r'TAMD $k$='+f'{args.kmax:.1f}'+r' $T^\prime/T$='+f'{args.Tprime/args.T:.0f}'+r' $\gamma^\prime/\gamma$='+f'{args.gamprime/args.gam:.0f}',alpha=0.6)
    ax.legend(loc='upper center',ncol=2)
    plt.savefig('pbd.png',bbox_inches='tight')

