#!/usr/bin/env python

"""
Jens Luebeck
UC San Diego, Bioinformatics & Systems Biology
BNFO 284 (Prof. Jeff Hasty)
12/9/16
jluebeck@ucsd.edu

Simulates Gray-Scott reaction diffusion model and can produce images for use in animations

Example usage:

    python GrayScottSimulation -o my_simulation -moviemode -p dots -n 10000

    This will output files with the prefix "OutputPrefix" into a directory of the same name, and as
the moviemode flag is set, it will store 250 images for animation. The -p command selects the preset
parameter group

Example usage 2:

    python GrayScottSimulation -o my_simulation2 -gm -n 5000

    Instead uses the Gierer-Meinhardt activator-inhibitor model (-gm). Does not support -p arguments.

Example usage 3:

    python GrayScottSimulation -o my_simulation3 -p manual

    Simulates Gray Scott using the 'myParams' values given immediately below. Alter as desired.

"""

#Below dictionary must maintain the same set of keys otherwise it will fail.
myParams = {"Du":0.19, "Dv":0.05, "F":0.06, "k":0.062, "ltype":"single","myCmap":plt.cm.cubehelix,"edgeMax":False,"dt":1,"dx":1}
#Du: Diffusivity of U
#Dv: Diffusivity of V
#F: "feed rate"
#k: dimensionless rate constant for V
#ltype: initialization [single | dual | random]
#myCmap: colormap object compatible with matplotlib
#edgeMax: Set a constant for heatmap scaling purposes (True | False)
#dt: timestep
#dx: spatial stepsize (same as dy)

import sys
import os
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import random
import copy
import argparse

#makes a heatmap of the given matrix (M)
def makeImg(M,fname, myCmap, colorbar=False, bg='black', setEdge=True):
    plt.figure()
    plt.rcParams['axes.facecolor'] = bg
    plt.rcParams['savefig.facecolor'] = bg
    plt.axis('off')
    #Hackish way to ensure constant color scale across images
    if setEdge:
        M[-1,-1] = 1
    plt.imshow(M, cmap=myCmap, extent=[-1,1,-1,1]);
    if colorbar:
        plt.colorbar()
    #reset value
    if setEdge:
        M[-1,-1] = 0

    plt.savefig(savPath + runName + "_" + fname + ".png",dpi=myDPI)
    plt.close()

#uses the five-point stencil method of finite differences to estimate the laplace operator
def laplacian_operator(U,V,dx):
    Lu = (U[0:-2,1:-1] + U[1:-1,0:-2] + U[1:-1,2:] + U[2:,1:-1] - 4*U[1:-1,1:-1])/dx**2
    Lv = (V[0:-2,1:-1] + V[1:-1,0:-2] + V[1:-1,2:] + V[2:,1:-1] - 4*V[1:-1,1:-1])/dx**2
    return Lu,Lv

#Numerical simulation function for Gray Scott
def GrayScott(params,initial_matrices):
    ltype,edgeMax,dx,F,Dv,k,Du,dt,myCmap = [params[x] for x in params.keys()]
    U,V = initial_matrices
    u,v = U[1:-1,1:-1], V[1:-1,1:-1]
    for i in range(n):
        Lu,Lv = laplacian_operator(U,V,dx)
        uvv = u*v*v
        su = Du*Lu - uvv + F *(1-u)
        sv = Dv*Lv + uvv - (F+k)*v
        u += dt*su
        v += dt*sv

        if movieOutput:
            #Some manually set initial frames to grab so we can see how the system evolves early on
            if i % frameMod == 0 or i in [1,2,3,4,5,10,15,20,30,40,50,60,70,80,90,100,110,120,150]:
                makeImg(v,"v_" + str(i),params["myCmap"],setEdge=params["edgeMax"])
                if i % 1000*frameMod == 0:
                    print str(i)

    return u,v

#Numerical simulation function for Gierer Meinhardt
def GM(params,initial_matrices):
    ltype,edgeMax,dx,F,Dv,k,Du,dt,myCmap = [params[x] for x in params.keys()]
    U,V = initial_matrices
    u,v = U[1:-1,1:-1], V[1:-1,1:-1]
    for i in range(n):
        Lu,Lv = laplacian_operator(U,V,dx)
        vv = v*v
        s = .01 * (.99 + .01 * np.random.rand(len(u),len(u)))
        sv = s*(vv/(u*(1+0.3*vv))) - 0.01*v + Dv*Lv
        su = s*vv - 0.015*u + Du*Lu
        u += dt*su
        v += dt*sv

        if movieOutput:
            #Some manually set initial frames to grab so we can see how the system evolves early on
            if i % frameMod == 0 or i in [1,2,3,4,5,10,40,80,150]:
                makeImg(v,"v_" + str(i),params["myCmap"],setEdge=params["edgeMax"])
                if i % 1000*frameMod == 0:
                    print str(i)

    return u,v

if __name__ == "__main__":
    #Parses the command line arguments
    parser = argparse.ArgumentParser(description="Gray-Scott simulation")
    parser.add_argument("-o", "--outname", help="simulation run output name",required=True)
    parser.add_argument("-m", "--moviemode", action='store_true',help="Run the script in \"movie mode\"\
        to store sequential images (default Off)",default=False)
    parser.add_argument("-p", "--preset",choices=['dots', 'lines', 'unstable','waves','manual'], \
        help="Run the script with preset parameters (default lines). Choose 'manual' if \
        you have set parameters in the script itself (see beginning of script)",default="lines")
    parser.add_argument("-n", "--timesteps", type=int, help="Number of timesteps to simulate (default 8000)",default=8000)
    parser.add_argument("-gm",action="store_true",help="Run this to simulate Gierer-Meinhardt\
        activator-inhibitor equation instead. Overrides \"-p\" argument (Default False).")


    args = parser.parse_args()
    print args

    #create output dir
    runName = args.outname
    savPath = args.outname + "_images/"
    if not os.path.exists(savPath):
        os.makedirs(savPath)

    #Selects initial conditions depending on pattern requested
    args.preset = args.preset.lower()
    if args.preset == "dots":
        params = {"Du":0.14, "Dv":0.06, "F":0.035, "k":0.065, "ltype":"single","myCmap":plt.cm.copper,"edgeMax":False,"dt":1,"dx":1}

    elif args.preset == "unstable":
        params = {"Du":0.16, "Dv":0.08, "F":0.02, "k":0.055, "ltype":"dual","myCmap":plt.cm.cubehelix,"edgeMax":True,"dt":.9,"dx":.9}

    elif args.preset == "waves":
        params = {"Du":0.12, "Dv":0.08, "F":0.02, "k":0.05, "ltype":"single","myCmap":cmocean.cm.ice,"edgeMax":True,"dt":.9,"dx":.9}

    elif args.preset == "lines":
        params = {"Du":0.19, "Dv":0.05, "F":0.06, "k":0.062, "ltype":"single","myCmap":plt.cm.cubehelix,"edgeMax":False,"dt":.9,"dx":.9}

    else:
        if args.gm:
            print "GIERER-MEINHARDT MODEL SELECTED"
            params = {"Du":0.2, "Dv":0.002, "F":0, "k":0.0, "ltype":"rand","myCmap":plt.cm.copper,"edgeMax":True,"dt":1,"dx":1}

        else:
            params = myParams


    #set up image saving
    totFrames = 200
    movieOutput = False
    n = args.timesteps
    print "Running simulation with " + str(n) + " timesteps."
    frameMod = n/totFrames
    size = 200

    myDPI = 300 #Image resolution DPI
    if args.moviemode:
        print "Movie mode set to ON"
        myDPI = 200 #reduce image DPI if movie mode
        movieOutput = True
        print str(totFrames) + " movie images will be produced"
    else:
        print "Movie mode set to OFF"


    #set initial conditions
    U = np.zeros((size, size))
    V = np.zeros((size, size))
    u,v = U[1:-1,1:-1], V[1:-1,1:-1]
    u+=1.0
    #sets initialization of single or double squares, or completely random
    if params["ltype"] == "single":
        r = 20
        U[size/2-r:size/2+r,size/2-r:size/2+r] = 0.50
        V[size/2-r:size/2+r,size/2-r:size/2+r] = 0.25
    elif params["ltype"] == "dual":
        r = 15
        U[size/4-r:size/4+r,size/4-r:size/4+r] = 0.50
        V[size/4-r:size/4+r,size/4-r:size/4+r] = 0.25
        U[3*size/4-r:3*size/4+r,3*size/4-r:3*size/4+r] = 0.50
        V[3*size/4-r:3*size/4+r,3*size/4-r:3*size/4+r] = 0.25

    else:
        u-=1
        u+=np.random.rand(len(u),len(u))
        v+=np.random.rand(len(u),len(u))

    #add small amount of noise everywhere
    u += (0.05 + 0.05*(np.random.random((size-2,size-2))*2-1))
    v += (0.05 + 0.05*(np.random.random((size-2,size-2))*2-1))

    initial_matrices = (U,V)

    makeImg(v,"initial_v",params["myCmap"],setEdge=params["edgeMax"])

    #RUN SIM
    if args.gm:
        u,v = GM(params,initial_matrices)
    else:
        u,v = GrayScott(params,initial_matrices)

    makeImg(v,"final_v",params["myCmap"],setEdge=params["edgeMax"])
