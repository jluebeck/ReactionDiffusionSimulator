#!/usr/bin/env python

"""
Jens Luebeck
UC San Diego, Bioinformatics & Systems Biology
jluebeck@ucsd.edu

Simulates Gray-Scott reaction diffusion model and can produce images for use in animations

Example usage:

    python ReactionDiffusion.py -o my_simulation --moviemode -n 10000

    This will output files with the prefix "OutputPrefix" into a directory of the same name, and as
the moviemode flag is set, it will store 250 images for animation. User will be prompted for model type

Example usage 2:

    python ReactionDiffusion.py -o my_simulation2 -m GM -n 5000

    Instead uses the Gierer-Meinhardt activator-inhibitor model (-gm).

"""

#Other notes:

#Du: Diffusivity of U
#Dv: Diffusivity of V
#F: "feed rate"
#k: dimensionless rate constant for V
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

#Numerical simulation function for Gray-Scott equations
def GS(params,initial_matrices):
    Du,Dv,k,F,dt,dx = params['Du'],params['Dv'],params['k'],params['F'],params['dt'],params['dx']
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

#Numerical simulation function for Gierer-Meinhardt equations
def GM(params,initial_matrices):
    Du,Dv,ku,kv,sv,dt,dx = params['Du'],params['Dv'],params['ku'],params['kv'],params['sv'],params['dt'],params['dx']
    U,V = initial_matrices
    u,v = U[1:-1,1:-1], V[1:-1,1:-1]
    for i in range(n):
        Lu,Lv = laplacian_operator(U,V,dx)
        vv = v*v
        s = kv * (.99 + .01 * np.random.rand(len(u),len(u)))
        sv = s*(vv/(u*(1+sv*vv))) - kv*v + Dv*Lv
        su = s*vv - ku*u + Du*Lu
        u += dt*su
        v += dt*sv

        if movieOutput:
            #Some manually set initial frames to grab so we can see how the system evolves early on
            if i % frameMod == 0 or i in [1,2,3,4,5,10,40,80,150]:
                makeImg(v,"v_" + str(i),params["myCmap"],setEdge=params["edgeMax"])
                if i % 1000*frameMod == 0:
                    print str(i)

    return u,v

#Numerical simulation function for FitzHugh-Nagumo equations
def FN(params,initial_matrices):
    U,V = initial_matrices
    u,v = U[1:-1,1:-1], V[1:-1,1:-1]
    Du,Dv,k,tau,dt,dx = params['Du'],params['Dv'],params['k'],params['tau'],params['dt'],params['dx']
    for i in range(n):
        Lu,Lv = laplacian_operator(U,V,dx)
        sv = Dv*Lv + v - v*v*v - u + k
        su = (Du*Lu + v - u)/tau
        u += dt*su
        v += dt*sv

        if movieOutput:
            #Some manually set initial frames to grab so we can see how the system evolves early on
            if i % frameMod == 0 or i in [1,2,3,4,5,10,40,80,150]:
                makeImg(v,"v_" + str(i),params["myCmap"],setEdge=params["edgeMax"])
                if i % 1000*frameMod == 0:
                    print str(i)

    return u,v

#Set model function and params
def setModelParams(model):
    global size
    if model == "FN":
        print "FitzHugh-Nagumo model selected"
        modelFunc = FN
        print "SETTING GRID SIZE TO 120 - STABILITY REASONS"
        print "WARNING: IF TIMESTEPS LESS THAN 80000, NO PATTERN MAY APPEAR"
        #FitzHugh-Nagumo requires fine-timestepping
        size = 120
        dx = 2./size
        dt = 0.9 * dx**2/2
        params = {"Du":5e-3, "Dv":2.8e-4, "tau":0.1, "k":-0.005,"myCmap":plt.cm.PRGn,"edgeMax":False,"dt":dt,"dx":dx,"seed":"noise"}
        params["seed"] = "noise"

    elif model == "GM":
        print "Gierer-Meinhardt model selected"
        modelFunc = GM
        params = {"Du":0.2, "Dv":0.002, "ku":0.015, "kv":0.01, "sv":0.3,"myCmap":plt.cm.copper,"edgeMax":False,"dt":1,"dx":1,"seed":"noise"}
        params["seed"] = "noise"

    elif model == "GS":
        print "Gray-Scott model selected"
        modelFunc = GS

    if model =="GS":
        pnames = ["solitons","coral","maze","waves","flicker","worms"]
        pvals = [{"Du":0.14, "Dv":0.06, "F":0.035, "k":0.065,"myCmap":plt.cm.copper,"edgeMax":False,"dt":1,"dx":1}, \
                {"Du":0.16, "Dv":0.08, "F":0.06, "k":0.062,"myCmap":plt.cm.cubehelix,"edgeMax":False,"dt":1,"dx":1}, \
                {"Du":0.19, "Dv":0.05, "F":0.06, "k":0.062,"myCmap":plt.cm.cubehelix,"edgeMax":False,"dt":1,"dx":1}, \
                {"Du":0.12, "Dv":0.08, "F":0.02, "k":0.05,"myCmap":plt.cm.cubehelix,"edgeMax":True,"dt":1,"dx":1},\
                {"Du":0.16, "Dv":0.08, "F":0.02, "k":0.055,"myCmap":plt.cm.cubehelix,"edgeMax":True,"dt":1,"dx":1},\
                {"Du":0.16, "Dv":0.08, "F":0.054, "k":0.064,"myCmap":plt.cm.cubehelix,"edgeMax":False,"dt":1,"dx":1}]

        pchoices = dict(zip(pnames,pvals))
        pattern = ""
        while not pattern:
            print "--- Pattern choices ---"
            print "pick 1: " + ", ".join(pnames)
            pattern = raw_input('Select reaction-diffusion model: ').rstrip().lower()
            if not pattern in pchoices:
                print "please select pattern name from list"
                pattern = ""

        params = pchoices[pattern]

        #Set initial seed pattern for GS
        seeds = ["single","dual","noise"]
        seed = ""
        while not seed:
            print "Initial seeding choices: " + ", ".join(seeds)
            seed = raw_input("initial seeding choice: ").rstrip().lower()
            if not seed in seeds:
                print "please select seed condition from list"
                seed = ""

        params["seed"] = seed

    return params,modelFunc


if __name__ == "__main__":
    #Parses the command line arguments
    parser = argparse.ArgumentParser(description="Gray-Scott simulation")
    parser.add_argument("-o", "--outname", help="simulation run output name",default="simulation_output")
    parser.add_argument("-mov", "--moviemode", action='store_true',help="Run the script in \"movie mode\"\
        to store sequential images (default off)",default=False)
    parser.add_argument("-n", "--timesteps", type=int, help="Number of timesteps to simulate (default 8000)",default=8000)
    parser.add_argument("-m","--model", choices=['FN','GM','GS'], help="Model simulation choice FN [FitzHugh-Nagumo]\
        \nGM [Gierer-Meinhardt]\nGS [Gray-Scott].")
    args = parser.parse_args()

    #Set simulation grid size
    size = 200

    model = args.model

    #Get input interactively if no command line args set
    while not model:
        print "--- Model choices --- \nFN [FitzHugh-Nagumo]\nGM [Gierer-Meinhardt]\nGS [Gray-Scott]"
        model = raw_input("(choose one): ").rstrip().upper()
        if model not in ["FN","GM","GS"]:
            print "please enter two-letter model name from list\n"
            model = ""


    #get params for model
    params,modelFunc = setModelParams(model)

    #create output dir
    runName = args.outname
    savPath = args.outname + "_images/"
    if not os.path.exists(savPath):
        os.makedirs(savPath)

    #set up image saving
    totFrames = 250
    movieOutput = False
    n = args.timesteps
    print "Running simulation with " + str(n) + " timesteps."
    frameMod = n/totFrames

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
    if params["seed"] == "single":
        r = 20
        U[size/2-r:size/2+r,size/2-r:size/2+r] = 0.50
        V[size/2-r:size/2+r,size/2-r:size/2+r] = 0.25
    elif params["seed"] == "dual":
        r = 15
        U[size/4-r:size/4+r,size/4-r:size/4+r] = 0.50
        V[size/4-r:size/4+r,size/4-r:size/4+r] = 0.25
        U[3*size/4-r:3*size/4+r,3*size/4-r:3*size/4+r] = 0.50
        V[3*size/4-r:3*size/4+r,3*size/4-r:3*size/4+r] = 0.25

    else: #seed with random noise
        u-=1
        u+=np.random.rand(len(u),len(u))
        v+=np.random.rand(len(u),len(u))

    #add small amount of noise everywhere
    #u += (0.01 + 0.01*(np.random.random((size-2,size-2))*2-1))
    #v += (0.01 + 0.01*(np.random.random((size-2,size-2))*2-1))

    initial_matrices = (U,V)

    #RUN SIM
    makeImg(v,"initial_v",params["myCmap"],setEdge=params["edgeMax"])
    u,v = modelFunc(params,initial_matrices)
    makeImg(v,"final_v",params["myCmap"],setEdge=params["edgeMax"])
