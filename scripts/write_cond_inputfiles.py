#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Feb 14 14:03:42 2022
modified 20221021 & 20240816
@author: dblaschke

PFDD -- Phase Field Dislocation Dynamics

© 2022. Triad National Security, LLC. All rights reserved.
This program was produced under U.S. Government contract 89233218CNA000001 for
Los Alamos National Laboratory (LANL), which is operated by Triad National
Security, LLC for the U.S. Department of Energy/National Nuclear Security
Administration. All rights in the program are reserved by Triad National
Security, LLC, and the U.S. Department of Energy/National Nuclear Security
Administration. The Government is granted for itself and others acting on its
behalf a nonexclusive, paid-up, irrevocable worldwide license in this material 
to reproduce, prepare derivative works, distribute copies to the public, perform
 publicly and display publicly, and to permit others to do so.
"""

from dataclasses import dataclass
import os
import sys
import time
# from shutil import copyfile
import numpy as np

lattice_resolution = 128

@dataclass
class material:
    '''this data class stores various material information, such as material name, layer thickness, number of grains across this thickness, dislocation density etc.'''
    name : str
    dlayer : float = 100e-9
    Ngrains: int = 1
    rho_disloc: float = 1e15
    Rgrain: float = -1 ## -1==take value from material.cpp
    
    

def write_inputfile(fname,material1,material2,temperature=300,El=(0.,1.0e-12,0.),Nlayers=1,pinterface=0.5):
    '''as its name indicates, this function writes an input file for the version of our pfdd c++ code that includes conductivity calculations'''
    # if not (isinstance(material1, material) and isinstance(material2, material)):
    #     raise ValueError("material1 and material2 must be instances of the material class")
    with open(fname,'w', encoding="utf8") as infile:
        infile.write("# PFDD, conductivity\nseed             576\napp_style        cond100B core_energy 1\n\n")
        infile.write("fft_style       fftw_slab pbc 1 1 1 mode 1 primitive x 1 0 0 y 0 1 0 z 0 0 1 # cubic\n")
        infile.write("dimension        3\nlattice         sc/6n 1\n")
        infile.write(f"region          box block 0 {lattice_resolution} 0 {lattice_resolution} 0 {lattice_resolution}\n")
        infile.write("create_box       box\ncreate_sites     box\n")
        ## material 1 also contains El and pinterface:
        infile.write(f"material       preset {material1.name} dlayer {material1.dlayer} El {El[0]} {El[1]} {El[2]} ")
        if material1.Rgrain != -1: infile.write(f"Rgrain {material1.Rgrain} ")
        if pinterface != 0.5: infile.write(f"pinterface {pinterface} ")
        infile.write(f"Nlayers {Nlayers} Ngrains {material1.Ngrains} rho_disloc {material1.rho_disloc}\n")
        ## material 2:
        infile.write(f"material       preset {material2.name} dlayer {material2.dlayer} ")
        if material2.Rgrain != -1: infile.write(f"Rgrain {material2.Rgrain} ")
        infile.write(f"Ngrains {material2.Ngrains} rho_disloc {material2.rho_disloc}\n")
        ##
        infile.write("solve_style      GLCond max_iter 300 tol 1e-4\n\n")
        infile.write("diag_style      strain\ndiag_style      stress\n\n")
        infile.write(f"temperature	     {temperature}\n")
        infile.write("scale           1.0 1.0 1.0\n")
        infile.write("sigma       initial 0.0 0.0 0.0 0.0 0.06 0.0 delta 0.0 0.0 0.0 0.0 0.0 0.0\n\n")
        infile.write("stats		1.0e-8\ndump            1 text 50 PFDD-composite.dump x y z xir1\nrun              1\n")


def volumefraction(dlayer1,dlayer2):
    '''calculates the colume fraction of material 2 given the layer thicknesses of materials 1 and 2'''
    return dlayer2/(dlayer1+dlayer2)


def approximateVf(target_Vf,ref_dlayer,lattice_resolution=lattice_resolution):
    '''To avoid quantization errors, this function calculates a volume fraction Vf which is as close as possible to
       target_Vf, but which can be achieved exactly given the lattice resolution.
       We return a tuple containing said Vf as well as dlayer for the second material.'''
    gridpts = int(round(target_Vf*lattice_resolution))
    Vf = gridpts/lattice_resolution
    #Vf = (materialB->dlayer) / ((materialA->dlayer) + (materialB->dlayer))
    dlayer = Vf*ref_dlayer / (1- Vf)
    return (Vf,dlayer)


def chdir(targetfolder):
    '''same as os.chdir, except 'targetfolder' is created first unless it already exists'''
    if not os.path.exists(targetfolder):
        os.mkdir(targetfolder)
    os.chdir(targetfolder)

## define command line options:
OPTIONS = {"threshold2g":int, "exefile":str, "Ncores":int, "dlayer":str,"temperatures":str,"Vf_target":str,"lattice_resolution":int,"pinterface":float}
## fct parse_options: copied from PyDislocDyn, see https://github.com/dblaschke-LANL/PyDislocDyn (bsd license)
def parse_options(arglist,optionlist=OPTIONS,globaldict=globals()):
    '''Search commandline arguments passed to this script for known options to set by comaring to a list of keyword strings "optionlist".
    These will then override default varibles set above in this script. This function also returns a copy of 'arglist' stripped of all 
    option calls for further processing (e.g. opening input files that were passed etc.).'''
    out = arglist
    setoptions = [i for i in out if "--" in i and "=" in i and i[:2]=="--"]
    for i in setoptions:
        out.remove(i)
        key,val = i[2:].split("=")
        if key in optionlist:
            globaldict[key] = optionlist[key](val)
            print(f"setting {key}={globaldict[key]}")
    time.sleep(1) ## avoid race conditions after changing global variables
    return out

if __name__ == '__main__':
    mater1 = 'Cu'
    mater2 = 'Nb'
    Ncores = 8
    pinterface=0.5
    exefile = "../src/pfdd_openmpi"
    threshold2g = 100 # [nm] *e-9
    ## layer thicknesses of material 1 (Cu by default)
    dlayer = '25,300,12' #np.linspace(25,300,12)*1e-9
    temperatures = '100,450,8' #np.linspace(100,450,8)
    Vf_target = '50,10,5' #np.linspace(0.5,0.1,5)
    
    if len(sys.argv) >= 3:
        args = parse_options(sys.argv[1:])
        mater1 = args[0]
        mater2 = args[1]
    else:
        print(f"Usage: {sys.argv[0]} <material1> <material2>\n"
            "each material must be one of three formats: 'material'\n"
            "or 'material,rhodisloc' or 'material,rhodisloc,Rgrain'\n"
            "additional options may be passed via '--keyword=value' anywhere on the commandline.\n"
            f"available options are: {OPTIONS}\n"
            "options 'dlayer', 'temperature', or 'Vf_target' require three comma-separated values\n"
            "to be passed to np.linspace().\n"
            "aborting.")
        sys.exit()
    
    if (lenmat:=len(mater1.split(",")))>1:
        mater1 = mater1.split(",")
        if lenmat>2:
            Rgrain1 = mater1[2]
        else:
            Rgrain1 = -1
        rhodisloc1 = mater1[1]
        mater1 = mater1[0]
    else:
        rhodisloc1 = 1e15
        Rgrain1 = -1
    
    if (lenmat:=len(mater2.split(",")))>1:
        mater2 = mater2.split(",")
        if lenmat>2:
            Rgrain2 = mater2[2]
        else:
            Rgrain2 = -1
        rhodisloc2 = mater2[1]
        mater2 = mater2[0]
    else:
        rhodisloc2 = 1e15
        Rgrain2 = -1
    
    title = f"{mater1}/{mater2} composite"
    print(f"generating input files for {title}")
    # print(f"setting exec = {exefile}")
    # print(f"using {Ncores=}")
    threshold2g = threshold2g*1e-9
    dlayer = np.asarray(dlayer.split(","),dtype=float)
    dlayer = np.linspace(dlayer[0],dlayer[1],int(dlayer[2]))*1e-9
    temperatures = np.linspace(*np.asarray(temperatures.split(","),dtype=int))
    Vf_target = np.linspace(*np.asarray(Vf_target.split(","),dtype=int))/100
    Vf_and_dl = [] ## contains tuples of volume fraction and layer thickness relative to ref_layer thickness
    for vftar in Vf_target:
        Vf_and_dl.append(approximateVf(vftar, 1))
    
    cwd = os.getcwd()
    njobs = 0
    with open(f"submitjobs_{mater1}{mater2}","w", encoding="utf8") as jobfile:
        jobfile.write("#!/bin/sh\n")
        jobfile.write(f"\necho \"{title}\"\n\n")
        jobfile.write('if [ "$1" != "" ]; then')
        jobfile.write("\n Ncores=$1 \nelse \n ")
        jobfile.write(f"{Ncores=}\nfi\n\n")
        jobfile.write(f'{exefile=}')
        for T in temperatures:
            Tfolder = f"T{T:.0f}"
            chdir(Tfolder)
            for iV,Vf in enumerate(Vf_and_dl):
                Vfolder = f"Vf{Vf[0]:.4f}"
                chdir(Vfolder)
                for i,dl in enumerate(dlayer):
                    strdl = f"{1e9*dl:.0f}"
                    fname = f"in{mater1}{mater2}_{strdl}"
                    fpath = f"{Tfolder}/{Vfolder}/"
                    M1 = material(mater1,dl,rho_disloc=rhodisloc1,Rgrain=Rgrain1)
                    M2 = material(mater2,Vf[1]*dl,rho_disloc=rhodisloc2,Rgrain=Rgrain2)
                    write_inputfile(fname, M1, M2,temperature=T,pinterface=pinterface)
                    jobfile.write(f"\necho \"calculating {fpath}{fname}\"\n")
                    jobfile.write(f"mpirun -np $Ncores $exefile -in {fpath}{fname} -log {fpath}{fname[2:]}.log | tee {fpath}{fname[2:]}.txt\n")
                    njobs += 1
                    if dl>=threshold2g:
                        chdir("2g")
                        fname = f"in{mater1}{mater2}_{strdl}_g2"
                        fpath += "2g/"
                        M1 = material(mater1,dl,Ngrains=2,rho_disloc=rhodisloc1,Rgrain=Rgrain1)
                        if Vf[1]*dl>=threshold2g:
                            M2 = material(mater2,Vf[1]*dl,Ngrains=2,rho_disloc=rhodisloc2,Rgrain=Rgrain2) # otherwise use previous M2 with only 1 grain
                        write_inputfile(fname, M1, M2,temperature=T,pinterface=pinterface)
                        jobfile.write(f"\necho \"calculating {fpath}{fname}\"\n")
                        jobfile.write(f"mpirun -np $Ncores $exefile -in {fpath}{fname} -log {fpath}{fname[2:]}.log | tee {fpath}{fname[2:]}.txt\n")
                        njobs += 1
                        os.chdir("..")
                os.chdir("..")
            
            os.chdir(cwd)
        jobfile.write(f"\necho \"completed {njobs} jobs\"\n")
    
    with open(f"{mater1}{mater2}_composite.dat","w", encoding="utf8") as logfile:
        logfile.write(f"#{title}\n")
        logfile.write(f"material1 = {mater1}\nmaterial2 = {mater2}\n")
        logfile.write(f"rhodisloc1 = {M1.rho_disloc}\nrhodisloc2 = {M2.rho_disloc}\n")
        logfile.write(f"Rgrain1 = {M1.Rgrain}\nRgrain2 = {M2.Rgrain}\n")
        logfile.write(f"{pinterface=}\n")
        logfile.write(f"dlayer = {(dlayer[0]*1e9):0f} {(dlayer[-1]*1e9):0f} {len(dlayer)}\n")
        logfile.write(f"threshold2g = {threshold2g*1e9:.0f}\n")
        logfile.write(f"temperatures = {temperatures[0]} {temperatures[-1]} {len(temperatures)}\n")
        logfile.write(f"Vf_target = {Vf_target[0]} {Vf_target[-1]} {len(Vf_target)}\n")
        logfile.write(f"lattice_resolution = {lattice_resolution}\n")
        logfile.write("\n##postprocessing parameters may be added/adapted below: \n")
        logfile.write("#stddev = 40\n#gperc = 0.1\n")
    