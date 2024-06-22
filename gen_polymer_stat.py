"""
Find the minimum and maximum z-coordinate of the polymer
"""

import os
import sys
import numpy as np
import pandas as pd
import re
import MDAnalysis as mda
import math
import statistics
import argparse

def create_parser():
    parser = argparse.ArgumentParser(prog = 'HeatMapUsingMDAnalysis', 
                                     usage = '%(prog)s [-h for help]', \
                                     description = 'Generate the minimym and maximum z-coordinate of the polymer.')
    parser.add_argument('-f', "--f", help = 'Input xtc file (Required).')
    parser.add_argument('-o',"--o")
    return parser


def main(args):
    GROFile=args.f
    OutFile=args.o
    det_polymer_extremes(GROFile,OutFile)

def det_polymer_extremes(filename,outfile):
    """
    Determine the minimum and maximum z-coordinate of the polymer

    
    Parameters
    ----------
    filename : str
        Inpit file name
    outfile : str
        Output file name

    Returns
    -------
    Tuple
        Minimum and maximum z-coordinate of the polymer
    """    
    
    try:
        print("Opening file...\n ",filename)
        f=open(filename)

    except IOError as e:
        print("Unable to open" + filename +". Please check the file")
   
    #Load coordinates in pandas dataframe
    coord=[]
    box=[] 
    n_atoms=0
    for idx,line in enumerate(f):
        if idx==1:
            n_atoms=int(line.split()[0])
        if idx>=2:
            parts=line.split()
            if len(parts)==6:
                coord.append([parts[0],parts[1],parts[2],float(parts[3]),float(parts[4]), float(parts[5])])
            elif len(parts)==5:
                 coord.append([parts[0],parts[1][0:-5],parts[1][-5:],float(parts[2]),float(parts[3]), float(parts[4])])
            elif len(parts)==3:
                box.append([parts[0],parts[1], parts[2]])

    box_x=float(box[0][0])
    box_y=float(box[0][1])
    print("Box dimensions:",box_x,box_y)
    labels=["resid","atom_type","index","x","y","z"]
    df=pd.DataFrame(coord,columns=labels)
    df.astype({'x': 'float64'}).dtypes
    df.astype({'y': 'float64'}).dtypes 
    df.astype({'z': 'float64'}).dtypes
    pd.set_option('display.max_rows',df.shape[0]+1)

    
    #Find the maximum of polymer atoms

    polymer=df.loc[df['resid'].str.contains("MPD|TMC", case=False)]
    max_polymer=polymer["z"].max()
    min_polymer=polymer["z"].min()
    print(polymer["z"].idxmin())
    print(min_polymer,max_polymer)

    return min_polymer,max_polymer

   
if __name__ == "__main__":
    args = create_parser().parse_args()
    main(args)


