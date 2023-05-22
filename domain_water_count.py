"""
Count the number of water molecules in feed and permeate domain
"""

import os
import sys
import numpy as np
import pandas as pd
import re
import math
import statistics
import argparse

class DomainWaterCount(object):
    def __init__(self):
        self.parser = argparse.ArgumentParser(prog = 'domain_water_count', usage = '%(prog)s [-h for help]', \
                                          description = 'Find polymer domain')
        self.parser.add_argument('-gro_file', "--gro_file", help = 'Input gro file (Required).')
        self.parser.add_argument('-o',"--o")
        self.parser.add_argument('-n_atoms_water',"--n_atoms_water", default=0, help='Number of atoms per water molecule')
        self.parser.add_argument('-dzp',"--dzp",help="Distance of feed inflection point from COG ")
        self.parser.add_argument('-dzn',"--dzn", help="Diatance of permeate inflection point from COG")
        self.parser.add_argument('-polymer_resnames',"--polymer_resnames")
    def main(self,bash=True,**kwargs):
        if bash==False:
            self.__dict__.update(kwargs)
        else:
            args=self.parser.parse_args()
            param=vars(args)
            for key in param:
                setattr(self,key,param[key])
        return self.domain_water_count()

    def domain_water_count(self):
        try:
            gro_file=open(os.path.join(self.current_path,self.gro_file))

        except IOError as e:
            print("Unable to open" + gro_file +". Please check the file")
    
        #Load atom coordinates in pandas dataframe
        coord=[]
        box=[] 
        n_atoms=0
        for idx,line in enumerate(f):
            if idx==1:
                n_atoms=int(line.split()[0])
            if idx>=2:
                parts=line.split()
                if len(parts)==6:
                    coord.append([parts[0],parts[1],int(parts[2]),float(parts[3]),float(parts[4]), float(parts[5])])
                elif len(parts)==5:
                     coord.append([parts[0],parts[1][0:-5],int(parts[1][-5:]),float(parts[2]),float(parts[3]), float(parts[4])])
                elif len(parts)==3:
                    box.append([parts[0],parts[1], parts[2]])

        box_x=float(box[0][0])
        box_y=float(box[0][1])
        box_z=float(box[0][2])
        labels=["resid","atom_type","index","x","y","z"]
        df=pd.DataFrame(coord,columns=labels)
        f.close()

        #Optional statement to print out complete dataframe
        #pd.set_option('display.max_rows',df.shape[0]+1)

        #print("Read the gro file\n",df)
        polymer_resnames="|".join(re.split)

        #Find the maximum and minimum z of polymer domain
        polymer=df.loc[df['resid'].str.contains(polymer_resnames, case=False)]
        #max_z_polymer=polymer["z"].max()
        #min_z_polymer=polymer["z"].min()

        max_z_polymer=polymer["z"].mean()+self.dzp
        min_z_polymer=polymer["z"].mean()-self.dzn
        

        #Finding the number of water molecules in the feed side and permeate side
        sol_frame=df.loc[df['resid'].str.contains("SOL",case=False)]
        keep_O_ind_feed=sol_frame.loc[((sol_frame['z']>max_z_polymer) & (sol_frame['atom_type']=='OW'))]
        keep_O_ind_permeate=sol_frame.loc[((sol_frame['z']<min_z_polymer) & (sol_frame['atom_type']=='OW'))]
        df_O_feed=sol_frame[sol_frame.index.isin(keep_O_ind_feed.index)]
        df_O_permeate=sol_frame[sol_frame.index.isin(keep_O_ind_permeate.index)]

        #print(df_O_feed)
        #print(df_O_permeate)

              
        o=open(self.o,"w+")
        o.write(";Water moleclues in the feed\n")
        o.write("%d\n"%(df_O_feed.shape[0]))
        o.write(";Water moleclues in the permeate\n")
        o.write("%d\n"%(df_O_permeate.shape[0]))
        #o.write(";Total water count\n")
        #o.write("%d\n"%(sol_frame.shape[0])/float(self.n_atoms_water))
        
        o.close()
        

        return df_O_feed.shape[0],df_O_permeate.shape[0],

    
if __name__ == "__main__":
    try:
        run=DomainWaterCount()
        run.main()
    except (RuntimeError, TypeError, NameError) as e:
       print(e)