"""
Retain polymer and water molecules within a certain distance from polymer
"""

import os
import sys
import numpy as np
import pandas as pd
import re
import math
import statistics
import argparse

class PolymerDomain(object):
    def __init__(self):
        self.parser = argparse.ArgumentParser(prog = 'obtain_polymer_domain', usage = '%(prog)s [-h for help]', \
                                          description = 'Find polymer domain')
        self.parser.add_argument('-gro_file', "--gro_file")
        self.parser.add_argument('-o',"--o")
        self.parser.add_argument('-dzp',"--dzp", default=0.0, help='z-distance from the top of polymer to retain solvent')
        self.parser.add_argument('-dzn',"--dzn", default=0.0, help='z-distance from the bottom of polymer to retain solvent')
        self.parser.add_argument('-polymer_resnames',"--polymer_resnames", default="SOL", help='z-distance from the bottom of polymer to retain solvent')
        self.parser.add_argument('-solvent_resnames',"--solvent_resnames", help='z-distance from the bottom of polymer to retain solvent')
        self.parser.add_argument('-n_atoms_water','--n_atoms_water',default=0,help='number of atoms per water molecules')
        self.current_path=os.getcwd()

    def main(self,bash=True,**kwargs):
        if bash==False:
            self.__dict__.update(kwargs)
        else:
            args=self.parser.parse_args()
            param=vars(args)
            for key in param:
                setattr(self,key,param[key])
        print(sys.argv)
        print(self.current_path)
        return self.obtain_polymer_domain()

    def obtain_polymer_domain(self):
        try:
            gro_file=open(self.gro_file)
            print(os.path.abspath(gro_file))
        except IOError as e:
            print("Unable to open" + self.gro_file +". Please check the file")
         
        try:
            #Load atom coordinates from gro file to pandas dataframe
            coord=[]
            box=[] 
            n_atoms=0
            for idx,line in enumerate(gro_file):
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
            labels=["res","atom_type","index","x","y","z"]
            df=pd.DataFrame(coord,columns=labels)
            gro_file.close()

            #Optional statement to print out complete dataframe
            #pd.set_option('display.max_rows',df.shape[0]+1)

            #Find the maximum and minimum z of polymer
            polymer_resnames='|'.join(self.polymer_resnames.split())
            polymer=df.loc[df['res'].str.contains(polymer_resnames, case=False)]
            max_polymer=polymer["z"].max()
            min_polymer=polymer["z"].min()
            max_z_polymer=max_polymer+float(self.dzp)
            min_z_polymer=min_polymer-float(self.dzn)

            solvent_resnames='|'.join(self.solvent_resnames.split())
            #Finding the water molecules in polymer domain
            sol_frame=df.loc[df['res'].str.contains(solvent_resnames,case=False)]
            p_strip=sol_frame[(sol_frame["z"] < max_z_polymer) & (sol_frame["z"] > max_polymer) & (sol_frame["atom_type"]=="OW")].shape[0]
            n_strip=sol_frame[(sol_frame["z"] > min_z_polymer) & (sol_frame["z"] < min_polymer) & (sol_frame["atom_type"]=="OW")].shape[0]
            keep_O_ind=sol_frame.loc[((sol_frame['z']<max_z_polymer) & (sol_frame['z']>min_z_polymer) & (sol_frame['atom_type']=='OW'))]
            df_O=sol_frame[sol_frame.index.isin(keep_O_ind.index)]
            polymer_water_molecules=df_O[(df_O["z"]>=min_polymer)& (df_O["z"]<=max_polymer)].shape[0]
            df_H1=sol_frame[sol_frame.index.isin(keep_O_ind.index+1)]
            df_H2=sol_frame[sol_frame.index.isin(keep_O_ind.index+2)]
            df_D=sol_frame[sol_frame.index.isin(keep_O_ind.index+3)]
            atom_dataframes=[df_O,df_H1,df_H2,df_D]
            sol_frame_domain=pd.concat(atom_dataframes).sort_index(kind='merge')
            

            #Combining the two dataframes
            polymer_sol_domain_list=[polymer,sol_frame_domain]
            polymer_sol_domain=pd.concat(polymer_sol_domain_list)

            res_split = polymer_sol_domain['res'].str.split('(\d+)([A-Za-z]+)',n=1,expand=True)
            res_split.rename(columns={2:'resname_string', 3:'resname_numeric'}, inplace=True)
            res_split['resname']=res_split['resname_string']+res_split['resname_numeric']
            polymer_sol_domain.insert(loc=1,column='resid',value=res_split[1])
            polymer_sol_domain['resid']=polymer_sol_domain['resid'].astype(int)
            polymer_sol_domain.insert(loc=2,column='resname',value=res_split['resname'])
            polymer_sol_domain.drop('res',axis=1,inplace=True)
            polymer_sol_domain=polymer_sol_domain.reset_index(drop=True)
            #print("Extracted the atoms present in the domain:\n",polymer_sol_domain)
            #n_water_molecules=polymer_sol_domain.loc[df['team1'].str[0]=='S', 'team1'].unique()
            #_water_molecules=len(pd.unique(sol_frame['resid']))
            
            self.n_water_molecules=sol_frame_domain.shape[0]/int(self.n_atoms_water)
        except Exception as e:
            print(e)
        try:     
            o=open(self.o,"w+")
        except IOError as e:
            print(e)

        o.write("Generated by obtain_polymer_domain.py")
        o.write("\n")
        o.write("%d"%(polymer_sol_domain.shape[0]))
        o.write("\n")
        print("Writing the gro file...")
        #Writing to gro file
        for i in range(0,len(polymer_sol_domain)):
            #o.write("%5d%-5s%5s%5d%8.3f%8.3f%8.3f\n" % (row["resid"],row["resname"],row["atom_type"],row["index"],row["x"],row["y"],row["z"]))
            o.write("%5d%-5s%5s%5d%8.3f%8.3f%8.3f\n" % (polymer_sol_domain.loc[i,"resid"],polymer_sol_domain.loc[i,"resname"],polymer_sol_domain.loc[i,"atom_type"],polymer_sol_domain.loc[i,"index"],polymer_sol_domain.loc[i,"x"],polymer_sol_domain.loc[i,"y"],polymer_sol_domain.loc[i,"z"]))
        o.write("%3s%8.3f%8.3f%8.3f"%(" ",box_x,box_y,box_z))
        o.write("\n")
        #Writing the minimum coordinates
        dir_name=os.path.dirname(os.path.abspath(self.o))    
        log=open(os.path.join(dir_name,"polymer_info.txt"),"w+")
        log.write(";Minimum polymer coordinate\n")
        log.write("%f"%(min_polymer))
        log.write("\n;Maximum polymer coordinate\n")
        log.write("%f"%(max_polymer))
        log.write("\n;Number of water molecules after cutting\n")
        log.write("%d"%(self.n_water_molecules))
        log.write("\n;Number of water molecules in feed strip\n")
        log.write("%d"%(p_strip))
        log.write("\n;Number of water molecules permeate strip\n")
        log.write("%d"%(n_strip))
        log.write("\n;Number of polymer water molecules\n")
        log.write("%d"%(polymer_water_molecules))
        
        output={"Box x(nm)":box_x,"Box y(nm)":box_y,"Min polymer coordinate":min_polymer,"Max polymer coordinate":max_polymer,"Domain water molecules":self.n_water_molecules,"Feed strip water molecules":p_strip,"Permeate strip water molecules":n_strip,"Polymer water molecules":polymer_water_molecules}
        return output

        
if __name__ == "__main__":
    try:
        run=PolymerDomain()
        run.main()
    except (RuntimeError, TypeError, NameError) as e:
       print(e)