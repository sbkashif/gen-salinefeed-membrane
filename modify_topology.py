"""
Count the number of water molecules in a domain -- feed, polymer, or permeate
"""

import os
import sys
import numpy as np
import pandas as pd
import re
import math
import statistics
import argparse

class ModifyTopology(object):
    def __init__(self):
        self.parser = argparse.ArgumentParser(prog = 'modify_topology', usage = '%(prog)s [-h for help]', \
                                          description = 'Find polymer domain')
        self.parser.add_argument('-top_file', "--top_file", help = 'Input xtc file (Required).')
        self.parser.add_argument('-o',"--o")
        self.parser.add_argument('-gro',"--gro", help='reference gro file')
        self.parser.add_argument('-target_resname','--target_resname',default="SOL",help="residue name of water")
        self.current_path=os.getcwd()
        
        
    def main(self,bash=True,**kwargs):
        if bash==False:
            self.__dict__.update(kwargs)
        else:
            args=self.parser.parse_args()
            param=vars(args)
            for key in param:
                setattr(self,key,param[key])
        print("Running",sys.argv)
        return self.modify_topology()

    def modify_topology(self):
        try:
            print("Opening gro file...\n ",self.gro)
            gro_file=open(os.path.join(self.current_path,self.gro))

        except IOError as e:
            print("Unable to open" + gro_file +". Please check the file")
    
        #Load atom coordinates in pandas dataframe
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
        res_split = df['res'].str.split('(\d+)([A-Za-z]+)',n=1,expand=True)
        res_split.rename(columns={2:'resname_string', 3:'resname_numeric'}, inplace=True)
        res_split['resname']=res_split['resname_string']+res_split['resname_numeric']
        df.insert(loc=1,column='resid',value=res_split[1])
        df['resid']=df['resid'].astype(int)
        df.insert(loc=2,column='resname',value=res_split['resname'])
        df.drop('res',axis=1,inplace=True)
        df=df.reset_index(drop=True)
        gro_file.close()

        #Optional statement to print out complete dataframe
        #pd.set_option('display.max_rows',df.shape[0]+1)

        try:
            print("Opening topology file...\n ",self.top_file)
            top_file=open(os.path.join(self.current_path,self.top_file))
            list_of_lines=top_file.readlines()

        except IOError as e:
            print("Unable to open" + top_file +". Please check the file")

        target_resname_lines=[]
        molecule_line=False
        for idx,line in enumerate(list_of_lines):
            if re.search("molecules",line):
                molecule_line=True
            if re.search(self.target_resname,line) and molecule_line==True:
                target_resname_lines.append(idx)
        if not target_resname_lines:
            target_resname_lines.append(len(list_of_lines))
            list_of_lines.append("")
        line_num=target_resname_lines[0]
        target_resname_lines.pop(0)

        #Finding the number of water molecules in the gro file
        sol_frame=df.loc[df['resname'].str.contains(self.target_resname,case=False)]
        n_target_molecules=len(pd.unique(sol_frame['resid']))
        #n_target_molecules=int(sol_frame.shape[0]/int(self.n_atoms_water))
        print("Number of molecules of the chosen residue:",n_target_molecules)
        
        list_of_lines[line_num]=str(self.target_resname)+"\t\t\t\t"+str(n_target_molecules)+"\n"
        top_file.close()      
        
        print("Writing new topology file:", self.o)
        o=open(os.path.join(self.current_path,self.o),"w+")
        for idx,line in enumerate(list_of_lines):
            if idx not in target_resname_lines:
                o.write(line)

        #o.writelines(list_of_lines)
        
        o.close()
        

        return None

    
if __name__ == "__main__":
    try:
        run=ModifyTopology()
        run.main()
    except (RuntimeError, TypeError, NameError) as e:
       print(e)
