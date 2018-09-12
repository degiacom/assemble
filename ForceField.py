# Copyright (c) 2014-2018 Matteo Degiacomi and Valentina Erastova
#
# Assemble is free software ;
# you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation ;
# either version 2 of the License, or (at your option) any later version.
# Assemble is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY ;
# without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
# See the GNU General Public License for more details.
# You should have received a copy of the GNU General Public License along with Assemble ;
# if not, write to the Free Software Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307 USA.
#
# Authors : Matteo Degiacomi, matteo.degiacomi@gmail.com, Valentina Erastova, valentina.erastova@gmail.com


import numpy as np
import logging

class ForceField(object):

    def __init__(self):
        self.bonded={}
        self.nonbonded={}
        self.combination=[]
        self.fftype=[]

        #default return values
        self.default_bond=1.5
        self.default_angle=114
        self.default_dihedral=120
        
        self.logger=logging.getLogger('assemble')
        
        
    def load(self,fffile):
        
        self.logger.info("\n> loading force field %s..."%fffile)
        f = open(fffile, 'r+')
        line = f.readline()
        
        while line:
            w=line.split()
 
            if "bondedtypes" in line:
                break

            line=f.readline()
        
        #extract equations type
        while line:
            w=line.split()
            if not line.isspace() and not ";" in w[0] and len(w)==4:
                self.fftype=np.array(w).astype(int)
                break
            line=f.readline()

        line=f.readline()

        #extract bonded potential constants 
        while line:
            w=line.split()
            
            if "atomtypes" in line:
                break
            
            if not line.isspace() and not ";" in w[0]:
                self.bonded[w[0]]=np.array(w[1:]).astype(float)
                        
            line=f.readline()

        line=f.readline()

        #non bonded potential
        while line:         
            w=line.split()

            if "defaults" in line:
                break
            
            if not line.isspace() and not ";" in w[0]:
                self.nonbonded[w[0]]=np.array(w[1:])

            
            line=f.readline()

        line=f.readline()

        #get combination rules
        while line:         
            w=line.split()
            if not line.isspace() and not ";" in w[0]:
                self.combination=np.array(w)
                break
        
            line=f.readline()
        
        f.close()
 
        if len(self.fftype)==0:
            raise IOError("bond types not found in force field %s!"%fffile)

        if len(self.bonded)==0:
            raise IOError("bonded parameters not found in force field %s!"%fffile)
        
        if len(self.nonbonded)==0:
            raise IOError("non-bonded parameters not found in force field %s!"%fffile)

        if len(self.combination)==0:
            raise IOError("combination rules not found in force field %s!"%fffile)


    def get_bond(self,name):
        if self.fftype[0]>=1 and self.fftype[0]<=7:
            return self.bonded[name][0]*10
        #tabulated potential, no way to know where the minimum is. Return default value
        else:
            return self.default_bond
    
    def get_angle(self,name):
        if self.fftype[1]>=1 and self.fftype[1]<=2:
            return self.bonded[name][0]
        #no analytical minimum exists, return default
        else:
            return self.default_angle
    
    def get_dihedral(self,name):
        if self.fftype[2]>=1 and self.fftype[2]<=2:
            return self.bonded[name][0]
        #no analytical minimum exists, return default
        else:
            return self.default_dihedral
        

if __name__=="__main__":
    
    FF=ForceField()
    FF.load("./database/forcefield/trappe.ff.txt")
    
    #FF.load("C:\Users\Matteo\workspace\polymer\database\forcefield\trappe.ff")
    
    print(FF.bonded)
    #print FF.nonbonded
    #print FF.combination
    print(FF.fftype)
