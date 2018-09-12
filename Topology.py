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

class Topology(object):

    def __init__(self):
        
        self.head=-1
        self.tail=-1
    
    def load(self,topfile):
        
        b=[]
        a=[]
        d=[]
        i=[]
        m=[]
        self.replace_cter=[] #cter substitutions
        self.replace_nter=[] #nter substitutions

        f = open(topfile, 'r+')
        
        target=""
        
        for line in f:
            w=line.split()
            if len(w)>0 and w[0]=="[" :
                target=w[1]

            elif len(w)>0 and ";" in w[0]:
                continue

            elif len(w)>0:
                if target=="bonds":
                    b.append(w)
                elif target=="angles":
                    a.append(w)
                elif target =="dihedrals" :       
                    d.append(w)
                elif target == "impropers":
                    i.append(w)
                elif target =="mapping":
                    m.append(w)
                elif target == "cterminal":
                    if len(w)==2:
                        self.tail=w
                    self.replace_cter.append(w)
                elif target == "nterminal":
                    if len(w)==2:
                        self.head=w
                    self.replace_nter.append(w)                
                else:
                    raise IOError("ERROR: unknown %s section in topology file %s"%(target, topfile))
                  
            
        f.close()
        
        self.bonds=np.array(b).astype(str)
        self.angles=np.array(a).astype(str)
        self.dihedrals=np.array(d).astype(str)
        self.impropers=np.array(i).astype(str)
        self.mapping=np.array(m).astype(str)
        #self.replace_nter=np.array(sn).astype("S10") #nter replacements
        #self.replace_cter=np.array(sc).astype("S10") #cter replacements

        if self.head==-1:
            raise IOError("head not found!")
        if self.tail==-1:
            raise IOError("tail not found!")


    #substitute bonds, angles and dihedrals according to terminal section
    #( call once to transform molecule into terminal molecule)
    # WARNING: number of chars of tailing element should not be greater than the greatest length in topology
    def make_terminal(self,target):
    
        if target=="nterminal":
            subst=self.replace_nter
        elif target=="cterminal":
            subst=self.replace_cter
        else:
            raise IOError("target should be equal to nterminal or cterminal")
    
        for i in range(0,len(subst),1):

            n=subst[i]

            if len(n)==2: #atom substitution
                test=self.mapping[:,0]==n[0]
                if np.any(test):
                    self.mapping[test,1]=n[-1]

            elif len(n)==3 and len(self.bonds)>0: #bond substitution
                test=np.logical_and(self.bonds[:,0]==n[0],self.bonds[:,1]==n[1])
                if np.any(test):
                    self.bonds[test,2]=n[-1]
                
            elif len(n)==4 and len(self.angles)>0: #angle substitution
                test=np.logical_and(np.logical_and(self.angles[:,0]==n[0],self.angles[:,1]==n[1]),self.angles[:,2]==n[2])
                if np.any(test): 
                    self.angles[test,3]=n[-1]

            elif len(n)==5: #dihedral or improper substitution
                if len(self.dihedrals)>0:
                    test=np.logical_and(np.logical_and(np.logical_and(self.dihedrals[:,0]==n[0],self.dihedrals[:,1]==n[1]),\
                                                self.dihedrals[:,2]==n[2]),self.dihedrals[:,3]==n[3])
                    if np.any(test):
                        self.dihedrals[test,4]=n[-1]
                        
                if len(self.impropers)>0:
                    test=np.logical_and(np.logical_and(np.logical_and(self.impropers[:,0]==n[0],self.impropers[:,1]==n[1]),\
												self.impropers[:,2]==n[2]),self.impropers[:,3]==n[3]) 
                    if np.any(test):            
                        self.impropers[test,4]=n[-1]

            else:
                raise IOError("Cannot understand line %s in %s section!"%(n,target))


    #search bond type connecting current molecule with next one
    def search_next_bond(self,a,b):
        #a is current molecule, b next one
        #case a +b
        test1=np.array(np.where(self.bonds[:,0:2]=="+%s"%b))[0]
        test2=np.array(np.where(self.bonds[:,0:2]==a))[0]
        r=self.bonds[np.intersect1d(test1,test2),:]

        return r
    
      
    #search bond type connecting current molecule with previous one                
    def search_prev_bond(self,a,b):
        #a is current molecule, b next one        
        #case -a b
        test1=np.array(np.where(self.bonds[:,0:2]==b))[0]
        test2=np.array(np.where(self.bonds[:,0:2]=="-%s"%a))[0]
        r=self.bonds[np.intersect1d(test1,test2),:]

        return r

        
    #search angle type connecting current molecule with next one     
    def search_next_angle(self,a,b,sign):
        #a (and c) is current molecule, b next one
        #c a +b
        if sign == "+": #caller is current molecule in polymer
            r=filter(lambda x: np.any(x=="+%s"%b) and np.any(x==a) and len([elem for elem in x[0:3] if "+" in elem])==1, self.angles)
            
        #-c -a b
        elif sign=="-": #caller is next molecule in polymer
            r=filter(lambda x: np.any(x=="-%s"%a) and np.any(x==b) and len([elem for elem in x[0:3] if "-" in elem])==2, self.angles)

        else:
            raise IOError("ERROR: sign must be - or +")
  
        return np.array(list(r))
        
    #search angle type connecting current molecule with previous one     
    def search_prev_angle(self,a,b,sign):
        #a is current molecule, c b next one        
        #c b -a
        if sign=="-": #caller is next molecule in polymer
            r=filter(lambda x: np.any(x=="-%s"%a) and np.any(x==b) and len([elem for elem in x[0:3] if "-" in elem])==1, self.angles)
            
        #+c +b a
        elif sign=="+": #caller is current molecule in polymer
            r=filter(lambda x: np.any(x=="+%s"%b) and np.any(x==a) and len([elem for elem in x[0:3] if "+" in elem])==2, self.angles)

        else:
            raise IOError("ERROR: sign must be - or +")
  
        return np.array(list(r))            
        
    #search dihedral type connecting current molecule with next one     
    def search_next_dihedral(self,a,b,sign):
        #a is current molecule, b next one
        #-d -c -a b
        if sign=="-": #caller is next molecule in polymer
            r=filter(lambda x: np.any(x=="-%s"%a) and np.any(x==b) and len([elem for elem in x[0:4] if "-" in elem])==3, self.dihedrals)
            
        #d c a +b
        elif sign=="+": #caller is current molecule in polymer
            r=filter(lambda x: np.any(x=="+%s"%b) and np.any(x==a) and len([elem for elem in x[0:4] if "+" in elem])==1, self.dihedrals)

        else:
            raise IOError("ERROR: sign must be - or +")

        return np.array(list(r))
        
    #search dihedral type connecting current molecule with previous one     
    def search_prev_dihedral(self,a,b,sign):
        #a is current molecule, b next one        
        #d c b -a        
        if sign=="-": #caller is next molecule in polymer
            r=filter(lambda x: np.any(x=="-%s"%a) and np.any(x==b) and len([elem for elem in x[0:4] if "-" in elem])==1, self.dihedrals)
            
        #+d +c +b a    
        elif sign=="+": #caller is current molecule in polymer
            r=filter(lambda x: np.any(x=="+%s"%b) and np.any(x==a) and len([elem for elem in x[0:4] if "+" in elem])==3, self.dihedrals)

        else:
            raise IOError("ERROR: sign must be - or +")

        return np.array(list(r))


    def check_existence(self,atom):
        if np.any(self.mapping[:,1]==atom):
            return True
        else:
            return False
        

if __name__=="__main__":
    
    Top=Topology()
    #Top.load("database/cis-PI-monomer.txt")
    #Top.load("database/vinyl-PBD-monomer.txt")
    Top.load("database/cis-PI-monomer.txt")    
       
    #print(Top.search_prev_dihedral('ZZ','AA','-'))
    #print(Top.search_next_dihedral('ZZ','AA','+'))
  
    Top.make_terminal("nterminal")
    Top.make_terminal("cterminal")
    
    print(Top.mapping)    
    print(Top.bonds)
    print(Top.angles)
    print(Top.dihedrals)
    