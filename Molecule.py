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


# Usage example:
#
# from Protein import Protein
# p=Molecule()
# p.import_pdb("filename.pdb")
# d=p.get_xyz()
# p.set_xyz(d+100)
# p.write_pdb("translated_protein.pdb")
# a=p.atomselect('A','*','CA')
# p.write_pdb("translated_protein_CA.pdb",a)

#packages
import numpy as np
import Topology as T

class Molecule:

    def __init__(self):

        #mhashtables (name to number)
        self.atom={}
        self.res={}
        self.chain={}
        self.atomtype={}
        
        #structure(whole pdb in numeric format, via hashtable)
        self.data=[]

        #atoms connectivity (hash calling bonds, angles, and dihedrals, corresponding a table)
        self.topology=T.Topology()
        
        #head, tail, hook head, hook tail
        self.limit={}

        self.pdbfile=""
        self.topfile=""
    
    #read a PDB as input
    def import_pdb(self,pdb,mode):

        try:
            f_in=open(pdb,"r")
        except IOError:
            raise IOError('ERROR: PDB file %s not found!'%pdb)

        self.pdbfile=pdb
        
        data_in=[]
        for line in f_in:
            record=line[0:6].strip()
   
            #extract info about boundaries
            if record=='LIMIT':
                w=line.split()
                if len(w)!=3:
                    raise IOError("ERROR: LIMIT keyword not correctly formatted")
                
                if w[1]=="HEAD":
                    self.limit['head']=w[2]
                elif w[1]=="TAIL":
                    self.limit['tail']=w[2]
                elif w[1]=="HEAD_HOOK":
                    self.limit['head_hook']=w[2]
                elif w[1]=="TAIL_HOOK":
                    self.limit['tail_hook']=w[2]
                else:
                    raise IOError("LIMIT keyword %s not found"%w[1])
            
            elif record=='ATOM':
                w=[]
                w.append(int(line[6:11]))

                a_name=line[12:17].strip()
                if not a_name in self.atom:
                    self.atom[a_name]=len(self.atom.values())+1
                w.append(self.atom[a_name])

                r_name=line[17:20].strip()
                if not r_name in self.res:
                    self.res[r_name]=len(self.res.values())+1
                w.append(self.res[r_name])

                ch=line[21].strip()
                if not ch in self.chain:
                    self.chain[ch]=len(self.chain.values())+1
                w.append(self.chain[ch])

                w.append(int(line[22:26]))	
                w.append(float(line[30:38]))
                w.append(float(line[38:46]))
                w.append(float(line[46:54]))

                try:
                    w.append(float(line[54:60]))
                except:
                    w.append(1.0)
                try:
                    w.append(float(line[60:66]))
                except:
                    w.append(0.0)
                try:
                    #w.append(line[76:78])
                    w.append(self.chain[ch])
                except:
                    w.append(self.chain[ch])

                data_in.append(w)
                               
                
        self.data=np.array(data_in).astype(float)
        
        if self.data.shape[0]<3:
            raise IOError("molecule %s should contain at least three atoms!")
            
    
        #check if termini have been provided
        if len(self.limit)<4 and mode=="pdb":
            raise IOError("LIMIT statement missing (HEAD, TAIL, HEAD_HOOK and TAIL_HOOK required)!")

        if len(self.limit)>0 and mode!="pdb":
            raise IOError("LIMIT statements should not be defined in pdb %s while using %s mode!"%(pdb,mode))
        
        f_in.close()


    def import_topology(self,top):
        
        #read the topology associated to the molecule
        try:
            f_in=open(top,"r")
        except IOError:
            raise IOError('Topology file %s not found!'%top)

        #load topology from file
        self.topology.load(top)
        self.topfile=top
        #verify also that every atom in topology corresponds to an atom in the structure
        for i in range(0,len(self.topology.mapping[:,0]),1):
            if not self.topology.mapping[i,0] in self.atom:
                raise IOError("Atom %s found in topology but not on coordinates!"%self.topology.mapping[i,0])
        
        #verify also that every atom in structure corresponds to an atom in topology            
        for i in range(0, len(self.atom), 1):
            if not np.any(self.topology.mapping[:,0]==list(self.atom)[i]):
                raise IOError("Atom %s found in topology but not on coordinates"%list(self.atom)[i])
            
        #compute index of termini (extract termini name, and look for its index in pdb file)
        self.limit['head']=self.index_from_name(self.topology.head[0])
        self.limit['tail']=self.index_from_name(self.topology.tail[0])
        
        return


    def index_from_name(self,name):
        return int(self.atomselect("*","*",name,True)[0][0,0])
        

    def get_xyz(self):
        return self.data[:,5:8]


    def set_xyz(self,coords):
        self.data[:,5:8]=coords


    def set_resid(self,resid):
        self.data[:,4]=resid


    def set_incr_index(self,start):
        for x in range(0,len(self.data),1):
            self.data[x,0]=start+x
            
        #if topology is loaded, increase too!


    def rotation(self,x,y,z):
        #angle in numpy need to be given in rad -> rad = deg * pi/180
        alpha = np.radians(x)
        beta = np.radians(y)
        gamma = np.radians(z)
        #ex.: rotation around axis x
        #|1 	0 		 0	   |
        #|0 	np.cos(alpha)      -np.sin(alpha)|
        #|0	np.sin(alpha)	np.cos(alpha) |
        Rx = np.array([[1,0,0], [0, np.cos(alpha), -np.sin(alpha)], [0, np.sin(alpha), np.cos(alpha)]])
        Ry = np.array([[np.cos(beta), 0, np.sin(beta)], [0, 1, 0], [-np.sin(beta), 0, np.cos(beta)]])
        Rz = np.array([[np.cos(gamma), -np.sin(gamma), 0], [np.sin(gamma), np.cos(gamma), 0], [0,0,1]])
        rotation = np.dot(Rx,np.dot(Ry,Rz))
        #multiply rotation matrix with each atom of the monomer
        self.data[:,5:8] = np.dot(self.data[:,5:8],rotation)


    def center(self):
        #compute protein center of geometry
        return np.mean(self.data[:,5:8],axis=0)


    def atomselect(self,chain,res,atom,get_index=False):

        #chain name boolean selector
        if chain=='*':
            chain_query=True
        else:
            chain_query=self.data[:,3]==self.chain[chain]

        #resid boolean selector
        if res=='*':
            res_query=True
        else:
            res_query=self.data[:,4]==res

        #atom name boolean selector
        if atom=='*':
            atom_query=True
        else:
            atom_query=self.data[:,1]==self.atom[atom]

        #mask needed in case chain, res and atom selection are all set to true
        select_all=self.data[:,0]!=-1

        #slice data array and return result (colums 5 to 7 contain xyz coords)
        query=np.logical_and(select_all,np.logical_and(np.logical_and(chain_query,res_query),atom_query))

        if get_index==True:
            return [self.data[query],np.where(query==True)[0]]
        else:
            #UPDATED! Was initially just return self.data[query]
            return self.data[query,5:8]

    
    def mapping(self,data):

        data_list=[]
        for i in range(0,len(data),1):
            #backmap to strings
            atom=[k for k, v in self.atom.items() if v == data[i,1]][0]
            res=[k for k, v in self.res.items() if v == data[i,2]][0]
            chain=[k for k, v in self.chain.items() if v == data[i,3]][0]
            #atomtype=[k for k, v in self.atomtype.items() if v == data[i,10]][0]
            atomtype=""
             
            l=(int(self.data[i,0]),atom,res,chain,int(self.data[i,4]),self.data[i,5],self.data[i,6],self.data[i,7],self.data[i,8],self.data[i,9],atomtype)
            data_list.append(l)
    
        return data_list


    def write_pdb(self,outname,data=[]):
        
        #map intergers to characters from input data (default: all the protein)
        if len(data)==0:
            data_list=self.mapping(self.data)
        else:
            data_list=self.mapping(data)

        f_out=open(outname,"w")

        for i in range(0,len(data_list),1):

            #create and write PDB line
            l=(data_list[i][0],data_list[i][1],data_list[i][2],data_list[i][3],data_list[i][4],data_list[i][5],data_list[i][6],data_list[i][7],data_list[i][8],data_list[i][9],data_list[i][10])
            L='ATOM  %5i  %-4s%-4s%1s%4i    %8.3f%8.3f%8.3f%6.2f%6.2f          %2s\n'%l
            f_out.write(L)            

        f_out.close()
        return
