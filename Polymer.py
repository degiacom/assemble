# Copyright (c) 2014 Matteo Degiacomi and Valentina Erastova
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
from copy import deepcopy
import logging

class Polymer(object):

    def __init__(self,db,ff,molname,mode):
        self.db=db
        self.poly=[]
        self.ff=ff
        self.molname=molname
        self.mode=mode        
        self.clash_thresh=0.9
        
        self.search_grid=self._make_search_grid()
        
        self.logger=logging.getLogger('assemble')
        
    def make(self,chain):
        
        #debug=0 #debug mode prints hooks positions in a separate file
        
        self.chain=chain
        
        self.logger.info("\n> generating polymer %s..."%self.molname)
        self.logger.info(">> sequence: %s"%self.chain)

        #add first element in the chain
        m=deepcopy(self.db.molecules[self.chain[0]])
        self.poly.append(m)

        #if debug: #DEBUG: print positions of HOOKS IN A SEPARATE FILE
        #    f_out = open("hooks.pdb", 'w')
                
        #iterate over chain string and build data structure
        for x in xrange(1,len(self.chain),1):
        
            #get new monomer
            m_new=deepcopy(self.db.molecules[self.chain[x]])

            #get head (of new monomer) and tail (of chain) coordinates
            tail=deepcopy(m.data[m.data[:,0]==int(m.limit['tail']),5:8][0])
            head=deepcopy(m_new.data[m_new.data[:,0]==int(m_new.limit['head']),5:8][0])

            #compute hooking point of existing polymer chain for new monomer
            if self.mode=="pdb":
                tail_hook=deepcopy(m.data[m.data[:,0]==int(m.limit['tail_hook']),5:8][0])
                head_hook=deepcopy(m_new.data[m_new.data[:,0]==int(m_new.limit['head_hook']),5:8][0])

            else:
                
                #get head and tail atomnames
                tailname=m.topology.tail[0]
                headname=m_new.topology.head[0]
                
                ###GET BOND###
                #get bond type within head and tail (must be one and one only!)
                b1=m.topology.search_next_bond(tailname,headname) #with plus
                b2=m_new.topology.search_prev_bond(tailname,headname) #with minus
                
                #check consistency between head and tail molecule topologies for bonds
                b=[]              
                if len(b1)>0 and len(b2)>0:             
                    b=np.unique(np.concatenate((b1[:,2],b2[:,2])))
                elif len(b1)>0 and len(b2)==0:
                    b=np.unique(b1[:,2])
                elif len(b1)==0 and len(b2)>0:
                    b=np.unique(b2[:,2])
                if len(b)>1:
                    #print "ERROR: inconsistency in bond descriptions in topologies of %s and %s"%(self.chain[x],self.chain[x-1])
                    raise IOError("inconsistency in bond descriptions in topologies of %s and %s"%(self.chain[x],self.chain[x-1]))
                if len(b)==0:
                    #print "ERROR: connection between %s and %s not found!"%(tailname,headname)
                    raise IOError("connection between %s and %s not found!"%(tailname,headname))
               
                #get bond distance in angstrom (same for both current and new molecules)
                try:
                    bond=self.ff.get_bond(b[0])
                except:
                    #print "ERROR: bond type %s not found in force field"%b[0]
                    raise IOError("bond type %s not found in force field"%b[0])

                ###GET DIHEDRAL CURRENT###
                keep=[0]
                a_tmp=m.topology.search_next_dihedral(tailname,headname,'+') #need 1 plus     
                a=self._remove_prev_to_next(a_tmp) 
                                                                                                                                                                                                                                                                                                       
                #if search failed on current molecule, look for parameters on next one
                if len(a)==0:  
                    keep=[]                  
                    a_tmp=m_new.topology.search_next_dihedral(tailname,headname,'-') #need double minus!
                    a=self._remove_prev_to_next(a_tmp) 
                    
                    #verify existence of all atoms in current molecule, if getting info from next one, and reformat
                    for i in xrange(0,len(a),1):
                                       
                        for j in xrange(0,4,1):
                            if "-" in  a[i,j]:
                                atomname=a[i,j].split("-")[1]
                                if atomname in m.atom:
                                    keep.append(i)
                                    a[i,j]=atomname
                            else:
                                    keep.append(i)
                                    atomname="+%s"%a[i,j]
                                    a[i,j]=atomname
                                    
                #check if a match was found in topologies of current or next molecule
                if len(keep)==0:
                    #print "no match found for dihedral!"
                    raise IOError("no match found for dihedral for forward hook involving atoms\n %s in %s and %s in %s!"%(tailname, m.topfile, headname, m_new.topfile))
  
  
                #get dihedral angle value
                try:
                    dihedral_val_tail=self.ff.get_dihedral(a[keep[0],4])
                except:
                    #print "ERROR: dihedral type %s not found in force field"%a[0,4]          
                    raise IOError("dihedral type %s not found in force field"%a[keep[0],4])  
                
                #extract names of atoms forming angle and dihedral with head and tail
                if "+" in a[keep[0],0] and a[keep[0],1]==tailname:
                    anglename=a[keep[0],2]
                    dihedralname=a[keep[0],3]
                elif "+" in a[keep[0],3] and a[keep[0],2]==tailname:
                    anglename=a[keep[0],1]
                    dihedralname=a[keep[0],0]
                else:
                    #print "umm... a dihedral potential looks weird..."
                    IOError("umm... a dihedral potential looks weird...")
                
                               
                #SEARCH FOR ANGLENAME AS WELL!
                ###GET ANGLE CURRENT###
                keep=[0]
                a=m.topology.search_next_angle(tailname,headname,'+') #need 1 plus   
                #if search failed on current molecule, look for parameters on next one
                if len(a)==0:
                    keep=[]                    
                    a=m_new.topology.search_next_angle(tailname,headname,'-') #need double minus!
                    #verify existence of all atoms in current molecule, if getting info from next one
                    for i in xrange(0,len(a),1):
                        for j in xrange(0,3,1):
                            if "-" in  a[i,j]:
                                atomname=a[i,j].split("-")[1]
                                if  atomname in m.atom:
                                    keep.append(i)
                                    a[i,j]=atomname
                            else:
                                    keep.append(i)
                                    atomname="+%s"%a[i,j]
                                    a[i,j]=atomname
                                              
                #check that one solution was found
                if len(keep)==0:
                    #print "no match found for angle!"
                    raise IOError("no match found for angle for forward hook involving atoms\n %s in %s and %s in %s!"%(tailname, m.topfile, headname, m_new.topfile))
                
                try:
                    angle_val_tail=self.ff.get_angle(a[keep[0],3])
                except:
                    #print "ERROR: angle type %s not found in force field"%a[0,3]                
                    raise IOError("angle type %s not found in force field"%a[keep[0],3])         

                
                #compute hooking point position for current molecule
                coor_bond_tail=m.atomselect("*","*",tailname)[0]
                coor_angle_tail=m.atomselect("*","*",anglename)[0]
                coor_dihedral_tail=m.atomselect("*","*",dihedralname)[0]
                tail_hook=self._place_pseudoatom(coor_bond_tail,coor_angle_tail,coor_dihedral_tail,bond,angle_val_tail,dihedral_val_tail)

                ###GET DIHEDRAL NEXT###
                keep=[0]
                a_tmp=m_new.topology.search_prev_dihedral(tailname,headname,'-') #need 1 minus 
                a=self._remove_prev_to_next(a_tmp) 
                
                #if search failed on current molecule, look for parameters on next one
                if len(a)==0:  
                    keep=[]                  
                    a_tmp=m.topology.search_prev_dihedral(tailname,headname,'+') #need double plus!
                    a=self._remove_prev_to_next(a_tmp) 
                    
                    #verify existence of all atoms in current molecule, if getting info from next one
                    for i in xrange(0,len(a),1):
                        for j in xrange(0,4,1):
                            if "+" in a[i,j]:
                                atomname=a[i,j].split("+")[1]
                                if  atomname in m.atom:
                                    keep.append(i)
                                    a[i,j]=atomname
                            else:
                                    keep.append(i)
                                    atomname="-%s"%a[i,j]
                                    a[i,j]=atomname
                                
                #check that one solution was found
                if len(keep)==0:
                    #print "no match found for angle!"
                    raise IOError("no match found for dihedral for backward hook involving atoms\n %s in %s and %s in %s!"%(headname, m_new.topfile,tailname, m.topfile))
                
                try:
                    dihedral_val_head=self.ff.get_dihedral(a[keep[0],4])
                except:
                    #print "ERROR: dihedral type %s not found in force field"%a[0,4]                
                    raise IOError("dihedral type %s not found in force field"%a[keep[0],4])
                
                #extract names of atoms forming angle and dihedral with head and tail
                if "-" in a[keep[0],0] and a[keep[0],1]==headname:
                    anglename=a[keep[0],2]
                    dihedralname=a[keep[0],3]
                elif "-" in a[keep[0],3] and a[keep[0],2]==headname:
                    anglename=a[keep[0],1]
                    dihedralname=a[keep[0],0]
                else:
                    #print a[0,:]
                    IOError("%s dihedral potential looks weird..."%a[keep[0],:])
                 
                ###GET ANGLE NEXT###
                keep=[0]
                a=m_new.topology.search_prev_angle(tailname,headname,'-') #need 1 minus 
                #if search failed on current molecule, look for parameters on next one
                if len(a)==0:                    
                    a=m.topology.search_prev_angle(tailname,headname,'+') #need double plus!
                    #verify existence of all atoms in current molecule, if getting info from next one
                    keep=[]
                    for i in xrange(0,len(a),1):
                        for j in xrange(0,3,1):
                            if "+" in  a[i,j]:
                                atomname=a[i,j].split("+")[1]
                                if atomname in m.atom:
                                    keep.append(i)
                                    a[i,j]=atomname
                            else:
                                    keep.append(i)
                                    atomname="-%s"%a[i,j]
                                    a[i,j]=atomname
                
                #check that one solution was found
                if len(keep)==0:
                    raise IOError("no match found for angle for backward hook involving atoms\n %s in %s and %s in %s!"%(headname, m_new.topfile,tailname, m.topfile))
               
                try:
                    angle_val_head=self.ff.get_angle(a[keep[0],3])
                except:
                    #print "ERROR: angle type %s not found in force field"%a[0,3]                
                    raise IOError("angle type %s not found in force field"%a[keep[0],3])

                #compute hooking point position for next molecule
                coor_bond_head=m_new.atomselect("*","*",headname)[0]
                coor_angle_head=m_new.atomselect("*","*",anglename)[0]
                coor_dihedral_head=m_new.atomselect("*","*",dihedralname)[0]
                head_hook=self._place_pseudoatom(coor_bond_head,coor_angle_head,coor_dihedral_head,bond,angle_val_head,dihedral_val_head)
            
            #ADD NEW CHAIN WITH CLASH DETECTION, IF IN GROMACS MODE###
            solved=False
            for i in xrange(0,len(self.search_grid),1):
                    
                ###compute superimposition###
                #prepare data to compute superimposition within previous tail and new head           
                t=np.array([tail,tail_hook])        
                h=np.array([head_hook,head])

                COM_h=np.sum(h,axis=0)/float(len(h))
                COM_t=np.sum(t,axis=0)/float(len(t))
            
                #compute the rotation matrix superimposing the current head with the previous tail
                trans=self._rmsd(t,h)[0]

                #bring new monomer to origin, rotate it, and translate it so that tail_hook superimposes with head
                crds_new=np.dot(m_new.get_xyz()-COM_h,trans)+COM_t
            
                if self.mode=="pdb": #no clash detection for pdb mode
                    solved=True
                    break
            
                if not self._clash_test(crds_new,self.get_xyz()): #if new molecule is clash free, add to chain
                    solved=True
                    break
                
                #perturb hooks position
                head_hook=self._place_pseudoatom(coor_bond_head,coor_angle_head,coor_dihedral_head,bond,angle_val_head,dihedral_val_head+self.search_grid[i,0])
                tail_hook=self._place_pseudoatom(coor_bond_tail,coor_angle_tail,coor_dihedral_tail,bond,angle_val_tail,dihedral_val_tail+self.search_grid[i,1])
                    
                if solved:
                    break

            if not solved:
                self.logger.info(">> WARNING: unsolved clash between %s and %s. Continuing..."%(m.topfile, m_new.topfile))
            
            m_new.set_xyz(crds_new)

            #coor_bond=m_new.atomselect("*","*",headname)[0]
            
            '''
            if debug: #DEBUG: print positions of HOOKS IN A SEPARATE FILE
                l=(x,"TAI","BND",'P',x,tail[0],tail[1],tail[2],1.0,1.0,"P")
                L='ATOM  %5i  %-4s%-4s%1s%4i    %8.3f%8.3f%8.3f%6.2f%6.2f          %2s\n'%l
                f_out.write(L)        
                
                l=(x+1,"THK","BND",'P',x,tail_hook[0],tail_hook[1],tail_hook[2],1.0,1.0,"P")
                L='ATOM  %5i  %-4s%-4s%1s%4i    %8.3f%8.3f%8.3f%6.2f%6.2f          %2s\n'%l
                f_out.write(L)        
                
                head_hook_new=np.dot(head_hook-COM_h,trans)+COM_t
                l=(x+2,"HHK","BND",'P',x+1,head_hook_new[0],head_hook_new[1],head_hook_new[2],1.0,1.0,"P")
                L='ATOM  %5i  %-4s%-4s%1s%4i    %8.3f%8.3f%8.3f%6.2f%6.2f          %2s\n'%l
                f_out.write(L)               
                
                head_new=np.dot(head-COM_h,trans)+COM_t
                l=(x+3,"HEA","BND",'P',x+1,head_new[0],head_new[1],head_new[2],1.0,1.0,"P")
                L='ATOM  %5i  %-4s%-4s%1s%4i    %8.3f%8.3f%8.3f%6.2f%6.2f          %2s\n'%l
                f_out.write(L)               
            '''    
                
            #push new monomer in chain
            self.poly.append(m_new)
            
            #newly added monomer becomes first of existing chain (and therefore reference for next monomer to hook to the chain)
            m=deepcopy(m_new)


        #if debug: #DEBUG: print positions of HOOKS IN A SEPARATE FILE       
        #    f_out.close()

        #extract all atom coordinates and align them along inertia tensor
        crds=self.get_xyz()
        crds=self._align_axes(crds)
        self.set_xyz(crds)
        
        
    def write_polymer(self,typef="pdb",mypath="."):
                
        #renumber atoms index and resid, set same chain name to all polymer
    
        self.logger.info(">> writing PDB file %s.pdb"%self.molname)
    
        f_out = open("%s/%s.pdb"%(mypath, self.molname), 'w')

        f_out.write("REMARK generated with Assemble.py, by Matteo Degiacomi and Valentina Erastova, 2014\n")
        f_out.write("REMARK sequence: %s\n"%self.chain)

        index=1    
        for j in xrange(0,len(self.poly),1):
            data_list=self.poly[j].mapping(self.poly[j].data)
            
            if typef=="pdb":
                tail_hook=int(self.poly[j].limit['tail_hook'])
                head_hook=int(self.poly[j].limit['head_hook'])
            else:
                tail_hook=""
                head_hook=""
            
            for i in xrange(0,len(data_list),1):
                #create and write PDB line of non pseudoatoms
                if data_list[i][0]!=head_hook and data_list[i][0]!=tail_hook:        
                #if True:       
                    l=(index,data_list[i][1],data_list[i][2],'P',j+1,data_list[i][5],data_list[i][6],data_list[i][7],data_list[i][8],data_list[i][9],data_list[i][10])
                    L='ATOM  %5i  %-4s%-4s%1s%4i    %8.3f%8.3f%8.3f%6.2f%6.2f          %2s\n'%l
                    f_out.write(L)        
                    index+=1

        f_out.close()


    #return coordinates of all the atoms in the system
    def get_xyz(self):

        #count amount of atoms and extract all atomic coordinates
        cnt=0
        pos=[] #atom positions (need for easily computing box size)
        for j in xrange(0,len(self.poly),1):
            data_list=self.poly[j].mapping(self.poly[j].data)
            cnt+=len(data_list)
            for i in xrange(0,len(data_list),1):
                pos.append([data_list[i][5],data_list[i][6],data_list[i][7]]) #store atom position

        return np.array(pos)

    #push atom coordinates in all molecules composing the polymer (suppose length matching)
    def set_xyz(self,crds):

        cnt=0
        for j in xrange(0,len(self.poly),1):        
            for i in xrange(0, len(self.poly[j].data),1):
                self.poly[j].data[i][5]=crds[cnt,0]
                self.poly[j].data[i][6]=crds[cnt,1]
                self.poly[j].data[i][7]=crds[cnt,2]
                cnt+=1
    
    
    def write_gromacs(self,mypath="."):

        ###WRITE .GRO FILE###        
 
        #count amount of atoms and extract all atomic coordinates
        cnt=0
        pos=[] #atom positions (need for easily computing box size)
        for j in xrange(0,len(self.poly),1):
            data_list=self.poly[j].mapping(self.poly[j].data)
            cnt+=len(data_list)
            for i in xrange(0,len(data_list),1):
                pos.append([data_list[i][5]/10.0,data_list[i][6]/10.0,data_list[i][7]/10.0]) #store atom position
        
        self.p=self.get_xyz()/10.0

        #compute box size and position of minimal value (needed to shift the protein in a region defined by the box)
        #NOTE: a 1A padding is added in every direction
        minpos=np.min(self.p,axis=0)-0.1
        self.box=np.max(self.p,axis=0)-minpos+0.2
        
        #prepare gromacs coordinates file
        self.logger.info(">> writing Gromacs coordinates file %s.gro"%self.molname)
        f_out = open("%s/%s.gro"%(mypath,self.molname), 'w')
        f_out.write("%s\n"%self.chain)
        f_out.write("%s\n"%len(self.p))
        
        #temporary storages for topology-related information
        at=[] #atom topology information       
        b=[] #bonds topology
        a=[] #angles topology
        d=[] #dihedrals topology
        imp=[] #impropers topology
        
        index=1 #atom counter
        for j in xrange(0,len(self.poly),1):

            data_list=self.poly[j].mapping(self.poly[j].data)

            #get topology information of current molecule
            top=self.poly[j].topology
            
            #if molecule is terminal, modify its topology accordingly
            if j==0:
                top.make_terminal("nterminal")
            if j==len(self.poly)-1:
                top.make_terminal("cterminal")                
            
            #gather bonds, angles and dihedrals information for topologies
            b.append(top.bonds)
            a.append(top.angles)
            d.append(top.dihedrals)
            imp.append(top.impropers)
            
            #@todo add helper lines in topology lines, defining what parameters are
            for i in xrange(0,len(data_list),1):
                #create and write line in gromacs format (.gro file)
                #NOTE: the molecule is moved so that its minimal position is at the origin
                #resname="%s%s"%(j+1,data_list[i][2])
                #L='%7s%7s%5i%8.3f%8.3f%8.3f\n'%(resname,data_list[i][1],index,data_list[i][5]/10.0-minpos[0],data_list[i][6]/10.0-minpos[1],data_list[i][7]/10.0-minpos[2])
                L='%5d%-5s%5s%5d%8.3f%8.3f%8.3f\n'%(j+1,data_list[i][2],data_list[i][1],index,data_list[i][5]/10.0-minpos[0],data_list[i][6]/10.0-minpos[1],data_list[i][7]/10.0-minpos[2])
                f_out.write(L)
                
                ##store topology information for later writing
                ##atomID, atomtype, resid, resname, atomname, number like resid (?), charge, mass
                #if j==0 and data_list[i][1]==top.head[0]:
                #    atomtype=top.head[1]
                #elif j==len(self.poly)-1 and data_list[i][1]==top.tail[0]:
                #    atomtype=top.tail[1]
                #else:
                atomtype=top.mapping[top.mapping[:,0]==data_list[i][1],1][0]
                    
                mass=self.ff.nonbonded[atomtype][1]
                charge=self.ff.nonbonded[atomtype][2]

                at.append([index,atomtype,j+1,data_list[i][2],data_list[i][1],j+1,charge,mass])
                
                index+=1
                if index>99999:
                    index=0

        #print box size at the end of gromacs coordinates file
        f_out.write("%10.5f%10.5f%10.5f\n"%(self.box[0],self.box[1],self.box[2]))
        f_out.close()

        ###GENERATE TOPOLOGY FILE### 
        self.logger.info(">> writing Gromacs topology file %s.itp"%self.molname) 
        #prepare data structure for indexing in numpy format
        atom_top=np.array(at).astype(str)
        
        f_out = open("%s/%s.itp"%(mypath,self.molname), 'w')
        f_out.write("; generated with Assemble.py, by Matteo Degiacomi and Valentina Erastova, 2014\n")
        f_out.write("; sequence: %s\n"%self.chain)
        
        #write header statements
        f_out.write("\n[ moleculetype ]\n%s       3\n"%self.molname)
        
        #write atoms lines
        f_out.write("\n [ atoms ]\n")
        for j in xrange(0,len(atom_top),1):
            f_out.write("%6s%11s%7s%7s%7s%7s%11s%11s\n"%(atom_top[j,0],atom_top[j,1],atom_top[j,2],atom_top[j,3],atom_top[j,4],atom_top[j,5],int(float(atom_top[j,6])),atom_top[j,7]))
    
        #INSERT PARAMETERS VALUES FROM FORCE FIELD INSTEAD OF S-BN...
    
        #write bond lines
        f_out.write("\n [ bonds ] \n")
        for j in xrange(0,len(b),1):
            for x in xrange(0,len(b[j]),1):
                b0=self._get_index(atom_top,j,b[j][x][0])
                b1=self._get_index(atom_top,j,b[j][x][1])
                if b0!=False and b1!=False:
                    vals_text='  '.join(self.ff.bonded[b[j][x][2]].astype(str))
                    f_out.write("%5s %6s %6s %8s\n"%(b0,b1,self.ff.fftype[0],vals_text))

        #write angles lines
        f_out.write("\n [ angles ] \n")
        for j in xrange(0,len(a),1):
            for x in xrange(0,len(a[j]),1):
                b0=self._get_index(atom_top,j,a[j][x][0])
                b1=self._get_index(atom_top,j,a[j][x][1])
                b2=self._get_index(atom_top,j,a[j][x][2])
                if b0!=False and b1!=False and b2!=False:
                    vals_text='  '.join(self.ff.bonded[a[j][x][3]].astype(str))
                    f_out.write("%5s %6s %6s %6s %8s\n"%(b0,b1,b2,self.ff.fftype[1],vals_text))

        #write dihedrals lines
        f_out.write("\n [ dihedrals ] \n")
        for j in xrange(0,len(d),1):
            for x in xrange(0,len(d[j]),1):
                b0=self._get_index(atom_top,j,d[j][x][0])
                b1=self._get_index(atom_top,j,d[j][x][1])
                b2=self._get_index(atom_top,j,d[j][x][2])
                b3=self._get_index(atom_top,j,d[j][x][3])
                           
                if b0!=False and b1!=False and b2!=False and b3!=0:
                    vals_text='  '.join(self.ff.bonded[d[j][x][4]].astype(str))
                    f_out.write("%5s %6s %6s %6s %6s %8s\n"%(b0,b1,b2,b3,self.ff.fftype[2],vals_text))
         
        #write impropers lines
        #f_out.write("\n [ impropers ] \n")
        for j in xrange(0,len(imp),1):
            for x in xrange(0,len(imp[j]),1):
                b0=self._get_index(atom_top,j,imp[j][x][0])
                b1=self._get_index(atom_top,j,imp[j][x][1])
                b2=self._get_index(atom_top,j,imp[j][x][2])
                b3=self._get_index(atom_top,j,imp[j][x][3])
                           
                if b0!=False and b1!=False and b2!=False and b3!=0:
                    vals_text='  '.join(self.ff.bonded[imp[j][x][4]].astype(str))
                    f_out.write("%5s %6s %6s %6s %6s %8s\n"%(b0,b1,b2,b3,self.ff.fftype[3],vals_text))


        f_out.write("\n#ifdef POSRES\n#include \"posre.itp\"\n#endif\n")
     
        f_out.close()
        
        return
  
    #generate and return list for clash avoidance scan
    def _make_search_grid(self,step=5):

        a=[]

        for i in xrange(0,180,step):
            for j in xrange(0,i+step,step):
                a.append([i,j])
                if i>0:
                    a.append([-i,j])
                if j>0:
                    a.append([i,-j])
                if i>0 and j>0:
                    a.append([-i,-j])
        
        return np.array(a)
  
    #check whether two ensembles of points have a couple at a distance less than a given threshold
    #return true if clash detected, false otherwise
    def _clash_test(self,points1,points2):
        
        for i in xrange(0,len(points1),1):
            dists=np.sqrt(np.sum((points2-points1[i])**2,axis=1))
            if np.any(dists<self.clash_thresh):
                return True
        
        return False
    
    ## compute matrix needed to rotate the system around an arbitrary axis (using Euler-Rodrigues formula).
    # @param axis 3d vector (numpy array), representing the axis around which to rotate
    # @param theta desired rotation angle
    # @retval 3x3 rotation matrix
    def _rotation_matrix(self, axis, theta):

        #if rotation angle is equal to zero, no rotation is needed
        if theta==0:
                return np.identity(3)

        #method taken from http://stackoverflow.com/questions/6802577/python-rotation-of-3d-vector
        axis = axis/np.sqrt(np.dot(axis,axis))
        a = np.cos(theta/2)
        b,c,d = -axis*np.sin(theta/2)
        return np.array([[a*a+b*b-c*c-d*d, 2*(b*c-a*d), 2*(b*d+a*c)],
                         [2*(b*c+a*d), a*a+c*c-b*b-d*d, 2*(c*d-a*b)],
                         [2*(b*d-a*c), 2*(c*d+a*b), a*a+d*d-b*b-c*c]])

  
  
    ## compute Structure's principal axes.
    # @retval 3x3 numpy array, containing the 3 principal axes ranked from smallest to biggest.
    def _get_principal_axes(self,points):
        #method taken from chempy source code, geometry.py, method getMomentOfInertiaTensor()

        #compute moment of inertia tensor
        I0 = np.zeros((3,3), np.float64)
        for i in xrange(0,len(points),1):
            mass = 1#self.mass[atom] / constants.Na
            I0[0,0] += mass * (points[i,1] * points[i,1] + points[i,2] * points[i,2])
            I0[1,1] += mass * (points[i,0] * points[i,0] + points[i,2] * points[i,2])
            I0[2,2] += mass * (points[i,0] * points[i,0] + points[i,1] * points[i,1])
            I0[0,1] -= mass * points[i,0] * points[i,1]
            I0[0,2] -= mass * points[i,0] * points[i,2]
            I0[1,2] -= mass * points[i,1] * points[i,2]

        I0[1,0] = I0[0,1]
        I0[2,0] = I0[0,2]
        I0[2,1] = I0[1,2]

        #Calculate and return the principal moments of inertia and corresponding
        #principal axes for the current geometry.
        e_values,e_vectors = np.linalg.eig(I0)

        indices = np.argsort(e_values)
        e_values = e_values[indices]
        e_vectors = e_vectors.T[indices]

        return e_vectors

     
    ##align structure on its principal axes.
    # first principal axis aligned along x, second along y and third along z.
    def _align_axes(self,points):
    
        #this method is inspired from the procedure followed in in VMD's orient package:
        ## set I [draw principalaxes $sel]           <--- show/calc the principal axes
        ## set A [orient $sel [lindex $I 2] {0 0 1}] <--- rotate axis 2 to match Z
        ## $sel move $A
        ## set I [draw principalaxes $sel]           <--- recalc principal axes to check
        ## set A [orient $sel [lindex $I 1] {0 1 0}] <--- rotate axis 1 to match Y
        ## $sel move $A
        ## set I [draw principalaxes $sel]           <--- recalc principal axes to check

        #center polymer to origin
        c=np.mean(points,axis=0)
        points[:,0]-=c[0]
        points[:,1]-=c[1]
        points[:,2]-=c[2]

        #get principal axes (ranked from smallest to biggest)
        axes=self._get_principal_axes(points)

        #align smallest principal axis against z axis
        rotvec=np.cross(axes[0],np.array([1,0,0])) #rotation axis
        sine=np.linalg.norm(rotvec)
        cosine=np.dot(axes[0],np.array([1,0,0]))
        angle=np.arctan2(sine,cosine) #angle to rotate around axis

        rotmatrix=self._rotation_matrix(rotvec,angle)
        points=np.dot(points, rotmatrix)

        #compute new principal axes (after previous rotation)
        axes=self._get_principal_axes(points)

        #align second principal axis against y axis
        rotvec=np.cross(axes[1],np.array([0,1,0])) #rotation axis
        sine=np.linalg.norm(rotvec)
        cosine=np.dot(axes[1],np.array([0,1,0]))
        angle=np.arctan2(sine,cosine) #angle to rotate around axis

        rotmatrix=self._rotation_matrix(rotvec,angle)
        points=np.dot(points, rotmatrix)
    
        return points

    #test for connection between previous and next molecule. If existing, remove from pool
    def _remove_prev_to_next(self, a):
        keep=[]
        for i in xrange(0,len(a),1):    
            prevmol=False
            nextmol=False
            for j in xrange(0,4,1):
                if "-" in  a[i,j]:
                    prevmol=True
                if "+" in  a[i,j]:
                    nextmol=True                                
                
            if not (prevmol and nextmol):
                keep.append(i)
        
        if len(keep)>0:
            return a[keep]
        else:
            return []


    def _get_index(self,top,resid,b0):
        if "+" in b0:
            name=top[np.logical_and(top[:,2]==str(resid+2), top[:,4]==b0.split("+")[1]),0]
        elif "-" in b0:
            name=top[np.logical_and(top[:,2]==str(resid), top[:,4]==b0.split("-")[1]),0]
        else:
            name=top[np.logical_and(top[:,2]==str(resid+1), top[:,4]==b0),0]
                
        if len(name)==0:
            return False
        elif len(name)==1:
            return name[0]
        else:
            raise IOError("ERROR: multiple instances of atom %s found in residue %s!"%(b0,resid))
    
    def _rmsd(self,m1,m2):
        
        L = len(m1)
        
        ##protein is already centered, don't need centering!
        COM1 = np.sum(m1,axis=0) / float(L)
        COM2 = np.sum(m2,axis=0) / float(L)
        m1 -= COM1
        m2 -= COM2
     
        E0 = np.sum( np.sum(m1*m1,axis=0),axis=0) + np.sum( np.sum(m2*m2,axis=0),axis=0)
     
        #This beautiful step provides the answer. V and Wt are the orthonormal
        # bases that when multiplied by each other give us the rotation matrix, U.
        # S, (Sigma, from SVD) provides us with the error!  Isn't SVD great!
        V, S, Wt = np.linalg.svd( np.dot( np.transpose(m2), m1))

        reflect = float(str(float(np.linalg.det(V) * np.linalg.det(Wt))))
        
        if reflect == -1.0:
            S[-1] = -S[-1]
            V[:,-1] = -V[:,-1]

        RMSD = E0 - (2.0 * sum(S))
        RMSD = np.sqrt(abs(RMSD / L))

        U = np.dot(V, Wt)
                
        return U, RMSD
        
    
    def _place_pseudoatom(self,coor_bond,coor_angle,coord_dihedral,bond,ang,di):
        
        angle=np.deg2rad(ang)
        dihed=np.deg2rad(di)

        ###inspired by https://github.com/molmod/molmod/blob/master/molmod/zmatrix.py###
        #coor_bond      #bond atom (tail/head)
        #coord_angle    #angle atom
        #coord_dihedral #furthest dihedral atom
        # define frame axes
        new_z = coor_angle - coor_bond #pos of angle atom - pos of tail/head
        norm_z = np.linalg.norm(new_z)
        if norm_z < 1e-15:
            new_z = np.array([0, 0, 1], float)
        else:
            new_z /= np.linalg.norm(new_z)
        new_x = coord_dihedral - coor_bond #bond atom - origin
        new_x -= np.dot(new_x, new_z)*new_z
        norm_x = np.linalg.norm(new_x)
        if norm_x < 1e-15:
            new_x = self._random_orthonormal(new_z)
        else:
            new_x /= np.linalg.norm(new_x)
        # we must make our axes frame left handed due to the poor IUPAC
        # definition of the sign of a dihedral angle.
        new_y = -np.cross(new_z, new_x)
        # coordinates of new atom:
        x = bond*np.cos(dihed)*np.sin(angle)
        y = bond*np.sin(dihed)*np.sin(angle)
        z = bond*np.cos(angle)
        coordinates = coor_bond + x*new_x + y*new_y + z*new_z

        return coordinates
        

    def _random_orthonormal(self,normal):
        #"""Return a random normalized vector orthogonal to the given vector"""
        normal_fns = [lambda a: np.array([0.0, -a[2], a[1]]),
                      lambda a: np.array([a[2], 0.0, -a[0]]),
                      lambda a: np.array([-a[1], a[0], 0.0])]
        
        u = normal_fns[np.argmin(np.fabs(normal))](normal)
        u /= np.linalg.norm(u)
        v = np.cross(normal, u)
        v /= np.linalg.norm(v)
        alpha = np.random.uniform(0.0, np.pi*2)
        return np.cos(alpha)*u + np.sin(alpha)*v