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

class System:
     
    def __init__(self,polymers,ff,params): 
        self.ff=ff
        self.polymers=polymers
        self.params=params
        
        self.logger=logging.getLogger('assemble')
        
    #generate a random distribution of molecules in a box, given their percentages
    def make_box(self, dim, data_in, use_fractional_mass):
        
        #prepare data
        names=[]
        v=[]
        for d in data_in:
            names.append(d[0])
            v.append(d[1])
        
        #if reference percentage is a fractional mass, tweak it using monomers mass
        if use_fractional_mass:
            percentage=self.convert_concentration(np.array(v), get_fractional=True)
        else:
            percentage=np.array(v)
        
        #normalize to 100 percent
        test=0           
        for l in percentage:
            test+=float(l)

        for i in range(len(percentage)):
            #v=float(percentage[i])
            percentage[i]=(percentage[i]/test)*100.0
        
        if len(dim)!=3:
            raise IOError("ERROR: expected 3 values for dimensions")
        
        length=dim[0]*dim[1]*dim[2]
        
        #make roulette array
        roulette=[percentage[0]]
        t=[percentage[0]]
        for r in range(1,len(percentage),1):
            t.append(percentage[r])
            roulette.append(roulette[r-1]+percentage[r])
    
        target=np.array(t)
        error=100000.0
        cbest=np.array([])
        #attempt distributing molecules in box according to desired concentration
        #keep best of 100 attemps
        for i in range(1,100):

            c=[]
    
            #pick a random monomer according to desired percentages
            counter=np.zeros(len(roulette))
            for x in range(0,length,1):
                rnd=np.random.rand(1)[0]*100
                index=0
                while True:
                    if roulette[index]>rnd:
                        break
                    else:
                        index+=1
                c.append(names[index])
                counter[index]+=1
    
            #chain completed, report produced percentages:
            counter/=float(length)
            counter*=100.0
        
            tmp_err=np.sum((counter-target)**2)
        
            if tmp_err<error:
                error=tmp_err
                cbest=np.array(c)
                counterbest=counter.copy()

        #print(counterbest)

        #if reference percentage is a fractional mass, tweak it using monomers mass
        if use_fractional_mass:
            perc=self.convert_concentration(counterbest, get_fractional=False) #convert into percentage
            perc/=np.sum(perc)
            perc*=100.0
            mw=counterbest.copy() #this is already a fractional mass
        else:
            perc=counterbest.copy()
            mw=self.convert_concentration(counterbest, get_fractional=True)
            mw/=np.sum(mw)
            mw*=100.0
        
        self.logger.info("> polymer units in box:")
        for x in range(0,len(counterbest),1):
            self.logger.info(">> %s: %s percent"%(names[x],perc[x]))
                      
                      
        self.logger.info("\n> polymers fractional weight:")
        for x in range(0,len(counterbest),1):
            self.logger.info(">> %s: %s percent"%(names[x],mw[x]))
                      
      
        return np.reshape(cbest,dim)

        
    #if get_fractional=True, percentage is edited using monomers mass
    #if get_fractional=False, fractional mass is converted into percentage
    def convert_concentration(self, c, get_fractional=True):
        m=[]
        for x in range(0,len(self.polymers),1):
            m.append(self.polymers[x].mass)
        allmasses=np.array(m)
        allmasses/=np.sum(allmasses)

        for x in range(0,len(allmasses),1):
            if get_fractional:   
                rescaled=c[x]/allmasses[x]
            else:   
                rescaled=c[x]*allmasses[x]
            
            c[x]=rescaled

        return c
        
        
    def create_system(self,mypath="."):

        ### CREATE BOX ###                         
        self.systembox=self.make_box(self.params.box_grid_shape, self.params.concentration, use_fractional_mass=False)

        #get maximal box between existing polymers, to define voxel size
        voxel_size=np.array([0.,0.,0.])
        for x in range(0,len(self.polymers),1):
            atoms=self.polymers[x].get_xyz()
            currentbox=np.max(atoms,axis=0)-np.min(atoms,axis=0)
            for k in range(0,3,1):
                if currentbox[k]>voxel_size[k]:
                    voxel_size[k]=currentbox[k]

        #use nanometers, increase voxel size by interchain_dist input
        voxel_size/=10.0
        voxel_size+=self.params.interchain_dist

        #create indicization (match name of polymer with index in list), and determine atomcount
        index_poly={}
        atomcount=0
        for i in range(0,len(self.polymers),1):
            name=self.polymers[i].molname
            index_poly[name]=i
            atomcount+=len(self.polymers[i].get_xyz())*len(self.systembox[self.systembox==name])


        ###STATISTICS###
        #iterate over every polymer in the system, and count the amount of individual elements
        monomer_count={}
        total_monomers=0
        average_length=0
        for i in np.unique(self.systembox):
            cnt=np.sum(self.systembox==i)
            chain_elements=np.array(list(self.polymers[index_poly[i]].chain))
            total_monomers+=cnt*len(chain_elements)
            
            for c in np.unique(chain_elements):
                element_increment=np.sum(chain_elements==c)*cnt
                if c not in list(monomer_count):
                    monomer_count[c]=element_increment
                else:
                    monomer_count[c]+=element_increment

        #print(statistics)
        self.logger.info("\n> monomers distribution in box:")
        for element in list(monomer_count):                    
            self.logger.info(">> %s: %s percent"%(element, 100.0*(float(monomer_count[element])/total_monomers)))

        avg_degree=float(total_monomers)/self.systembox.size
        self.logger.info("\n> number average degree of polymerization: %s"%avg_degree)


        ### GENERATE SYSTEM GRO FILE ###
        
        fout_index=open("%s/index_%s.ndx"%(mypath,self.params.output),'w')
        
        f_out=open("%s/%s.gro"%(mypath,self.params.output),'w')
        f_out.write("system\n")
        f_out.write("%s\n"%atomcount)

        #iterate over all box
        index=1 #atom counter
        minpos=[]
        maxpos=[]
        cntres=0
        
        #must sort the molecules according to their polymer name
        flatsystem=self.systembox.flatten()
        indices=flatsystem.argsort() #unique names of molecules, sorted alphabetically

        #count how many molecules per kind there are
        sortedflat=flatsystem[indices]
        uniquesorted=np.unique(sortedflat)
        contmols=[]
        for i in range(0,len(uniquesorted),1):
            name=uniquesorted[i]
            cnt=len(self.systembox[self.systembox==name])
            contmols.append([name,cnt])
   
        #counters needed for index file creation
        index_full=1 #atom counter without wrapping
        index_per_mol=1 #atom counter per molecule type (for pretty formatting)
        cont_mol=0 #count how many molecules of a specific kind passed
        cont_mol_kind=0 #count how many molecule kinds passed

        fout_index.write("[ %s ]\n"%contmols[0][0])
        for ndx in indices:

            #insert new molecule kind in index file
            cont_mol+=1
            if cont_mol>contmols[cont_mol_kind][1] and cont_mol_kind+1<len(contmols):
                cont_mol=1
                cont_mol_kind+=1
                fout_index.write("\n[ %s ]\n"%contmols[cont_mol_kind][0])
                index_per_mol=1 #atom counter per molecule type (for pretty formatting)
            
            pos=np.unravel_index(ndx,self.systembox.shape)
            x=pos[0]
            y=pos[1]
            z=pos[2]
                
            p=self.polymers[index_poly[self.systembox[x,y,z]]]
                          
            #compute polyhedron center
            pos=p.get_xyz()/10.0
            cntr=np.mean(pos,axis=0)/10.0
            
            #store min and max atom positions, for box size definition
            minpos.append(np.min(pos-cntr+voxel_size*np.array([x,y,z]),axis=0))
            maxpos.append(np.max(pos-cntr+voxel_size*np.array([x,y,z]),axis=0))
            
            #iterate over elements
            for j in range(0,len(p.poly),1):
                cntres+=1
                #wrap counting in case too many residues are present
                if cntres>99999:
                    cntres=0
                data_list=p.poly[j].mapping(p.poly[j].data)

                for i in range(0,len(data_list),1):
                    #create and write line in gromacs format (.gro file)
                    L='%5d%-5s%5s%5d%8.3f%8.3f%8.3f\n'%(cntres,data_list[i][2],data_list[i][1],index,\
													data_list[i][5]/10.0-cntr[0]+voxel_size[0]*x+voxel_size[0]/2.0,\
													data_list[i][6]/10.0-cntr[1]+voxel_size[1]*y+voxel_size[1]/2.0,\
													data_list[i][7]/10.0-cntr[2]+voxel_size[2]*z+voxel_size[2]/2.0)
                    f_out.write(L)

                    fout_index.write("%s "%index_full)
                    if np.mod(index_per_mol,15)==0:
                        fout_index.write("\n")

                    index+=1
                    index_full+=1
                    index_per_mol+=1
                    
                    #wrapping!
                    if index>99999:
                        index=0

        minbox=np.min(np.array(minpos),axis=0)
        maxbox=np.max(np.array(maxpos),axis=0)
        box=maxbox-minbox
        # account for intermolecular separation beyond periodic boundary condition
        box+=self.params.interchain_dist

        self.logger.info("\n> box size: %10.5f x %10.5f x %10.5f nm^3"%(box[0],box[1],box[2]))

        ### BOX INFORMATION TO ADD! ###
        f_out.write("%10.5f%10.5f%10.5f\n"%(box[0],box[1],box[2]))
        f_out.close()

        fout_index.write("\n")
        fout_index.close()

        ### CREATE TOP FILE ###
    
        f_out=open("%s/%s.top"%(mypath,self.params.output),'w')

        f_out.write("; generated with Assemble.py, by Matteo Degiacomi and Valentina Erastova, 2014-2018\n")

        f_out.write("\n[ defaults ]\n")
        f_out.write("; nbfunc        comb-rule       gen-pairs       fudgeLJ fudgeQQ\n")
        #f_out.write("{:<16}{:<16}{:<16}{:<8}{:<8}\n".format(self.ff.combination[0],self.ff.combination[1],self.ff.combination[2],self.ff.combination[3],self.ff.combination[4]))
        f_out.write("%16s%16s%16s%8s%8s\n"%(self.ff.combination[0],self.ff.combination[1],self.ff.combination[2],self.ff.combination[3],self.ff.combination[4]))


        f_out.write("\n[ atomtypes ]\n")
        f_out.write(";name   at.num  mass     charge  ptype  sigma     epsilon\n")
        for l in list(self.ff.nonbonded):
            line=self.ff.nonbonded[l]
            #f_out.write("{:<8}{:<8}{:<9}{:<8}{:<7}{:<10}{:<15}\n".format(l,line[0],line[1],line[2],line[3],line[4],line[5]))
            f_out.write("%8s%8s%9s%8s%7s%10s%15s\n"%(l,line[0],line[1],line[2],line[3],line[4],line[5]))                  

        for x in range(0,len(self.polymers),1):
            f_out.write("\n#include \"%s.itp\""%self.polymers[x].molname)
 
        f_out.write("\n\n[ system ]\n%s\n"%self.params.output)
        f_out.write("\n[ molecules ]\n")
        

        for i in range(0,len(contmols),1):
            f_out.write("%s %s\n"%(contmols[i][0],contmols[i][1]))
            
        f_out.close()
