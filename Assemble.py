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

from Database import *
from Polymer import *
from Parser import *
from ForceField import *
from System import *

import sys, os

import logging


def make_chain(length, percentage):
        
        logger=logging.getLogger('assemble')
        
        c=""
        #make roulette array
        roulette=[percentage[0][1]]
        for r in range(1,len(percentage),1):
            roulette.append(roulette[r-1]+percentage[r][1])
    
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
            c+=percentage[index][0]
            counter[index]+=1
    
        #chain completed, report produced percentages:
        counter/=float(length)
        counter*=100.0
        logger.info(">> chain: %s"%c)
        logger.info(">> monomer percentages:")
        for x in range(0,len(counter),1):
            logger.info(">> %s: %s percent"%(percentage[x][0],counter[x]))
            
        return c


def run(infile):
     
    if os.path.isfile(infile)!=1 :
        print("ERROR: setup file not found!")
        sys.exit(1)
    
    #parse input file
    params=Parser()
    params.set_default_values()
    params.parse(infile) #parse input file
    params.check_variables() #check consistency of defined variables
   
    #folder name equal path + system name
    folder="%s/%s"%(params.output_folder,params.output)

    #if output folder does not exist, create it
    if os.path.isdir(folder)!=1 :
        try:
            os.makedirs(folder)
        except:
            print("\n> could not create output folder %s!"%folder)
            return

    else:
        print("\n> WARNING: found folder %s. Existing files will be overwritten..."%folder)
   
    #prepare logger
    logname='%s/%s.log'%(folder,params.output)
    if os.path.isfile(logname):
        os.remove(logname)
    
    #@todo kill all previous handles and close logfiles before proceeding (in case an unexpected crash happened before)
    logger = logging.getLogger('assemble')
    logger.setLevel(logging.INFO)
    fh = logging.FileHandler(logname)
    ch = logging.StreamHandler()
    logger.addHandler(fh)
    logger.addHandler(ch)
    logger.info(">                 Assemble v1.0                     <")
    logger.info("> (c) 2014, Matteo Degiacomi and Valentina Erastova <")
    
    if params.mode=="gromacs":
        ff=ForceField()
        ff.load(params.ff)
        ff.default_bond=params.default_bond
        ff.default_angle=params.default_angle
        ff.default_dihedral=params.default_dihedral
    else:
        ff=""
    
    #load molecules objects into hash table from database list
    db=Database()
    if params.db!="":
        try:
            db.load(params.db,params.mode)
        except Exception as e:
            logger.exception(e)
            fh.close()
            logger.removeHandler(fh)
            logger.removeHandler(ch)
            sys.exit(1)
    
    if len(params.residue.keys())>0:
        logger.info("\n> adding residues in database...")
        
    for r in params.residue.keys():
        try:
            if params.mode=="pdb":
                db.add(r,params.residue[r])
            if params.mode=="gromacs":
                db.add(r,params.residue[r][0],params.residue[r][1])
        except:
            logger.exception(e)
            fh.close()
            logger.removeHandler(fh)
            logger.removeHandler(ch)

    
    #create random chains for all the molecules not having a chain explicitly defined
    for m in params.molecule:
        if m not in params.chain:
            logger.info("\n> randomizing polymer chain for molecule %s..."%m)
            params.chain[m]=make_chain(params.length[m],params.percentage[m])
   
    
    polymers=[]
    for m in params.molecule:        

        try:

            #generate polymer
            poly=Polymer(db,ff,m,params.mode,params.gromacs_nrxl)
            poly.clash_thresh=params.clash_thresh
            poly.make(params.chain[m])
        
            if params.mode=="pdb":
                poly.write_polymer(typef="pdb",mypath=folder)
            elif params.mode=="gromacs":
                poly.write_polymer(typef="gromacs",mypath=folder)
                poly.write_gromacs(mypath=folder)
            
            polymers.append(poly)
            
        except Exception as e:
            logger.exception(e)
            fh.close()
            logger.removeHandler(fh)
            logger.removeHandler(ch)
    
    if np.any(params.box_grid_shape==0):
        logger.warning("\n> no box size information provided, skipping system generation...")
        fh.close()
        logger.removeHandler(fh)
        logger.removeHandler(ch)
        return
    
    if params.mode=="gromacs":
        logger.info("\n> generating system...\n")
        system=System(polymers,ff,params)
        system.create_system(mypath=folder)
    
    
    logger.info("\n> All done! Thank you for using Assemble! <")
    fh.close()
    logger.removeHandler(fh)
    logger.removeHandler(ch)
    
if __name__=="__main__":

    cwd=os.getcwd()
    assembled=os.path.abspath(os.path.dirname(str(sys.argv[0])))
    os.environ["ASSEMBLEPATH"]="%s;%s"%(cwd,assembled)
    
    infile=str(sys.argv[1])
    run(infile)