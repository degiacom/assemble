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

from Molecule import *
import logging
import os

class Database(object):

    def __init__(self):
        self.molecules={}
        self.logger=logging.getLogger('assemble')


    def findfile(self, infile):

        #search path for relative paths:
        #1. working directory of database file
        #2. current working directory
        #3. assemble home directory
        mypath=os.environ["ASSEMBLEPATH"].split(";")

        #find pdbfile 
        fname = ""
        if os.path.isabs(infile): #absolute path, just check existence
            if os.path.isfile(infile):
                fname=infile
        else: #relative path, look for first hit in different folders
            for p in mypath:
                test="%s%s%s"%(p,os.sep,infile)
                if os.path.isfile(test):
                    fname=test
                    break

        return fname
    

    def load(self, infile, mode):
             
        files=[]
        
        self.logger.info("\n> Preparing molecules database...")
        try:
            f = open(infile, 'r+')
        except:
            raise IOError("database file %s not found!"%infile)
    
        #append folder of database file at the beginning of search path
        thisdir=os.path.dirname(os.path.abspath(infile))
        os.environ["ASSEMBLEPATH"]="%s;%s"%(thisdir, os.environ["ASSEMBLEPATH"])
        
        for line in f:  
            
            w=line.split()
       
            if len(w)==0: #to new line if line is empty
                continue
            
            if len(w)>0 and "#" in w[0]: #skip commented line
                continue
           
            if len(w[0])>1:
                raise IOError("found %s identifier in database file %s.\nOne letter code expected!"%(w[0],infile))
           
            fname=self.findfile(w[1])            
            if fname=="":
                raise IOError("PDB file %s not found for molecule %s"%(w[1],w[0]))
            
            try:
                self.logger.info(">> loading PDB %s"%fname)
                m=Molecule()
                m.import_pdb(fname, mode)
                
            except IOError:
                raise IOError("loading of PDB file %s failed for molecule %s"%(fname,w[0]))

            try:
                if mode=="gromacs":
                    if len(w)==3:
 
                        fname=self.findfile(w[2])            
                        if fname=="":
                            raise IOError("topology file %s not found for molecule %s"%(w[2],w[0]))

                        self.logger.info(">> loading topology %s"%fname)
                        m.import_topology(fname)
                    else:
                        raise IOError("in gromacs mode, a topology file is expected for molecule %s!"%w[0])

            except Exception as e:
                raise IOError("Could not load topology file %s for molecule %s!\n%s"%(w[2],w[0],e))

            if w[0] in list(self.molecules):
                self.logger.warning("\n> WARNING: duplicate key %s in database %s. Overwriting."%(w[0], infile))
                            
            self.molecules[w[0]]=m
            
        f.close()

            
    def add(self,code,pdb,topology=-1):
        
        if code in list(self.molecules):
            self.logger.warning("\n> WARNING: residue %s already present in database. Overwriting."%code)
        
        try:
            self.logger.info(">> loading PDB %s"%pdb)
            m=Molecule()
            m.import_pdb(pdb, "gromacs")
        except Exception as e:
            raise IOError(e)
        
        try:
            if topology!=-1:
                self.logger.info(">> loading topology %s"%topology)
                m.import_topology(topology)
        except Exception as e:
            raise IOError(e)
        
        self.molecules[code]=m
    
        
    def remove(self,code):
        try:
            del self.molecules[code]
        except IOError:
            raise IOError("ERROR: molecule %s not found, cannot remove!"%code)
        
        
    def save(self,filename):
        
        fout=open(filename,"w")
        
        for x in list(self.molecules):
            line="%s %s %s\n"%(x,self.molecules[x].pdbfile,self.molecules[x].topfile)
            fout.write(line)
            
        fout.close()