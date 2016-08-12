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


import os, sys
import numpy as np

class Parser:

	parameters={}
	
	def __init__(self):

		#keyword_name, variable_name, type, default value
		self.add('mode','mode','str',"pdb")
		
		self.add('ForceField','ff','str',"NA")
		self.add('default_bond','default_bond','float',1.5)
		self.add('default_angle','default_angle','float',114)
		self.add('default_dihedral','default_dihedral','float',120)
		self.add('clash_threshold','clash_thresh','float',0.9)
		
		self.add('database','db','str',"")
		self.add('residue','residue','dictionary',{})
		
		self.add('molecule','molecule','array str',np.array([]))
		self.add('chain','chain','dictionary',{})
		self.add('composition','percentage','dictionary',{})
		self.add('length','length','dictionary',{})
		
		self.add('concentration','concentration','dictionary',{})		
		self.add('box_grid_shape','box_grid_shape','array int',np.array([0.0,0.0,0.0]))
				
		self.add('system_name','output','str',"system")
		self.add('output_folder','output_folder','str',".")

		self.set_default_values()


	#insert a new keyword entry in parameters dictionary
	def add(self,key,variable,vartype,default):
		self.parameters[key]=[variable,vartype,default]

	#set default values for all defined keywords
	def set_default_values(self):
		for k,v in self.parameters.iteritems():
			exec 'self.%s=v[2]'%v[0]

	#parse input file
	def parse(self,infile):

		f = open(infile, 'r+')
		line = f.readline()
		while line:
			w = line.split()

			if len(w) > 0 and str(w[0][0])!='#':

				#val=[variable_name,variable_type,default_value]
				try:
					val=self.parameters[w[0]]
				except KeyError:
					print "unrecognised keyword %s"%w[0]
					sys.exit(1)

				#if type is string
				if val[1].split()[0]=='str':
					exec 'self.%s=%s("%s")'%(val[0],val[1],w[1])
	
				#if type is int or float
				elif val[1].split()[0]=='int' or val[1].split()[0]=='float':
					exec 'self.%s=%s(%s)'%(val[0],val[1],w[1])
	
				#if type is an array of int, float, or str
				elif val[1].split()[0]=='array' and (val[1].split()[1]=='int' or val[1].split()[1]=='float' or val[1].split()[1]=='str'):
					exec 'self.%s=np.array(%s).astype(%s)'%(val[0],w[1:len(w)],val[1].split()[1])				

				elif val[1].split()[0]=='dictionary':
					exec 'self.%s[\"%s\"]=%s'%(val[0],w[1],w[2:])
				
				else:	
					print "unrecognised type for keyword %s: %s"%(w[0],val[1])
					sys.exit(1)

			line = f.readline()
		f.close()


	#verify standard variables consistency
	def check_variables(self):
		
		if self.db!="" and os.path.isfile(self.db)!=1 :
			print "ERROR: database file not found!"
			sys.exit(1)

		if self.mode!="gromacs" and self.mode!="pdb":
			print "ERROR: mode should be equal to pdb or gromacs!"
			sys.exit(1)

		if self.mode=="gromacs" and self.ff=="NA":
			print "ERROR: gromacs mode requires ForceField keyword!"
			sys.exit(1)
						
		if self.ff!="NA" and os.path.isfile(self.ff)!=1 :
			print "ERROR: forcefield file not found!"
			sys.exit(1)
		
		if len(self.molecule)==0:
			print "ERROR: no molecule name has been provided (keyword \"molecule\")!"
			sys.exit(1)
		
		for m in self.molecule:
			
			if m in self.chain:
			
				v=self.chain[m]
				self.chain[m]=v[0]
				
				if m in self.percentage:
					print "ERROR: chain undefined for molecule %s, must define monomer percentages!"%m
					sys.exit(1)			
				if m in self.length:
					print "ERROR: chain undefined for molecule %s, must define polymer length!"%m
					sys.exit(1)			
			else:
				if m in self.percentage:
					print "WARNING: chain defined for molecule %s, percentage keyword will be ignored!"%m
				if m in self.length:
					print "WARNING: chain defined for molecule %s, length keyword will be ignored!"	%m	

			if m in self.length:
				v=self.length[m]
				self.length[m]=int(v[0])


			#reformat percentages statement
			if m in self.percentage:
				if len(self.percentage[m])%2!=0:
					print "ERROR: incorrect number of arguments provided in percentage statement!"
	
				p=[]
				cnt=0
				for x in xrange(0,len(self.percentage[m])/2,1):
					try:
						p_f=float(self.percentage[m][x*2+1])
					except ValueError:
						print "ERROR: percentage keyword should be a list of one letter code followed by a float"
						return -1
					
					p.append([self.percentage[m][x*2],p_f])
					cnt+=p_f
				
				if cnt!=100:
					print "ERROR: sum of provided percentages equal to %s. Should be 100!"%cnt
					return -1	
				
				self.percentage[m]=p
				
			'''
			#test concentration statement
			if m in self.concentration:
				
				if len(self.concentration[m])%2!=0:
					print "ERROR: incorrect number of arguments provided in concentration statement!"
					return -1
			'''
		
		p2=[]
		cnt=0
		for m in self.concentration.keys():
			val=float(self.concentration[m][0])
			cnt+=val
			p2.append([m,val])
		
		#if cnt!=100:	
		#	print "ERROR: sum of provided percentages equal to %s. Should be 100!"%cnt
		#	return -1
			
		self.concentration=p2
		