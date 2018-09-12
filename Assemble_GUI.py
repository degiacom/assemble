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


import os
import sys
import wx
from Database import Database
from ForceField import ForceField
import Assemble as A
from Parser import Parser
#from Polymer import Polymer


class MainWindow(wx.Frame):
    
    def __init__(self, parent, title):
        self.dirname=''
        # A "-1" in the size parameter instructs wxWidgets to use the default size.
        # In this case, we select 200px width and the default height.
        wx.Frame.__init__(self, parent, title=title, size=(700,600))
        self.CreateStatusBar() # A Status bar in the bottom of the window
        self.SetBackgroundColour("white")
        
        cwd=os.getcwd()
        assembled=os.path.abspath(os.path.dirname(str(sys.argv[0])))
        os.environ["ASSEMBLEPATH"]="%s;%s"%(cwd,assembled)
        
        #initialize datastructures
        self.db=Database()
        self.ff=ForceField()

        # Setting up the menu bar
        filemenu= wx.Menu()
        menuLoad = filemenu.Append(wx.ID_OPEN, "&Load..."," Load a setup file")
        menuSave = filemenu.Append(wx.ID_SAVEAS, "&Save as..."," Save the current setup")
        menuHelp = filemenu.Append(wx.ID_HELP, "&Help"," Optimization user manual")
        menuAbout= filemenu.Append(wx.ID_ABOUT, "&About"," Information about this program")
        filemenu.AppendSeparator()
        menuExit = filemenu.Append(wx.ID_EXIT,"E&xit"," Terminate the program")

        # Creating the menu bar
        menuBar = wx.MenuBar()
        menuBar.Append(filemenu,"&File") # Adding the "filemenu" to the MenuBar
        self.SetMenuBar(menuBar)  # Adding the MenuBar to the Frame content.
        

        #panels
        self.I=InputData(self)
        self.F=InputForceField(self)
        self.P=Polymers(self)
        self.S=SystemOut(self)

        self.make=wx.Button(self, label="MAKE!")

        box = wx.BoxSizer(wx.VERTICAL)
        box.Add(self.I, 1, wx.EXPAND | wx.ALL, 3 )
        box.Add(self.F, 0, wx.EXPAND | wx.ALL, 3 )
        box.Add(self.P, 1, wx.EXPAND | wx.ALL, 3 )
        box.Add(self.S, 0, wx.EXPAND | wx.ALL, 3 )
        box.Add(self.make, 0, wx.ALIGN_CENTER | wx.ALL, 3 )
        self.SetSizer(box)
        self.SetSize((800,700))
        self.Centre()

        # bind events to menu items
        self.Bind(wx.EVT_MENU, self.OnLoad, menuLoad)
        self.Bind(wx.EVT_MENU, self.OnSave, menuSave)
        self.Bind(wx.EVT_MENU, self.OnHelp, menuHelp)
        self.Bind(wx.EVT_MENU, self.OnAbout, menuAbout)
        self.Bind(wx.EVT_MENU, self.OnExit, menuExit)
        self.Bind(wx.EVT_BUTTON, self.OnMake,self.make)

        self.Show()

    def OnAbout(self,e):
        # Create a message dialog box
        dlg = wx.MessageDialog(self, " Assemble v1.1! \n (c) Matteo Degiacomi and Valentina Erastova\n 2014-2018", "About Assemble", wx.OK)
        dlg.ShowModal() # Shows it
        dlg.Destroy() # finally destroy it when finished.

    def OnExit(self,e):
        """ Close the interface"""
        self.Close(True)  # Close the frame.

    def OnSave(self,e):
        """ Save a file"""
        dlg = wx.FileDialog(self, "Chose a name for your setup file!", self.dirname, "", "*", wx.SAVE)
        if dlg.ShowModal() == wx.ID_OK:
            self.filename = dlg.GetFilename()
            self.dirname = dlg.GetDirectory()
            self.makeFile(os.path.join(self.dirname, self.filename))
        dlg.Destroy()

    def OnHelp(self,e):
        """ Open webpage """
        dlg = wx.MessageDialog(self, "Please see user manual for usage information. \n For bug reports, please contact the authors", "Assemble help", wx.OK)
        dlg.ShowModal() # Shows it
        dlg.Destroy() # finally destroy it when finished.

    def OnLoad(self,e):
        """ Open a file"""
        infile=-1
        dlg = wx.FileDialog(self, "Choose a file", self.dirname, "", "*", wx.FD_OPEN)
        if dlg.ShowModal() == wx.ID_OK:
            self.filename = dlg.GetFilename()
            self.dirname = dlg.GetDirectory()
            infile =os.path.join(self.dirname, self.filename)
                
        dlg.Destroy()

        #call file parsing
        if infile!=-1:
            if os.path.isfile(infile):
                self.parse(infile)
            else:
                self.errorPopup("File not found!", "ERROR!")

    def parse(self,infile):

        #WARNING: length and composition statements are currently ignored!

        try:
            P=Parser()
            P.set_default_values()
            P.parse(infile)
        except:
            self.errorPopup("Cannot load %s!"%(infile), "ERROR!")
            return

        if os.path.isfile(P.ff):
            self.F.setff.SetValue(P.ff)
        else:
            self.errorPopup("%s not found!"%P.ff, "WARNING!")
        
        try:
            if os.path.isfile(P.db):            
                self.I.setdatabase.SetValue(P.db)
                self.I.OnBLoad(None)
        except Exception as e:
            return
          
        pos=self.I.lc.GetItemCount()
        for r in list(P.residue):
            try:
                self.db.add(r,P.residue[r][0],P.residue[r][1])
            except Exception as e:
                return
 
            p=self.I.lc.InsertStringItem(pos,r)
            self.I.lc.SetStringItem(p,1,P.residue[r][0])
            self.I.lc.SetStringItem(p,2,P.residue[r][1])
            pos+=1



        self.S.setgridx.SetValue(str(P.box_grid_shape[0]))
        self.S.setgridy.SetValue(str(P.box_grid_shape[1]))
        self.S.setgridz.SetValue(str(P.box_grid_shape[2]))
        self.S.setsystem.SetValue(P.output)
        self.S.setpath.SetValue(P.output_folder)

        cnt=0
        for name in P.molecule:   
                        
            if name in list(P.chain):
                
                textperc=""
                cnt2={}
                for x in range(0,len(P.chain[name][0]),1):
                    if P.chain[name][0][x] in list(cnt2):
                        cnt2[P.chain[name][0][x]]=cnt2[P.chain[name][0][x]]+1
                    else:
                        cnt2[P.chain[name][0][x]]=1
                
                for x in list(cnt2):
                
                    if x not in list(self.db.molecules):
                        self.errorPopup("Chain in polymer %s features %s molecule,\n that is not contained in loaded database!"%(name,x), "ERROR!")
                
                    val=cnt2[x]/float(len(P.chain[name][0]))*100.0
                    textperc+="%s=%6.2f, "%(x,val)

                self.P.lc.InsertStringItem(cnt, str(cnt+1))
                self.P.lc.SetStringItem(cnt, 0, name) 
                self.P.lc.SetStringItem(cnt, 1, P.chain[name][0])
                self.P.lc.SetStringItem(cnt, 2, str(len(P.chain[name][0])))
                self.P.lc.SetStringItem(cnt, 3, textperc)
            
                if name in list(P.concentration):
                    try: #in case concentration value is not given
                        self.P.lc.SetStringItem(cnt, 4, P.concentration[name][0])
                    except:
                        continue
                        
        cnt+=1
        

    def OnMake(self,event):

        #dbname="db_tmp"
        #try:
        #    self.db.save(dbname)
        #except:
        #    self.errorPopup("could not prepare database!", "ERROR!")
        #    return
          
        #check that concentrations are 100%. If one element only, force 100%, else check sum
        if self.P.lc.GetItemCount()==1:
            self.P.lc.SetStringItem(0,4,str(100.0))  
        
        
        else:
            conc=0
            for i in range(0,self.P.lc.GetItemCount(),1):
                c=self.P.lc.GetItem(itemIdx=i,col=4).GetText()
                try:
                    conc+=float(c)
                except:
                    self.errorPopup("concentration of polymer %s must be\na number between 0 and 100!"%self.P.lc.GetItem(itemIdx=i,col=0).GetText(), "ERROR!")
                    return  
            
            #if conc!=100:
                #self.errorPopup("sum of polymer concentrations is equal to %s (expected 100)\nTreating concentrations as ratios!"%conc, "WARNING!")                             
                ##self.errorPopup("sum of polymer concentrations\n should be equal to 100 percent!", "ERROR!")
                ##return            
             
        #check box size information consistency                      
        try:
            int(self.S.setgridx.GetValue())
        except:
            self.errorPopup("box x size must be an integer number!", "ERROR!")      
            return          

        try:
            int(self.S.setgridy.GetValue())
        except:
            self.errorPopup("box y size must be an integer number!", "ERROR!")
            return    
            
        try:
            int(self.S.setgridz.GetValue())
        except:
            self.errorPopup("box z size must be an integer number!", "ERROR!")
            return
                    
        if int(self.S.setgridx.GetValue())==0 or int(self.S.setgridy.GetValue())==0 or int(self.S.setgridx.GetValue())==0:
            self.errorPopup("box information contains zero values. System generation will be skipped...", "WARNING!")     


        fname_in = "input_tmp"
        dirname = os.path.dirname(os.path.realpath(__file__))
        fname = os.path.join(dirname, fname_in)
 
        try:
            self.makeFile(fname)
        except Exception as e:
            self.errorPopup("problem in input file!\n"%e, "ERROR!")
            #os.remove(dbname)
            os.remove(fname)
            return

        try:
            A.run(fname)
        except Exception as e:
            self.errorPopup("polymer generation failed!\n%s"%e, "ERROR!")
            #os.remove(dbname)
            os.remove(fname)
            return
            
        #os.remove(dbname)
        os.remove(fname)
        
        self.errorPopup("Files generation complete!", "SUCCESS!")

            
    def makeFile(self,fname):
        
        fout=open(fname,"w")
        
        fout.write("mode gromacs\n")
        
        #try:
        #    fout.write("database %s\n"%dbname)
        #except:
        #    self.errorPopup("could not prepare database!", "ERROR!")

        for i in range(0,self.I.lc.GetItemCount(),1):
            item=self.I.lc.GetItem(itemIdx=i,col=0)
            pdb=self.I.lc.GetItem(itemIdx=i,col=1)
            top=self.I.lc.GetItem(itemIdx=i,col=2)
            
            fout.write("residue %s %s %s\n"%(item.GetText(),pdb.GetText(),top.GetText()))


        ff=self.F.setff.GetValue().replace('\\','/')
        if os.path.isfile(ff):
            fout.write("ForceField %s\n"%ff)
        else:
            self.errorPopup("force field file not found!", "ERROR!")


        #print(self.S.setpath.GetValue())
        outfolder=self.S.setpath.GetValue().replace('\\','/')
        if len(outfolder)==0:
            outfolder="."
        fout.write("output_folder %s\n"%outfolder)
        
        fout.write("system_name %s\n"%self.S.setsystem.GetValue())
        fout.write("box_grid_shape %s %s %s\n"%(self.S.setgridx.GetValue(),self.S.setgridy.GetValue(),self.S.setgridz.GetValue()))
        
        molecules=[]
        for i in range(0,self.P.lc.GetItemCount(),1):
            item=self.P.lc.GetItem(itemIdx=i,col=0)
            chain=self.P.lc.GetItem(itemIdx=i,col=1)
            conc=self.P.lc.GetItem(itemIdx=i,col=4)
            
            fout.write("chain %s %s\n"%(item.GetText(),chain.GetText()))
            fout.write("concentration %s %s\n"%(item.GetText(),conc.GetText()))
            molecules.append(str(item.GetText()))

        fout.write("molecule %s"%' '.join(molecules))

        fout.close()
    
    def errorPopup(self,msg,title):
        """ Popup error message"""
        print(msg)
        dlg = wx.MessageDialog(self, msg, title, wx.OK)
        dlg.ShowModal() # Shows it
        dlg.Destroy() # finally destroy it when finished.


class InputData(wx.Panel):
    def __init__(self, parent):

        wx.Panel.__init__(self, parent,style=wx.SIMPLE_BORDER|wx.TAB_TRAVERSAL)

        self.dirname=''
        self.parent=parent

        #structure line
        #topology line
        self.lbldatabase=wx.StaticText(self, label="database:")
        self.setdatabase=wx.TextCtrl(self, value="")
        self.databasebutton=wx.Button(self, label="Select")

        #list of database entries
        self.lc = wx.ListCtrl(self, -1, style=wx.LC_REPORT)
        self.lc.InsertColumn(0, 'name',)
        self.lc.InsertColumn(1, 'PDB')
        self.lc.InsertColumn(2, 'topology')
        self.lc.SetColumnWidth(0, 50)
        self.lc.SetColumnWidth(1, 240)
        self.lc.SetColumnWidth(2, 240)

        #boundary conditions buttons box
        self.buttonsbox = wx.BoxSizer(wx.VERTICAL)
        self.buttonsbox.Add(wx.Button(self, 1, "Load"), 0, wx.EXPAND )
        self.buttonsbox.Add(wx.Button(self, 2, "Add"), 0, wx.EXPAND )
        self.buttonsbox.Add(wx.Button(self, 3, "Edit"), 0, wx.EXPAND )
        self.buttonsbox.Add(wx.Button(self, 4, "Remove"), 0, wx.EXPAND )
        self.buttonsbox.Add(wx.Button(self, 5, "Save"), 0, wx.EXPAND )


        gs1 = wx.FlexGridSizer(2, 3, 3, 5)
        #gs1.SetFlexibleDirection(wx.HORIZONTAL)
        gs1.AddGrowableCol(1,1)
        gs1.AddGrowableRow(1,1)
        gs1.AddMany([(self.lbldatabase, 0, wx.ALIGN_CENTER_VERTICAL | wx.ALIGN_RIGHT,3),
                    (self.setdatabase, 0, wx.ALIGN_LEFT | wx.EXPAND, 3),
                    (self.databasebutton, 0, wx.ALIGN_LEFT, 3),
                    (wx.Size(1,1), 0, wx.ALIGN_LEFT, 3),
                    (self.lc, 1,wx.ALIGN_LEFT | wx.EXPAND, 3),
                    (self.buttonsbox, 0, wx.ALIGN_LEFT, 3)])

        self.Bind(wx.EVT_BUTTON, self.OnBSelect,self.databasebutton)
        self.Bind(wx.EVT_BUTTON, self.OnBLoad,id=1)
        self.Bind(wx.EVT_BUTTON, self.OnBAdd,id=2)
        self.Bind(wx.EVT_BUTTON, self.OnBEdit,id=3)
        self.Bind(wx.EVT_BUTTON, self.OnBRemove,id=4)
        self.Bind(wx.EVT_BUTTON, self.OnBSave,id=5)

        self.Bind(wx.EVT_LIST_ITEM_ACTIVATED, self.OnBEdit,self.lc)
        
        box = wx.BoxSizer(wx.VERTICAL)
        box.Add(gs1, 1, wx.EXPAND | wx.ALL, 3 )
        self.SetSizer(box,1)
        self.SetSize((400,400))
        self.Centre()
        
    def OnBSelect(self,event):
        dlg = wx.FileDialog(self, "Load a database file", self.dirname, "", "*.*", wx.FD_OPEN)
        if dlg.ShowModal() == wx.ID_OK:
            self.databasefile = dlg.GetFilename()
            self.databasedir = dlg.GetDirectory()
            self.setdatabase.SetValue(os.path.join(self.databasedir, self.databasefile))
        dlg.Destroy()

    def OnBLoad(self,event):

        #if hasattr(self, 'databasefile') and hasattr(self, 'databasedir'):
        if self.setdatabase.GetValue()!="":    
            #clear the current database, and load a new one
            if self.lc.GetItemCount()>0:
                self.lc.DeleteAllItems()
                for x in list(self.parent.db.molecules):
                    self.parent.db.remove(x) 
            
            try:
                self.parent.db.load(self.setdatabase.GetValue(), "gromacs")
            except IOError as e:
                error="Error in database file loading!\n%s"%e
                self.parent.errorPopup(error,"ERROR!")
                raise IOError(error)
            
            pos=0    
            lines=sorted(list(self.parent.db.molecules))
            for x in lines:      
                filename=self.parent.db.molecules[x].pdbfile
                topologyfile=self.parent.db.molecules[x].topfile
                p=self.lc.InsertStringItem(pos,x)
                self.lc.SetStringItem(p,1,filename)
                self.lc.SetStringItem(p,2,topologyfile)
                pos+=1

        else:
            self.parent.errorPopup("No database file has been specified!","ERROR!")
            raise IOError("")
            

    def OnBSave(self,event):
        dlg = wx.FileDialog(self, "Save a database file", self.dirname, "", "*.*",  wx.SAVE)
        if dlg.ShowModal() == wx.ID_OK:
            self.filename = dlg.GetFilename()
            self.dirname = dlg.GetDirectory()
            self.setdatabase.SetValue(os.path.join(self.dirname, self.filename))
            try:
                self.parent.db.save(os.path.join(self.dirname, self.filename))
            except:
                self.parent.errorPopup("Could not save file!","ERROR!")
        dlg.Destroy()

    def OnBAdd(self,event):

        dia=DatabaseEditor(self)
        dia.ShowModal()
        dia.Destroy()

    def OnBEdit(self,event):
        
        index = self.lc.GetFocusedItem()

        params=[]
        params.append(index)
        if index != -1 :
            params.append(self.lc.GetItem(index, 0).GetText())
            params.append(self.lc.GetItem(index, 1).GetText())
            params.append(self.lc.GetItem(index, 2).GetText())

            dia=DatabaseEditor(self,params)
            dia.ShowModal()
            dia.Destroy()

    def OnBRemove(self,event):
        index = self.lc.GetFocusedItem()
        if index!=-1:
            name=str(self.lc.GetItem(index, 0).GetText())
            try:
                self.parent.db.remove(name)
            except:
                self.parent.errorPopup("Cannot remove molecule... umm... weird...","ERROR!")
                #print("Cannot remove molecule... umm... weird...")
                return
            
            self.lc.DeleteItem(index)


class InputForceField(wx.Panel):
    def __init__(self, parent):

        wx.Panel.__init__(self, parent,style=wx.SIMPLE_BORDER|wx.TAB_TRAVERSAL)

        self.dirname=''
        self.parent=parent

        #topology line
        self.lblff=wx.StaticText(self, label="force field:")
        self.setff=wx.TextCtrl(self, value="")
        self.ffbutton=wx.Button(self, label="Select")
        self.Bind(wx.EVT_BUTTON, self.OnClickSelectff,self.ffbutton)

        gs1 = wx.FlexGridSizer(1, 3, 3, 5)
        #gs1.SetFlexibleDirection(wx.HORIZONTAL)
        gs1.AddGrowableCol(1,1)
        #gs1.AddGrowableRow(1,1)
        gs1.AddMany([(self.lblff, 0, wx.ALIGN_CENTER_VERTICAL | wx.ALIGN_RIGHT, 3),
                    (self.setff, 0, wx.ALIGN_LEFT | wx.EXPAND, 3),
                    (self.ffbutton, 0, wx.ALIGN_LEFT, 3)])
    
        box = wx.BoxSizer(wx.VERTICAL)
        box.Add(gs1, 1, wx.EXPAND | wx.ALL, 3 )
        self.SetSizer(box,1)
        self.SetSize((400,400))
        self.Centre()


    def OnClickSelectff(self,event):
        dlg = wx.FileDialog(self, "Choose Force Field file", self.dirname, "", "*.*", wx.FD_OPEN)
        if dlg.ShowModal() == wx.ID_OK:
            self.filename = dlg.GetFilename()
            self.dirname = dlg.GetDirectory()
            self.setff.SetValue(os.path.join(self.dirname, self.filename))
        dlg.Destroy()


class DatabaseEditor(wx.Dialog):
    def __init__(self, parent, params=-1):
        wx.Dialog.__init__(self, parent,title="Database Editor",style=wx.RESIZE_BORDER|wx.DEFAULT_DIALOG_STYLE,size=(300,150))

        self.SetSizeHints(300,200, maxH = 300)

        self.dirname=parent.dirname
        self.params=params
        self.parent=parent

        self.lblname = wx.StaticText(self, label="name:")
        self.editname = wx.TextCtrl(self, value="", size=(40,-1))

        #pdb file
        self.lblpdb=wx.StaticText(self, label="PDB:")
        self.setpdb=wx.TextCtrl(self, value="")
        self.pdbbutton=wx.Button(self, label="Select")
        self.Bind(wx.EVT_BUTTON, self.OnClickSelectPDB,self.pdbbutton)

        #topology file
        self.lbltop=wx.StaticText(self, label="topology:")
        self.settop=wx.TextCtrl(self, value="")
        self.topbutton=wx.Button(self, label="Select")
        self.Bind(wx.EVT_BUTTON, self.OnClickSelectTop,self.topbutton)

        gs = wx.FlexGridSizer(3, 3, 3, 3)
        gs.AddGrowableCol(1,1)
        gs.AddMany([(self.lblname, 0,  wx.ALIGN_LEFT | wx.ALIGN_CENTER_VERTICAL, 3),
                    (self.editname, 1, wx.ALIGN_LEFT | wx.ALIGN_CENTER_VERTICAL, 3),
                    (wx.Size(1,1), 0, wx.ALIGN_LEFT, 3),
                    (self.lblpdb, 0,  wx.ALIGN_LEFT | wx.ALIGN_CENTER_VERTICAL, 3),
                    (self.setpdb, 1, wx.EXPAND | wx.ALIGN_LEFT | wx.ALIGN_CENTER_VERTICAL, 3),
                    (self.pdbbutton, 0,  wx.ALIGN_RIGHT | wx.ALIGN_CENTER_VERTICAL, 3),
                    (self.lbltop, 0,  wx.ALIGN_LEFT | wx.ALIGN_CENTER_VERTICAL, 3),
                    (self.settop, 1, wx.EXPAND | wx.ALIGN_LEFT | wx.ALIGN_CENTER_VERTICAL, 3),
                    (self.topbutton, 0,  wx.ALIGN_RIGHT | wx.ALIGN_CENTER_VERTICAL, 3)])

        if params!=-1 and len(params)==4:
            self.editname.SetValue(self.params[1])
            self.setpdb.SetValue(self.params[2])
            self.settop.SetValue(self.params[3])

        self.okbutton=wx.Button(self, label="OK!")
        self.Bind(wx.EVT_BUTTON, self.OnOk, self.okbutton,id=wx.ID_OK)
        self.cancel =wx.Button(self, label="Cancel")
        self.Bind(wx.EVT_BUTTON, self.OnCancel,self.cancel,id=wx.ID_CANCEL)
        self.buttonbox = wx.BoxSizer(wx.HORIZONTAL)
        self.buttonbox.Add(self.okbutton, 0, wx.RIGHT, 3 )
        self.buttonbox.Add(self.cancel, 0, wx.LEFT, 3 )

        #draw title, and grid in a vertical box
        box = wx.BoxSizer(wx.VERTICAL)
        box.Add(gs, 0, wx.ALL | wx.EXPAND, 3 )
        box.Add(self.buttonbox, 1, wx.CENTER, 3 )
        self.SetSizer(box)
        self.Centre()


    def OnClickSelectPDB(self,event):
        dlg = wx.FileDialog(self, "Choose PDB file", self.dirname, "", "*.*", wx.FD_OPEN)
        if dlg.ShowModal() == wx.ID_OK:
            self.filename = dlg.GetFilename()
            self.dirname = dlg.GetDirectory()
            self.setpdb.SetValue(os.path.join(self.dirname, self.filename))
        dlg.Destroy()


    def OnClickSelectTop(self,event):
        dlg = wx.FileDialog(self, "Choose topology file", self.dirname, "", "*.*", wx.FD_OPEN)
        if dlg.ShowModal() == wx.ID_OK:
            self.filename = dlg.GetFilename()
            self.dirname = dlg.GetDirectory()
            self.settop.SetValue(os.path.join(self.dirname, self.filename))
        dlg.Destroy()


    def OnOk(self,e):

        #if editing remove first from database the past entry
        if self.params!=-1:
            num_items=self.params[0]
            name=str(self.parent.lc.GetItem(num_items, 0).GetText())
            self.parent.parent.db.remove(name)
            self.parent.lc.DeleteItem(num_items)
        
        name=str(self.editname.GetValue())
        pdb=str(self.setpdb.GetValue())
        top=str(self.settop.GetValue())
        
        if os.path.isfile(pdb)!=1:
            self.parent.parent.errorPopup("PDB file not found!","ERROR!")
            return            
        
        if os.path.isfile(top)!=1:
            self.parent.parent.errorPopup("topology file not found!","ERROR!")
            return             
        
        if name=="":
            self.parent.parent.errorPopup("a residue name must be provided!","ERROR!")
            return             
        
        if len(name)>1:
            self.parent.parent.errorPopup("The residue name should be a unique character\nNote: the database is case sensitive!","ERROR!")
            return             
        
        
        residues=[]
        for i in range(0,self.parent.lc.GetItemCount(),1):
            item=self.parent.lc.GetItem(itemIdx=i,col=0)
            residues.append(str(item.GetText()))
        if name in residues:
            self.parent.parent.errorPopup("The desired residue name already exists in the database!\nPlease choose a unique name!","ERROR!")
            return            
        
        try:
            self.parent.parent.db.add(name,pdb,top)
        except Exception as e:
            self.parent.parent.errorPopup("Cannot add residue:\n%s!"%e,"ERROR!")
            return
            
        num_items = self.parent.lc.GetItemCount()
        self.parent.lc.InsertItem(num_items,name)    
        self.parent.lc.SetItem(num_items, 1, pdb)
        self.parent.lc.SetItem(num_items, 2, top) 

        self.Close(True)  # Close the frame.


    def OnCancel(self,e):
        self.Close(True)  # Close the frame.


class Polymers(wx.Panel):
    def __init__(self, parent):

        self.dirname=''
        self.parent=parent

        wx.Panel.__init__(self, parent,style=wx.SIMPLE_BORDER|wx.TAB_TRAVERSAL)

        #self.font1 = wx.Font(10, wx.NORMAL, wx.ITALIC, wx.BOLD)
        self.font1 = wx.Font(10, wx.FONTFAMILY_DEFAULT, wx.FONTSTYLE_NORMAL, wx.BOLD)

        self.title = wx.StaticText(self, label="POLYMERS")
        self.title.SetFont(self.font1)

        #list of polymers
        self.lc = wx.ListCtrl(self, -1, style=wx.LC_REPORT)
        self.lc.InsertColumn(0, 'name')
        self.lc.InsertColumn(1, 'chain')
        self.lc.InsertColumn(2, 'length')
        self.lc.InsertColumn(3, 'components (%)')
        self.lc.InsertColumn(4, 'concentr.')
        self.lc.SetColumnWidth(0, 65)
        self.lc.SetColumnWidth(1, 230)
        self.lc.SetColumnWidth(2, 60)
        self.lc.SetColumnWidth(3, 150)
        self.lc.SetColumnWidth(4, 80)

        #boundary conditions buttons box
        self.buttonsbox = wx.BoxSizer(wx.VERTICAL)
        self.buttonsbox.Add(wx.Button(self, 1, "Add"), 0 , wx.EXPAND)
        self.buttonsbox.Add(wx.Button(self, 2, "Edit"), 0 , wx.EXPAND)
        self.buttonsbox.Add(wx.Button(self, 3, "Remove"), 0 , wx.EXPAND)

        #box the group
        gs = wx.FlexGridSizer(2, 3, 3, 5)
        #gs.SetFlexibleDirection(wx.HORIZONTAL)
        gs.AddGrowableCol(0,1)
        gs.AddGrowableRow(0,1)
        gs.AddMany([(self.lc, 1, wx.ALIGN_LEFT | wx.EXPAND, 3),
                    (self.buttonsbox, 0, wx.ALIGN_TOP | wx.ALIGN_LEFT, 3)])

        #draw title, boundaries and buttons list
        boundbox = wx.BoxSizer(wx.VERTICAL)
        boundbox.Add(self.title, 0, wx.EXPAND | wx.ALL, 3 )
        boundbox.Add(gs, 1, wx.EXPAND | wx.ALL, 3 )
        self.SetSizer(boundbox,1)
        self.SetSize((400,400))
        self.Centre()

        #buttons bindings
        self.Bind(wx.EVT_BUTTON, self.OnBadd,id=1)
        self.Bind(wx.EVT_BUTTON, self.OnBedit,id=2)
        self.Bind(wx.EVT_BUTTON, self.OnBremove,id=3)

        self.Bind(wx.EVT_LIST_ITEM_ACTIVATED, self.OnBedit,self.lc)
        

    def OnBadd(self,event):
        
        #add error message
        if self.parent.I.lc.GetItemCount()>0:
            dia=PolymerEditor(self)
            dia.ShowModal()
            dia.Destroy()
        else:
            self.parent.errorPopup("No residue is loaded\nin the database!","ERROR!")
            
            
    def OnBedit(self,event):

        index = self.lc.GetFocusedItem()

        params=[]
        params.append(index)
        if index != -1 :
            params.append(str(self.lc.GetItem(index, 0).GetText()))
            params.append(str(self.lc.GetItem(index, 1).GetText()))
            params.append(str(self.lc.GetItem(index, 2).GetText()))
            params.append(str(self.lc.GetItem(index, 3).GetText()))
            params.append(str(self.lc.GetItem(index, 4).GetText()))

            dia=PolymerEditor(self,params)
            dia.ShowModal()
            dia.Destroy()

    def OnBremove(self,event):
        
        index = self.lc.GetFocusedItem()
        if index!=-1:
            self.lc.DeleteItem(index)

            #renumber items after deletion
            num_items = self.lc.GetItemCount()
            for x in range(index,num_items,1):
                self.lc.SetStringItem(x, 0, str(x+1))


class SystemOut(wx.Panel):
    def __init__(self, parent):

        self.dirname=''
        self.parent=parent

        wx.Panel.__init__(self, parent,style=wx.SIMPLE_BORDER|wx.TAB_TRAVERSAL)

        self.font1 = wx.Font(10, wx.FONTFAMILY_DEFAULT, wx.FONTSTYLE_NORMAL, wx.BOLD)

        self.lblgrid=wx.StaticText(self, label=" grid x,y,z: ")
        self.setgridx=wx.TextCtrl(self, value="0", size=(70,-1)) # set size here!
        self.setgridy=wx.TextCtrl(self, value="0", size=(70,-1)) # set size here!
        self.setgridz=wx.TextCtrl(self, value="0", size=(70,-1)) # set size here!
        self.lblsystem=wx.StaticText(self, label="system name: ")
        self.setsystem=wx.TextCtrl(self, value="system")

        #output path
        self.lblpath=wx.StaticText(self, label="output folder:")
        self.setpath=wx.TextCtrl(self, value=os.getcwd())
        self.pathbutton=wx.Button(self, label="Select")
        self.Bind(wx.EVT_BUTTON, self.OnClickSelectpath,self.pathbutton)

        buttons = wx.BoxSizer(wx.HORIZONTAL)
        buttons.Add(self.setgridx, 0, wx.ALIGN_CENTER_VERTICAL |wx.ALIGN_LEFT, 1 )
        buttons.Add(self.setgridy, 0, wx.ALIGN_CENTER_VERTICAL |wx.ALIGN_LEFT, 1 )
        buttons.Add(self.setgridz, 0, wx.ALIGN_CENTER_VERTICAL |wx.ALIGN_LEFT, 1 )
        
        buttons.Add(wx.Size(25,25), 0, wx.ALIGN_CENTER_VERTICAL |wx.ALIGN_LEFT, 3)
        buttons.Add(self.lblsystem, 0, wx.ALIGN_CENTER_VERTICAL | wx.ALIGN_RIGHT, 3  )
        buttons.Add(self.setsystem, 1, wx.ALIGN_CENTER_VERTICAL |wx.ALIGN_LEFT, 3 )

        gs1 = wx.FlexGridSizer(2, 3, 3, 5)
        gs1.AddGrowableCol(1,1)
        gs1.AddMany([(self.lblgrid, 0, wx.ALIGN_CENTER_VERTICAL | wx.ALIGN_RIGHT, 3 ),
                    (buttons, 0, wx.ALIGN_LEFT | wx.EXPAND, 3),
                    (wx.Size(1,1), 0, wx.ALIGN_CENTER_VERTICAL |wx.ALIGN_LEFT, 3),
                    (self.lblpath, 0, wx.ALIGN_CENTER_VERTICAL | wx.ALIGN_RIGHT, 3),
                    (self.setpath, 0, wx.ALIGN_LEFT | wx.EXPAND, 3),
                    (self.pathbutton, 0, wx.ALIGN_LEFT, 3)])
    
  
        box = wx.BoxSizer(wx.VERTICAL)
        box.Add(gs1, 1, wx.EXPAND | wx.ALL, 3 )
        #box.Add(buttons, 1, wx.EXPAND | wx.ALL, 3 )
        self.SetSizer(box,1)
        #self.SetSizer(gs1)
        self.SetSize((700,600))
        self.Centre()



    def OnClickSelectpath(self,event):
        dlg = wx.DirDialog(self, "Choose path for output folder", self.dirname, style=wx.DD_DEFAULT_STYLE | wx.DD_NEW_DIR_BUTTON)
        if dlg.ShowModal() == wx.ID_OK:
            self.dirname = dlg.GetPath()
            self.setpath.SetValue(self.dirname)
        dlg.Destroy()


class PolymerEditor(wx.Dialog):
    def __init__(self, parent,params=-1):
        wx.Dialog.__init__(self, parent,title="Polymer Editor",style=wx.RESIZE_BORDER|wx.DEFAULT_DIALOG_STYLE|wx.TAB_TRAVERSAL,size=(400,500))

        self.SetSizeHints(400,500)

        self.params=params
        self.parent=parent

        #extract names of available molecules from database
        self.molecules=[]
        for i in range(0,self.parent.parent.I.lc.GetItemCount(),1):
            item=self.parent.parent.I.lc.GetItem(itemIdx=i, col=0)
            self.molecules.append(str(item.GetText()))

        self.lblname = wx.StaticText(self, label="name:")
        self.editname = wx.TextCtrl(self, value="", size=(40,-1))
        self.lblsource = wx.StaticText(self, label="source:")
        self.sourceList = ['chain','percentage']
        self.editsource = wx.ComboBox(self, value=self.sourceList[0], size=(80, -1), choices=self.sourceList, style=wx.CB_DROPDOWN)
        self.Bind(wx.EVT_COMBOBOX, self.OnSourceCombobox,self.editsource)

        gs0 = wx.FlexGridSizer(1, 4, 3, 3)
        gs0.AddGrowableCol(1,1)
        gs0.AddGrowableCol(2,1)
        gs0.AddMany([(self.lblname, 0,  wx.ALIGN_LEFT | wx.ALIGN_CENTER_VERTICAL, 3),
                    (self.editname, 1, wx.EXPAND | wx.ALIGN_LEFT | wx.ALIGN_CENTER_VERTICAL, 3),
                    #(wx.Size(1,1), 0, wx.ALIGN_LEFT, 3),    
                    (self.lblsource, 0,  wx.ALIGN_RIGHT | wx.ALIGN_CENTER_VERTICAL, 3),
                    (self.editsource, 0, wx.ALIGN_RIGHT | wx.ALIGN_CENTER_VERTICAL, 3)])

        #list of polymers
        self.lblconcentration = wx.StaticText(self, label="percent.:")
        self.lc = wx.ListCtrl(self, -1, style=wx.LC_REPORT)
        self.lc.InsertColumn(0, 'res.')
        self.lc.InsertColumn(1, '%')
        self.lc.SetColumnWidth(0, 50)
        self.lc.SetColumnWidth(1, 50)
        self.lbllength = wx.StaticText(self, label="length:")
        self.editlength = wx.TextCtrl(self, value="", size=(40,-1))
     
        #boundary conditions buttons box
        self.add=wx.Button(self, 1, "Add")
        self.edit=wx.Button(self, 2, "Edit")
        self.remove=wx.Button(self, 3, "Remove")
        self.generate=wx.Button(self, 4, "Generate")
        self.buttonsbox = wx.BoxSizer(wx.VERTICAL) 
        self.buttonsbox.Add(self.add, 0 , wx.EXPAND)
        self.buttonsbox.Add(self.edit, 0 , wx.EXPAND)
        self.buttonsbox.Add(self.remove, 0 , wx.EXPAND)
        self.buttonsbox.Add(self.generate, 0 , wx.EXPAND)

        #buttons bindings for chain generator
        self.Bind(wx.EVT_BUTTON, self.OnBadd,self.add)
        self.Bind(wx.EVT_BUTTON, self.OnBedit,self.edit)
        self.Bind(wx.EVT_BUTTON, self.OnBremove,self.remove)
        self.Bind(wx.EVT_BUTTON, self.OnGenerate,self.generate)

        self.Bind(wx.EVT_LIST_ITEM_ACTIVATED, self.OnBedit,self.lc)

        gs = wx.FlexGridSizer(2, 3, 3, 3)
        gs.AddGrowableCol(1,1)
        gs.AddGrowableRow(1,1)
        gs.AddMany([(self.lbllength, 0,  wx.ALIGN_LEFT | wx.ALIGN_CENTER_VERTICAL, 3),
                    (self.editlength, 0, wx.ALIGN_LEFT | wx.ALIGN_CENTER_VERTICAL, 3),
                    (wx.Size(1,1), 0, wx.ALIGN_LEFT, 3),
                    (self.lblconcentration, 0,  wx.ALIGN_LEFT | wx.ALIGN_CENTER_VERTICAL, 3),  
                    (self.lc, 0, wx.EXPAND | wx.ALIGN_LEFT | wx.ALIGN_CENTER_VERTICAL, 3),
                    (self.buttonsbox, 1,  wx.ALIGN_RIGHT | wx.ALIGN_TOP, 3)])

        self.lblchain = wx.StaticText(self, label="chain:")
        self.editchain = wx.TextCtrl(self, value="")

        self.lblconc = wx.StaticText(self, label="concentr.:")
        self.editconc = wx.TextCtrl(self, value="", size=(40,-1))
     

        gs1 = wx.FlexGridSizer(2, 2, 3, 3)
        gs1.AddGrowableCol(1,1)
        gs1.AddMany([(self.lblchain, 0,  wx.ALIGN_LEFT | wx.ALIGN_CENTER_VERTICAL, 3),
                    (self.editchain, 0, wx.EXPAND | wx.ALIGN_LEFT | wx.ALIGN_CENTER_VERTICAL, 3),
                    (self.lblconc, 0,  wx.ALIGN_LEFT | wx.ALIGN_CENTER_VERTICAL, 3),
                    (self.editconc, 0, wx.EXPAND | wx.ALIGN_LEFT | wx.ALIGN_CENTER_VERTICAL, 3)])

        self.okbutton=wx.Button(self, label="OK!")
        self.Bind(wx.EVT_BUTTON, self.OnOk, self.okbutton,id=wx.ID_OK)
        self.cancel =wx.Button(self, label="Cancel")
        self.Bind(wx.EVT_BUTTON, self.OnCancel,self.cancel,id=wx.ID_CANCEL)
        self.buttonbox = wx.BoxSizer(wx.HORIZONTAL)
        self.buttonbox.Add(self.okbutton, 0, wx.RIGHT, 3 )
        self.buttonbox.Add(self.cancel, 0, wx.LEFT, 3 )

        #draw title, and grid in a vertical box
        box = wx.BoxSizer(wx.VERTICAL)
        box.Add(gs0, 0, wx.EXPAND | wx.ALL, 3 )
        box.AddSpacer(10)
        box.Add(gs, 1, wx.EXPAND | wx.ALL, 3 )
        box.AddSpacer(10)         
        box.Add(gs1, 0, wx.EXPAND | wx.ALL, 3 )
        box.Add(self.buttonbox, 0, wx.CENTER, 3 )
        self.SetSizer(box)
        self.Centre()

        self.editchain.Enable(True)
        self.lblchain.Enable(True)
        
        self.lc.Enable(False)
        self.add.Enable(False)
        self.edit.Enable(False)
        self.remove.Enable(False)
        self.generate.Enable(False)
        self.editlength.Enable(False)
        self.lbllength.Enable(False)
        self.lblconcentration.Enable(False)

        if params!=-1:
            self.editname.SetValue(self.params[1])
            self.editchain.SetValue(self.params[2])
            self.editlength.SetValue(self.params[3])

            c=self.params[4].split(",")
            step=0
            for x in c:
                data=x.split("=")
                if len(data)==2:
                    self.lc.InsertStringItem(step, str(step))
                    if data[0]!='':
                        self.lc.SetStringItem(step, 0, str(data[0]).strip())
                    if data[1]!='':
                        self.lc.SetStringItem(step, 1, str(data[1]).strip())
                        
                    step+=1

            self.editconc.SetValue(self.params[5])


    def OnSourceCombobox(self,event):
        if self.editsource.Value=="percentage":
            self.activator=False
        else:
            self.activator=True

        self.editchain.Enable(self.activator)
        self.lblchain.Enable(self.activator)
        
        self.lc.Enable(not self.activator)
        self.add.Enable(not self.activator)
        self.edit.Enable(not self.activator)
        self.remove.Enable(not self.activator)
        self.generate.Enable(not self.activator)
        self.editlength.Enable(not self.activator)
        self.lbllength.Enable(not self.activator)
        self.lblconcentration.Enable(not self.activator)



    def OnOk(self,e):    
        
        name=str(self.editname.GetValue())
        
        if name=="":
            self.parent.parent.errorPopup("A name for this polymer must be provided!","ERROR!")
            return
        
        polymers=[]
        for i in range(0,self.parent.lc.GetItemCount(),1):
            item=self.parent.lc.GetItem(itemIdx=i,col=0)
            polymers.append(str(item.GetText()))
        if name in polymers:
            if self.params==-1:
                self.parent.parent.errorPopup("The desired polymer name already exists in the database!\nPlease choose a unique name!","ERROR!")
                return            
            elif self.params[1]!=name:            
                self.parent.parent.errorPopup("The desired polymer name already exists in the database!\nPlease choose a unique name!","ERROR!")
                return
                    
        #recalculate length and percentages values on the available chain!
        chain=str(self.editchain.GetValue()) 
        if len(chain)==0:
            self.parent.parent.errorPopup("A chain must be defined!","ERROR!")
            return
        
        for x in range(0,len(chain),1):                
            if chain[x] not in self.molecules:
                self.parent.parent.errorPopup("Residue %s is not included in the\nresidues database"%chain[x],"ERROR!")
                return
            
        self.editlength.SetValue(str(len(chain)))

        textperc=""
        #if self.editsource.GetValue()=="chain":
        cnt={}
        for x in range(0,len(chain),1):
            if chain[x] in list(cnt):
                cnt[chain[x]]=cnt[chain[x]]+1
            else:
                cnt[chain[x]]=1
        
        i=0
        for x in list(cnt):
            val=cnt[x]/float(self.editlength.GetValue())*100.0
            textperc+="%s=%6.2f, "%(x,val)
   
            self.lc.InsertStringItem(i, str(i+1))
            self.lc.SetStringItem(i, 0, str(x))
            self.lc.SetStringItem(i, 1, str(val))
            i+=1       

        #else:
        #    for index in range(0,self.lc.GetItemCount(),1):
        #        name=str(self.lc.GetItem(index,0).GetText())
        #        percentage=str(self.lc.GetItem(index,1).GetText())
        #        textperc+="%s=%s, "%(name,percentage)
        
        if self.editsource=="percentage":
            try:
                self.OnGenerate(e)
            except:
                return
                
        if self.params!=-1:
            num_items=self.params[0]
        else:
            num_items = self.parent.lc.GetItemCount()
            self.parent.lc.InsertStringItem(num_items, str(num_items+1))

        self.parent.lc.SetStringItem(num_items, 0, str(self.editname.GetValue()))
        self.parent.lc.SetStringItem(num_items, 1, str(self.editchain.GetValue()))
        self.parent.lc.SetStringItem(num_items, 2, str(self.editlength.GetValue()))
        self.parent.lc.SetStringItem(num_items, 3, str(textperc))
        self.parent.lc.SetStringItem(num_items, 4, str(self.editconc.GetValue()))

        self.Close(True)  # Close the frame.


    def OnBadd(self,event):
        #add error message
        dia=ChainEditor(self)
        dia.ShowModal()
        dia.Destroy()
            
    def OnBedit(self,event):

        index = self.lc.GetFocusedItem()

        params=[]
        params.append(index)
        if index != -1 :
            params.append(str(self.lc.GetItem(index, 0).GetText()))
            params.append(str(self.lc.GetItem(index, 1).GetText()))
            dia=ChainEditor(self,params)
            dia.ShowModal()
            dia.Destroy()


    def OnGenerate(self,e):
        
        try:
            editlength_f=int(self.editlength.GetValue())
        except ValueError:
            self.parent.parent.errorPopup("length should be an integer number\n greater than 0!","ERROR!")
            return -1
        
        if editlength_f<=0:
            self.parent.parent.errorPopup("length should be an integer number\n greater than 0!","ERROR!")
            return -1
    
        perc=[]
        test=0.0
        for index in range(0,self.lc.GetItemCount(),1):
            percentage=float(self.lc.GetItem(index,1).GetText())
            test+=percentage
            perc.append([self.lc.GetItem(index,0).GetText(),percentage])
                                    
        if test!=100.0:
            self.parent.parent.errorPopup("Sum of percentages equal to %s (100 is expected)!\nTreating given values as ratios..."%test,"WARNING!")
            for i in range(0,len(perc),1):
                perc[i][1]/=float(test)
                perc[i][1]*=100.0
            #self.parent.parent.errorPopup("Sum of percentages is equal to %s! 100 is expected!"%test,"ERROR!")
            #return -1
        
        try:
            chain=A.make_chain(editlength_f, perc)
        except:
            print("ERROR: something unexpected happened during chain generation, how weird...")
            return -1
        
        self.editchain.SetValue(chain)


    def OnBremove(self,event):
        
        index = self.lc.GetFocusedItem()
        if index!=-1:
            self.lc.DeleteItem(index)
    
            #renumber items after deletion
            #num_items = self.lc.GetItemCount()
            #for x in range(index,num_items,1):
            #    self.lc.SetStringItem(x, 0, str(x+1))

    def OnCancel(self,e):
        self.Close(True)  # Close the frame.
        

class ChainEditor(wx.Dialog):
    def __init__(self, parent,params=-1):
        wx.Dialog.__init__(self, parent,title="Random Chain Generator",style=wx.RESIZE_BORDER|wx.DEFAULT_DIALOG_STYLE|wx.TAB_TRAVERSAL,size=(250,90))

        self.params=params
        self.parent=parent

        self.editmol = wx.ComboBox(self, value="", size=(60, -1), choices=self.parent.molecules, style=wx.CB_DROPDOWN)
        #self.Bind(wx.EVT_COMBOBOX, self.OnSourceCombobox,self.editsource)

        self.lblperc = wx.StaticText(self, label="%:")
        self.editperc = wx.TextCtrl(self, value="", size=(60,-1))

        if self.params!=-1:
            self.editmol.SetValue(self.params[1])
            self.editperc.SetValue(self.params[2])

        gs = wx.GridSizer(1, 3, 3, 3)
        gs.AddMany([(self.editmol, 0, wx.ALIGN_RIGHT | wx.ALIGN_CENTER_VERTICAL, 3),
                    (self.lblperc, 0,  wx.ALIGN_RIGHT | wx.ALIGN_CENTER_VERTICAL, 3),
                    (self.editperc, 0, wx.ALIGN_LEFT | wx.ALIGN_CENTER_VERTICAL, 3)])

        self.okbutton=wx.Button(self, label="OK!")
        self.Bind(wx.EVT_BUTTON, self.OnOk, self.okbutton,id=wx.ID_OK)
        self.cancel =wx.Button(self, label="Cancel")
        self.Bind(wx.EVT_BUTTON, self.OnCancel,self.cancel,id=wx.ID_CANCEL)
        self.buttonbox = wx.BoxSizer(wx.HORIZONTAL)
        self.buttonbox.Add(self.okbutton, 0, wx.RIGHT, 3 )
        self.buttonbox.Add(self.cancel, 0, wx.LEFT, 3 )

        #draw title, and grid in a vertical box
        box = wx.BoxSizer(wx.VERTICAL)
        box.Add(gs, 0, wx.EXPAND | wx.ALL, 3 )
        box.Add(self.buttonbox, 1, wx.CENTER, 3 )
        self.SetSizer(box)
        self.Centre()

    def OnOk(self,e):

        name=str(self.editmol.GetValue())
        perc=str(self.editperc.GetValue())
        
        if name=="":
            self.parent.parent.parent.errorPopup("a residue must be selected!","ERROR!")
            return            

        if name not in self.parent.molecules:
            self.parent.parent.parent.errorPopup("The desired residue does not exist in the database!","ERROR!")
            return

        if perc=="":
            self.parent.parent.parent.errorPopup("a percentage value must be selected!","ERROR!")
            return  
        
        try:
            perc_f=float(perc)
        except:
            self.parent.parent.parent.errorPopup("percentage value must be a number between 0 and 100!","ERROR!")        
            return
          
        if perc_f>100.0 or perc_f<0.0:
            self.parent.parent.parent.errorPopup("percentage value must be a number between 0 and 100!","ERROR!")                
            return

        residues=[]
        for i in range(0,self.parent.lc.GetItemCount(),1):
            item=self.parent.lc.GetItem(itemIdx=i,col=0)
            residues.append(str(item.GetText()))
        if name in residues:
            if self.params==-1:
                    self.parent.parent.parent.errorPopup("A percentage for has been already set for this residue!","ERROR!")
                    return
            elif self.params[1]!=name:
                    self.parent.parent.parent.errorPopup("A percentage for has been already set for this residue!","ERROR!")
                    return
        
        if self.params!=-1:
            num_items=self.params[0]
        else:
            num_items = self.parent.lc.GetItemCount()
            self.parent.lc.InsertStringItem(num_items, str(num_items+1))

        self.parent.lc.SetStringItem(num_items, 0, name)
        self.parent.lc.SetStringItem(num_items, 1, perc)

        self.Close(True)  # Close the frame.


    def OnCancel(self,e):
        self.Close(True)  # Close the frame.
        
    def OnClickMakeChain(self,e):
        pass


app = wx.App(False)
frame = MainWindow(None, "Assemble!")
frame.Show()
app.MainLoop()