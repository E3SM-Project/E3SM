import wx, wx.html
import os, sys, csv, glob
from wx.lib.expando import ExpandoTextCtrl
import getpass, socket


aboutText = """<p>ACME land model (ALM) version 1 point model,  
running on version %(wxpy)s of <b>wxPython</b> and %(python)s of <b>Python</b>.
See <a href="http://wiki.wxpython.org">wxPython Wiki</a></p>""" 

def Contains(str,set):
    return 0 not in [c in str for c in set]


class HtmlWindow(wx.html.HtmlWindow):
    def __init__(self, parent, id, size=(640,480)):
        wx.html.HtmlWindow.__init__(self,parent, id, size=size)
        if "gtk2" in wx.PlatformInfo:
            self.SetStandardFonts()

    def OnLinkClicked(self, link):
        wx.LaunchDefaultBrowser(link.GetHref())
        
class AboutBox(wx.Dialog):
    def __init__(self):
        wx.Dialog.__init__(self, None, -1, "About <<project>>",
            style=wx.DEFAULT_DIALOG_STYLE|wx.THICK_FRAME|wx.RESIZE_BORDER|
                wx.TAB_TRAVERSAL)
        hwin = HtmlWindow(self, -1, size=(480,240))
        vers = {}
        vers["python"] = sys.version.split()[0]
        vers["wxpy"] = wx.VERSION_STRING
        hwin.SetPage(aboutText % vers)
        btn = hwin.FindWindowById(wx.ID_OK)
        irep = hwin.GetInternalRepresentation()
        hwin.SetSize((irep.GetWidth()+25, irep.GetHeight()+10))
        self.SetClientSize(hwin.GetSize())
        self.CentreOnParent(wx.BOTH)
        self.SetFocus()

class Frame(wx.Frame):
    global file_opened
    global filename
    def __init__(self, title):

        wx.Frame.__init__(self, None, title=title, pos=(100,100), size=(900,700))
        self.Bind(wx.EVT_CLOSE, self.OnClose)

        file_opened=0
        #menu bar
        menuBar = wx.MenuBar()
        menu = wx.Menu()
        #  open function on menu bar
        #m_open = menu.Append(wx.ID_OPEN, "Open", "Open a file")
        #self.Bind(wx.EVT_MENU, self.OnOpen, m_open)
        #  exit function on menu bar
        m_exit = menu.Append(wx.ID_EXIT, "E&xit\tAlt-X", "Close window and exit program.")
        self.Bind(wx.EVT_MENU, self.OnClose, m_exit)
        #  file function on menu bar
        menuBar.Append(menu, "&File")
        menu = wx.Menu()
        m_about = menu.Append(wx.ID_ABOUT, "&About", "Information about this program")
        self.Bind(wx.EVT_MENU, self.OnAbout, m_about)
        menuBar.Append(menu, "&Help")
        self.SetMenuBar(menuBar)
        print(file_opened)

        #optpanel = wx.Panel(self)
        #self.sizer = wxGridBagSizer(5,5)
        #self.abc = wxCheckBox(self.panel,-1,'abc')
        
        self.statusbar = self.CreateStatusBar()
        panel = wx.Panel(self)

        hbox = wx.BoxSizer(wx.HORIZONTAL)
        box  = wx.BoxSizer(wx.VERTICAL)        
        #box1 = wx.BoxSixer(wx.HORIZONTAL)
        
        m_text = wx.StaticText(panel, -1, "Point ALM options", style=wx.ALIGN_LEFT)
        m_text.SetFont(wx.Font(12, wx.SWISS, wx.NORMAL, wx.BOLD))
        m_text.SetSize(m_text.GetBestSize())
        box.Add(m_text, 0, wx.ALL, 3)
        
        #CCSM input data directory
        m_text9 = wx.StaticText(panel, -1, "Input data directory", style=wx.ALIGN_LEFT)
        m_text9.SetFont(wx.Font(10, wx.SWISS, wx.NORMAL, wx.BOLD))
        m_text9.SetSize(m_text9.GetBestSize())
        box.Add(m_text9, 0, wx.ALL, 3)
        self.inputdirtxt = wx.TextCtrl(panel, -1, ccsm_input)
        self.inputdirtxt.Bind(wx.EVT_TEXT, self.OnInputdirText)
        box.Add(self.inputdirtxt, 0, wx.EXPAND)
        browse_inputdir = wx.Button(panel, -1, "Change", (30,30))
        browse_inputdir.Bind(wx.EVT_BUTTON, self.OnInputdirOpen)
        box.Add(browse_inputdir, 0, wx.ALL, 3)

        #Run directory
        m_text9a = wx.StaticText(panel, -1, "Model run directory", style=wx.ALIGN_LEFT)
        m_text9a.SetFont(wx.Font(10, wx.SWISS, wx.NORMAL, wx.BOLD))
        m_text9a.SetSize(m_text9a.GetBestSize())
        box.Add(m_text9a, 0, wx.ALL, 3)
        self.rundirtxt = wx.TextCtrl(panel, -1, rundir)
        self.rundirtxt.Bind(wx.EVT_TEXT, self.OnRundirText)
        box.Add(self.rundirtxt, 0, wx.EXPAND)
        browse_rundir = wx.Button(panel, -1, "Change", (30,30))
        browse_rundir.Bind(wx.EVT_BUTTON, self.OnRundirOpen)
        box.Add(browse_rundir, 0, wx.ALL, 3)

        #site group selection
        m_text2a = wx.StaticText(panel, -1, "site group", style=wx.ALIGN_LEFT)
        m_text2a.SetFont(wx.Font(10, wx.SWISS, wx.NORMAL, wx.BOLD))
        m_text2a.SetSize(m_text2a.GetBestSize())
        box.Add(m_text2a, 0, wx.ALL, 3)
        self.mysitegroup = wx.Choice(panel, -1, choices = mysitegroups)
        self.mysitegroup.Bind(wx.EVT_CHOICE, self.OnSiteGroupSelect)
        box.Add(self.mysitegroup, 0, wx.ALL, 3)

        #site selection
        m_text3 = wx.StaticText(panel, -1, "site", style=wx.ALIGN_LEFT)
        m_text3.SetFont(wx.Font(10, wx.SWISS, wx.NORMAL, wx.BOLD))
        m_text3.SetSize(m_text3.GetBestSize())
        box.Add(m_text3, 0, wx.ALL, 3)
        self.m_site = wx.Choice(panel, -1, choices = mysites)
        self.m_site.Bind(wx.EVT_CHOICE, self.OnSiteSelect)
        box.Add(self.m_site, 0, wx.ALL, 3)

        #case name entry
        m_text7 = wx.StaticText(panel, -1, "Case prefix", style=wx.ALIGN_LEFT)
        m_text7.SetFont(wx.Font(12, wx.SWISS, wx.NORMAL, wx.BOLD))
        m_text7.SetSize(m_text7.GetBestSize())
        box.Add(m_text7, 0, wx.ALL, 3)
        self.m_caseprefix = wx.TextCtrl(panel, -1, mycaseprefix)
        self.m_caseprefix.Bind(wx.EVT_TEXT, self.OnCaseText)
        box.Add(self.m_caseprefix,0,wx.ALL,3)

        #spinup selection
        m_text4 = wx.StaticText(panel, -1, "Run type", style=wx.ALIGN_LEFT)
        m_text4.SetFont(wx.Font(10, wx.SWISS, wx.NORMAL, wx.BOLD))
        m_text4.SetSize(m_text4.GetBestSize())
        box.Add(m_text4, 0, wx.ALL, 3)
        self.m_spinup = wx.Choice(panel, -1, choices=['ad spinup','final spinup','transient','full simulation'])
        self.m_spinup.Bind(wx.EVT_CHOICE, self.OnSpinupButton)
        box.Add(self.m_spinup,0,wx.ALL,3)

        #number of years entry
        m_text5 = wx.StaticText(panel, -1, "# of years", style=wx.ALIGN_LEFT)
        m_text5.SetFont(wx.Font(10, wx.SWISS, wx.NORMAL, wx.BOLD))
        m_text5.SetSize(m_text4.GetBestSize())
        box.Add(m_text5, 0, wx.ALL, 3)
        self.nyears = wx.TextCtrl(panel, -1, mynyears)
        box.Add(self.nyears,0,wx.ALL,3)

        #restart file entry (option to browse)
        m_text6 = wx.StaticText(panel, -1, "Finidat file", style=wx.ALIGN_LEFT)
        m_text6.SetFont(wx.Font(10, wx.SWISS, wx.NORMAL, wx.BOLD))
        m_text6.SetSize(m_text4.GetBestSize())
        box.Add(m_text6, 0, wx.ALL, 3)
        self.m_finidat = wx.TextCtrl(panel, txtID, myfinidat)
        self.m_finidat.Bind(wx.EVT_TEXT, self.OnFinidatText)
        box.Add(self.m_finidat,0,wx.EXPAND)
        browse_finidat = wx.Button(panel, -1, "Browse", (30,30))
        browse_finidat.Bind(wx.EVT_BUTTON, self.OnFinidatOpen)
        box.Add(browse_finidat, 0, wx.ALL, 3)

        #pft physiology file entry
        m_text8 = wx.StaticText(panel, -1, "Custom parameter file", style=wx.ALIGN_LEFT)
        m_text8.SetFont(wx.Font(10, wx.SWISS, wx.NORMAL, wx.BOLD))
        m_text8.SetSize(m_text4.GetBestSize())
        box.Add(m_text8, 0, wx.ALL, 3)
        self.m_pftfile = wx.TextCtrl(panel, -1, mypftfile)
        self.m_pftfile.Bind(wx.EVT_TEXT, self.OnPftfileText)
        box.Add(self.m_pftfile,0,wx.EXPAND)
        browse_finidat = wx.Button(panel, -1, "Browse", (20,20))
        browse_finidat.Bind(wx.EVT_BUTTON, self.OnPftfileOpen)
        box.Add(browse_finidat, 0, wx.ALL, 3)
 

        #-----------Box 2------------------------
        box2 = wx.BoxSizer(wx.VERTICAL)
        
        #site information
        self.m_siteinfo = wx.TextCtrl(panel, -1, 'Site information for: \n\n Longitude: \n Latitude: \n Primary PFTS:\n PFT1: \n PFT2: ', style=wx.TE_READONLY | wx.TE_MULTILINE)
        self.m_siteinfo.SetFont(wx.Font(10, wx.SWISS, wx.NORMAL, wx.BOLD))
        self.m_siteinfo.SetSize(m_text2a.GetBestSize())
        #self.m_siteinfo.Bind(wx.EVT_TEXT, self.OnFinidatText)
        box2.Add(self.m_siteinfo,0,wx.EXPAND)
    
        m_text10= wx.StaticText(panel, -1, "Site-specific options", style=wx.ALIGN_LEFT)
        m_text10.SetFont(wx.Font(10, wx.SWISS, wx.NORMAL, wx.BOLD))
        m_text10.SetSize(m_text4.GetBestSize())
        box2.Add(m_text10, 0, wx.ALL, 3)
        self.mksrfdat = wx.CheckBox(panel, -1, 'Make surface data files')
        self.mksrfdat.SetValue(True)
        box2.Add(self.mksrfdat,0,wx.EXPAND)
        self.sitedata = wx.CheckBox(panel, -1,'Use site-level data in surface data files')
        self.sitedata.SetValue(True)
        box2.Add(self.sitedata,0,wx.EXPAND)
        self.siteparms = wx.CheckBox(panel, -1,'Use site parameters')
        self.siteparms.SetValue(True)
        box2.Add(self.siteparms,0,wx.EXPAND)

        m_text11= wx.StaticText(panel, -1, "Model configuration", style=wx.ALIGN_LEFT)
        m_text11.SetFont(wx.Font(10, wx.SWISS, wx.NORMAL, wx.BOLD))
        m_text11.SetSize(m_text4.GetBestSize())
        box2.Add(m_text11, 0, wx.ALL, 3)
        self.spmode =  wx.CheckBox(panel, -1,'Use satellite phenology')
        self.spmode.SetValue(False)
        box2.Add(self.spmode,0,wx.EXPAND)
        self.ecamode =  wx.CheckBox(panel, -1,'Use ECA (RD if not selected)')
        self.ecamode.SetValue(False)
        box2.Add(self.ecamode,0,wx.EXPAND)
        self.conly =  wx.CheckBox(panel, -1,'Carbon only mode')
        self.conly.SetValue(False)
        box2.Add(self.conly,0,wx.EXPAND)
        self.cnonly =  wx.CheckBox(panel, -1,'Carbon-nitrogen mode')
        self.cnonly.SetValue(False)
        box2.Add(self.cnonly,0,wx.EXPAND)
        self.fire = wx.CheckBox(panel, -1,'Use wildfire submodel')
        self.fire.SetValue(True)
        box2.Add(self.fire,0,wx.EXPAND)
        self.dynroot = wx.CheckBox(panel, -1,'Use dynamic rooting')
        self.dynroot.SetValue(False)
        box2.Add(self.dynroot,0,wx.EXPAND)
        self.cpl_bypass = wx.CheckBox(panel, -1,'No datm (coupler bypass)')
        self.cpl_bypass.SetValue(True)
        box2.Add(self.cpl_bypass,0,wx.EXPAND)
        self.onehour = wx.CheckBox(panel, -1,'Use hourly timestep')
        self.onehour.SetValue(True)
        box2.Add(self.onehour,0,wx.EXPAND)
        self.cruncep = wx.CheckBox(panel, -1,'Use cru-ncep forcing data')
        self.cruncep.SetValue(False)
        box2.Add(self.cruncep,0,wx.EXPAND)

        #Run PTCLM button
        m_text12= wx.StaticText(panel, -1, "\nConfigure, build and run", style=wx.ALIGN_LEFT)
        m_text12.SetFont(wx.Font(12, wx.SWISS, wx.NORMAL, wx.BOLD))
        m_text12.SetSize(m_text4.GetBestSize())
        box2.Add(m_text12, 0, wx.ALL, 3)
        m_run = wx.Button(panel, -1, "   Begin simulation   ", style=wx.ALIGN_CENTER)
        m_run.Bind(wx.EVT_BUTTON, self.OnRun)
        box2.Add(m_run, 0, wx.ALL, 10)

        #box 3 - UQ 
        box3 = wx.BoxSizer(wx.VERTICAL)
        m_text30= wx.StaticText(panel, -1, "Uncertainty Quantification", style=wx.ALIGN_LEFT)
        m_text30.SetFont(wx.Font(10, wx.SWISS, wx.NORMAL, wx.BOLD))
        m_text30.SetSize(m_text30.GetBestSize())
        box3.Add(m_text30, 0, wx.ALL, 3)

        m_text31 = wx.StaticText(panel, -1, "Ensemble File", style=wx.ALIGN_LEFT)
        m_text31.SetFont(wx.Font(10, wx.SWISS, wx.NORMAL, wx.BOLD))
        m_text31.SetSize(m_text31.GetBestSize())
        box3.Add(m_text31, 0, wx.ALL, 3)
        self.ensemblefiletxt = wx.TextCtrl(panel, -1, 'none')
        self.ensemblefiletxt.Bind(wx.EVT_TEXT, self.OnEnsembleFileText)
        box3.Add(self.ensemblefiletxt, 0, wx.EXPAND)
        browse_ensembledir = wx.Button(panel, -1, "Change", (30,30))
        browse_ensembledir.Bind(wx.EVT_BUTTON, self.OnEnsembleFileOpen)
        box3.Add(browse_ensembledir, 0, wx.ALL, 3)

        self.mc_ensemble = wx.CheckBox(panel, -1,'Run Monte Carlo ensemble')
        self.mc_ensemble.SetValue(False)
        box3.Add(self.mc_ensemble,0,wx.EXPAND)
        m_text32 = wx.StaticText(panel, -1, "Number of ensemble members", \
                                 style=wx.ALIGN_LEFT)
        m_text32.SetFont(wx.Font(10, wx.SWISS, wx.NORMAL, wx.BOLD))
        m_text32.SetSize(m_text31.GetBestSize())
        box3.Add(m_text32, 0, wx.ALL, 3)
        self.n_ensemble = wx.TextCtrl(panel, -1, myn_ensemble)
        box3.Add(self.n_ensemble,0,wx.ALL,3)
        m_text33 = wx.StaticText(panel, -1, "Number of processors to use", \
                                 style=wx.ALIGN_LEFT)
        m_text33.SetFont(wx.Font(10, wx.SWISS, wx.NORMAL, wx.BOLD))
        m_text33.SetSize(m_text31.GetBestSize())
        box3.Add(m_text33, 0, wx.ALL, 3)
        self.n_proc = wx.TextCtrl(panel, -1, myn_proc)
        self.n_proc.SetValue('256')
        box3.Add(self.n_proc,0,wx.ALL,3)

        #Text editor
        edit_file = wx.Button(panel, -1, "View/Edit Files", (20,20))
        edit_file.Bind(wx.EVT_BUTTON, self.OnPftfileView)
        box3.Add(edit_file, 0, wx.ALL, 3)

        hbox.Add(box,1,wx.EXPAND)
        hbox.Add(box2,1,wx.EXPAND)
        hbox.Add(box3,1,wx.EXPAND)
        panel.SetSizer(hbox)
        panel.Layout()

    def OnCaseText(self, event):
        mycaseprefix=event.GetString()
        
    def OnClose(self, event):
        dlg = wx.MessageDialog(self, 
            "Do you really want to close this application?",
            "Confirm Exit", wx.OK|wx.CANCEL|wx.ICON_QUESTION)
        result = dlg.ShowModal()
        dlg.Destroy()
        if result == wx.ID_OK:
            self.Destroy()

    def OnEnsembleFileOpen(self, event):
        file_opened = 0
        dialog = wx.FileDialog(None, wildcard="*.dat;*.txt", defaultDir="./", style = wx.OPEN)
        #dialog = wx.FileDialog(None, defaultDir="./", style = wx.OPEN)
        result = dialog.ShowModal()
        if result == wx.ID_OK:
            ensemble_file = dialog.GetPath()
            self.ensemblefiletxt.SetValue(ensemble_file)
            #count number of lines in the file
            with open(ensemble_file) as f:
                myn_ensemble = sum(1 for _ in f)
            self.n_ensemble.SetValue(str(myn_ensemble))
            myn_proc = self.n_proc.GetValue()
            if (int(myn_proc) > int(myn_ensemble+1)):
                self.n_proc.SetValue(str(int(myn_ensemble+1)))
            self.mc_ensemble.SetValue(False)
            dialog.Destroy()
        else:
            dialog.Destroy()
  
    def OnEnsembleFileText(self, event):
        ensemble_file=event.GetString()

    def OnInputdirOpen(self, event):
        file_opened = 0
        dialog = wx.DirDialog(None, style = wx.OPEN)
        result = dialog.ShowModal()
        if result == wx.ID_OK:
            ccsm_input = dialog.GetPath()
            self.inputdirtxt.SetValue(ccsm_input)
            dialog.Destroy()
        else:
            dialog.Destroy()
  
    def OnInputdirText(self, event):
        ccsm_input=event.GetString()

    def OnRundirOpen(self, event):
        file_opened = 0
        dialog = wx.DirDialog(None, style = wx.OPEN)
        result = dialog.ShowModal()
        if result == wx.ID_OK:
            rundir = dialog.GetPath()
            self.inputdirtxt.SetValue(rundir)
            dialog.Destroy()
        else:
            dialog.Destroy()
  
    def OnRundirText(self, event):
        rundir=event.GetString()

    def OnFinidatOpen(self, event):
        rundir=self.rundirtxt.GetValue()
        dialog = wx.FileDialog(None, wildcard = "Restart files|*clm2*r*.nc", defaultDir=rundir, style = wx.OPEN)
        result = dialog.ShowModal()
        if result == wx.ID_OK:
            myfinidat = dialog.GetPath()
            self.m_finidat.SetValue(myfinidat)
            print(myfinidat)
            dialog.Destroy()
        else:
            dialog.Destroy()

    def OnFinidatText(self, event):
        myfinidat=event.GetString()
        print(myfinidat)

    def OnSiteGroupSelect(self,event):
        index=event.GetSelection()
        mysitegroup_current=mysitegroups[index]
        #load site information for selected group
        ccsm_input=self.inputdirtxt.GetValue()
        fname=ccsm_input+'/lnd/clm2/PTCLM/'+mysitegroup_current+"_sitedata.txt"
        AFdatareader = csv.reader(open(fname,"rb"))
        nsites=1
        mysites=['all']
        mylons=['all sites']
        mylats=['all sites']
        for row in AFdatareader:
            if nsites > 1:
                mysites.append(row[0])
                mylons.append(row[3])
                mylats.append(row[4])
            nsites=nsites+1
        mysite=mysites[0]   #set new default site
        #change the site selection drop-down box
        self.m_site.Clear()
        for site in mysites:
            self.m_site.Append(str(site))
        self.m_site.SetSelection(0)
        
    def OnSiteSelect(self,event):
        #get the site group
        indexsg = self.mysitegroup.GetSelection()
        mysitegroup_current=mysitegroups[indexsg]
        #load site information for selected group
        ccsm_input=self.inputdirtxt.GetValue()
        fname=ccsm_input+'/lnd/clm2/PTCLM/'+mysitegroup_current+"_sitedata.txt"
        AFdatareader = csv.reader(open(fname,"rb"))
        nsites=1
        mysites=['all']
        mynames=['all sites']
        mylons=['all sites']
        mylats=['all sites']
        for row in AFdatareader:
            if nsites > 1:
                mysites.append(row[0])
                mynames.append(row[1])
                mylons.append(row[3])
                mylats.append(row[4])
            nsites=nsites+1
        #get the selected site and get appropriate info
        fname = ccsm_input+'/lnd/clm2/PTCLM/'+mysitegroup_current+"_pftdata.txt"
        AFdatareader = csv.reader(open(fname,"rb"))
        index=event.GetSelection()
        mysite=mysites[index]
        pft1_type = 'N/A'
        pft1_pct = '0'
        pft2_type = 'N/A'
        pft2_pct = '0'
        for row in AFdatareader:
            if (row[0] == mysite):
                pft1_type = row[2]
                pft1_pct  = row[1]
                pft2_type = row[4]
                pft2_pct  = row[3]
        print pft2_type
        if (pft1_type.strip() != 'N/A' and pft2_pct.strip() == '0.0'):
            self.m_siteinfo.SetValue('Site information for '+mysite+'\n Site name: '+mynames[index]+ \
                                         '\n Longitude: '+str(mylons[index])+'\n Latitude: '+str(mylats[index])+ \
                                         '\n primary PFTS:\n '+mypfts[int(pft1_type)]+': '+pft1_pct+'%')
        elif (pft2_type.strip() != 'N/A'):
            self.m_siteinfo.SetValue('Site information for '+mysite+'\n Site name: '+mynames[index]+ \
                                         '\n Longitude: '+str(mylons[index])+'\n Latitude: '+str(mylats[index])+ \
                                         '\n primary PFTS:\n '+mypfts[int(pft1_type)]+': '+pft1_pct+'%\n '+ \
                                         mypfts[int(pft2_type)]+': '+pft2_pct+'%')
        else:
            self.m_siteinfo.SetValue('Site information for '+mysite+'\n Site name: '+mynames[index]+ \
                                         '\n Longitude: '+str(mylons[index])+'\n Latitude: '+str(mylats[index])+ \
                                         '\n primary PFTS:\n PFT1: N/A\n PFT2: N/A')

    def OnSpinupButton(self,event):
        myindex=event.GetSelection()
        index=self.m_site.GetSelection()
        mysite=mysites[index]
        if myindex == 0:
            myspinup='--ad_spinup'
            self.nyears.SetValue('250')
        if myindex == 1:
            myspinup='--final_spinup'
            self.nyears.SetValue('250')
        if myindex == 2:
            myspinup='--transient'
            self.nyears.SetValue('150')
        if myindex == 3:
            myspinup='--fullrun'
            self.nyears.SetValue('250')
        print(myspinup)
            
    def OnPftfileOpen(self, event):
        ccsm_input=self.inputdirtxt.GetValue()
        dialog = wx.FileDialog(self, defaultDir=ccsm_input+"/lnd/clm2/pftdata", style = wx.OPEN)
        result = dialog.ShowModal()
        if result == wx.ID_OK:
            mypftfile = dialog.GetPath()
            self.m_pftfile.SetValue(mypftfile)
            print(mypftfile)
            dialog.Destroy()
        else:
            dialog.Destroy()

    def OnPftfileText(self, event):
        mypftfile=event.GetString()
        print(mypftfile)
        
    def OnPftfileView(self, event):
        mypftfile=self.m_pftfile.GetValue()
        print(mypftfile)
        os.system('emacs '+mypftfile+' &')
            
    def OnRun(self, event):
        dorun = True
        #get all of the information for PTCLM
        ccsm_input=self.inputdirtxt.GetValue()
        rundir=self.rundirtxt.GetValue()
        ensemble_file=self.ensemblefiletxt.GetValue()
        myn_ensemble=self.n_ensemble.GetValue()
        myn_proc=self.n_proc.GetValue()
        myindex=self.mysitegroup.GetSelection()
        mysitegroup_current=mysitegroups[myindex]
        ccsm_input=self.inputdirtxt.GetValue()
        fname=ccsm_input+'/lnd/clm2/PTCLM/'+mysitegroup_current+"_sitedata.txt"
        AFdatareader = csv.reader(open(fname,"rb"))
        nsites=1
        mysites=['all']
        mylons=['all sites']
        mylats=['all sites']
        for row in AFdatareader:
            if nsites > 1:
                mysites.append(row[0])
                mylons.append(row[3])
                mylats.append(row[4])
            nsites=nsites+1
        mynyears=self.nyears.GetValue()
        mycaseprefix=self.m_caseprefix.GetValue()
        mypftfile=self.m_pftfile.GetValue()
        myfinidat=self.m_finidat.GetValue()
        myindex=self.m_site.GetSelection()
        mysite=mysites[myindex]
        myindex=self.m_spinup.GetSelection()
        mymksrfdat=self.mksrfdat.GetValue()
        mysitedata=self.sitedata.GetValue()
        myspmode=self.spmode.GetValue()
        myfire=self.fire.GetValue()
        mydynroot=self.dynroot.GetValue()
        mymcensemble=self.mc_ensemble.GetValue()
        myconly=self.conly.GetValue()
        mycnonly=self.cnonly.GetValue()
        mycpl_bypass=self.cpl_bypass.GetValue()
	myonehour=self.onehour.GetValue()     
        mycruncep=self.cruncep.GetValue()

        cdate="1850"
        if (myindex == 2):
            cdate="20TR"
        if (myspmode):
            compset="CLM45"
            if (mycpl_bypass):
                compset="CLM45CB"
        else:
            compset="CLM45CN"
            if (mycpl_bypass):
                compset="CLM45CBCN"
        compset = "I"+cdate+compset
        
        cmd='python pointCLM.py --site '+mysite+' --machine '+machine+' --rmold --ccsm_input '+ccsm_input+' --sitegroup '+ \
            mysitegroup_current+' --compset '+compset+' --rmold'
        if (myindex == 0):
            cmd = cmd+' --ad_spinup --nyears_ad_spinup '+str(mynyears)
        elif (myindex < 3):
            cmd = cmd+' --run_n '+str(mynyears)
        elif (myindex == 3):
            cmd = 'python site_fullrun.py --site '+mysite+' --machine '+machine+' --ccsm_input '+ccsm_input+' --sitegroup '+ \
                mysitegroup_current+' --nyears_ad_spinup '+str(mynyears)+' --nyears_final_spinup '+str(mynyears)+ \
                ' --spinup_vars'
            if (mycpl_bypass):
                cmd = cmd+' --cpl_bypass'
        if (mycaseprefix != 'none'):
            cmd = cmd+' --caseidprefix '+mycaseprefix
        if mymksrfdat == False:
            cmd=cmd+' --nopointdata'
        if mysitedata == False:
            cmd=cmd+' --pftgrid --soilgrid'
        if mypftfile != '':
            cmd=cmd+' --pftfile '+mypftfile
        if myfinidat != '':
            cmd=cmd+' --finidat '+myfinidat
        elif (myindex == 1):
            cmd=cmd+' --coldstart'
        if (myfire == False):
            cmd = cmd+' --nofire'
        if (mydynroot == False):
            cmd = cmd+' --no_dynroot'
        if (myconly == True):
            cmd = cmd+' --c_only'
        if (mycnonly == True):
            cmd = cmd+' --cn_only'
        if (ensemble_file != 'none'):
            cmd = cmd+' --ensemble_file '+ensemble_file+' --ng '+myn_proc
        if (mymcensemble):
            cmd = cmd+' --mc_ensemble '+myn_ensemble+' --ng '+myn_proc
        if (myonehour):
            cmd = cmd+' --tstep 1'
        if (mycruncep):
            cmd = cmd+' --cruncep'
        cmd = cmd+' --runroot '+rundir
        
        print cmd
        os.system(cmd)
          
    def OnAbout(self, event):
        dlg = AboutBox()
        dlg.ShowModal()
        dlg.Destroy()  

#global variables
#global txtID
#global mycompsets, mysites
#global mycasename, myfinidat, mypftfile, mynyears, ccsm_input,
global PTCLMdir, mypfts, mysite, mysites, machine

username = getpass.getuser()
hostname = socket.gethostname()

PTCLMdir = os.path.abspath('./')
rundir = os.path.abspath(PTCLMdir+'../../../../../../run')
if ('or-condo' in hostname):
  ccsm_input = '/lustre/or-hydra/cades-ccsi/proj-shared/project_acme/ACME_inputdata/'
  rundir = '/lustre/or-hydra/cades-ccsi/scratch/'+username
  machine = 'cades'
elif ('eos' in hostname or 'titan' in hostname):
  ccsm_input = '/lustre/atlas/world-shared/cli900/cesm/inputdata'
  rundir = '/lustre/atlas/scratch/cli112/'+username
else:
  ccsm_input='/home/'+username+'/models/ccsm_inputdata'
  machine=oic2

os.chdir(ccsm_input+'/lnd/clm2/PTCLM')


mysiteindex=1
mysites=['all']
mylats=['all sites']
mylons=['all sites']
mysitegroups=[]
mysitegroup_current="AmeriFlux"
mypath=os.path.abspath("./")

#get available site groups
for filename in os.listdir("./"):
    if Contains(filename,'sitedata.txt'):
        groupname, type = filename.split('_')
        isrepeat=False
        for g in mysitegroups:
            if groupname == g:
                isrepeat=True
        if isrepeat == False:
            mysitegroups.append(groupname)

#load site information for default group
fname="./"+mysitegroup_current+"_sitedata.txt"
AFdatareader = csv.reader(open(fname,"rb"))
nsites=0
for row in AFdatareader:
    if nsites > 0:
        mysites.append(row[0])
        mylons.append(row[3])
        mylats.append(row[4])
    nsites=nsites+1

txtID=1
mysite=mysites[mysiteindex]
mycaseprefix='none'     #default case name
myfinidat='<none>'    #default finidat
myspinup=''
myfinidat=''
mypftfile=''
mynyears=''
myn_ensemble=''
os.chdir(PTCLMdir)
mypftfile=''
myn_proc=''

myoutvars=['NEE','GPP','NPP', 'ER', 'GR', 'AR', 'MR', 'TLAI', 'LEAFC', 'FROOTC', 'LIVESTEMC', 'DEADSTEMC','PFT_FIRE_CLOSS', 'LEAFC_ALLOC']
mypfts=['Bare ground','Temperate ENF','Boreal ENF','DNF','Tropical EBF','Temperate EBF','Tropical DBF','Temperate DBF','Boreal DBF','Broadleaf evergreen shrub','Temperate broadleaf deciduous shrub', 'Broadleaf deciduous boreal shrub', 'C3 arctic grass', 'C3 non-arctic grass', 'C4 grass', 'C3 crop', 'C3 irrigated', 'Corn'] 

app = wx.App(redirect=False)   # Error messages go to popup window
top = Frame("ALM point model GUI")
top.Show()
app.MainLoop()
