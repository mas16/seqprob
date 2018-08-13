#SeqProb GUI

#Statistical analysis of chemical shift distribution and protein amino acid sequence
#Chemical shift statistics last updated from BMRB 02/2018

#Since all chemical shift distributions are normal
#The statistical score is based on the probability density function for noraml variables

#by Matthew Stetz
#University of Pennsylvania, Perelman School of Medicine
#09/2009

import wx
import wx.grid as gridlib
import wx.lib.scrolledpanel as scrolled
import math
import os

class CustomDataTable(gridlib.PyGridTableBase):
    def __init__(self, log):
        gridlib.PyGridTableBase.__init__(self)
        #dictionary where all user CA and CB shifts stored
        #by default the dictionary is empty so the user can input anything
        #the program is set-up to also read in data from a txt file
        self.data = {}      
        #print self.data
        self.colLabels = ["CA (ppm)","CB (ppm)"]        
        self.odd=wx.grid.GridCellAttr()
        self.odd.SetBackgroundColour("gray")
        self.odd.SetFont(wx.Font(10, wx.SWISS, wx.NORMAL, wx.BOLD))
        self.even=wx.grid.GridCellAttr()
        self.even.SetBackgroundColour("white")
        self.even.SetFont(wx.Font(10, wx.SWISS, wx.NORMAL, wx.BOLD))

    #Define dimensions of input grid
    def GetNumberRows(self):
        return 500

    def GetNumberCols(self):
        return 2

    #Recognize and store input data
    def IsEmptyCell(self, row, col):
        return self.data.get((row, col)) is not None

    def GetValue(self, row, col):
        value = self.data.get((row, col))
        if value is not None:
            return value
        else:
            return ''

    def SetValue(self, row, col, value):
        self.data[(row, col)] = value

    def GetAttr(self, row, col, kind):
        attr = [self.even, self.odd] [row % 2]
        attr.IncRef()
        return attr

    def GetColLabelValue(self, col):
        return self.colLabels[col]

###FROM wxpython DEMO -- This gives you some more control over the grid parameters
class CustTableGrid(gridlib.Grid):
    def __init__(self, parent, log):
        gridlib.Grid.__init__(self, parent, -1)
        table = CustomDataTable(log)
        self.shifttabledata=table
        # The second parameter means that the grid is to take ownership of the
        # table and will destroy it when done.  Otherwise you would need to keep
        # a reference to it and call its Destroy method later.
        self.SetTable(table, True)
        self.SetRowLabelSize(35)
        self.SetMargins(0,0)
        self.AutoSizeColumns(False)
        gridlib.EVT_GRID_CELL_LEFT_DCLICK(self, self.OnLeftDClick)
        # I do this because I don't like the default behaviour of not starting the
        # cell editor on double clicks, but only a second click.

    def OnLeftDClick(self, evt):
        if self.CanEnableCellControl():
            self.EnableCellEditControl()
###End wxpyton DEMO code

#Generates a Frame to display all over your input and input controls
class TestFrame(wx.Frame):
    def __init__(self, parent, log):

        wx.Frame.__init__(
                    self, parent, -1, "SeqProb", size=(400,600)
                    )
        self.log=log

        #All widgets must be placed on a 'panel'
        #so you have control over their placement
        p = wx.Panel(self, -1, style=0)  
        grid = CustTableGrid(p, log)
        grid.SetColSize(0, 150)     
        grid.SetColSize(1, 150)     
        self.shifttable = grid
        
        #Define your menubar with functional menu options
        menuBar = wx.MenuBar()  
        menu = wx.Menu()
        menuBar.Append(menu, "File")        
        
        openf = menu.Append(-1, "Import CA and CB Shifts")
        self.Bind(wx.EVT_MENU, self.OnOpen, openf)

        menu.AppendSeparator()
        
        exit = menu.Append(-1, "Exit")
        self.Bind(wx.EVT_MENU, self.OnExit, exit)

        menu2 = wx.Menu()
        
        menu3 = wx.Menu()
        menuBar.Append(menu3, "About")
        
        usermanual=menu3.Append(-1, "User Manual")
        self.Bind(wx.EVT_MENU, self.OnUserManual, usermanual)       
                
        about=menu3.Append(-1, "About") 
        self.Bind(wx.EVT_MENU, self.OnAbout, about)     

        self.SetMenuBar(menuBar)

        #Define a button that will clear your sequence input--
        #Note: this code just defines the button, not its functionality
        b2 = wx.Button(p, -1, "Click to Clear Sequence", size=(335, 35))
        b2.SetDefault()
        self.Bind(wx.EVT_BUTTON, self.ClearSequence, b2)    

        #Define a button that will do your probability score calculations and display results--
        #Again, this is just code to make the button, not its function
        b4 = wx.Button(p, -1, "Click to Calculate", size=(335, 35))
        b4.SetDefault()
        self.Bind(wx.EVT_BUTTON, self.Calculate, b4)

        #Create and label your sequence input text box.
        self.sequencestring=''      
        basicLabel = wx.StaticText(p, -1, "Enter Protein Sequence:", pos=(50,30))   
        basicText = wx.TextCtrl(p, -1, self.sequencestring, size=(335, 200), pos=(50,30), style=wx.TE_MULTILINE)      
        self.bt=basicText  

        #Create a sizer to add your widgets to the panel and organize them. 
        bs = wx.BoxSizer(wx.VERTICAL)   
        bs.AddMany([basicLabel, basicText])
        bs.Add(b2)  
        bs.Add(grid, 1)  
        bs.Add(b4)        
        p.SetSizer(bs)

    #Function to trigger the File Open Dialog box (OS non-specific) for importing chemical shifts
    def OnOpen(self, evt):
        filetype = "Text file (*.txt)|*.txt|"
        dialog = wx.FileDialog(None, "Choose a file", os.getcwd(), "", filetype, wx.OPEN)
        if dialog.ShowModal()==wx.ID_OK:
            path = dialog.GetPath()
            dialog.Destroy()
        shifts = open(path)
        shifts = shifts.readlines()
        i = 0
        split=[]
        ca = []
        cb = []
    #Data reads in as a list of strings where each string corresponds to a single CA/CB pair. The CA and CB chemical shifts are delimited by a whitespace (\t). Therefore strings are split whenever a \t is encountered.
        while i < len(shifts):
            split.append(shifts[i].rpartition('\t'))
            i = i + 1
        i = 0   
        while i < len(split):
            ca.append(split[i][0])
            cb.append(split[i][2])
            i = i + 1
        i = 0
        dca = []
        dcb = []
        while (i < len(ca)):
            dca.append((i,0))
            dcb.append((i,1))   
            i=i+1
        i = 0
        dictionary={}
        while (i < len(ca)):
            dictionary.update({(dca[i]):(ca[i])})
            dictionary.update({(dcb[i]):(cb[i])})   
            i = i + 1
        self.shifttable.shifttabledata.data = dictionary
        
        #For some reason (on linux) the GUI needs to refresh
        #before it displays the newly imported data. These commands take care of that.
        self.shifttable.SelectAll()
        self.shifttable.ClearSelection()
        
    #Function to close the program when you select "Exit" from the menubar
    def OnExit(self, evt):
        self.Close()

    #Function to make the "About" information pop up as a separate mini-frame
    #when you select "About" from the menubar
    def OnAbout(self, evt):
        text="This program written by:\n\nMatt Stetz\n\nFall 2009\nUniversity of Pennsylvania"
        self.text=text
        win = Help(self, "About", self.text, style=wx.DEFAULT_FRAME_STYLE | wx.THICK_FRAME)
        win.Show(True)

    #Function to make the user manual pop up as a separate mini-frame when you select "User Manual" from the menubar
    def OnUserManual(self,evt):
        text="***SeqProb***\n\n1. Copy and paste your sequence in single-letter amino acid codes. The input is not case sensitive. Try to avoid whitespaces between letters. The program is designed to ignore whitespaces in the sequence input but it's better to avoid them when possible. The string should begin with the residue closest to the N-terminus.\n\n 2. Input however many CA and CB shifts you have determined to be sequential or import shifts from a '.txt' file. The program is set-up to read chemical shifts tabulated using SPARKY then exported to EXCEL meaning they should be arranged in delimited columns. If you open a SPARKY '.list' file in EXCEL then save the data as a '.txt. file only with two consecutive columns (CA, then CB in that order only), the program should work. The delimiter is assumed to be whitespace. For clarity all imported data is displayed before used in calculations. Please double check your data. \n\n 3. Click calculate. Your sequence should appear in a pop-up window with the most probable portion corresponding to your chemical shift input highlighted in blue. Another pop-up window should appear displaying a grid. The grid lists your input sequence along the horizontal axis and all calculated probability scores along the verticle axis. Assignments for input CA and CB shifts are displayed.\n\n Single amino acid probability scores are a product of normal distribution probability densities for each input CA and CB shift relative to average chemical shifts and standard deviations calculated by the BMRB (http://www.bmrb.wisc.edu/ref_info/statful.htm) current as of 10/01/2009.\n\n You can clear your sequence input by clicking the 'Sequence Clear' button.\n\n"
        self.text=text
        win = Help(self, "User Manual", text, style=wx.DEFAULT_FRAME_STYLE | wx.THICK_FRAME)
        win.Show(True)

    #Function to clear sequence input text box when you click the button "Clear Sequence"
    def ClearSequence(self, evt):
        self.bt.Clear()

    #Function to do all of your calculations and display your results when you click the button "Calculate"
    def Calculate(self, evt):
        sequence=self.bt.GetValue() 
        sequence=sequence.replace(' ',"") #Remove whitespace if any present in input sequence       
        self.inputsequence=sequence.upper() #Convert sequence input to uppercase            
    
    #This collects all of your CA and CB shift inputs into ordered lists        
        i = 0
        k = 0
        CAshifts = []
        CBshifts = []
                
    #The following conditionals basically make sure the user input is interpretable...error windows are designed to pop up if they aren't
        if (len((self.shifttable.shifttabledata.data))%2) == 0:
            while i < (len((self.shifttable.shifttabledata.data))/2):
                if (self.shifttable.shifttabledata.data[(i,0)]==''):    #Check for whitespaces      
                    del self.shifttable.shifttabledata.data[(i,0)]
                    break   
                try:
                    float(self.shifttable.shifttabledata.data[(i,0)])   #Check for letters
                except (TypeError, ValueError):  
                    text="ERROR: CA shifts must be numbers!\n\n"
                    self.text=text          
                    win = PopWarning(self, "Error Message", self.text, style=wx.DEFAULT_FRAME_STYLE)
                    win.Show(True)                      
                    break           
                    del self.shifttable.shifttabledata.data[(i,0)]
                else:
                    CAshifts.append(float(self.shifttable.shifttabledata.data[(i,0)]))  
                if (self.shifttable.shifttabledata.data[(i,1)]==''):
                    del self.shifttable.shifttabledata.data[(i,1)]
                    break
                try:
                    float(self.shifttable.shifttabledata.data[(i,1)])   #Check for letters
                except (TypeError, ValueError):                    
                    text="ERROR: CB shifts must be numbers!\n\n"
                    self.text=text          
                    win = PopWarning(self, "Error Message", self.text, style=wx.DEFAULT_FRAME_STYLE)
                    win.Show(True)  
                    break                       
                    del self.shifttable.shifttabledata.data[(i,1)]
                else:
                    CBshifts.append(float(self.shifttable.shifttabledata.data[(i,1)]))
                i = i + 1       
                self.CAshiftlist = CAshifts
                self.CBshiftlist = CBshifts
        else: 
            text="ERROR: The number of CA shifts does not equal the number of CB shifts!\n\n"
            self.text=text          
            win = PopWarning(self, "Error Message", self.text, style=wx.DEFAULT_FRAME_STYLE)
            win.Show(True)
            self.CAshiftlist = []
            self.CBshiftlist = []           

        #This calls a function that returns a list of all probability scores calculated     
        scorel=ScoreL(self, self.inputsequence, self.CAshiftlist, self.CBshiftlist)     
        self.scorel=scorel
    
        #This calls a function that returns a list of all indices for all probability scores calculated     
        scoreindex=ScoreIndex(self, scorel)
        self.scoreindex=scoreindex
    
        #This displays a simple frame that prints your input sequence
        #and highlights the most probable location within your sequence that your input CA and CB shifts could correspond to       
        win = ResultsFrame(self, "Your Protein Sequence",self.inputsequence, self.scoreindex[0], style=wx.DEFAULT_FRAME_STYLE | wx.THICK_FRAME)
        win.Show(True)      
    
        #This displays a detailed results grid showing all locations within the sequence corresponding to all probability scores        
        win2 = RGridFrame(self, self.log, self.inputsequence, self.scorel, self.scoreindex, self.CAshiftlist, self.CBshiftlist)     
        win2.Show(True)

#This defines the simple frame used to display your input sequence
#and highlights the most probable residues that your CA and CB shifts could correspond to		
class ResultsFrame(wx.Frame):
    def __init__(self, parent, title, tseq, maxindex, pos=wx.DefaultPosition, size=(500,300),style=wx.DEFAULT_FRAME_STYLE ):
        self.tseq=tseq
        self.maxindex=maxindex
        wx.Frame.__init__(self, parent, -1, title, pos, size, style)
        panel = wx.Panel(self, -1)
        multitext=wx.TextCtrl(panel, -1, self.tseq, size=(500,300), style=wx.TE_MULTILINE|wx.TE_RICH2|wx.TE_READONLY)
        multitext.SetStyle(maxindex[0], maxindex[-1]+1, wx.TextAttr('white', 'blue'))

class Help(wx.Frame):
    def __init__(self, parent, title, text, pos=wx.DefaultPosition, size=(300,300) ,style=wx.DEFAULT_FRAME_STYLE ):
        wx.Frame.__init__(self, parent, -1, title, pos, size, style)
        panel = wx.ScrolledWindow(self, -1)
        panel.SetScrollbars(1, 1, 500, 500)
        panel.Scroll(1,1)
        multitext=wx.TextCtrl(panel, -1, text, size=(500, 500), style=wx.TE_MULTILINE|wx.TE_RICH2|wx.TE_WORDWRAP|wx.TE_READONLY)    
        bs = wx.BoxSizer(wx.HORIZONTAL) 
        bs.Add(multitext, wx.EXPAND)       
        panel.SetSizer(bs)      
        panel.Fit()     

#This creates a mini-frame that will pop-up when there is a problem with the input 
class PopWarning(wx.Frame):
    def __init__(self, parent, title, text, pos=wx.DefaultPosition, size=(400, 200), style=wx.DEFAULT_FRAME_STYLE):
        wx.Frame.__init__(self, parent, -1, title, pos, size, style)
        panel = wx.Panel(self, -1)
        multitext=wx.TextCtrl(panel, -1, text, size=(400,140), style=wx.TE_MULTILINE|wx.TE_WORDWRAP|wx.TE_READONLY)
        returnb = wx.Button(panel, -1, "Return")
        returnb.SetDefault()
        self.Bind(wx.EVT_BUTTON, self.Return, returnb)
        bs = wx.BoxSizer(wx.VERTICAL)   
        bs.Add(multitext, wx.GROW)
        bs.Add(returnb, 0, wx.ALIGN_CENTER)        
        panel.SetSizer(bs)

    def Return(self, evt):
        self.Close()

###################################################################################################################

#This defines another grid that is used to display all of the results.
class ResultsDataTable(gridlib.PyGridTableBase):
    def __init__(self, log, xseq, xscorel, xscoreindex, xcashiftlist, xcbshiftlist):
        gridlib.PyGridTableBase.__init__(self)

        #Define all of your passed values locally       
        self.xseq=xseq
        self.xscorel=xscorel
        self.xscoreindex=xscoreindex
        self.xcashiftlist=xcashiftlist
        self.xcbshiftlist=xcbshiftlist

        #This stores your sequence such that each amino acid can serve as the label for a column
        list1=[]
        list2=[]
        i = 0   
        self.data={}

        while i < len(xseq):
            list1.append((0, i+1))
            list2.append(xseq[i])           
            self.data.update({list1[i]:list2[i]})
            i = i + 1
            
        #This assigns some constant labels for clarity
        self.data.update({(0,0):'Protein Sequence:'})
        self.data.update({(1,0):'Probability Score'})
        
        #This prints the probability scores onto the grid such that they are double spaced for clarity  
        i = 0       
        while i < len(xscorel):
            self.data.update({(((2*i)+2),0): (self.xscorel[i][0])})
            i = i + 1   
    
        #This prints your amino acid column labels
        i = 1
        self.colLabels = ['']
        colLab=self.colLabels
        while i < (len(xseq)+1):
            colLab.append(i)
            i = i + 1       

        #This loop is probably more complicated than it needs to be but it works.
        #This is used to display your chemical shift inputs, in order, for each probaility score at the corresponding sequence location
        m = 0
        i = 0
        j = int(xscoreindex[i][m])
        k = 0
        while i < len(xscoreindex):     
            k = 0
            m = 0
            while k < len(xcashiftlist):                
                j = int(xscoreindex[i][m])
                while j <= (int(xscoreindex[i][-1])):
                    self.data.update({(((2*i)+2), (j+1)):xcashiftlist[k]})
                    self.data.update({(((2*i)+3), (j+1)):xcbshiftlist[k]})      
                    j = j + 1           
                m = m + 1                   
                k = k + 1                   
            i = i + 1 
            
        #This defines your grid atttributes
        self.odd=wx.grid.GridCellAttr()
        self.odd.SetBackgroundColour("white")
        self.odd.SetFont(wx.Font(10, wx.SWISS, wx.NORMAL, wx.BOLD))
    
        self.even=wx.grid.GridCellAttr()
        self.even.SetBackgroundColour("white")
        self.even.SetFont(wx.Font(10, wx.SWISS, wx.NORMAL, wx.BOLD))

    def GetNumberRows(self):
        return ((len(self.xscorel))*2)+2

    def GetNumberCols(self):
        return (len(self.xseq)+1)

    def IsEmptyCell(self, row, col):
        return self.data.get((row, col)) is not None

    def GetValue(self, row, col):
        value = self.data.get((row, col))
        if value is not None:
            return value
        else:
            return ''

    def SetValue(self, row, col, value):
        self.data[(row, col)] = value

    def GetAttr(self, row, col, kind):
        attr = [self.even, self.odd] [row % 2]
        attr.IncRef()
        attr.SetReadOnly(1)
        return attr
    
    def GetColLabelValue(self, col):
        return self.colLabels[col]

###FROM wxpython DEMO...same as above
class ResultsGrid(gridlib.Grid):
        def __init__(self, parent, log, xseq, xscorel, xscoreindex, xcashiftlist, xcbshiftlist):
            gridlib.Grid.__init__(self, parent, -1)
    
            table = ResultsDataTable(log, xseq, xscorel, xscoreindex, xcashiftlist, xcbshiftlist)
            self.shiftresultstable=table

            # The second parameter means that the grid is to take ownership of the
            # table and will destroy it when done.  Otherwise you would need to keep
            # a reference to it and call it's Destroy method later.
            self.SetTable(table, True)

            self.SetRowLabelSize(0)
            self.SetMargins(0,0)
            self.AutoSizeColumns(False)
 
            gridlib.EVT_GRID_CELL_LEFT_DCLICK(self, self.OnLeftDClick)


            # I do this because I don't like the default behaviour of not starting the
            # cell editor on double clicks, but only a second click.
        def OnLeftDClick(self, evt):
            if self.CanEnableCellControl():
                self.EnableCellEditControl()
###End wxpython DEMO code

#This defines your results frame that will be used to display your results grid
class RGridFrame(wx.Frame):
    def __init__(self, parent, log, xseq, xscorel, xscoreindex, xcashiftlist, xcbshiftlist):

        wx.Frame.__init__(
                self, parent, -1, "Your Probabilty Scores", size=(500, 500)
                )

        #Again, define all passed values locally
        self.log=log
        self.xseq=xseq
        self.xscorel=xscorel
        self.xscoreindex=xscoreindex
        self.xcashiftlist=xcashiftlist
        self.xcbshiftlist=xcbshiftlist              
    
        #Make a panel with scroll bars to store the widget...wx.ScrolledWindow=wx.Panel 
        p = wx.ScrolledWindow(self, -1)
        p.SetScrollbars(1, 1, 500, 500)
        p.Scroll(1,1)
            
        #Add the results grid to the panel
        grid = ResultsGrid(p, log, xseq, xscorel, xscoreindex, xcashiftlist, xcbshiftlist)
        self.resultstable = grid

        #Define a sizer to control the placement of the grid onto the panel (doesn't matter which kind you use since there is only one widget)
        bs = wx.GridSizer(rows=1, cols=1, hgap=5, vgap=5)

        #Add the results grid to the sizer      
        bs.Add(grid, 1, wx.GROW)    
        p.SetSizer(bs)

###Functions
#Basically one function is used to calculate the probablilty densities
#and store them as a list of values in descending order
#All other values returned are calculated from this initial list

def ScoreL(self, sequence, CA, CB):

#chemshift is a dictionary of average chemical shifts and standard deviations from the BMRB (calculated from a total of 3615400 chemical shifts current as of 10/01/2009)

    chemshift = {'A':[[53.168, 2.892], [19.066, 3.235]],         
                 'R':[[56.791, 3.507], [30.711, 2.612]],
                 'D':[[54.689, 2.831], [40.896, 2.691]],
                 'N':[[53.560, 3.732], [38.738, 3.827]],
                 'C':[[58.058, 3.460], [33.232, 6.442]],
                 'E':[[57.337, 3.432], [30.030, 3.352]],
                 'Q':[[56.566, 2.716], [29.208, 2.627]],
                 'H':[[56.530, 3.522], [30.320, 3.296]],
                 'I':[[61.624, 3.449], [38.607, 3.051]],
                 'L':[[55.647, 2.250], [42.268, 2.038]],
                 'K':[[56.945, 3.354], [32.811, 3.083]],
                 'M':[[56.138, 2.340], [32.973, 3.376]],
                 'F':[[58.113, 3.992], [39.980, 3.804]],
                 'P':[[63.345, 4.030], [31.895, 3.364]],
                 'S':[[58.707, 2.900], [63.714, 5.345]],
                 'T':[[62.211, 2.785], [69.586, 6.055]],
                 'W':[[57.747, 5.210], [30.123, 5.065]],
                 'Y':[[58.145, 3.178], [39.321, 3.200]],
                 'V':[[62.505, 3.236], [32.737, 2.430]],
                 'G':[[45.378, 2.289], [1.00, 1.00]]
                }

    #Stops program if a non-amino acid code is entered into the input sequence
    for i in sequence.upper():
        if not chemshift.has_key(i):
            text="ERROR: You have entered an invalid single-letter amino acid code!\n\n"
            self.text=text
            win = PopWarning(self, "Error Message", text, style=wx.DEFAULT_FRAME_STYLE)
            win.Show(True)
            return None 

    incr = len(CA)
    i = 0
    k = 0
    seqstore=[]
    seqincr=[]
    list1=[]
    list2=[]
    self.list2=list2
    sequence.upper()

    while i <= (len(sequence)-int(incr)):   
        while k < int(incr):
            seqstore.append(sequence[k+i])              
            list1.append(k+i)
            k = k + 1
        seqincr.append(seqstore)
        list2.append(list1)
        list1 = []          
        seqstore = []
        k = 0           
        i = i + 1

    i = 0
    k = 0
    j = 0
    prodlist=[]
    scorelist=[]
    testlist=[] 
    
    while (k < len(seqincr)):
        while (i < int(incr)):
            comCA = (1/(chemshift[(seqincr[k][i])][0][1]*math.sqrt(2*math.pi)))*math.exp(-1*(abs((CA[int(i)] - chemshift[(seqincr[k][i])][0][0]))**2)/(2*((float(chemshift[(seqincr[k][i])][0][1]))**2)))
            comCB = (1/(chemshift[(seqincr[k][i])][1][1]*math.sqrt(2*math.pi)))*math.exp(-1*(abs((CB[int(i)] - chemshift[(seqincr[k][i])][1][0]))**2)/(2*((float(chemshift[(seqincr[k][i])][1][1]))**2)))
            prod = comCA*comCB
            prodlist.append(prod)
            i = i + 1
        scorelist.append([sum(prodlist)])
        prodlist=[]
        i = 0
        k = k+1
    
#This creates an ordered list of scores (highest to lowest) and each score's respective index.
#Once stored in the list, scores were set to a very very small number instead of removed.
#This is because identical scores are given the SAME index by default.
#By doing this you can separate identical scores by index.
    
    scolist=[]
    k=0
    while k < len(scorelist):
        scolist.append([(max(scorelist)), ((scorelist.index(max(scorelist))))])
        scorelist[scorelist.index(max(scorelist))]=-10000000000000000000
        k = k + 1
    return scolist

#This creates an ordered list of the locations of each probaility score within the sequence. In the same order as 'scolist' above.

def ScoreIndex(self, scolist):
	i = 0
	scolistseq=[]
	while i < len(scolist):
	    scolistseq.append(self.list2[scolist[i][1]])
	    i = i + 1
	return scolistseq

#############################################################################

if __name__ == '__main__':
    import sys
    app = wx.PySimpleApp()
    frame = TestFrame(None, sys.stdout)
    frame.Show(True)
    app.MainLoop()

#############################################################################
