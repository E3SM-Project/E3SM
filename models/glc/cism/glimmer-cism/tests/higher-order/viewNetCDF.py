#!/usr/bin/env python
# A netCDF file viewer for Glimmer I/O files
# Written by Glen Granzow on March 28, 2010

######## Check to see that there is a command line argument ########
from sys import argv, exit
if len(argv) == 1 or argv[1][0] == '-':
  print 'USAGE:\n'
  print '  python viewNetCDF.py filename [-m]\n'
  print '    -m   Place additional options in a menubar at the top of the window.'
  print '         The default is to have these options in a pop-up menu which is'
  print '         opened using the right mouse button.\n'
  exit(0)

######## Import the required Python modules ########
import numpy, thread
from netCDF import *
from Tkinter import *
from math import log, exp
from matplotlib import pyplot
import matplotlib.cm

######## Get and sort the variable names in the netCDF file ########
netCDFfile = NetCDFFile(argv[1],'r')
variables = netCDFfile.variables.keys()
variables.sort()

######## Get and sort the available color map names ########
names = matplotlib.cm.datad.keys()
lowercasenames = list()
mixedcasenames = list()
for name in names:
  if name.islower():
    lowercasenames.append(name)
  else:
    mixedcasenames.append(name)
lowercasenames.sort()
mixedcasenames.sort(key=(lambda x: x.lower()))
i = j = 0
sortednames = list()
while i < len(lowercasenames) and j < len(mixedcasenames):
  if lowercasenames[i] < mixedcasenames[j].lower():
    sortednames.append(lowercasenames[i])
    i += 1
    if i < len(lowercasenames) and lowercasenames[i] == lowercasenames[i-1]+'_r':
      sortednames.append(lowercasenames[i])
      i += 1
  else:
    sortednames.append(mixedcasenames[j])
    j += 1
    if j < len(mixedcasenames) and mixedcasenames[j] == mixedcasenames[j-1]+'_r':
      sortednames.append(mixedcasenames[j])
      j += 1
colormaps = sortednames

######## Code to return a discretized version of a colormap ########
# From http://www.scipy.org/Cookbook/Matplotlib/ColormapTransformations 
def cmap_discretize(cmap, N):
  import scipy.interpolate
  cdict = cmap._segmentdata.copy()
  colors_i = numpy.linspace(0,1,N)
  indices  = numpy.linspace(0,1,N+1)
  for key in ('red','green','blue'):
    D = numpy.array(cdict[key])
    I = scipy.interpolate.interp1d(D[:,0], D[:,1])
    colors = I(colors_i)
    A = numpy.zeros((N+1,3), float)
    A[:,0] = indices
    A[1:,1] = colors
    A[:-1,2] = colors
    L = []
    for l in A:
      L.append(tuple(l))
      cdict[key] = tuple(L)
  return matplotlib.colors.LinearSegmentedColormap('colormap',cdict,1024)

######## A parsing tool for entry widgets ########
def parse(entrywidget,dtype,default):
  try:
    value = dtype(entrywidget.get())
  except:
    value = default
  entrywidget.delete(0,END)
  entrywidget.insert(0,str(value))
  return value

######## A function to reset the plotting parameters ########
def reset(*args):
  cmaplist.selection_clear(0,END)
  cmaplist.selection_set(colormaps.index('jet'))
  cmaplist.see(colormaps.index('jet'))
  for widget in entryWidgetsList:
    widget.delete(0,END)
  clearfigure.set(1)
  colorbar.set(1)
  logarithm.set(0)
  animation.set(0)
  delayentry.delete(0,END)
    
######## A function to clear the plot window ########
def clf():
  pyplot.clf()
  pyplot.show()

######## A function to plot in a new window ########
def new():
  pyplot.figure()
  pyplot.gcf().canvas.mpl_connect('button_press_event',getValue)
  animationloop()

######## A function to get a value from a plot ########
def getValue(event):
  x, y = map(int,map(round,[event.xdata, event.ydata]))
  if len(numpy.shape(data)) == 1:
    print 'The value of',v,'at',x,'is',data[x]
  else:
    print 'The value of',v,'at ('+str(x)+','+str(y)+') is',data[y,x]

######## A function to loop for animation ########
ANIMATION_IS_DISABLED = True
def animationloop():
# If animation is not selected just do one plot
  if ANIMATION_IS_DISABLED or animation.get() == 0:
    plot()
    return
# Otherwise loop over time or level
# This never worked; I think matplotlib's use of Tkinter caused problems
  import time
  use_thread = False
  v = variables[int(variablelist.curselection()[0])]
  time_or_level = animation.get()-1
  if len(netCDFfile.variables[v].dimensions)<(3+time_or_level):
    d = ('time','level')[time_or_level]
    print 'WARNING:',v,'does not have a',d,'dimension;'
    print '         Animation was aborted.'
    plot()
    return
  entryWidget = (timeentry,levelentry)[time_or_level]
  start  = parse(entryWidget,int,0)
  finish = netCDFfile.variables[v].shape[time_or_level]
  if start >= finish-1:
    start = 0
  delay = parse(delayentry,float,0.0)
  plotbutton.config(text='STOP',fg='red',command=stopanimation)
  pyplot.clf()
  print 'entering actual loop',start,finish
  if use_thread:
    thread.start_new_thread(actualloop,(start,finish,entryWidget,delay))
  else:
    actualloop(start,finish,entryWidget,delay)
  print 'Hello!'

def actualloop(start,finish,entryWidget,delay):
  global stopflag
  stopflag = False
  for i in xrange(start,finish):
    entryWidget.delete(0,END)
    entryWidget.insert(0,str(i))
    plot()
    if stopflag: break
    time.sleep(delay)
  plotbutton.config(text='Plot',fg='black',command=animationloop)

def stopanimation():
  global stopflag
  stopflag = True
    
######## The function that does the plotting ########
plotCounter = 0
def plot():
  global plotCounter, v, data
  v = variables[int(variablelist.curselection()[0])]
  dimensions = len(numpy.array(netCDFfile.variables[v].dimensions))
  if clearfigure.get():
    pyplot.clf()
  if dimensions == 1:
    data = numpy.array(netCDFfile.variables[v][:])
    plotfunction = (pyplot.plot,pyplot.semilogy)[logarithm.get()]
    plotfunction(data)
  else:
    indices = [':',':']
    if dimensions > 2: 
      time = parse(timeentry,int,0)
      if animation.get() == 1: # Increment time
        time += 1
        timeentry.delete(0,END)
        timeentry.insert(0,str(time))
      indices.insert(0,str(time))
    if dimensions > 3:
      level = parse(levelentry,int,0)
      if animation.get() == 2: # Increment level
        level += 1
        levelentry.delete(0,END)
        levelentry.insert(0,str(level))
      indices.insert(1,str(level))
    indices = ','.join(indices)
    data = numpy.array(eval('netCDFfile.variables[v]['+indices+']'))
    mask = (data != data) # look for NaN
    if numpy.any(mask):
      print 'Masking',numpy.sum(mask),'NaN values in',v
      data = numpy.ma.masked_array(data,mask)
    if logarithm.get() == 1:
      mask0 = (data <= 0)
      if numpy.any(mask0):
        print 'Masking',numpy.sum(mask0),'values in',v,'that are less than or equal to zero.'
        data = numpy.ma.masked_array(data,mask0)
      if parse(vminentry,float,0) <= 0: vminentry.delete(0,END)
      if parse(vmaxentry,float,0) <= 0: vmaxentry.delete(0,END)
    vmin = parse(vminentry,float,numpy.amin(data))
    vmax = parse(vmaxentry,float,numpy.amax(data))
    colormap = eval('pyplot.cm.'+colormaps[int(cmaplist.curselection()[0])])
    norm = (None,matplotlib.colors.LogNorm())[logarithm.get()]
    ticks = None
    n = parse(contoursentry,int,0)
    if n > 1:
      colormap = cmap_discretize(colormap,n)
      ticks = vmin + numpy.arange(n+1) * (vmax-vmin)/n
      if logarithm.get() == 1:
        ticks = map(exp,log(vmin) + numpy.arange(n+1) * (log(vmax)-log(vmin))/n)
    pyplot.subplot(*subplot)
    pyplot.imshow(data,origin=origin.get(),interpolation='nearest',vmin=vmin,vmax=vmax,cmap=colormap,norm=norm)
    if colorbar.get():
      cb = pyplot.colorbar(ticks=ticks)
      if n > 1 and logarithm.get() == 1:
        cb.ax.set_yticklabels(['%g' % t for t in ticks])
  if hasattr(netCDFfile.variables[v],'long_name'):
    title = netCDFfile.variables[v].long_name
  elif hasattr(netCDFfile.variables[v],'standard_name'):
    title = netCDFfile.variables[v].standard_name
  else:
    title = v
  if hasattr(netCDFfile.variables[v],'units'):
    title += ' ('+netCDFfile.variables[v].units+')'
  pyplot.title(title)
  plotCounter += 1
  if suptitle != None:
    pyplot.suptitle(suptitle)
  elif plotCounter == 1:
    pyplot.suptitle('Some systems may require you to close this first window to continue...')
  elif plotCounter == 2:
    pyplot.suptitle('This second window may not close properly...')
  if plotCounter == 1:
    pyplot.gcf().canvas.mpl_connect('button_press_event',getValue)
  pyplot.show()

######## A function for directional (1D) plotting ########
def directionalplot():
  DirectionalPlotWindow()
  
class DirectionalPlotWindow:
  def __init__(self):
    window = Toplevel()
    window.title('1D plots')
    self.entryList = list()
    self.direction = IntVar()
    for (i,label) in enumerate(['time','level','y','x']):
      Radiobutton(window,text=label,variable=self.direction,value=i).grid(row=i,column=0,sticky=W)
      self.entryList.append(Entry(window,width=12))
      self.entryList[i].grid(row=i,column=1,columnspan=2,padx=5,sticky=W)
    Button(window,text='Clear',width=5,command=self.clear ).grid(row=4,column=0,sticky=E,padx=2)
    Button(window,text='Plot', width=5,command=self.plot1D).grid(row=4,column=1,sticky=W,padx=2,pady=5)
    Button(window,text='New',  width=5,command=self.new1D ).grid(row=4,column=2,sticky=W,padx=2)
    
  def clear(self):
    clf()
    for entry in self.entryList:
      entry.delete(0,END)

  def new1D(self):
    pyplot.figure()
    self.plot1D()

  def plot1D(self):
    indices = [str(parse(entry,int,0)) for entry in self.entryList]
    indices[self.direction.get()] = ':'
    v = variables[int(variablelist.curselection()[0])]
    dimensions = len(numpy.array(netCDFfile.variables[v].dimensions))
    if dimensions == 1:
      indices = ':'
    if dimensions == 2:
      indices = indices[2:]
    if dimensions == 3:
      indices.pop(1)
    indices = ','.join(indices)
    data = numpy.array(eval('netCDFfile.variables[v]['+indices+']'))
    plotfunction = (pyplot.plot,pyplot.semilogy)[logarithm.get()]
    plotfunction(data)
    pyplot.xlabel(('time','level','y','x')[self.direction.get()])
    pyplot.ylabel(v)
    pyplot.show()
    
######## A function to print variable attributes ########
ignore = None
def printattributes():
  global ignore
  if ignore == None:
#   Initialization (occurs the first time printattributes is called)
    from glob import glob
    from os import remove
    verbose = False
#   Find a filename that does not already exist
    for i in xrange(1000000):
      filename = 'viewNetCDF'+str(i)+'.trash'
      if filename not in glob('*.trash'): break
    if verbose:
      print 'Creating temporary file:',filename
    trash = NetCDFFile(filename,'w')
    print '\nGLOBAL ATTRIBUTES OF',argv[1]
    for attribute in dir(netCDFfile):
      if attribute not in dir(trash):
        print attribute +":", getattr(netCDFfile,attribute )
    variable = trash.createVariable('trash','f',tuple())
    ignore = dir(variable)
    trash.close()
    if verbose:
      print 'Deleting:',filename
    remove(filename)

  v = variables[int(variablelist.curselection()[0])]
  print '\nATTRIBUTES OF',v
  v = netCDFfile.variables[v]
  for attribute in dir(v):
    if attribute not in ignore:
      print attribute+':',getattr(v,attribute)

######## A function to set the super title ########
suptitle = None
def setSuptitle():
  thread.start_new_thread(threadedSetSuptitle,tuple())

def threadedSetSuptitle():
  global suptitle
  suptitle = raw_input('Enter the title: ')

######## A function to set subplot parameters ########
subplot = (1,1,1)
def setSubplot():
  thread.start_new_thread(threadedSetSubplot,tuple())

def threadedSetSubplot():
  global subplot
  subplot = raw_input('Enter nrows, ncolumns, row, column: ').strip()
  if subplot.find(',') > 0: 
    subplot = subplot.split(',')
  elif subplot.find(' ') > 0: 
    subplot = subplot.split(' ')
  subplot = map(int,subplot)
  if len(subplot) > 3:
    subplot = subplot[:2] + [(subplot[2]-1)*subplot[1]+subplot[3]]

######## Create the GUI window ########
root = Tk()
root.title(argv[1])
leftframe  = Frame(root)
rightframe = Frame(root)

######## Create a GUI frame containing the list of netCDF variables ########
variablesframe = LabelFrame(leftframe,text='Variable',labelanchor=N)
variablelist = Listbox(variablesframe,selectmode=SINGLE,height=10,exportselection=0)
for v in variables:
  variablelist.insert(END,v)
variablelist.selection_set(0)
scrollbar = Scrollbar(variablesframe,command=variablelist.yview)
variablelist.config(yscrollcommand=scrollbar.set)
scrollbar.pack(side=RIGHT, fill=Y)
variablelist.pack()

######## Create a GUI frame containing the list of available color maps ########
cmapframe = LabelFrame(leftframe,text='Color Map',labelanchor=N)
cmaplist = Listbox(cmapframe,selectmode=SINGLE,height=5,exportselection=0)
for c in colormaps:
  cmaplist.insert(END,c)
cmaplist.selection_set(colormaps.index('jet'))
cmapscrollbar = Scrollbar(cmapframe,command=cmaplist.yview)
cmaplist.config(yscrollcommand=cmapscrollbar.set)
cmaplist.see(colormaps.index('jet'))
cmapscrollbar.pack(side=RIGHT, fill=Y)
cmaplist.pack()

######## Create a GUI frame containing entry widgets ########
entryframe = Frame(rightframe)
entryWidgetsList = list()
def createEntryWidget(label):
  row    = len(entryWidgetsList)
  label  = Label(entryframe,text=label)
  widget = Entry(entryframe,width=12)
  label.grid(column=0,row=row,sticky=E)
  widget.grid(column=1,row=row)
  entryWidgetsList.append(widget)
  return widget
timeentry  = createEntryWidget('time:')
levelentry = createEntryWidget('level:')
vminentry  = createEntryWidget('vmin:')
vmaxentry  = createEntryWidget('vmax:')
contoursentry = createEntryWidget('contours:')

######## Create a GUI frame containing buttons ########
buttonframe = Frame(rightframe)
resetbutton = Button(buttonframe,text='Reset',command=reset)
plotbutton  = Button(buttonframe,text='Plot',command=animationloop)
newbutton   = Button(buttonframe,text='New',command=new)
resetbutton.pack(side=LEFT,padx=1,pady=4)
plotbutton.pack(side=LEFT,padx=1,pady=4)
newbutton.pack(side=LEFT,padx=1,pady=4)
resetbutton.config(width=4)
plotbutton.config(width=4)
newbutton.config(width=4)

######## Create a GIU frame containing check buttons ########
checkbuttonframe = Frame(rightframe)
clearfigure = IntVar()
colorbar    = IntVar()
logarithm   = IntVar()
Checkbutton(checkbuttonframe,text='clear before plotting',variable=clearfigure).pack(anchor=W)
Checkbutton(checkbuttonframe,text='draw colorbar',variable=colorbar).pack(anchor=W)
Checkbutton(checkbuttonframe,text='log scale',variable=logarithm).pack(anchor=W)
clearfigure.set(1)
colorbar.set(1)

######## Create a GUI frame for animation ########
animationframe = LabelFrame(rightframe,text='Increment',labelanchor=N)
animation = IntVar()
nonebutton  = Radiobutton(animationframe,text='none',variable=animation,value=0)
timebutton  = Radiobutton(animationframe,text='time',variable=animation,value=1)
levelbutton = Radiobutton(animationframe,text='level',variable=animation,value=2)
label = Label(animationframe,text='delay:')
delayentry = Entry(animationframe,width=10)
nonebutton.grid(row=0,columnspan=2,sticky=W)
timebutton.grid(row=1,columnspan=2,sticky=W)
levelbutton.grid(row=2,columnspan=2,sticky=W)
if not ANIMATION_IS_DISABLED:
  label.grid(row=3,column=0,sticky=E)
  delayentry.grid(row=3,column=1,padx=5,pady=2,sticky=W)

######## Create a menubar or pop up menu ########
menu = Menu(root,tearoff=0)
menu.add_command(label='clear figure',command=clf)
menu.add_command(label='directional plots',command=directionalplot)
menu.add_command(label='print attributes',command=printattributes)
menu.add_command(label='set title',command=setSuptitle)
menu.add_command(label='set subplot',command=setSubplot)
#menu.add_command(label='get value',command=getValue)
origin = StringVar()
menu.add_checkbutton(label='invert image',var=origin,onvalue='upper',offvalue='lower')

reuseparameters = IntVar()
# If the value of reuseparameters changes, toggleBinding is called
def toggleBinding(*args):
  if reuseparameters.get() == 0:
    variablelist.bind('<Button-1>',reset)
  else:
    variablelist.unbind('<Button-1>')
reuseparameters.trace('w',toggleBinding)
reuseparameters.set(0)
# Clicking on the checkbutton will change the value of reuseparameters
menu.add_checkbutton(label='reuse plot parameters',var=reuseparameters)

menubar = IntVar()
# If the value of menubar changes, toggleMenubar is called
def toggleMenubar(*args):
  if menubar.get() == 0:
    root.config(menu=Menu(root)) # Make the menubar an empty menu
    menu.add_command(label='close this menu')
    root.bind('<Button-3>',popup)
  else:
    menu.delete(8) # Remove the 'close this menu' option
    root.config(menu=menu) # Put the menu in a menubar
    root.unbind('<Button-3>') # Disable the popup
def popup(event):
  menu.post(event.x_root, event.y_root)
menubar.trace('w',toggleMenubar)
# Clicking on the checkbutton changes the value of menubar
menu.add_checkbutton(label='menu bar',var=menubar)
menubar.set(0)
if '-m' in argv:  # Check if -m was on the command line
  menubar.set(1)

#menu.add_command(label='quit',command=root.quit)

######## Layout the GUI ########
layout = ('pack','grid')[1]
if layout=='pack':
  variablesframe.pack(side=TOP,padx=5,pady=0,anchor=N)
  cmapframe.pack(side=TOP,padx=5,pady=2)
  leftframe.pack(side=LEFT)
  entryframe.pack(side=TOP,padx=5,pady=2)
  checkbuttonframe.pack(side=TOP,pady=0)
  animationframe.pack(side=TOP)
  buttonframe.pack(side=TOP,pady=2)
  rightframe.pack(side=LEFT)
else:
  variablesframe.grid(row=0,column=0,padx=5,pady=0)
  cmapframe.grid(row=1,column=0,padx=5,pady=4)
  leftframe.grid(row=0,column=0)
  entryframe.grid(row=0,column=0,padx=5,pady=5)
  checkbuttonframe.grid(row=1,column=0)
  animationframe.grid(row=2,column=0)
  buttonframe.grid(row=3,column=0)
  rightframe.grid(row=0,column=1)

######## Display the GUI window ########
mainloop()
netCDFfile.close()
print plotCounter,'plots were produced'
