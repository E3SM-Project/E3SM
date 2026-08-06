import sys
import os
import getopt
import copy

class toleranceClass:
  def __init__(self):
    self.tol = 0.0

  def checkTol(self,num1,num2):
    denom = max([abs(num1),abs(num2)])
    num = abs(num1-num2)
    relDiff = num/denom
    if relDiff > self.tol:
      self.tol = relDiff

usage="Usage: Difference two ascii files \n\
  python diff.py file1 file2\n\
  python diff.py --maxTol=1.e3 file1 file2"

def openFile(fileName):
  try:
    fid = open(fileName,'r')
  except IOError as e:
     print "File ", fileName, " cannot be found"
     print usage
     sys.exit(-1)
  return fid,fileName.strip()

def is_number(s):
  """
  Take from: 
  http://stackoverflow.com/questions/354038/how-do-i-check-if-a-string-is-a-number-in-python
  """
  try:
    float(s)
    return True
  except ValueError:
    return False

if __name__=="__main__":

  if (len(sys.argv) < 3):
    print usage
    sys.exit(1)

  try:
    opts, args = getopt.getopt(sys.argv[1:], '', ['maxTol='])
  except getopt.GetoptError as err:
    # print help information and exit:
    print str(err) # will print something like "option -a not recognized"
    print usage
    sys.exit(2)

  maxTol = 0.0
  for option, arg in opts:
    if option == "--maxTol":
      if is_number(arg):
        maxTol = float(arg)
      else:
        print "Tolerance is not a number"
        sys.exit(2) 
        
    else:
      print usage
      sys.exit(3)

  # open the files
  file1,fileName1 = openFile(sys.argv[-2])
  file2,fileName2 = openFile(sys.argv[-1])

  #print "Running difference on", fileName1, ",", fileName2

  # Now loop through the files
  file1EOF = False
  file2EOF = False

  tolerance = toleranceClass()

  while ( not file1EOF and not file2EOF):
   
    # Read a line from file1
    line1 = file1.readline()
    line2 = file2.readline()

    if not line1:
      file1EOF = True

    if not line2:
      file2EOF = True

    # 
    if (file1EOF and not file2EOF) or (not file1EOF and file2EOF):
      print "Files are different"
      sys.exit(2) 

    # Now loop through both of the lines
    strings1 = line1.split()
    strings2 = line2.split()
    
    if len(strings1) != len(strings2):
      print "Files are different"
      sys.exit(3) 
    else:
      # loop over the items in the range
      for index in range(len(strings1)):
        if strings1[index] == strings2[index]:
          continue
        elif is_number(strings1[index]) and is_number(strings2[index]):
          tolerance.checkTol(float(strings1[index]),float(strings2[index]))
        else:    
          print "Files are different"
          sys.exit(4) 
  
  # Test to determine the level of similarity
  if tolerance.tol == 0.0:
    print "Files are identical"
    sys.exit(0) 
  elif tolerance.tol < maxTol:
    print "Files are similar with differences", tolerance.tol, "less than", maxTol
    sys.exit(0) 
  else:
    print "Files are different at tolerance", maxTol
    sys.exit(-1) 
