#!/usr/bin/env python

import sys
import re

def IsContinued(string):
  """Indicate if a line of Fortran ends with a continuation character."""
  string = string.rstrip()
  if string[len(string)-1] == '&':
    return True
  else:
    return False

def GetContinuation(string):
  """Concatenate lines as long as continuation characters encountered."""
  if IsContinued(string):
    nextline = sys.stdin.readline()
    while IsContinued(nextline):
      nextline += sys.stdin.readline()
    return string + nextline
  else:
    return string

def IsComment(string):
  """Indicate if a string is a Fortran comment line or not."""
  if len(string.lstrip()) == 0:
    return False
  elif (string[0] == 'C' or string[0] == 'c' or string.lstrip()[0] == '!'):
    return True
  else:
    return False

def StripComments(string):
  if IsComment(string):
    # The entire line is a comment, so return the empty string.
    return ""
  else:
    # Return everything before the comment character, if there is one.
    # Note that, as written, this will break if there is a '!' inside 
    # of a text string.
    return string.split('!')[0]

def ExtractSymbols(lines):
  
  # Concatenate everything into one big string.
  bigline = ""
  for i in range(len(lines)):
    bigline = bigline + StripComments(lines[i])
  
  # Take out everything that isn't a variable name.
  bigline = bigline.replace("pointer", "")
  bigline = bigline.replace("real*8", "")
  bigline = bigline.replace("PetscScalar", "")
  bigline = bigline.replace("Scalar", "")
  bigline = bigline.replace("&", " ")
  bigline = bigline.replace("(:)", " ")
  bigline = bigline.replace(",", " ")
  bigline = bigline.replace("::", " ")

  symbols = bigline.split()

  return symbols

# We need the !PETSC_ARRAYF90_MACROS functionality when pointers defined 
# in a module are used outside of that module with ArrayF90 statements.
# In that case, we don't need to define any variables in the code outside 
# the module, but that code does need to have the proper macros defined.
def PrintMacros(line):
  line = line.replace("!PETSC_ARRAYF90_MACROS","")
  symbols = line.split()
  for symbol in symbols:
    line = "#define "
    line += symbol + "(ib) "
    line += symbol+"_v("+symbol+"_i + (ib))"
    print line


def PrintDummyArrays(lines):
  symbols = ExtractSymbols(lines)

  # Print the _v 1-D real*8 arrays that will be used to access the values.
  for line in lines:
    line = line.replace("pointer", "dimension(1), target")
    line = line.replace("(:)", "")
    for symbol in symbols:
      # We use regular expressions here to make sure that we don't do 
      # replacement on symbols inside other words.
      line = re.sub(r"\b" + symbol + r"\b", symbol + "_v", line)
    print line,

  # Print the _i arrays that will be used to index the _v arrays.
  regexp = re.compile(r"Scalar,|PetscScalar,|Scalar|PetscScalar")
  for line in lines:
    line = line.replace("pointer", "")
    line = line.replace("(:)", "")
    line = line.replace("real*8,", "PetscOffset")
    line = line.replace("real*8", "PetscOffset")
    line = regexp.sub("PetscOffset", line)
      # Had to use regular expressions above, since the code might 
      # use either "Scalar" or "PetscScalar", which causes problems 
      # if using the simple 'replace' method.
    for symbol in symbols:
      # We use regular expressions here to make sure that we don't do 
      # replacement on symbols inside other words.
      line = re.sub(r"\b" + symbol + r"\b", symbol + "_i", line)
    print line,

  # Print the macros that will be used to access the arrays.
  for symbol in symbols:
    line = "#define "
    line += symbol + "(ib) "
    line += symbol+"_v("+symbol+"_i + (ib))"
    print line


# Gets the argument list from a GetArrayF90() or RestoreArrayF90() call.
def GetArgList(line):
  argstring = line.split("ArrayF90(")[1]
  argstring = argstring.split(")")[0]
  argstring = argstring.replace(")", "")
  argstring = argstring.replace(",", " ")
  return argstring.split()
  

# Given a GetArrayF90 statement (concatenated into one line), convert 
# it into the F77-style statement.
def TranslateGetArrayF90(line):
  workline = line.replace("&", "")
  arglist = GetArgList(workline)
  workline = workline.replace("ArrayF90", "Array")
  newline = workline.split('GetArray(')[0] + 'GetArray('
  remainderline = workline.split('GetArray(')[1]
  remainderline = remainderline.split(')',1)[1]

  newline += arglist[0] + ','
  newline += arglist[1] + '_v,'
  newline += arglist[1] + '_i,'
  newline += arglist[2] + ')'
  newline += remainderline
  
  # If the length of the new line exceeds the 132 characters allowed 
  # by Fortran in one line, then print a warning.
  # If I get the time, I should make this split the line using 
  # '&' continuation.  It's probably not necessary for my purposes, though.
  if len(newline) > 132:
    sys.stderr.write("WARNING: A GetArray line exceeds 132 characters")

  print newline


# Given a RestoreArrayF90 statement (concatenated into one line), 
# convert it into the F77-style statement and undefine the macros 
# that were used to access the arrays.
def TranslateRestoreArrayF90(line):
  workline = line.replace("&", "")
  arglist = GetArgList(workline)
  workline = workline.replace("ArrayF90", "Array")
  newline = workline.split('RestoreArray(')[0] + 'RestoreArray('
  remainderline = workline.split('RestoreArray(')[1]
  remainderline = remainderline.split(')',1)[1]

  newline += arglist[0] + ','
  newline += arglist[1] + '_v,'
  newline += arglist[1] + '_i,'
  newline += arglist[2] + ')'
  newline += remainderline
  
  # If the length of the new line exceeds the 132 characters allowed 
  # by Fortran in one line, then print a warning.
  # If I get the time, I should make this split the line using 
  # '&' continuation.  It's probably not necessary for my purposes, though.
  if len(newline) > 132:
    sys.stderr.write("WARNING: A RestoreArray line exceeds 132 characters")

  print newline
  if no_undefs != True:
    print "#undef", arglist[1]


def main(argv=None):
  if argv is None:
    argv = sys.argv
  
  # We need an option to specify that appropriate macros NOT be undefined 
  # after a VecRestoreArrayF90().  Ideally we could undefine the macro, 
  # but this causes problems if preprocessor #ifdef's are being used:
  # A macro might be undefined, and then used in code that will actually 
  # be ignored by the preprocessor.  However, even though that code would 
  # actually not be compiled, the compiler will complain at us.
  global no_undefs
  no_undefs = False
  
  cmd_offset = 0
  if(sys.argv[1].find("--") > -1):
    cmd_offset = 1
    no_undefs = True

  if len(argv) > 3 + cmd_offset:
    sys.stderr.write("Usage: remove_ArrayF90.py [--no_undefs] infile outfile")
    return 1
  if len(argv) > 1 + cmd_offset:
    infile = open(sys.argv[1 + cmd_offset], 'r')
    sys.stdin = infile
  if len(argv) == 3 + cmd_offset:
    outfile = open(sys.argv[2 + cmd_offset], 'w')
    sys.stdout = outfile

  rgxp_getarray = re.compile("GetArrayF90", re.IGNORECASE)
  rgxp_restorearray = re.compile("RestoreArrayF90", re.IGNORECASE)
  rgxp_beginF90 = re.compile("!BEGIN_PETSC_ARRAYF90", re.IGNORECASE)
  rgxp_endF90 = re.compile("!END_PETSC_ARRAYF90", re.IGNORECASE)
  rgxp_macroF90 = re.compile("!PETSC_ARRAYF90_MACROS", re.IGNORECASE)
  
  while True:
    line = sys.stdin.readline()
    if line == "":
      return 0

    if rgxp_getarray.search(StripComments(line)):
      line = GetContinuation(line)
      TranslateGetArrayF90(line)
    elif rgxp_restorearray.search(StripComments(line)):
      line = GetContinuation(line)
      TranslateRestoreArrayF90(line)
    elif rgxp_macroF90.search(line):
      PrintMacros(line)
    elif rgxp_beginF90.search(line):
      lines = []
      while True:
        line = sys.stdin.readline()
        if rgxp_endF90.search(line):
          break
        lines.append(line)
      PrintDummyArrays(lines)
    else:
      print line,
    

if __name__ == "__main__":
  sys.exit(main())
