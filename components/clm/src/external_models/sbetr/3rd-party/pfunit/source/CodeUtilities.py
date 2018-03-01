#!/usr/bin/env python
# For python 2.6-2.7
from __future__ import print_function
# For python2.5
# from __future__ import with_statement

from Utilities import *
import textwrap

class module:
    def __init__(self, name):
        # what do we need?
        self.name = ''
        self.declarations = []
        self.implementations = []
        self.generation = []
        #
        self.name = name
        self.fileName = name+'.F90'
        return 
    def generate(self):
        generation = [ 'module '+self.name]
        generation.extend( [ i.generate() for i in self.declarations ] )
        generation.extend([ 'contains' ])
        generation.extend( [ i.generate() for i in self.implementations ] )
        generation.extend([ 'end module '+self.name])
        return generation
    def addDeclaration(self,declaration):
        # print('adding declaration: ',declaration)
        if type(declaration) is list :
            self.declarations.extend(declaration)
        else:
            self.declarations.append(declaration)
        return self
    def addImplementation(self,implementation):
        self.implementations.append(implementation)
        return self
    def addRoutineUnit(self, rUnit,expose=False):
        # Might need to add more than one decl.
        self.addDeclaration(rUnit.getDeclarations(expose=expose))
        self.addImplementation(rUnit.getImplementation())
        return self
    def addInterfaceBlock(self, interface,expose=False):
        self.addDeclaration(interface.getDeclaration(expose=expose))
        self.addImplementation(interface.getImplementation())
        return self
    def getName(self):
        return self.name
    def setFileName(self,fName):
        self.fileName = fName
        return self
    def getFileName(self):
        return self.fileName

class declaration:
    def __init__(self,name,simpleDeclaration):
        self.simpleDeclaration = simpleDeclaration
        self.fullDeclaration = ''
        self.name = name
        return
    def generate(self):
        return self.simpleDeclaration

class implementation:
    def __init__(self,name,source):
        self.name = name
        self.source = source
    def generate(self):
        return self.source

class routineUnit:
    def __init__(self,name,implementSource):
        self.name = name
        self.declaration = declaration(self.name, self.name) # get better later
        self.declarations = []
        self.declarations.append(self.declaration)
        self.implementation = implementation(self.name,implementSource)
        return
    def setName(self,name):
        self.name = name
        return
    def getName(self):
        return self.name
    def setDeclaration(self,declaration):
        self.declaration = declaration
        self.declarations = [self.declaration]
        return
    def addDeclaration(self,declaration):
        self.declarations.append(declaration)
        return
    def setImplementation(self,implementationSource):
        self.implementation = implementation(self.name, implementationSource)
        return
    def getDeclaration(self,expose=False):
        return self.declaration
    def getDeclarations(self,expose=False):
        return self.declarations
    def getImplementation(self):
        return self.implementation
    def clearDeclarations(self):
        self.declarations = []
        self.declaration = ''
        return self

class interfaceBlock:
#    name = ''
#    moduleProcedureAlternatives = []
#    moduleProcedureImplementations = []
    def __init__(self,name):
        self.name = name
        self.moduleProcedureAlternatives = []
        self.moduleProcedureImplementations = []
    def generateDeclaration(self):
        retStr = '\ninterface ' + self.name + '\n'
        if self.moduleProcedureAlternatives != [] :
            # Note that we need to treat the first line as a special case.
            retStr += \
                '\n   module procedure ' + \
                '\n   module procedure '.join(self.moduleProcedureAlternatives)
        retStr += '\n\nend interface ' + self.name + '\n'
        return declaration(self.name,retStr)

    def generateItemizedDeclarations(self):
        retStr = '\n!interface ' + self.name + '\n'
        if self.moduleProcedureAlternatives != [] :
            # Note that we need to treat the first line as a special case.
            retStr += \
                '\n   public :: ' + \
                '\n   public :: '.join(self.moduleProcedureAlternatives)
        retStr += '\n\n!end interface ' + self.name + '\n'
        return declaration(self.name,retStr)

    def generateImplementation(self):
        retStr = '! interface ' + self.name + ' implementations\n'
        if self.moduleProcedureImplementations != [] :
            retStr += \
                '\n   '.join(self.moduleProcedureImplementations)
        retStr += '\n! end interface ' + self.name + ' implementations'
        return implementation(self.name,retStr)
    def addModuleProcedureAlternative(self,newName):
        self.moduleProcedureAlternatives.append(newName)
        return self
    def addRoutineUnit(self,routineUnit,expose=False):
        for d in routineUnit.getDeclarations(expose=expose):
            self.addModuleProcedureAlternative(d.generate())
        # self.addModuleProcedureAlternative(routineUnit.getDeclaration().generate())
        self.moduleProcedureImplementations.append(routineUnit.getImplementation().generate())
        return self
    def getDeclaration(self,expose=False):
        if expose :
            return self.generateItemizedDeclarations()
        else :
            return self.generateDeclaration()
    def getImplementation(self):
        return self.generateImplementation()
    
class fortranSubroutineSignature:
    def __init__(self,name):
        self.name = name
        self.ArgumentToFType = {}
        self.ReturnFType = ''
        self.SubroutineType = 'subroutine'
    def setReturnFType(self,ReturnFType):
        self.ReturnFType = ReturnFType
        if not ReturnFType :
            self.SubroutineType = 'function'
        return self
    def addArg(self,arg,fType):
        self.ArgumentToFType[arg] = fType
        return self
    def generateInterfaceEntry(self):
        print('generateInterfaceEntryNotImplemented_'+this.name)
        return
    def generateImplementationSignature(self):
        print('generateImplementationSignatureNotImplemented_'+this.name)
        return
    def generateImplementationClose(self):
        print('generateImplementationCloseNotImplemented_'+this.name)
        return
    
def indentKluge(indentString,txt):
    wrapper = textwrap.TextWrapper(initial_indent=indentString, subsequent_indent=indentString);
    txtList=map(wrapper.fill,str.splitlines(txt))
    return "\n".join(txtList)

def iterateOverMultiRank(nr,variableName,shapeName,centralText):
    indent= str(' '*(3*(nr+1)));
    txt = indentKluge(indent,centralText)
    r = range(nr); rrev = range(nr); rrev.reverse();
    codeSnippet = ''.join([' '*(3*(nr-i))+'do '+variableName+str(i+1)+'= 1,'+shapeName+'('+str(i+1)+')\n' for i in rrev])+txt+'\n'+''.join([' '*(3*(nr-i))+'end do\n' for i in r])
    #    print(codeSnippet)
    return codeSnippet

# Text formatting functions for Fortran from the m4-based code generator.
    
def DIMS(rank):
    if rank > 0:
        return '('+','.join([':' for i in range(rank) ])+')'
    else:
        return ''

def DIMS_SET(dims):
    "Return a comma separated list of dimensions, delineated by parentheses."
    retStr = ''
    if len(dims) > 0:
        retStr = '('+','.join([str(i) for i in dims])+')'
    return retStr

def DIMS_RANDOM_INTS(rank,maxDim):
    return [random.randint(1,maxDim) for i in range(rank)]

def RANDOM_INDEX(dims):
    return [random.randint(1,dims[i]) for i in range(len(dims))]

def DIMS_RANDOM(rank,maxDim):
    return DIMS_SET(DIMS_RANDOM_INTS(rank,maxDim))

def DIMS_IncrementRandomElement(dims):
    newDims = copy.copy(dims)
    i = random.randint(0,len(dims)-1)
    newDims[i] = newDims[i] + 1
    return DIMS_SET(newDims)

def FULLTYPE(fType):
    fTypes = { 'int' : 'integer',
              'char' : 'character' }
    if fType in fTypes:
        ret = fTypes[fType]
    else:
        ret = fType
    return ret

typeTower = {
    'integer' : 0,
    'real'    : 1,
    'complex' : 2 }

def maxType(type1,type2) :
    retType = type1
    if typeTower[type1] < typeTower[type2] :
        retType = type2
    return retType

def maxPrecision(prec1,prec2) :
    retPrec = prec1
    if 'default' in [prec1,prec2] :
        if prec2 != 'default' :
            retPrec = prec2
    elif prec1 < prec2 :
        retPrec = prec2
    return retPrec

def KINDATTRIBUTE0(fType,precision):
    ret = ''
    if fType.lower() == 'real' :
        ret = 'kind=r'+str(precision)+''
    elif fType.lower() == 'complex' :
        ret = 'kind=r'+str(precision)+''
    elif fType.lower() == 'integer' :
        ret = 'kind=i'+str(precision)+''
        #ret = ''
    return ret

def KINDATTRIBUTE(fType,precision):
    ret = ''
    if fType.lower() == 'real' :
        ret = '(kind=r'+str(precision)+')'
    elif fType.lower() == 'complex' :
        ret = '(kind=r'+str(precision)+')'
    elif fType.lower() == 'integer' :
        ret = '(kind=i'+str(precision)+')'
        #ret = ''
    return ret

def testKINDATTRIBUTE():
    print('COMPLEX,32 -> ' + KINDATTRIBUTE('COMPLEX',32))
    print('REAL,32 -> ' + KINDATTRIBUTE('REAL',32))
    print('real,32 -> ' + KINDATTRIBUTE('real',32))
    print('integer,32 -> ' + KINDATTRIBUTE('integer',32))

def DECLARE(variableName,fType,precision,rank,opts=', intent(in)'):
    return FULLTYPE(fType)+KINDATTRIBUTE(fType,precision)+opts+' :: '+variableName+DIMS(rank)

def DECLARESCALAR(variableName,fType,precision,rank):
    return FULLTYPE(fType)+KINDATTRIBUTE(fType,precision)+' :: '+variableName

def testDECLARE():
    print('xVar,real,32,3 -> ' + DECLARE('xVar','real',32,3))
    print('iVar,int,32,3 -> ' + DECLARE('iVar','int',32,3))
    print('iVar,int,32,0 -> ' + DECLARE('iVar','int',32,0))
    print('iVar,int,default,1 -> ' + DECLARE('iVar','int','default',1))

def OVERLOAD(routineName,fType,precision,rank):
    routineNameModifier=str(fType)+'_'+str(precision)+'_'+str(rank)
    return routineName+'_'+routineNameModifier.lower()+'D'

def testOVERLOAD():
    print('testRoutine,real,32,2 -> '+OVERLOAD('testRoutine','real',32,2))
    print('testRoutine,integer,64,0 -> '+OVERLOAD('testRoutine','integer',64,0))
    print('testRoutine,int,32,4 -> '+OVERLOAD('testRoutine','int',32,4))

def DECLAREPOINTER(pointerName,fType,precision,rank):
    return FULLTYPE(fType)+KINDATTRIBUTE(fType,precision)+', pointer :: '+ \
        OVERLOAD(pointerName,fType,precision,rank)+DIMS(rank)+' = null()'

def testDECLAREPOINTER():
    print('d-pointer: p,real,32,0 -> '+DECLAREPOINTER('p','real',32,0))
    print('d-pointer: p,integer, 64, 1 -> '+DECLAREPOINTER('p','integer',64,1))
    print('d-pointer: p,integer, 64, 3 -> '+DECLAREPOINTER('p','integer',64,3))

def NAME(fType,kind,rank):
    if fType == 'real' :
        fTypeToken = 'r'
    elif fType == 'complex' :
        fTypeToken = 'c'
    else:
        fTypeToken = 'int'
    if kind == 'default':
        kindToken = ''
    else:
        kindToken = str(kind)
    return fTypeToken+kindToken+'_'+str(rank)+'D'

def testNAME():
    print('real,32,2 -> '+NAME('real',32,2))
    print('integer,64,0 -> '+NAME('integer',64,0))

def EXPANDSHAPE(rank, variableName):
    if rank == 0:
        return ''
    elif rank == 1:
        return '(size('+variableName+'))'
    else:
        return '('+','.join(['size('+variableName+','+str(i)+')' for i in range(1,rank+1)])+')'

def testEXPANDSHAPE():
    print('0,test -> ' + str(EXPANDSHAPE(0,'test')))
    print('1,test -> ' + str(EXPANDSHAPE(1,'test')))
    print('3,test -> ' + str(EXPANDSHAPE(3,'test')))
    print('5,test -> ' + str(EXPANDSHAPE(5,'test')))

class ArrayDescription:
    def __init__(self,fType,kind,rank):
        self.fType = fType
        self.kind = kind
        self.rank = rank
        #        if rank == 0:
        #            print('ArrayDescription:Warning: rank == 0!!!')
    def NAME(self):
        return NAME(self.fType, self.kind, self.rank)
    def DECLARE(self,variableName):
        return DECLARE(variableName,self.fType, self.kind, self.rank)
    def DECLARESCALAR(self,variableName):
        return DECLARESCALAR(variableName,self.fType, self.kind, 0)
    def KIND(self):
        return self.kind
    def RANK(self):
        return self.rank
    def FTYPE(self):
        return self.fType
    def EXPANDSHAPE(self,variableName):
        return EXPANDSHAPE(self.rank,variableName)
    def FailureMessageFork(self,messageForRank1,messageOtherwise):
        if self.rank == 1 :
            return messageForRank1
        else:
            return messageOtherwise

def AddBlockSymbols(predicate, blockSymbols, inStr):
    # Note a big change in how this works -- this will actually remove str...
    retStr = '' 
    if predicate :
        if blockSymbols :
            retStr = blockSymbols[0] + inStr + blockSymbols[1]
    return retStr
        
def MakeNamesWithRank(variableName, rank):
    retStr = ''
    if rank != 0:
        retStr = ','.join([variableName+str(i+1) for i in range(rank)])
    else:
        retStr = ''.join([variableName+str(rank+1)])
    return retStr

# def compareELEMENTS(varName1,varName2,itername,shape,rank):
#    if rank == 0:
#        retStr = ''

def reportKind(t,p):
    k = ''
    if t == 'real' :
        k = '(kind=r'+str(p)+')'
    elif t == 'complex' :
        k = '(kind=c'+str(p)+')'
    elif t == 'integer' :
        k = '(kind=i'+str(p)+')'
        # default integer
        # k = ''
    else :
        k = '<Generate...py-reportKind-error:  unsupported type>'
    return k

def CONJG(x,fType='complex'):
    # If x complex return conjg, else return as is.
    retStr = str(x)
    if fType == 'complex':
        retStr = 'conjg('+str(x)+')'
    return retStr
    
