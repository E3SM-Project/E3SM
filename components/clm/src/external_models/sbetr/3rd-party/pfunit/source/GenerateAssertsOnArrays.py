#!/usr/bin/env python
# For python 2.6-2.7
from __future__ import print_function
# For python2.5
from __future__ import with_statement
# 
# Generate Assert<Type>.F90, which provides assertEqual and others for arrays.
# 
# Usage:  ./GenerateAssertsOnArrays.py
#
# Outputs:
#    A large number of library files implementing assertRELATION routines.
#    generated.inc
#    AssertArrays.fh
#
# M. Rilee
#
#    2014-1215:  Minor revisions. Updated integer arrays to handle same number of ranks as others.
#
#    2014-0418:  Moved effective assert routines to own file.  Segregated other routines into files by rank.
#
#    2014-0324:  Fully implemented relational operators beyond "=".
#
#    2013-0814:  Added default r64 to call from assertEqual_w/o_tol to internal proc.
#                Added logical to makeExpectedFTypes - but not for prime time.
#
#    Initial: 2013-0304
# 

##  Abbreviations
#   WOTol =>  WithoutTolerance
#   def   =>  default

##### system code #####

# python2 - Deprecated in python 2.7+
import imp
try:
    imp.find_module('argparse')
    found = True
except ImportError:
    found = False

# Preferred for python 2.7+, python 3
# import importlib
# argparse_loader = importlib.find_loader('argparse')
# found = argparse_loader is not None

if found:    
    import argparse
else:
    print('GenerateAssertOnArrays.py::Error. pFUnit requires argparse module provided by python version >= 2.7.')
    print('Quitting!'); quit()

##### utility code #####

from Utilities import *
from CodeUtilities import *
import textwrap
import random
import copy

##### preliminaries #####

parser = argparse.ArgumentParser( \
                                  description='Generate assertions with relational operators on arrays.', \
                                  usage='%(prog) --maxRank MaximumRankOfArrays' \
                                  )
parser.add_argument('--maxRank', help='The maximum rank of the arrays for which to generate code.')
args = parser.parse_args()

                                  
##### begin generation code #####


### Restrictions on types and type combinations.

def dr_TolAllowedPrecisions(t,pFound='64') :
    "returns a list of strings corresponding to the precisions 'tolerance' may take on. \
    dr_ refers to 'Difference Report.'  Please see the DifferenceReport routines."
    allowed = []
    if t == 'logical' :
        allowed = []
    elif t == 'integer' :
        allowed = ['32','64']
    elif t == 'real' or 'complex' :
        if pFound == '32' :
            allowed = ['32']
        elif pFound == '64' :
            allowed = ['32','64']
        else :
            raise ValueError("dr_TolAllowedPrecisions: Bad value of pFound.")
    return allowed

def dr_TolAllowedPrecisions_orig(t,pFound='64') :
    "returns a list of strings corresponding to the precisions 'tolerance' may take on. \
    dr_ refers to 'Difference Report.'  Please see the DifferenceReport routines."
    allowed = []
    if t == 'logical' :
        allowed = []
    elif t == 'integer' :
        allowed = ['64']
    else:
        allowed = ['32','64']
    return allowed
# Sorry the following is confusing.
# It's partially a risky mess because default might be overloaded with 32 and 64.
def allowedPrecisions(t,tFound='',pFound='64') :
    "returns a list of strings corresponding to the precisions 'expected' may take on."
    allowed = []
    if tFound == 'integer' :
        # The new case
        if t == 'logical' :
            allowed = []
        else:
            # If found is integer allow t=tExpected to be either 32 or 64
            # without regard to pFound.
            allowed = ['32','64']
    else :
        # The old case, with 'def' below, which ladders the precisions appropriately.
        if t == 'logical' :
            allowed = []
        elif t == 'integer' :
            allowed =['32','64']
            # allowed = ['def']
        elif t == 'real' or 'complex' :
            if pFound == '32' :
                allowed = ['32']
            elif pFound == '64' : 
                allowed = ['32','64']
    return allowed

# 2014-1208-1646-25-UTC MLR 
def allowedExpected(tFound) :
    # allowed = []
    if tFound in 'integer' :
        allowed = ['integer','real']
    # mlr: old and good    allowed = ['integer']
    #    allowed = []
    #elif tFound in 'integer' :
    #    allowed = ['integer']
    elif tFound == 'real' :
        allowed = ['integer','real']
    elif tFound == 'complex' :
        allowed = ['integer','real','complex']
    else :
        allowed = []
    return allowed

#### Type coercions

def coerceReal(x,kind='r32') :
    if kind == 'ckDefault' :
        kind = 'r32'
    return 'real('+x+',kind='+kind+')'

def coerceComplex(x,kind='c32') :
    if kind == 'ckDefault' :
        kind = 'c32'
    return 'cmplx('+x+',kind='+kind+')'

def coerceKind(x,kind='ckDefault',t='real'):
    coerceStr = x
    if t == 'real' :
        coerceStr = coerceReal(x,kind=kind)
    elif t == 'complex' :
        coerceStr = coerceComplex(x,kind=kind)
    elif t == 'integer':
        coerceStr = coerceReal(x,kind=kind)
    else:
        coerceStr = 'coerceKind: ERROR - t = '+t+', kind = '+kind+', x = '+x
    return coerceStr

####

def tolDECLARE(tolerance,descr,opts=', optional, intent(in)',name='tolerance'):
    retStr = ''
    if tolerance == 0 :
        retStr = DECLARE(name,descr.FTYPE(),descr.KIND(),0,opts=opts)
    else:
        retStr = """real(kind=r"""+str(tolerance)+""")"""+opts+""" :: """+name
    return retStr

def makeSubroutineName(assertionName,expectedName,foundName,tolerance):
    return \
       """assert"""+assertionName+"""_""" + \
       expectedName + """_""" + foundName + """_tol""" + tolerance

       #
       # What does it mean to compare a 0D with a 1D array? MLR ***
       #

comparisonCase = {
    "NotEqual":"NEQP",
    "Equal":"EQP",
    "GreaterThan":"GTP",
    "GreaterThanOrEqual":"GEP",
    "LessThan":"LTP",
    "LessThanOrEqual":"LEP",
    "RelativelyEqual":"RELEQP"
    }

def generateASSERT(assertionName,expectedDescr, foundDescr, tolerance):
    subroutineName = makeSubroutineName(assertionName, \
                                        expectedDescr.NAME(), \
                                        foundDescr.NAME(), \
                                        str(tolerance))
    
    # Maybe set up an object where comments have some extra meaning.
    commentPreambleString = \
"""
  !---------------------------------------------------------------------------
  !> Asserts that two real numbers are equal.  If they are not, an
  !! Exception is thrown with the given message.
  !!
  !! @param expected - expected real numbers
  !! @param found -  found real numbers
  !! @param message - the identifying message for the Exception
  !!
  !! @throw Exception - when two real numbers are not equal.
  !---------------------------------------------------------------------------
"""

    declareExpected = \
"     " + expectedDescr.DECLARE('expected') + "\n"
    declareFound = \
"     " + foundDescr.DECLARE('found') + "\n"
    declareTolerance = \
"     " + tolDECLARE(tolerance,foundDescr,opts=', intent(in)')
    declareTolerance_orig = \
"     " + tolDECLARE(tolerance,foundDescr,opts=', optional, intent(in)')

# foundDescr.FTYPE() == integer breaks the tolerance cast near line 337.
#    toleranceKind = KINDATTRIBUTE0(foundDescr.FTYPE(),foundDescr.KIND())
# If found is an integer, then recast the kind parameter i32 -> r64.
# Note:  tolerances are set at r64 by default.

# Okay -- we need to recall the logic in tolDECLARE, which accounts for when
# tolerance=0, in which case we default to found->kind().
    if tolerance != 0:
        toleranceKind = KINDATTRIBUTE0('real',tolerance)
    else:
        fft0 = foundDescr.FTYPE()
        if fft0.lower() == 'integer':
            toleranceKind = KINDATTRIBUTE0('real',foundDescr.KIND())
        else:
            toleranceKind = KINDATTRIBUTE0(foundDescr.FTYPE(),foundDescr.KIND())

    declareExpectedScalar_expected = \
"      " + expectedDescr.DECLARESCALAR('expected') + "\n"
    declareFoundScalar_found = \
"      " + foundDescr.DECLARESCALAR('found') + "\n"

    declareExpectedScalar_expected0 = \
"      " + expectedDescr.DECLARESCALAR('expected0') + "\n"
    declareFoundScalar_found0 = \
"      " + foundDescr.DECLARESCALAR('found0') + "\n"

    # eType = expectedDescr.FTYPE(); fType=foundDescr.FTYPE()
    # ePrec = expectedDescr.KIND();  fPrec=foundDescr.KIND()
    internalSubroutineName = makeAssertInternalName(assertionName,
                                                    expectedDescr,foundDescr,tolerance)

    # 2014-0415 Need to go to an older style of reference since we're fooling
    # the TKR system of the compilers.
    # foundArray = 'found'+DIMS_SET([1]*foundDescr.RANK())
    foundArray = 'found'

# Subroutine name
    callInternalSubroutineName = internalSubroutineName
# Interface name
#    callInternalSubroutineName = 'Assert'+assertionName+'_internal'
    

    if 'complex' in [foundDescr.FTYPE(),expectedDescr.FTYPE()] :
        declareDelta = \
"      " + "complex(kind=kind(found)) :: delta"+foundDescr.EXPANDSHAPE('found') +"\n"+ \
"      " + "complex(kind=kind(found)) :: delta1" +"\n"
    else:
        declareDelta = \
"      " + "real(kind=kind(found)) :: delta"+foundDescr.EXPANDSHAPE('found') +"\n"+ \
"      " + "real(kind=kind(found)) :: delta1" +"\n"

# Need to handle scalar case...
    foundFirstElt = DIMS_SET([1]*foundDescr.RANK())
    expectedFirstElt = DIMS_SET([1]*expectedDescr.RANK())
    
    retString = \
        commentPreambleString + \
"""
   subroutine """+subroutineName+"""( &
   &  expected, found, tolerance, message, location )
! was tolerance, message -- need to propagate changes... e.g. to test files
     implicit none\n""" + \
     declareExpected + \
     declareFound + \
     declareTolerance + """
     character(len=*), optional, intent(in) :: message
     type (SourceLocation), optional, intent(in) :: location

     real(kind=kind(tolerance)) :: tolerance_
     character(len=:), allocatable :: message_
     type (SourceLocation) :: location_

! Tolerance is now not optional.
!     if(present(tolerance)) then
        tolerance_ = tolerance
!     else
!        tolerance_ = real(0."""+ifElseString(toleranceKind,', '+toleranceKind,'')+""")
!     end if

     if(present(location)) then 
        location_ = location
     else
        location_ = UNKNOWN_SOURCE_LOCATION
     end if

     if(present(message)) then
        message_ = message
     else
        message_ = NULL_MESSAGE
     end if

    call assertSameShape(shape(expected),shape(found), message=message_, location=location_)
    if (anyExceptions()) return

! Next allow call to here...
!mlr-NextStep-Begin
! 2014-0414 Call interface here instead of internal subroutine name, if possible. mlr
! """ + internalSubroutineName + """
     call """+callInternalSubroutineName+"""(&
     &  expected, shape(expected), """+foundArray+""", shape(found), &
     &  tolerance_, message_, location_, """+ comparisonCase[assertionName] +""" )
!mlr-NextStep-End
     
   end subroutine
"""

    declareTolerance_ = \
"     " + tolDECLARE(tolerance,foundDescr,opts='',name='tolerance_')

    retString += \
        commentPreambleString + \
"""
   subroutine """+subroutineName+'_WOTol'+"""( &
   &  expected, found, message, location )
     implicit none\n""" + \
     declareExpected + \
     declareFound + """
     character(len=*), optional, intent(in) :: message 
     type (SourceLocation), optional, intent(in) :: location

     call """+subroutineName+"""(&
   &   expected, found, &
   &   tolerance=real(0."""+ifElseString(toleranceKind,', '+toleranceKind,', '+KINDATTRIBUTE0('real',64))+"""), &
   &   message=message, location=location )
     
   end subroutine
"""

    return retString

def makeAssertInternalName(assertionName,eDescr,fDescr,tolerance):
    eType=eDescr.FTYPE(); fType=fDescr.FTYPE()
    ePrec=eDescr.KIND();  fPrec=fDescr.KIND()
    eRank=min(eDescr.RANK(),1)
    fRank=min(fDescr.RANK(),1)
    subroutineName = \
        'assert' + assertionName + \
        '_e'+ str(eRank) +'_'+eType+str(ePrec)+ \
        '_f'+ str(fRank) +'_'+fType+str(fPrec)+ \
        '_tol'+str(tolerance)+'_'
    return subroutineName

def makeAssertInternal_type(assertionName,eDescr,fDescr,tolerance):
    eType=eDescr.FTYPE(); fType=fDescr.FTYPE()
    ePrec=eDescr.KIND();  fPrec=fDescr.KIND()
    eRank=min(eDescr.RANK(),1)
    fRank=min(fDescr.RANK(),1)

    subroutineName = makeAssertInternalName(assertionName,eDescr,fDescr,tolerance)

    if(eRank != 0):
        expected_i = 'expected(i)'
        eOpts = ", dimension(product(eShape))"
    else:
        expected_i = 'expected'
        eOpts = ""
        
    if(fRank != 0):
        found_i = 'found(i)'
        fOpts = ", dimension(product(fShape))"
    else:
        found_i = 'found'
        fOpts = ""

    if(eRank != 0 or fRank != 0):
        all_i = "all"
    else:
        all_i = ""

    expectedDeclaration = \
"    " + DECLARE('expected',eType,ePrec,0,\
                                  opts=eOpts+', intent(in)') + "\n" + \
"    " + DECLARE('expected_',eType,ePrec,0,\
                                  opts='') + "\n" 

    foundDeclaration = \
"    " + DECLARE('found',fType,fPrec,0,\
                               opts=fOpts+', intent(in)') + "\n" + \
"    " + DECLARE('found_',fType,fPrec,0,\
                               opts= '' ) + "\n" + \
"    " + DECLARE('delta1',fType,fPrec,0,\
                               opts= '' ) + "\n"
    
    toleranceDeclaration = \
ifElseString(tolerance == 0,\
"""
    real(kind=kind(found)) :: tolerance
""", \
"""
    real(kind=r"""+str(tolerance)+"""), intent(in) :: tolerance
""" ) 

# Start to define routine...
    retStr = """
    subroutine """+subroutineName+"""( &
    & expected,eShape,found,fShape,tolerance,message,location, &
    & comparison )""" + \
"""
    use Params_mod
    use Exception_mod
    use StringConversionUtilities_mod
    use ThrowFundamentalTypes_mod, only : locationFormat
    implicit none
    integer, intent(in), dimension(:) :: eShape, fShape
    character(len=*), intent(in) :: message
    type (SourceLocation), intent(in) :: location
    integer, intent(in) :: comparison
""" + \
    expectedDeclaration + \
    foundDeclaration + \
    toleranceDeclaration + """
    
! mlr 2013-0908 Note:  Perhaps have tolerance_ with a type depending on found... incl. logical or int.
    real(kind=kind(tolerance)) :: tolerance_
!---    real(kind=kind(expected)) :: expected_
!---    real(kind=kind(found)) :: found_
    integer :: i,m,ir
    logical OK
    integer, dimension(size(fShape)) :: iLocation
    character(len=MAXLEN_SHAPE) :: locationInArray
    real :: denominator

    ! Return immediately if the two are precisely equal.
    ! This is necessary to deal with identical infinities, which cannot be
    ! subtracted.
    if (comparison .eq. EQP) then
       if (""" + all_i + """(expected == found)) return
    end if
    ! Note:  The above begs the question about how to handle Inf for non-.eq. cases...

! MLR: The following just might work...
    tolerance_ = tolerance

! fType != 'complex' = """ + str(fType != 'complex') + """

! Note:  Could assert size(expected) = size(found) and fShape = eShape...

!    print *,'0800 ',product(fShape),fShape
    m = product(fShape)
    i = 0
!    
    OK = .true.
    
! Note:  Comparison occurs here.  Could use isWithinTolerance or other comparison function.
! mlr 2013-0908 Other comparisons:  tolerance-less integer comparison... logical...
    if( m > 0 )then
       do while ( i < m .and. OK )
         i = i + 1

         delta1 = """ + expected_i + """-""" + found_i + """

         select case (comparison)
            case (EQP)
              OK = &
              &  isWithinTolerance( &
              &    delta1, &
              &    real(tolerance_,kind=r64), &
              &    L_INFINITY_NORM )
            case (NEQP)
              OK = &
              &  .not. &
              &  isWithinTolerance( &
              &    delta1, &
              &    real(tolerance_,kind=r64), &
              &    L_INFINITY_NORM ) """ + \
ifElseString(fType != 'complex', \
"""
            case (GTP)
              OK = delta1 .gt. 0
            case (GEP)
              OK = delta1 .ge. 0
            case (LTP)
              OK = delta1 .lt. 0
            case (LEP)
              OK = delta1 .le. 0 """,'') + \
"""
            case (RELEQP)
              if ( abs("""+expected_i+""") > 0 ) then
                 OK = &
              &  isWithinTolerance( &
              &    delta1 / """+expected_i+""", &
              &    real(tolerance_,kind=r64), &
              &    L_INFINITY_NORM )
              else
                 OK = &
              &  isWithinTolerance( &
              &    delta1, &
              &    real(tolerance_,kind=r64), &
              &    L_INFINITY_NORM )
              end if
            case default
              ! This case should not occur for this type-kind-rank.
              print *,'internal: """+subroutineName +""" select-error-1'
              OK = .false.
         end select
         
!         OK = .not. ( expected(i) /= found(i) )
!         OK = .not. ( """+expected_i+""" /= """+found_i+""" )
       end do
    else
!         i = 1
         delta1 = """ + expected_i + """-""" + found_i + """    

         select case (comparison)
            case (EQP)
              OK = &
              &  isWithinTolerance( &
              &    delta1, &
              &    real(tolerance_,kind=r64), &
              &    L_INFINITY_NORM )
            case (NEQP)
              OK = &
              &  .not. &
              &  isWithinTolerance( &
              &    delta1, &
              &    real(tolerance_,kind=r64), &
              &    L_INFINITY_NORM ) """ + \
ifElseString(fType != 'complex', \
"""
            case (GTP)
              OK = delta1 .gt. 0
            case (GEP)
              OK = delta1 .ge. 0
            case (LTP)
              OK = delta1 .lt. 0
            case (LEP)
              OK = delta1 .le. 0 """,'') + \
"""
            case (RELEQP)
              if ( abs("""+expected_i+""") > 0 ) then
                 OK = &
              &  isWithinTolerance( &
              &    delta1 / """+expected_i+""", &
              &    real(tolerance_,kind=r64), &
              &    L_INFINITY_NORM )
              else
                 OK = &
              &  isWithinTolerance( &
              &    delta1, &
              &    real(tolerance_,kind=r64), &
              &    L_INFINITY_NORM )
              end if
            case default
              ! This case should not occur for this type-kind-rank.
              print *,'internal: """+subroutineName +""" select-error-2'
              OK = .false.
         end select

!         OK = &
!         &  isWithinTolerance( &
!         &    delta1, &
!         &    real(tolerance_,kind=r64), &
!         &    L_INFINITY_NORM )
         
!         OK = .not. ( """+expected_i+""" /= """+found_i+""" )
    end if

    if( .not. OK )then

    ! Save the FirstBad...
    expected_ = """+expected_i+"""
    found_    = """+found_i+"""

!    if( m > 0 )then
    if( size(fshape) > 0 ) then

    i = i - 1
    do ir = 1,size(fShape)
      iLocation(ir) = mod(i,fShape(ir)) + 1
      i = i / fShape(ir)
    end do

!    print *,'0998 ',m
!    print *,'0999 ',size(fShape)
!    print *,'1000 ',iLocation
    write(locationInArray,locationFormat(iLocation)) iLocation

    else

    write(locationInArray,*) '[1]'

    end if

! Scalar
! Note use of abs

    select case (comparison)
    case (EQP)
       call throw( &
       & appendWithSpace(message, &
       & trim(valuesReport(expected_,found_)) // &
       & '; '//trim(differenceReport(abs(found_ - expected_), tolerance_)) // &
       & unlessScalar(fShape,';  first difference at element '//trim(locationInArray))//'.'), &
       & location = location &
       )
    case (NEQP)
       call throw( &
       & appendWithSpace(message, &
       & 'NOT: '//trim(valuesReport(expected_,found_)) // &
       & '; '//trim(differenceReport(abs(found_ - expected_), tolerance_)) // &
       & unlessScalar(fShape,';  first equality at element '//trim(locationInArray))//'.'), &
       & location = location &
       ) """ + \
ifElseString(fType != 'complex', \
"""
    case (GTP)
       call throw( &
       & appendWithSpace(message, &
       & trim(valuesReport(expected_,found_, &
       &   ePrefix='expected', &
       &   fPrefix='to be greater than:')) // &       
       & unlessScalar(fShape,';  first difference at element '//trim(locationInArray))//'.'), &
       & location = location &
       )
    case (GEP)
       call throw( &
       & appendWithSpace(message, &
       & trim(valuesReport(expected_,found_, &
       &   ePrefix='expected', &
       &   fPrefix='to be greater than or equal to:')) // &       
       & unlessScalar(fShape,';  first difference at element '//trim(locationInArray))//'.'), &
       & location = location &
       )
    case (LTP)
       call throw( &
       & appendWithSpace(message, &
       & trim(valuesReport(expected_,found_, &
       &   ePrefix='expected', &
       &   fPrefix='to be less than:')) // &       
       & unlessScalar(fShape,';  first difference at element '//trim(locationInArray))//'.'), &
       & location = location &
       )
    case (LEP)
       call throw( &
       & appendWithSpace(message, &
       & trim(valuesReport(expected_,found_, &
       &   ePrefix='expected', &
       &   fPrefix='to be less than or equal to:')) // &       
       & unlessScalar(fShape,';  first difference at element '//trim(locationInArray))//'.'), &
       & location = location &
       ) """,'') + \
"""
    case (RELEQP)
       if (expected_ .eq. 0) then
          denominator = expected_
       else
          denominator = 1.0
       end if 
       call throw( &
       & appendWithSpace(message, &
       & trim(valuesReport(expected_,found_, &
       &   ePrefix='RELEQ: expected', &
       &   fPrefix='to be near:')) // &
       & '; '//trim(differenceReport( &
       &   abs(found_ - expected_)/denominator, tolerance_)) // &
       & unlessScalar(fShape,';  first difference at element '//trim(locationInArray))//'.'), &
       & location = location &
       )
    case default
       ! This case should not occur for this type-kind-rank.
       print *,'internal: """+subroutineName +""" select-error-3'
       call throw( &
       & appendWithSpace(message, &
       & 'pFUnit internal error:  unexpected comparison given type-kind-rank'), &
       & location = location &
       & )
    end select

    end if

    end subroutine """+subroutineName+"""

       

"""

    runit = routineUnit(subroutineName,retStr)
    
    
    return runit
    

class constraintASSERT(routineUnit):
    "Defines the comparison code as a routineUnit so that it can be \
    used by the module generation code.  These declarations are used \
    to construct interface blocks as well as the routines themeselves."
    def __init__(self, assertionName, expectedDescr, foundDescr, tolerance):
        self.expectedDescr = expectedDescr
        self.foundDescr = foundDescr
        self.name = makeSubroutineName( assertionName, \
                                        expectedDescr.NAME(), \
                                        foundDescr.NAME(), \
                                        str(tolerance) )
        #! bad... too simple...
        ## Add in the extra module procedures... If needed...
        self.setDeclaration(declaration(self.name,self.name))
        ## Kluge.  Need to make makeSubroutineNames and load the extra interface entries there.
        self.name1 = self.name+'_WOTol'
        self.addDeclaration(declaration(self.name1,self.name1))

        ## If you need another kind of code generator, perhaps
        ## conditioned on eDesc., fDesc., or tol, then that logic
        ## would go here... E.g. to implement assertEqual(Logical(...))
        ##
        ## Dependency injection.  Will generate "assert"+assertionName
        ## assertionName="Equal"
        ## This next line actually generates the text of the code.
        self.setImplementation(generateASSERT(assertionName, \
                                              expectedDescr, \
                                              foundDescr, \
                                              tolerance ))
        self.tolerance = tolerance
        return

def constructAssertInterfaceBlock(assertionName,foundFTypes=['real'],foundRanks=[],patchIntXX=False):
    AssertInterfaceBlockName='assert'+assertionName
    AssertInterfaceBlock = interfaceBlock(AssertInterfaceBlockName)
    # Construct asserts generates the combinations based on what is passed in here.
    [AssertInterfaceBlock.addRoutineUnit(r) for r in constructASSERTS(assertionName,foundFTypes=foundFTypes,foundRanks=foundRanks,patchIntXX=patchIntXX)]
    return AssertInterfaceBlock

def VECTOR_NORM_NAME(rank,fType='real',precision=64):
    return """vectorNorm_"""+str(rank)+"D"+"_"+fType+str(precision)

def vnSQRT(x,fType,precision):
    retStr = ''
    if fType == 'integer' :
        retStr = 'sqrt(real('+x+',kind=r64))'
    else :
        retStr = 'sqrt('+x+')'
    return retStr

def generateVECTOR_NORM(rank,fType='real',precision=64):
    subroutineName = VECTOR_NORM_NAME(rank,fType=fType,precision=precision)
    dimStr = DIMS(rank)
    retstr = \
"""
  !---------------------------------------------------------------------------
  !> Returns the independent of norm in vector by the given diminsional
  !! double-precission real numbers and given integer norm
  !!
  !! The following is for rank = """+str(rank)+""".
  !!
  !! @param x - given dimensional double-precision real numbers
  !! @param norm - given norm
  !!
  !! @return independent of norm
  !---------------------------------------------------------------------------
  function """+subroutineName+"""(x, norm) result(y)
    """+DECLARE('x',fType,precision,rank,opts=', intent(in)')+"""
    integer :: norm
! mlr 2013-0908 Maybe we change the range of VECTOR_NORM to include integer & logical.
    real (kind=r64) :: y
""" + ifElseString(rank == 0, \
"""
    y = abs(x) ! independent of norm for rank=0 (scalar) case.
""", \
"""
! Note that abs(complex) is like the L2_NORM unless care is taken *here*.  Fix later...
    select case (norm)  ! code to support rank /= 0 cases.
    case (L_INFINITY_NORM)
       y = maxval(abs(x))
    case (L1_NORM)
       y = sum(abs(x))
    case (L2_NORM)
!       y = sqrt(sum(x**2))
!       y = sqrt(sum(x*conjg(x)))
!       y = sqrt(sum(x*"""+CONJG('x',fType)+"""))
       y = """ + vnSQRT("""sum(x*"""+CONJG('x',fType)+""")""",fType,precision) + """
    end select
""") + \
"""
  end function """ + subroutineName + """
"""
    return retstr

class VECTOR_NORM(routineUnit):
    def __init__(self,rank,fType='real',precision=32):
        self.rank = rank
        self.fType = fType
        self.precision = precision
        self.name = VECTOR_NORM_NAME(rank,fType=self.fType,precision=self.precision)
        self.declaration = declaration(self.name,self.name)
        self.declarations = [self.declaration]
        self.implementation \
            = implementation( \
                self.name, \
                generateVECTOR_NORM(rank,fType=self.fType,precision=self.precision))
        return

def constructVectorNormInterfaceBlock(foundRanks=range(6)):
    VectorNormInterface = interfaceBlock('vectorNorm')
    list(map(VectorNormInterface.addRoutineUnit,
             flattened( [[VECTOR_NORM(i,fType=t,precision=p) for i in foundRanks \
                          for p in allowedPrecisions(t) ] \
                         for t in ['real','complex','integer']
                     ])))
    return VectorNormInterface

def constructDifferenceReportInterfaceBlock():
    DifferenceReportInterface = interfaceBlock('differenceReport')
    list(map(DifferenceReportInterface.addRoutineUnit,
             flattened( [[[makeDifferenceReport_type(t=t,p=p,tol=tol) \
                           for tol in dr_TolAllowedPrecisions(t) ]
                          for p in allowedPrecisions(t) ]
                         for t in ['integer','real','complex'] \
                     ])))
    return DifferenceReportInterface

def constructValuesReportInterfaceBlock():
    ValuesReportInterface = interfaceBlock('valuesReport')
    list(map(ValuesReportInterface.addRoutineUnit,
             flattened( [[[[makeValuesReport_type(te=te,tf=tf,pe=pe,pf=pf) \
                            for pe in allowedPrecisions(te,tFound=tf,pFound=pf) ] \
                            #sigh for pe in allowedPrecisions(te,pFound=pf) ] \
                           for pf in allowedPrecisions(tf) ] \
                          for te in allowedExpected(tf) ] \
                         for tf in ['integer','real','complex'] \
                     ])))
    return ValuesReportInterface

# Scalar args?
# Need assertionName...
# def constructAssertInternalInterfaceBlock(assertionName="Equal"):
def constructAssertInternalInterfaceBlock(assertionName,expose=False):
    # foundRanks... mlr... how would foundRanks work here?
    AssertInternalInterfaceBlockName='assert'+assertionName+'_internal'
    AssertInternalInterface = interfaceBlock(AssertInternalInterfaceBlockName)
    # 2014-0415 expose may not be necessary here! MLR
    list(map((lambda x: AssertInternalInterface.addRoutineUnit(x,expose=expose)),
             [makeAssertInternal_type(assertionName,
                                      a.getExpectedDescription(),
                                           a.getFoundDescription(),
                                           a.getTolerance()
                                       )
              for a in
              flattened(
                  [[[[[[[[AssertRealArrayArgument(assertionName,te,pe,re,tf,pf,rf,tol)
                          for re in [0,1] ]
                         for rf in [0,1] ]
                        for tol in makeTolerances(pe,pf) ]
                       for pe in allowedPrecisions(te,tFound=tf,pFound=pf) ]
                       #sigh for pe in allowedPrecisions(te,pFound=pf) ]
                      for pf in allowedPrecisions(tf) ]
                     for te in allowedExpected(tf) ]
                    for tf in ['integer','real','complex']
                ]])]))
    return AssertInternalInterface

def isWithinToleranceName(rank,fType='real',precision=64):
    return """isWithinTolerance_"""+str(rank)+"""D"""+"_"+fType+str(precision)

def generateIsWithinTolerance(rank,fType='real',precision=64):
    "Generate the code for the comparison function. Calls \
    vectorNorm..."
    subroutineName = isWithinToleranceName(rank,fType=fType,precision=precision)
    dimStr = DIMS(rank)

    declareKind = ''
    if fType == 'integer' :
        declareKind = '(kind=i'+str(precision)+')'
        #declareKind = ''
    elif fType == 'real' :
        declareKind = '(kind=r'+str(precision)+')'
    elif fType == 'complex' :
        declareKind = '(kind=c'+str(precision)+')'
    else:
        print('isWithinToleranceTypeError')
        
    retstr = \
"""
   logical function """+subroutineName+"""(x, tolerance, norm)
!     """+fType+""" (kind=r"""+str(precision)+"""), intent(in) :: x"""+dimStr+"""
     """+fType+declareKind+""", intent(in) :: x"""+dimStr+"""
     real (kind=r64), intent(in) :: tolerance
     integer,         intent(in) :: norm

     """+subroutineName+""" = ( vectorNorm(x, norm) <= tolerance )

   end function """+subroutineName+"""
"""
    return retstr

class IsWithinTolerance(routineUnit):
    "A routineUnit specialized to the isWithinTolerance comparison function."
    def __init__(self,rank,fType='real',precision=64):
        self.rank = rank
        self.precision = precision
        self.name = isWithinToleranceName(rank,fType=fType,precision=precision)
        self.fType = fType
        self.declaration = declaration(self.name, self.name)
        self.declarations = [self.declaration]
        self.implementation \
            = implementation(self.name, \
                                           generateIsWithinTolerance(self.rank, \
                                                                     fType=self.fType, \
                                                                     precision=self.precision))
        return

def constructIsWithinToleranceInterfaceBlock(foundRanks=range(6)):
    "For the comparison function, make an interface block and \
    implementation for inclusion into a module."
    iwt_InterfaceBlock = interfaceBlock('isWithinTolerance')
    list(map(iwt_InterfaceBlock.addRoutineUnit,
             flattened(
                 [[IsWithinTolerance(i,fType=t,precision=p)
                   for i in foundRanks
                   for p in allowedPrecisions(t) ]
                  for t in ['real','complex','integer']
              ])))
    return iwt_InterfaceBlock


### Currently Active ###
def makeValuesReport_type(te='real',tf='real',pe='64',pf='64'):
    expectedKind = reportKind(te,pe)
    foundKind =    reportKind(tf,pf)
    mxType = maxType(te,tf) 
    mxPrec = maxPrecision(pe,pf)
    
    if te == 'integer':
        coercedExpected = 'expected'
    else:
        coercedExpected = coerceKind('expected',t=mxType)
    if tf == 'integer':
        coercedFound    = 'found'
    else:
        coercedFound    = coerceKind('found',t=mxType)
        
    runit = routineUnit('valuesReport_'+te+tf+pe+pf, \
"""
      character(len=MAXLEN_MESSAGE) &
      & function valuesReport_"""+te+tf+pe+pf+""" &
      & (expected,found,ePrefix,ePostfix,fPrefix,fPostfix) &
      & result(valuesReport)
        """+te+expectedKind+""", intent(in) :: expected
        """+tf+foundKind+""", intent(in) :: found
        character(len=*), optional, intent(in) :: &
      &   ePrefix, ePostfix, fPrefix, fPostfix
        character(len=MAXLEN_MESSAGE) :: &
      &   ePrefix_, ePostfix_, fPrefix_, fPostfix_

      if( .not.present(ePrefix) ) then
         ePrefix_ = 'expected'
      else
         ePrefix_ = ePrefix
      end if
      if( .not.present(ePostfix) ) then
         ePostfix_ = ''
      else
         ePostfix_ = ePostfix
      end if
      if( .not.present(fPrefix) ) then
         fPrefix_ = 'but found:'
      else
         fPrefix_ = fPrefix
      end if
      if( .not.present(fPostfix) ) then
         fPostfix_ = ''
      else
         fPostfix_ = fPostfix
      end if

! Note: removed '<.>'
        valuesReport = &
      & trim(ePrefix_)//' '// trim(toString("""+coercedExpected+""")) // &
      & trim(ePostfix_)//' '// &
      & trim(fPrefix_)//' '//trim(toString("""+coercedFound+""")) // &
      & trim(fPostfix_)// &
      & ''
      
      end function
""")
    # runit.setDeclaration(declaration(runit.getName(),'public '+runit.getName()))
    return runit

### Currently Active ###
def makeDifferenceReport_type(t='real',p='64',tol='64'):
    """Report on the difference found between expected and found.  This should be
    redesigned to support difference between real & integer, e.g. treatment of tolerances."""
    # TODO:  Fix integer for non-default values, e.g. INT64...
    # TODO:  Fix the papering-over of integer/real/tolerance treatment.
    if t == 'integer' :
        # expectedDeclaration = " integer, intent(in) :: difference"
        expectedDeclaration = t+"""(kind=i"""+p+"""), intent(in) :: difference"""
        coercedDifference = 'difference'
    else :
        expectedDeclaration = t+"""(kind=r"""+p+"""), intent(in) :: difference"""
        coercedDifference = coerceKind('difference',t=t)
    #    coercedTolerance =  coerceKind('tolerance',t=t)
    
    # runit = routineUnit(
    
    differenceReportSource=\
"""
    character(len=MAXLEN_MESSAGE) &
    & function differenceReport_"""+t+p+tol+"""(difference, tolerance) result(differenceReport)
     """+expectedDeclaration+"""
     real(kind=r"""+tol+"""), intent(in) :: tolerance
!     real(kind=r"""+tol+"""), optional, intent(in) :: tolerance
     character(len=2) rel
     if (abs("""+coercedDifference+""") .gt. tolerance) then
        rel = '> '
     else
        rel = '<='
     end if"""
    if t != 'integer':
        differenceReportSource=differenceReportSource+\
"""
     differenceReport = ' difference: |' // trim(toString("""+coercedDifference+""")) // &
      & '| '// trim(rel) //' tolerance:' // trim(toString("""+'tolerance'+"""))
    end function 
"""
    else:
        differenceReportSource=differenceReportSource+\
"""
     differenceReport = ' difference: |' // trim(toString("""+coercedDifference+""")) // &
      & '| '
    end function 
"""
    runit = routineUnit('differenceReport_'+t+p+tol,differenceReportSource)

# Don't need the following because we'll add to an interface block.
#    runit.setDeclaration(declaration(runit.getName(),'public '+runit.getName()))
    return runit

def declareUSES(internalRoutines=[]):
    retStr=\
"""
   use Params_mod
   use AssertBasic_mod
   use Exception_mod
   use SourceLocation_mod
   use StringConversionUtilities_mod
   use AssertArraysSupport_mod
"""
    for i in internalRoutines:
        retStr=retStr+\
"""   use AssertArraysInternal"""+i+"""_mod
"""
    return retStr

def declareDISCIPLINE():
    return \
"""
   implicit none
   private
"""

# Flag AssertEqual

def declareEXPORTS(assertionRoutines,basename='AssertReal'):
    retPublic = ""
    for aRoutineName in assertionRoutines:
        retPublic = retPublic + \
"""   public :: """+aRoutineName+"""
"""
    if basename == 'AssertReal' :
        retPublic = retPublic + """

   public :: vectorNorm
   public :: isWithinTolerance
 
   public :: L_INFINITY_NORM
   public :: L1_NORM
   public :: L2_NORM

"""
    return retPublic

def declareEXPORTS_PARAMETERS():
    return \
"""
   integer, parameter :: L_INFINITY_NORM = 0
   integer, parameter :: L1_NORM         = 1
   integer, parameter :: L2_NORM         = 2

   integer, parameter :: MAXLEN_SHAPE = 80
"""

#### Helper functions for constructASSERTS -- a main workhorse for generating the specific routines.

#
# In the following: expectedPrecision and foundFType are strings.
# 
def makeExpectedFTypes(expectedPrecision,foundFType,foundFTypes=['real']):
    """A very application-specific mapping to construct an fType list
    for expected.  Make sure that if we're looking at complex that we
    do not replicate real-real comparisons."""
    retTypes = ['makeExpectedFType::ERROR']
    if 'logical' in foundFTypes :
        if foundFType == 'logical' :
            retTypes=['logical']        
    elif expectedPrecision == 'def':
        # 
        if not 'complex' in foundFTypes :
            retTypes=['integer']
        else :
            # If we're in AssertComplex and we're not duplicating reals...
            if foundFType == 'real' :
                retTypes=[]
            else :
                retTypes=['integer']
    elif expectedPrecision == 32 or expectedPrecision == 64:
        # This logic is probably not correct.
        if foundFType == 'integer' :
            # Processing integer-found.
            # mlr - fingers crossed.
            retTypes=['integer','real']
        elif not 'complex' in foundFTypes :
            # Processing case where we're not combining with complex.
            # mlr 2 - ???
            retTypes=['integer','real']
            # retTypes=['real']
        else :
            if foundFType == 'real' :
                # Tom asserts that finding a real when expecting complex should be an error.
                # retTypes=['complex']
                retTypes=[]
            else :
                retTypes=['integer','real','complex']
                #? retTypes=['real','complex']
    return retTypes

# Compare with allowedExpected
def makeExpectedFTypesWithoutExpectedPrecision( \
        foundFType,foundFTypes=['ERROR']):
    if 'ERROR' in foundFTypes:
        print('makeExpectedFTypesWithoutExpectedPrecision::ERROR')
    retTypes = []
    # if found is something, then expected is in...
    if 'logical' == foundFType :
        retTypes = ['logical']
    elif 'integer' == foundFType :
        # Is Finding an integer when expecting complex an error too?
        retTypes = ['integer','real']
    elif 'real' == foundFType :
        # Finding a real when expecting complex is an error. TLC
        retTypes = ['integer','real']
    elif 'complex' :
        # We're finding complexes.
        retTypes = ['integer','real','complex']
    return retTypes

def makeExpectedRanks(foundRank):
    ranks = [0]
    if foundRank != 0:
        ranks.append(foundRank)
    return ranks

def makeTolerances(expectedP, foundP) :
    "unless default (int) is found, collect all of the tolerances \
    found in eP and fP and return the maximum"
    tol = -1
    if type(expectedP) is list :
        ep = expectedP
    else :
        ep = [expectedP]
    if type(foundP) is list :
        fp = foundP
    else:
        fp = [foundP]
    lp = []
    if not 'def' in ep :
        lp = lp + ep 
    if not 'def' in fp :
        lp = lp + fp
    if lp == [] :
        # 2013-1022 MLR Fix this!!!
        # print('tolerance error! setting lp to 64.')
        lp = [64]
    tol = max(lp)
    return [tol]

class AssertRealArrayArgument:
    def __init__(self,aName,eft,ep,er,fft,fp,fr,tol):
        # print(' ',eft,ep,er,fft,fp,fr,tol)
        self.assertionName = aName
        self.expectedFType = eft
        self.expectedPrecision = ep
        self.expectedRank = er
        self.foundFType = fft
        self.foundPrecision = fp
        self.foundRank = fr
        self.tolerance = tol
        # ArrayDescriptions
        self.expectedDescription = None
        self.foundDescription = None
        # Now set them...
        self.updateDescriptions()

    def updateDescriptions(self):
        self.expectedDescription = ArrayDescription( \
                                                     self.expectedFType, \
                                                     self.expectedPrecision, \
                                                     self.expectedRank )
        self.foundDescription = ArrayDescription( \
                                                  self.foundFType, \
                                                  self.foundPrecision, \
                                                  self.foundRank )
        return

    def getAssertionName(self):
        return self.assertionName
    def getExpectedDescription(self):
        return self.expectedDescription
    def getFoundDescription(self):
        return self.foundDescription
    def getTolerance(self):
        return self.tolerance

# See allowedPrecisions.
def makeExpectedPrecisions(foundPrecision,foundFType='real',expectedFType=''):
    if expectedFType == 'integer' :
        # New style, eFT known.
        expectedPrecisions = [32,64]
    elif expectedFType == 'real' :
        expectedPrecisions = [32]
        if foundPrecision > 32 :
            expectedPrecisions.append(foundPrecision)
    elif expectedFType == 'complex' :
        expectedPrecisions = [32]
        if foundPrecision > 32 :
            expectedPrecisions.append(foundPrecision)
    elif expectedFType == 'logical' :
        expectedPrecisions = ['']
    else:
        # Old style, eFT not known.
       if foundFType == 'logical' :
           expectedPrecisions = ['']
       elif foundFType == 'integer' :
           expectedPrecisions = [32,64]
           # expectedPrecisions = ['def']
       else :
           expectedPrecisions = [32]
           #expectedPrecisions = ['def',32]
           #expectedPrecisions = ['defXXX',32]
           if foundPrecision > 32 :
              expectedPrecisions.append(foundPrecision)
    return expectedPrecisions

def ca_MakeAllowedPrecisions(foundFType) :
    precs = []
    ap = allowedPrecisions(foundFType,pFound='64')
    for i in ap :
        if i == 'def' :
            precs = precs + [i]
        else :
            precs = precs + [int(i)]
    return precs

def constructASSERTS(assertionName,foundFTypes=['real','complex'],foundRanks = [0,1,2,3,4,5],patchIntXX=False):

    AssertList = []

    # Note:  expectedPrecision <= foundPrecision
    # Note:  Need to eliminate redundancy of real asserts that can arise in AssertComplex.
    #        I.e. remove real-real comparisons when complex is available.
 
    # expectedFTypes -> 'integer' if expectedPrecision 'def' else 'real'
    # was tolerances = [32,64], but replaced by the following:
    # tolerances = max expectedPrecisions & foundPrecisions
    # expectedPrecisions = ['def',32,64] that are < foundPrecision
    # expectedRanks(foundRank) -> [0,foundRank]
    # + passed in foundFTypes = ['real','complex']
    # + passed in foundFTypes = ['real']
    # foundPrecisions = [32,64] replaced with ca_MakeAllowedPrecisions
    # + passed in foundRanks = [0,1,2,3,4,5]

# -> foundFTypes --> adding 'complex'

# THE MAIN LOOP.
# May need a special case if we don't want to construct a real-real in AssertComplex...
# Many type-kind-rank combinations are not allowed.  The allowed combinations result from 
# some options depending on others.  These are implemented in the "make..." functions listed below.
# The argument foundFType (found Fortran Type) is a key independent variable, which drives the types
# chosen for other arguments.
#
# The variable a contains the arguments for the specialized AssertEqual being generated.
# The constraintASSERTEQUAL object is the specialized routine, which is then used to
# construct the module.
#
# To change the list of asserts constructed, one can either change the logic implemented in the
# network of "make..." functions below, or one could add other routineUnits to AssertList, as long
# as it make sense to include them in the list (and interface block).
#
#    print ('1000: ',foundRanks)

# The following provides basic functionality, originally thought through for
# the default integer type.  Moving to multiple integer kinds has a gap in 
# combinatorial type coverage. So we'll handle the gaps below.
#
    AssertList = \
    [ \
      #test      a \
      constraintASSERT(a.getAssertionName(),\
                            a.getExpectedDescription(),\
                            a.getFoundDescription(),\
                            a.getTolerance()\
                            ) \
    for a in \
    flattened( \
               [[[[[[[ \
                       AssertRealArrayArgument(assertionName,eft,ep,er,fft,fp,fr,tol) \
                       for tol in makeTolerances(ep,fp) ] \
                       for eft in makeExpectedFTypes(ep,fft,foundFTypes=foundFTypes) ]  \
                       for ep in makeExpectedPrecisions(fp,foundFType=fft)  ] \
                       for er in makeExpectedRanks(fr) ] \
                       for fp in ca_MakeAllowedPrecisions(fft) ] \
                       for fft in foundFTypes ] \
                       for fr in foundRanks ] \
                       )]

    ## To insert by hand, one might try the following (sketch...)...
    ## a = AssertRealArrayArgument('integer','def','1','integer','def','1',0)
    ## AssertList += [MyConstraintAssertEqual(a.getExpectedDescription(),a.getFoundDescription())]
    ## Any specialization of routineUnit should work here...
    ## The code is generated when the routineUnit is instantiated, so that support code would
    ## need to be available.

    #
    # For real and complex r32/r64 and c32/c64 precisions are the same kind of thing,
    # while i32/i64 for integers have completely different meanings. Therefore our
    # default logic, which does not allow the combination [expected-x64, found-x32],
    # is incorrect for integer.
    #
    # if len(foundFTypes) == 1 :
    if patchIntXX:
       if foundFTypes[0] != 'integer' :
          nameList = [ assertionName ]
          eftList  = [ 'integer' ]
          epList   = [ 64 ]
          fftList  = foundFTypes
          fpList   = [ 32 ]
          frList   = foundRanks
          # erList   = []
          # tolList  = [ 32 ]
          
          for name in nameList:
              for eft in eftList:
                  for ep in epList:
                      for fft in fftList:
                          for fp in fpList:
                              for fr in frList:
                                  for er in makeExpectedRanks(fr):
                                      for tol in makeTolerances(ep,fp):
                                          a = AssertRealArrayArgument( name, eft, ep, er, fft, fp, fr, tol )
                                          AssertList.append(\
                                              constraintASSERT(a.getAssertionName(),\
                                                  a.getExpectedDescription(),\
                                                  a.getFoundDescription(),\
                                                  a.getTolerance()\
                                                  ))
    
    return AssertList

def constructDeclarations(assertionRoutines,basename='',foundRanks=[]):
    "Construct declarations to be used at the beginning of the Module."
    # foundRanks works into this... how? mlr 2014-0407
    declarations = \
        [ \
          declaration('uses',declareUSES(internalRoutines=assertionRoutines)), \
          declaration('discipline',declareDISCIPLINE()), \
          declaration('exports',declareEXPORTS(assertionRoutines,basename)), \
          declaration('exportsParameters',declareEXPORTS_PARAMETERS()) \
          ]
    return declarations



def declareUSES_() :
    "Set up use statements for SupportModule."
    return \
"""
   use Params_mod
   use AssertBasic_mod
   use Exception_mod
   use SourceLocation_mod
   use StringConversionUtilities_mod
"""


# Use "global" declareDISCIPLINE()
def declareEXPORTS_(exportedRoutines) :
    retPublic=""
    for aRoutineName in exportedRoutines :
         retPublic = retPublic + \
"""   public :: """+aRoutineName+"""
"""
    retPublic = retPublic + \
"""
   public :: valuesReport
   public :: differenceReport

   public :: vectorNorm
   public :: isWithinTolerance
   
"""
#    print ('3000: retPublic: ',retPublic)
    return retPublic

def declareEXPORTS_INTERNAL_(exportedRoutines) :
    retPublic=""
    for aRoutineName in exportedRoutines :
         retPublic = retPublic + \
"""   public :: """+aRoutineName+"""
"""
#-    retPublic = retPublic + \
#-"""
#-   public :: valuesReport
#-   public :: differenceReport
#-
#-   public :: vectorNorm
#-   public :: isWithinTolerance
#-   
#-"""
#    print ('3000: retPublic: ',retPublic)
    return retPublic


# Use "global" declareEXPORTS_PARAMETERS():

def constructSupportModuleDeclarations(exportedRoutineNames=[]):
#    print ('2000: eRN: ',exportedRoutineNames)
# Start main of constructSupportModuleDeclarations
    declarations = \
        [ \
          declaration('uses',declareUSES_()), \
          declaration('discipline',declareDISCIPLINE()), \
          declaration('exports',declareEXPORTS_(exportedRoutineNames)), \
          declaration('exportsParameters',declareEXPORTS_PARAMETERS()) \
          ]
    return declarations

def constructInternalModuleDeclarations(exportedRoutineNames=[]):
#    print ('2000: eRN: ',exportedRoutineNames)
# Start main of constructSupportModuleDeclarations
    declarations = \
        [ \
          declaration('uses',declareUSES_()), \
          declaration('discipline',declareDISCIPLINE()), \
          declaration('exports',declareEXPORTS_INTERNAL_(exportedRoutineNames)), \
          declaration('exportsParameters',declareEXPORTS_PARAMETERS()) \
          ]
    return declarations


def constructSupportModule(baseName='AssertArraysSupport',assertionShortNames=[],foundRanks=[],maxRank=5):
    # Just a rename to capture an idea.  Will fix later. MLR
#    exportedRoutineNames = ['assert'+i+'_internal' for i in assertionShortNames]
    exportedRoutineNames = []
#    print ('1500: exportedRoutineNames: ',exportedRoutineNames)
    m1 = module(baseName+'_mod')
    m1.setFileName(baseName+'.F90')
    # Note exportedRoutineNames may be empty but still makes a bunch of declarations.
    [m1.addDeclaration(d) for d in constructSupportModuleDeclarations(exportedRoutineNames)]
    m1.addInterfaceBlock(constructDifferenceReportInterfaceBlock())
    m1.addInterfaceBlock(constructValuesReportInterfaceBlock())
    # Generate internals for all ranks.
    m1.addInterfaceBlock(constructVectorNormInterfaceBlock(foundRanks=range(maxRank+1)))
    m1.addInterfaceBlock(constructIsWithinToleranceInterfaceBlock(foundRanks=range(maxRank+1)))
    
    # The following will add "assert" to the shortname.
#+    for assertionShortName in assertionShortNames:
#+        m1.addInterfaceBlock(constructAssertInternalInterfaceBlock(assertionShortName,foundRanks=foundRanks,expose=True),expose=True)
#?        m1.addInterfaceBlock(constructAssertInternalInterfaceBlock(assertionShortName,expose=True),expose=True)

    return m1

def constructInternalModule(baseName='AssertArraysInternal',assertionShortName="",exportedRoutineNames=[],maxRank=[]):
# oops -- hardwired longname here
    baseName_ = baseName+"assert"+str(assertionShortName)
    m1 = module(baseName_+'_mod')
    m1.setFileName(baseName_+'.F90')
    
    [m1.addDeclaration(d) for d in constructInternalModuleDeclarations(exportedRoutineNames)]

    m1.addInterfaceBlock(constructDifferenceReportInterfaceBlock())
    m1.addInterfaceBlock(constructValuesReportInterfaceBlock())
    # Generate internals for all ranks.
    m1.addInterfaceBlock(constructVectorNormInterfaceBlock(foundRanks=range(maxRank+1)))
    m1.addInterfaceBlock(constructIsWithinToleranceInterfaceBlock(foundRanks=range(maxRank+1)))
    
    
    # The following will add "assert" to the shortname.
    m1.addInterfaceBlock(constructAssertInternalInterfaceBlock(assertionShortName,expose=True),expose=True)

    return m1


def constructModule(baseName='AssertReal',assertionShortNames=['NOP'],foundFTypes=['real'],foundRanks=[],patchIntXX=False):
#    assertionShortNames=['Equal']
#    assertionShortNames=['GreaterThan']
    assertionNames=['assert'+Suffix for Suffix in assertionShortNames]
    rankName = ''
    for iRank in foundRanks : rankName = rankName + str(iRank)
    "A main test of how to construct the module."
    m1 = module(baseName+rankName+'_mod')
    m1.setFileName(baseName+rankName+'.F90')
    [m1.addDeclaration(d) for d in constructDeclarations(assertionNames,basename=baseName,foundRanks=foundRanks)]
    # add interface blocks (and the implementations)
    foundRanks1=foundRanks
    if (not 0 in foundRanks) : foundRanks1 = [0] + foundRanks
#! Moved to support module.        
#!    m1.addInterfaceBlock(constructVectorNormInterfaceBlock(foundRanks=foundRanks1))
#!    m1.addInterfaceBlock(constructIsWithinToleranceInterfaceBlock(foundRanks=foundRanks1))
    for assertionName in assertionShortNames:
        m1.addInterfaceBlock(constructAssertInterfaceBlock(assertionName,foundFTypes=foundFTypes,foundRanks=foundRanks,patchIntXX=patchIntXX))
        
# The following lines have been moved to the support module.
#    m1.addInterfaceBlock(constructDifferenceReportInterfaceBlock())
#    m1.addInterfaceBlock(constructValuesReportInterfaceBlock())
#

# The following generates internal routines.  Moved to support module.
    # *in test*
    # This is where the "list of asserts" is generated.
#???    for assertionName in assertionShortNames:
#???        m1.addInterfaceBlock(constructAssertInternalInterfaceBlock(assertionName,foundRanks=foundRanks))

# Old removal. 2014-0414
    # add individual routine units
    #    m1.addRoutineUnit(makeThrowDifferentValues()) # arg. provides a routine unit.
    return m1

def filePreamble(filename):
    return """

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!  This file '"""+filename+"""' is automatically generated by
!  'GenerateAssertsOnArrays.py'.  Changes made here will be
!  overwritten the next time that script is run.
!
!  2013-0722 MLR Michael.L.Rilee-1@nasa.gov
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

"""

def addFileToMakefile(fileName,includeFile="generated.inc"):
    with open(includeFile,'a') as f:
        f.write('GENERATED_CODE += ' + fileName + '\n')
    return

def addModToF90Include(modName,includeFile="AssertArrays.fh",postFix=""):
    with open(includeFile,'a') as f:
        f.write('   use ' + modName + postFix + '\n')
    return

def makeGeneratedInclude(includeFile="generated.inc"):
    with open(includeFile,'w') as f:
        f.write('# This file automatically generated. Do not modify.\n')
    return

relationalOperatorShortNames=["NotEqual",\
                            "Equal",\
                            "GreaterThan",\
                            "GreaterThanOrEqual",\
                            "LessThan",\
                            "LessThanOrEqual",\
                            "RelativelyEqual"]
def makeModuleReal(maxRank=5):
    assertionShortNames=relationalOperatorShortNames
#    assertionShortNames=["NotEqual",\
#                        "Equal",\
#                        "GreaterThan",\
#                        "GreaterThanOrEqual",\
#                        "LessThan",\
#                        "LessThanOrEqual",\
#                        "RelativelyEqual"]
    for iRank in range(0,maxRank+1):
        foundRanks = [iRank]
        mod = constructModule( \
                            assertionShortNames=assertionShortNames,\
                            foundRanks=foundRanks)
        with open(mod.getFileName(),'w') as f:
            f.write(filePreamble(mod.getFileName()))
            f.write('\n'.join(mod.generate()))
        addFileToMakefile(mod.getFileName())
        for iName in assertionShortNames : 
            addModToF90Include(mod.getName(),postFix=", only : " + 'Assert'+iName)
    print('makeModuleReal: done')
    return

def makeModuleComplex(maxRank=5):
    assertionShortNames=['NotEqual','Equal','RelativelyEqual']
    for iRank in range(0,maxRank+1):
        foundRanks = [iRank]
        mod = constructModule(baseName='AssertComplex',\
                          assertionShortNames=assertionShortNames,\
                          foundFTypes=['real','complex'],\
                          foundRanks=foundRanks,\
                          patchIntXX=True)
        with open(mod.getFileName(),'w') as f:
            f.write(filePreamble(mod.getFileName()))
            f.write('\n'.join(mod.generate()))
        addFileToMakefile(mod.getFileName())
        for iName in assertionShortNames : 
           addModToF90Include(mod.getName(),postFix=", only : " + 'assert'+iName)
    print('makeModuleComplex: done')
    return

def makeModuleInteger(maxRank=5):
    # assertionShortNames=['Equal']
    assertionShortNames=\
        ["NotEqual",\
        "Equal",\
        "GreaterThan",\
        "GreaterThanOrEqual",\
        "LessThan",\
        "LessThanOrEqual"]

    for iRank in range(0,maxRank+1):
        foundRanks=[iRank]
        mod = constructModule(baseName='AssertInteger',\
                          assertionShortNames=assertionShortNames,\
                          foundFTypes=['integer'],\
                          foundRanks=foundRanks)
        with open(mod.getFileName(),'w') as f:
            f.write(filePreamble(mod.getFileName()))
            f.write('\n'.join(mod.generate()))
        addFileToMakefile(mod.getFileName())
# Not ready to add integers to AssertArrays...
        for iName in assertionShortNames : 
           addModToF90Include(mod.getName(),postFix=", only : " + 'assert'+iName)       
    print('makeModuleInteger: done')
    return

# def makeModuleLogical():
#     mod = constructModule(baseName='AssertLogical1',foundFTypes=['logical'])
#     with open(mod.getFileName(),'w') as f:
#         f.write(filePreamble(mod.getFileName()))
#         f.write('\n'.join(mod.generate()))
#     print('makeModuleInteger: done')
#     return

def makeSupportModule(assertionShortNames=[],maxRank=5):
#        print ('1000: aSN: ',assertionShortNames)
        mod = constructSupportModule(assertionShortNames=assertionShortNames,maxRank=maxRank)
        with open(mod.getFileName(),'w') as f:
            f.write(filePreamble(mod.getFileName()))
            f.write('\n'.join(mod.generate()))
        addFileToMakefile(mod.getFileName())
#        for iName in assertionShortNames : 
#           addModToF90Include(mod.getName(),postFix=", only : " + 'assert'+iName)
        print('makeSupportModule: done')
        return

def makeInternalModule(assertionShortNames=[],maxRank=5):
        for iAssertion in assertionShortNames :
                mod = constructInternalModule(assertionShortName=iAssertion,maxRank=maxRank)
                with open(mod.getFileName(),'w') as f:
                    f.write(filePreamble(mod.getFileName()))
                    f.write('\n'.join(mod.generate()))
                addFileToMakefile(mod.getFileName())
#???                
#                Do I need to addModToF90Include?

#????
#        for iName in assertionShortNames : 
#           addModToF90Include(mod.getName(),postFix=", only : " + 'assert'+iName)
        print('makeInternalModule: done')
        return

    
def main(maxRank=5):
    # Make the modules for the different types...
    #
    makeGeneratedInclude(includeFile="generated.inc")
    #++
    makeModuleReal(maxRank=maxRank)
    #++
    makeModuleComplex(maxRank=maxRank)
    #? The following requires testing.
    makeModuleInteger(maxRank=maxRank)
    #? Just started...
    #- makeModuleLogical()

    # Make the routines that do the work.
    makeInternalModule(assertionShortNames=relationalOperatorShortNames,maxRank=maxRank)
    
    # Make the intermediate routines.
    makeSupportModule(assertionShortNames=relationalOperatorShortNames,maxRank=maxRank)
    
    return

if __name__ == "__main__":
    main(maxRank=int(args.maxRank or 5))
