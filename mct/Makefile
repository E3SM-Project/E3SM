#!/bin/make

MODULE		= mct

SRCS_mks	= Makefile			\
		  Makefile.conf.in		\
	          Makefile.dependency		\
		  Makefile.rules		\
		  configure			\
		  README			

SRCS_C		=

SRCS_f		=

SRCS_F		=

SRCS_f90	=

SRCS_F90	= m_MCTWorld.F90		\
		  m_AttrVect.F90		\
		  m_GlobalMap.F90	        \
		  m_GlobalSegMap.F90		\
		  m_GlobalSegMapComms.F90	\
		  m_Accumulator.F90		\
		  m_SparseMatrix.F90      	\
		  m_Navigator.F90		\
		  m_AttrVectComms.F90		\
		  m_AttrVectReduce.F90		\
		  m_AccumulatorComms.F90  	\
		  m_GeneralGrid.F90		\
		  m_GeneralGridComms.F90	\
		  m_SpatialIntegral.F90		\
		  m_MatAttrVectMul.F90    	\
		  m_Merge.F90			\
		  m_GlobalToLocal.F90		\
		  m_ExchangeMaps.F90		\
		  m_ConvertMaps.F90		\
		  m_SparseMatrixDecomp.F90	\
		  m_SparseMatrixToMaps.F90	\
		  m_SparseMatrixComms.F90 	\
		  m_SparseMatrixPlus.F90	\
		  m_Router.F90			\
		  m_Rearranger.F90        	\
		  m_Transfer.F90          

TEXFILES	= m_AttrVect.tex		\
		  m_GlobalMap.tex         	\
		  m_GlobalSegMap.tex      	\
		  m_Accumulator.tex       	\
		  m_SparseMatrix.tex      	\
		  m_Navigator.tex         	\
		  m_AttrVectComms.tex     	\
		  m_AccumulatorComms.tex  	\
		  m_GeneralGrid.tex		\
		  m_MatAttrVectMul.tex    	\
		  m_Merge.tex 			\
		  m_Router.tex

SRCS_ALL	= $(SRCS_f)		\
		  $(SRCS_F)		\
		  $(SRCS_f90)		\
		  $(SRCS_F90)		\
		  $(SRCS_C)

OBJS_ALL	= $(SRCS_f:.f=.o)	\
		  $(SRCS_F:.F=.o)	\
		  $(SRCS_f90:.f90=.o)	\
		  $(SRCS_F90:.F90=.o)	\
		  $(SRCS_C:.c=.o)

ARCH_TARGETS	= $(SRCS_mks) $(SRCS_ALL)
ALL_TARGETS	= lib$(MODULE).a

include Makefile.rules
include Makefile.conf

# targets

lib$(MODULE).a:	$(OBJS_ALL)
	$(RM) $@
	@ echo $(AR) $@ $(OBJS_ALL)
	@ $(AR) $@ $(OBJS_ALL)

doc:	$(TEXFILES)

do-it-anyway:

#________________________________________
# Regular targets

#________________________________________
# Just dependencies

include Makefile.dependency

#.
