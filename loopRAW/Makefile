APP_NAME      := analysisDST

#SOURCE_FILES  := analysisDSTsim_new.cc
SOURCE_FILES  := analysisDST_2step.cc
USES_RFIO     := no
USES_ORACLE   := yes
USES_GFORTRAN := yes

include $(HADDIR)/hades.def.mk

HYDRA_LIBS    += -lDst


.PHONY:  default
default: clean build install

include $(HADDIR)/hades.app.mk

