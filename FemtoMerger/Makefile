APP_NAME      := femtoMerge

SOURCE_FILES  := femtoMerge.cc
USES_RFIO     := no
USES_ORACLE   := yes
USES_GFORTRAN := yes

include $(HADDIR)/hades.def.mk

.PHONY:  default
default: clean build install

include $(HADDIR)/hades.app.mk
