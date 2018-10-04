
APP_NAME      := res_ana
SOURCE_FILES  := analysis.cc angles_ana.cc 
INSTALL_DIR   := ..

USES_RFIO     := no
USES_ORACLE   := no
USES_GFORTRAN := yes
HOME:=/lustre/nyx/hades/user/shower

include $(HADDIR)/hades.def.mk

#LIB_DIRS += $(PLUTODIR) ${HOME}/usr/lib64
#INC_DIRS += ${HOME}/usr/include
#HYDRA_LIBS    += -lParticleEvent


.PHONY:  default
#default: clean build install
default: build install

include $(HADDIR)/hades.app.mk

