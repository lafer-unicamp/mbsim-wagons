PACKAGES=mbsim

SRCDIR:=$(dir $(lastword $(MAKEFILE_LIST)))
include $(SRCDIR)/default_build.mk
