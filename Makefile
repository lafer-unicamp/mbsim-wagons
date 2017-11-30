# Name of the executable
Target = src/main
exec_name = bin/test

# MBSim installation directory
mbsimdir = /usr/local/mbsim-env/

# Defining sources
SRCEXT = cpp
sources = $(shell find $(SRCDIR) -name "*."$(SRCEXT))

# Custom includes
incl = -I include

# Directory paths
BINDIR = bin
BUILDDIR = build
SRCDIR = src

# Do not edit the following lines
CXX = g++-4.8
libsources = 
objects = $(patsubst $(SRCDIR)/%,$(BUILDDIR)/%,$(sources:.$(SRCEXT)=.o))
CPPFLAGS= -g3 -m64 -std=c++11 -Wall -Wfatal-errors -Werror -Wno-unknown-pragmas -fopenmp \
`pkg-config --cflags mbsim`
#CPPFLAGS= -g3 -m64 -std=c++11 -Wall -D_GLIBCXX_USE_CXX11_ABI=0 -Wfatal-errors -Werror -Wno-unknown-pragmas -fopenmp \
`pkg-config --cflags mbsim`
#`$(mbsimdir)/bin/mbsim-config --cflags`
	#CPPFLAGS= -m32 -g3 -Wall -U_GLIBCXX_USE_CXX11_ABI -Wabi-tag -Werror -Wno-unknown-pragmas `pkg-config --cflags mbsim`
$(Target) : $(objects)
	@echo " Linking..."
	$(CXX) $^  --trace -o $(exec_name) `pkg-config --libs mbsim` -lpthread -lm -Wl,-rpath,$(mbsimdir)lib

# Build objects
$(BUILDDIR)/%.o: $(SRCDIR)/%.$(SRCEXT)
	@echo " Building..."
	@mkdir -p $(BUILDDIR)
	$(CXX) $(CPPFLAGS) $(incl) -c -o $@ $<
	
%.d: %.cc
	set -e; $(CXX) -MM $(CPPFLAGS) $< \
	  | sed 's/\(.*\)\.o[ :]*/\1.o \1.d : \g' > $@; \
	  [ -s $@ ] || rm -f $@
	
include $(sources:.cpp:.d)
	
.PHONY : clean
clean :
	@echo " Cleaning..."
	$(RM) -r $(BUILDDIR) $(exec_name)
