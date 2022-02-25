PACKAGES=mbsim

SRCDIR:=$(dir $(lastword $(MAKEFILE_LIST)))


# This file builds a "EXECNAME" executable from all sources found under $(SRCDIR).
# SRCDIR must be set .
# PACKAGES must also be set .
# Moreover the LDFLAGS, CPPFLAGS and CXXFLAGS are honored if set externally.

BUILDDIR = $(SRCDIR)build
EXECNAME = $(BUILDDIR)/teste

# Enable VPATH builds with sources at SRCDIR
VPATH=$(SRCDIR)

# use a sources all *.cc file under SRCDIR
SOURCES:=$(shell (cd $(SRCDIR)src; find -name "*.cpp"))
# object and dependency files (derived from SOURCES)
# OBJECTS=$(SOURCES:.cpp=.o)
OBJECTS=$(addprefix $(BUILDDIR)/,$(SOURCES:.cpp=.o))
# OBJECTS=$(addprefix $(BUILDDIR),$(notdir $(SOURCES:.cpp=.o)))
DEPFILES=$(addprefix $(BUILDDIR)/,$(SOURCES:.cpp=.o.d))

# enable C++17
CXXFLAGS += -ggdb -std=c++17 -D_USE_MATH_DEFINES -I$(SRCDIR)include

# platform specific settings
ifeq ($(PLATFORM),Windows)
  SHEXT=.dll
  EXEEXT=.exe
  PIC=
  LDFLAGSRPATH=
  WIN=1
else
  SHEXT=.so
  EXEEXT=
  PIC=-fpic
  LDFLAGSRPATH=-Wl,--disable-new-dtags
  WIN=0
endif

# default target
all: $(EXECNAME)$(EXEEXT)

# FMI export target
fmiexport: mbsimfmi_model$(SHEXT)

# mingw -Wl,-Map generates some dependencies named rtr0*.o which needs to be removed

# link main executable with pkg-config options from PACKAGES (runexamples.py executes always ./main)
$(EXECNAME)$(EXEEXT): $(OBJECTS)
	$(CXX) -Wl,-Map=$@.linkmap -o $@ $(OBJECTS) $(LDFLAGS) $(shell pkg-config --libs $(PACKAGES))
	@sed -rne "/^LOAD /s/^LOAD (.*)$$/ \1 \\\/p" $@.linkmap | grep -Ev rtr[0-9]+\.o > $@.d2 || true
	@test $(WIN) -eq 0 && (echo "$@: \\" > $@.d && cat $@.d2 >> $@.d && rm -f $@.linkmap $@.d2) || (rm -f $@.d2)

# FMI export target
mbsimfmi_model$(SHEXT): $(OBJECTS)
	$(CXX) -shared $(LDFLAGSRPATH) -Wl,-rpath,\$$ORIGIN,-Map=$@.linkmap -o $@ $(OBJECTS) $(LDFLAGS) $(shell pkg-config --libs $(PACKAGES))
	@sed -rne "/^LOAD /s/^LOAD (.*)$$/ \1 \\\/p" $@.linkmap | grep -Ev rtr[0-9]+\.o > $@.d2 || true
	@test $(WIN) -eq 0 && (echo "$@: \\" > $@.d && cat $@.d2 >> $@.d && rm -f $@.linkmap $@.d2) || (rm -f $@.d2)

rpath: $(OBJECTS)
	$(CXX) -o $@ $^ $(LDFLAGS) $(shell pkg-config --libs $(PACKAGES))  $(shell pkg-config --libs-only-L $(PACKAGES) | sed 's/-L/-Wl,-rpath,/g')

# compile source with pkg-config options from PACKAGES (and generate dependency file)
$(OBJECTS): $(BUILDDIR)%.o: src%.cpp | build
	@$(CXX) -MM -MP $(PIC) $(CPPFLAGS) $(CXXFLAGS) $(shell pkg-config --cflags $(PACKAGES)) $< > $@.d
	$(CXX) -c $(PIC) -o $@ $(CPPFLAGS) $(CXXFLAGS) $(shell pkg-config --cflags $(PACKAGES)) $<

# clean target: remove all generated files
clean:
	rm -f main$(EXEEXT) mbsimfmi_model$(SHEXT) $(OBJECTS) $(DEPFILES) main$(EXEEXT).d mbsimfmi_model$(SHEXT).d

# include the generated make rules (without print a warning about missing include files (at first run))
-include $(DEPFILES)
-include main$(EXEEXT).d
-include mbsimfmi_model$(SHEXT).d
