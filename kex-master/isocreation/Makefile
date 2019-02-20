# Notes: 
# * using -std=c++0x for std::remove_const<>
# * defining __ACLIC__ to signal that PHAST stuff (filling of objects)
#   is not required

SHELL=/bin/bash

CPPFLAGS=$(shell root-config --cflags) -D__ACLIC__
CXXFLAGS=-fpic -Wall -Wextra -Woverloaded-virtual -Wno-unused-parameter -Wno-unused-variable -g -pipe -I. 
LDFLAGS= -L /lib64  -fpic -pipe  $(shell root-config --ldflags --libs) -lpthread  -lz



# additional source files which don't need a dictionary
SRC=isotropic_cluster.cc 

DEPS=$(patsubst %.cc,%.d,$(DSRC) $(SRC))
DICTS=$(patsubst %.cc,%Dict.cc,$(DSRC))
DICTOBJS=$(patsubst %.cc,%Dict.o,$(DSRC))
CLASSOBJS=$(patsubst %.cc,%.o,$(DSRC))
OBJS=$(patsubst %.cc,%.o,$(SRC))


isotropic_cluster: $(DICTOBJS) $(CLASSOBJS) $(OBJS)
	$(CXX) -o $@ $^ $(LDFLAGS) 



.PHONY: clean
clean:
	rm -f *.so *.par isotropic_cluster  {,user/}*{.o,.gch,.d,.d.*,Dict.cc,Dict.h}


.PHONY: distclean
distclean: clean
	rm -f {,user/}*.d


.PHONY: proofclean
proofclean:
	rm -rf $$HOME/.proof /tmp/$$USER/.getproof xpdtut
	pkill xrootd || true


.PHONY: debug
debug:
	@echo "LDLIBS: $(LDLIBS)"
	@echo "DEPS:   $(DEPS)"


# adapted from https://www.gnu.org/software/make/manual/make.html#Automatic-Prerequisites
$(DEPS): %.d: %.cc
	@set -e; rm -f $@; \
              $(CXX) -MM -MF $@.$$$$ $(CPPFLAGS) $<; \
              sed 's,\($*\)\.o[ :]*,\1.o $@: ,g' < $@.$$$$ > $@; \
              rm -f $@.$$$$


ifneq ($(MAKECMDGOALS),clean)
ifneq ($(MAKECMDGOALS),distclean)
-include $(DEPS)
endif
endif


# First attempt at using precompiled headers, but further research is necessary.
#
# "Only one precompiled header can be used in a particular compilation."
# cf. http://gcc.gnu.org/onlinedocs/gcc/Precompiled-Headers.html
#
# see also http://stackoverflow.com/questions/58841/precompiled-headers-with-gcc
#
# gcc -H
#
#GCHS=$(patsubst %.cc,%.h.gch,$(DSRC) $(SRC))
#$(DICTOBJS) $(CLASSOBJS) $(OBJS): $(GCHS)
#
#$(GCHS): %.h.gch: %.h
#	$(CXX) $(CPPFLAGS) $(CXXFLAGS) $^

