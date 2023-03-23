# Make sure to set up GENIE before building
ROOT_CXX_FLAGS := $(shell root-config --cflags) -fPIC
ROOT_LIBS := $(shell root-config --libs --evelibs) -lMathMore

GENIE_CXX_FLAGS := -I$(shell genie-config --topsrcdir) -I${LOG4CPP_INC} -O3
GENIE_LIBS := $(shell genie-config --libs) -L${LOG4CPP_LIB} -lxml2 -llog4cpp
GENIE_LIBS += -lEGPythia6 $(gsl-config --libs) -lGenVector
GENIE_LIBS += -L${LHAPDF_LIB} -lLHAPDF

# When including the Fortran code, you will likely need to link the executable
# to -lgfortran as well
#test_xsec: test_xsec.o libFortranWrapper.so
#	$(CXX) $(ROOT_LIBS) $(GENIE_LIBS) -L. -lFortranWrapper \
#	  -O3 -o $@ $^ -lgfortran

test_tensor: test_tensor.o libFortranWrapper.so
	$(CXX) $(ROOT_LIBS) $(GENIE_LIBS) -L. -lFortranWrapper \
          -O3 -o $@ $^ -lgfortran

libFortranWrapper.so: root_dictionary.o carray.o
	$(CXX) $(ROOT_LIBS) $(GENIE_LIBS) -O3 --shared -o $@ $^ -lgfortran

%.o: %.cxx
	$(CXX) $(ROOT_CXX_FLAGS) $(GENIE_CXX_FLAGS) -o $@ -c $^

%.o: %.f90
	gfortran -fPIC -O3 -o $@ -c $^

root_dictionary.o:
	$(RM) root_dictionary*.*
	rootcling -f root_dictionary.cxx $(GENIE_CXX_FLAGS) \
	  -c LinkDef.h
	$(CXX) $(ROOT_CXX_FLAGS) $(GENIE_CXX_FLAGS) \
	-o $@ -c root_dictionary.cxx
	$(RM) root_dictionary.cxx

.PHONY: clean

clean:
	$(RM) *.o test_tensor libFortranWrapper.so root_dictionary_rdict.pcm

