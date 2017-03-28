
CXXFLAGS += -I${CLHEP_INCLUDE_DIR}

LDFLAGS += -L${CLHEP_LIB_DIR}

ROOTCFLAGS    = $(shell root-config --cflags --glibs )
ROOTLIBS      = $(shell root-config --libs)

%:%.cc
	$(CXX) $(CXXFLAGS) $(LDFLAGS) $(ROOTCFLAGS) $(ROOTLIBS) -lRooFit -lHtml -lMinuit $< -o $@  

%TARGETS = Sterile_fit Sterile_fit_scan simple_t2k simple_t2k_scan t2k_simple_scan2 t2k_simple_scan3 simple_t2k_2D t2k_simple_2D2 t2k_simple_2D3 simple_t2k_2D0
TARGETS = simple_DUNEBinning simple_t2kBinning

all: $(TARGETS) 

clean:
	rm -rf $(TARGETS)

 NugenDeuteriumGen :  NugenDeuteriumGen.o
	g++ -o  NugenDeuteriumGen  NugenDeuteriumGen.o -L${ROOTSYS}/lib $(ROOTLIBS) $(ROOTCFLAGS) -lm -lc
 NugenDeuteriumGen.o :  NugenDeuteriumGen.cc
	g++ -c ${ROOTCFLAGS} -I${DOGS_PATH}/DCNuGen2 -I${DOGS_PATH}/DCBase -lRooFit -lHtml -lMinuit -I${DOGS_PATH}/DCGeo -I${DOGS_PATH}/DCValidity -I${DOGS_PATH}/DCEvent  NugenDeuteriumGen.cc

