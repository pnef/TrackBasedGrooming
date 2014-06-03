# --------------------------------------------- #
# Makefile for FCNC code                        #
# Pascal Nef, March 6th 2014                    #
#                                               #
# Note: source setup.sh before make             #
# --------------------------------------------- #

CXXFLAGS =   -O2 -Wall 
LINKLIBS := $(FASTJETLOCATION)/lib/libScJet.a

.PHONY: clean debug all

all: FCNC

FCNC:  FCNC.so FCNCTools.so FCNCAnalysis.so $(LINKLIBS)
	$(CXX) FCNC.so FCNCTools.so FCNCAnalysis.so -o $@.exe \
	$(CXXFLAGS) -Wno-shadow  \
	`root-config --glibs`  \
	-L$(FASTJETLOCATION)/lib `$(FASTJETLOCATION)/bin/fastjet-config --libs --plugins ` -lSubjetJVF -lVertexJets -lJetCleanser \
	-L$(PYTHIA8LOCATION)/lib -lpythia8 -llhapdfdummy \
	-L$(BOOSTLIBLOCATION) -lboost_program_options 

FCNC.so: FCNC.C    FCNCTools.so FCNCAnalysis.so $(LINKLIBS)
	$(CXX) -o $@ -c $<  \
	$(CXXFLAGS) -Wno-shadow -fPIC -shared \
	`$(FASTJETLOCATION)/bin/fastjet-config --cxxflags --plugins` -lSubjetJVF -lVertexJets -lJetCleanser \
	-I$(PYTHIA8LOCATION)/include \
	-I $(BOOSTINCDIR) \
	`root-config --cflags` 

FCNCTools.so : FCNCTools.cc FCNCTools.h $(LINKLIBS)
	$(CXX) -o $@ -c $<  \
	$(CXXFLAGS) -Wno-shadow -fPIC -shared \
	`$(FASTJETLOCATION)/bin/fastjet-config --cxxflags --plugins` -lSubjetJVF -lVertexJets -lJetCleanser \
	-I$(PYTHIA8LOCATION)/include \
	`root-config --cflags --libs`

FCNCAnalysis.so : FCNCAnalysis.cc FCNCAnalysis.h $(LINKLIBS)
	$(CXX) -o $@ -c $<  \
	$(CXXFLAGS) -Wno-shadow -fPIC -shared \
	`$(FASTJETLOCATION)/bin/fastjet-config --cxxflags --plugins` -lSubjetJVF -lVertexJets  -lJetCleanser \
	-I$(PYTHIA8LOCATION)/include \
	`root-config --cflags --libs`   

clean:
	rm -rf *.exe
	rm -rf *.o
	rm -rf *.so
	rm -f *~
