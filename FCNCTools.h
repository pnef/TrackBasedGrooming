#ifndef FCNCTOOLS_H
#define FCNCTOOLS_H 

#include <vector>
#include <math.h>
#include <string>

#include "fastjet/ClusterSequence.hh"
#include "fastjet/PseudoJet.hh"  

#include "Pythia8/Pythia.h"

#include "myFastJetBase.h"

using namespace std;

class FCNCTools {
    private:
        int m_test;



    public:
        FCNCTools();
        
        // methods
    float JetCharge(fastjet::PseudoJet jet,float kappa);
    bool   BosonMatch(fastjet::PseudoJet jet, vector<fastjet::PseudoJet> Bosons, float jetrad, int BosonID );
    bool   TruthMatchDR(fastjet::PseudoJet jet, vector<fastjet::PseudoJet> TruthJets, float dR, float &truthmatchpt);
    float  dRTruth(fastjet::PseudoJet jet, vector<fastjet::PseudoJet> TruthJets, float ptThreshold);

    fastjet::PseudoJet SubJetJVFGroomer(fastjet::PseudoJet jet, vector<fastjet::PseudoJet> HStracks, vector<fastjet::PseudoJet> PUtracks);
};

#endif

