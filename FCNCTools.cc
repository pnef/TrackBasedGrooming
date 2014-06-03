#include <math.h>
#include <vector>
#include <string>
#include <sstream>
#include <set>

#include "fastjet/ClusterSequence.hh"
#include "fastjet/PseudoJet.hh"  
#include "fastjet/JetDefinition.hh"
#include "fastjet/contrib/SubjetJVF.hh"

#include "FCNCTools.h"
#include "myFastJetBase.h"

#include "TRandom3.h"
#include "TError.h"
#include "TVector3.h"

using namespace std;

// Constructor 
FCNCTools::FCNCTools(){
    m_test = 0;
}

fastjet::PseudoJet FCNCTools::SubJetJVFGroomer(fastjet::PseudoJet jet, vector<fastjet::PseudoJet> HStracks, vector<fastjet::PseudoJet> PUtracks){

    fastjet::JetDefinition subjet_def(fastjet::kt_algorithm, 0.3);
    fastjet::contrib::SubjetJVF subjetJVF(subjet_def, -1);
    fastjet::PseudoJet result = subjetJVF(jet, HStracks, PUtracks);
//    cout << result.has_user_info<fastjet::contrib::SubjetJVFUserInfo>() << endl;

    return result;
    
}



float FCNCTools::JetCharge(fastjet::PseudoJet jet,float kappa){
  //Returns the jet charge with weighting factor kappa
  float charge=0.;
  for (unsigned int i=0; i<jet.constituents().size(); i++){
      if(! jet.constituents()[i].has_user_info<MyUserInfo>()) continue;
      charge+=jet.constituents()[i].user_info<MyUserInfo>().charge()*pow(jet.constituents()[i].pt(),kappa);
  }
  return charge/pow(jet.pt(),kappa);
}

bool FCNCTools::BosonMatch(fastjet::PseudoJet jet, vector<fastjet::PseudoJet> Bosons, float jetrad, int BosonID ){
    for (unsigned int i=0; i<Bosons.size(); i++){
        if (Bosons[i].user_info<MyUserInfo>().pdg_id() != BosonID) continue;
        if (Bosons[i].delta_R(jet)<jetrad){
            return true;
        }
    }
    return false;
}

// matching to truth jets: return true if dR to closest truth jet < dR
bool FCNCTools::TruthMatchDR(fastjet::PseudoJet jet, vector<fastjet::PseudoJet> TruthJets, float dR, float& truthmatchpt){
    float highestmatchpt = -999;
    for(unsigned int i=0; i<TruthJets.size(); ++i){
        if(jet.delta_R(TruthJets[i]) < dR) {
            if (TruthJets[i].pt() > highestmatchpt) highestmatchpt = TruthJets[i].pt();
        }
    }
    if(highestmatchpt  > 0) {truthmatchpt = highestmatchpt; return true;}
    else                    {truthmatchpt = -1;             return false;} 
}

// dR to closest truth jet above pT threshold
float FCNCTools::dRTruth(fastjet::PseudoJet jet, vector<fastjet::PseudoJet> TruthJets, float ptThreshold){
    float dR = 999.;
    for(unsigned int i=0; i<TruthJets.size(); ++i){
        if (TruthJets[i].pt() < ptThreshold ) continue;
        if (jet.delta_R(TruthJets[i]) < dR  ) dR = jet.delta_R(TruthJets[i]);
    }
    return dR;
}

