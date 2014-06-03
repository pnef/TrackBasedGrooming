#include <math.h>
#include <vector>
#include <string>
#include <sstream>
#include <set>


#include "TFile.h"
#include "TTree.h"
#include "TClonesArray.h"
#include "TParticle.h"
#include "TDatabasePDG.h"
#include "TMath.h"
#include "TVector3.h"

#include "FCNCAnalysis.h"
#include "FCNCTools.h"

#include "myFastJetBase.h"
#include "fastjet/ClusterSequence.hh"
#include "fastjet/PseudoJet.hh"  
#include "fastjet/tools/Filter.hh"
#include "fastjet/Selector.hh"
#include "fastjet/tools/Subtractor.hh"
#include "fastjet/ClusterSequenceArea.hh"
#include "fastjet/contrib/SubjetJVF.hh"
#include "fastjet/contrib/VertexJets.hh"
#include "fastjet/contrib/JetCleanser.hh"
#include "fastjet/tools/JetMedianBackgroundEstimator.hh"

#include "Pythia8/Pythia.h"

using namespace std;

// ---------------- Selectors ----------------------
// ----------------------------------------------------------------------
fastjet::Selector SelectorPileupGhostTrack() {
return new SelectorWorkerPileupGhostTrack();
}

fastjet::Selector SelectorHardScatterGhostTrack() {
return new SelectorWorkerHardScatterGhostTrack();
}

fastjet::Selector SelectorPileupTrack() {
return new SelectorWorkerPileupTrack();
}

fastjet::Selector SelectorHardScatterTrack() {
return new SelectorWorkerHardScatterTrack();
}
// ----------------------------------------------------------------------




// Constructor 
FCNCAnalysis::FCNCAnalysis(){
    if(fDebug) cout << "FCNCAnalysis::FCNCAnalysis Start " << endl;
    ftest = 0;
    fDebug = false;
    fOutName = "test.root";
    tool = new FCNCTools();

    // jet def 
    m_jet_def               = new fastjet::JetDefinition(fastjet::antikt_algorithm, 0.4);
    m_jet_def_largeR_ALTAS  = new fastjet::JetDefinition(fastjet::antikt_algorithm, 1.0);

    if(fDebug) cout << "FCNCAnalysis::FCNCAnalysis End " << endl;
}

// Destructor 
FCNCAnalysis::~FCNCAnalysis(){
    delete tool;
    delete m_jet_def;
    delete m_jet_def_largeR_ALTAS;
}

// Begin method
void FCNCAnalysis::Begin(){
   // Declare TTree
   tF = new TFile(fOutName.c_str(), "RECREATE");
   tT = new TTree("EventTree", "Event Tree for FCNC");
    
   DeclareBranches();
   ResetBranches();
   

   return;
}

// End
void FCNCAnalysis::End(){
    
    tT->Write();
    tF->Close();
    return;
}

// GetParticles from pythia
bool FCNCAnalysis::GetPythiaParticles(Pythia8::Pythia* pythia8, event_type etype, int nPU){

    int n_runs = 1;
    if(etype == pileup){ n_runs = nPU;}

//    if(etype == pileup) cout << "generating nPU " << nPU << endl;


    for(int irun=0; irun<n_runs; ++irun){
        // next event 
        if (!pythia8->next()) return false;

        for (unsigned int ip=0; ip<pythia8->event.size(); ++ip){

            int origin = -1;
            if      (etype == hardscatter) origin = 0;
            else if (etype == pileup     ) origin = irun +1;

            fastjet::PseudoJet p(pythia8->event[ip].px(), pythia8->event[ip].py(), pythia8->event[ip].pz(),pythia8->event[ip].e() );
            p.reset_PtYPhiM(p.pt(), p.rapidity(), p.phi(), 0.);
            p.set_user_info(new MyUserInfo(pythia8->event[ip].id(), ip, pythia8->event[ip].charge(), origin));


            // tracks: charged particles with pt>0.5 GeV, |eta|<2.4 
            if(pythia8->event[ip].isFinal()    && 
               fabs(pythia8->event[ip].id())  !=11  && 
               fabs(pythia8->event[ip].id())  !=12  && 
               fabs(pythia8->event[ip].id())  !=13  && 
               fabs(pythia8->event[ip].id())  !=14  && 
               fabs(pythia8->event[ip].id())  !=16  && 
               pythia8->event[ip].isCharged()  && 
               pythia8->event[ip].pT() > 0.5   && 
               fabs(pythia8->event[ip].eta())<2.4){

                if      (origin == 0 ) {HStracks.push_back(p); p.set_user_info(new GhostUserInfo(false, false, true));}
                else if (origin >= 1 ) {PUtracks.push_back(p); p.set_user_info(new GhostUserInfo(false, true,  true));}

            }

            if(pythia8->event[ip].isFinal()         && 
               fabs(pythia8->event[ip].id())  !=11  && 
               fabs(pythia8->event[ip].id())  !=12  && 
               fabs(pythia8->event[ip].id())  !=13  && 
               fabs(pythia8->event[ip].id())  !=14  && 
               fabs(pythia8->event[ip].id())  !=16  && 
               pythia8->event[ip].pT()        > 0.5  ) {

               particlesForJets.push_back(p);

               if(origin == 0){
               particlesForTruthJets.push_back(p);
               }

            }

            // truth bosons: only looking at these that decay to two particles, this excluded intermediate boson that are their own children
            //               or in final state
            int absID = abs(pythia8->event[ip].id());
            if( (  (absID == 23 || absID == 24 || absID == 25 || absID == 32 || absID ==  99999999) && pythia8->event.daughterList(ip).size() ==2 ) || 
                (  (absID == 23 || absID == 24 || absID == 25 || absID == 32 || absID ==  99999999) && pythia8->event[ip].isFinal())) { 
                    if(fDebug){
                        cout << "FCNCAnalysis::AnalyzeEvent " ;
                        vector<int> daugthers = pythia8->event.daughterList(ip);
                        cout << "boson " << pythia8->event[ip].id();
                        if(daugthers.size() >0 ) cout << " with daughters " << pythia8->event[daugthers[0]].id() ;
                        if(daugthers.size() >1 ) cout << " and " << pythia8->event[daugthers[1]].id();
                        cout << endl;
                    }
                    Bosons.push_back(p);
                    if(fTNBosonsFilled == MaxNBosons) {cout << "Warning: More than " << MaxNBosons << " bosons found!" << endl; continue;}
                    fTBosonID [fTNBosonsFilled] = pythia8->event[ip].id();
                    fTBosonPt [fTNBosonsFilled] = pythia8->event[ip].pT();
                    fTBosonEta[fTNBosonsFilled] = pythia8->event[ip].eta();
                    fTBosonPhi[fTNBosonsFilled] = pythia8->event[ip].phi();
                    fTBosonM  [fTNBosonsFilled] = pythia8->event[ip].m();
                    fTNBosonsFilled++;
            }
        }
    } // end irun loop

    return true;
}

// Analyze
void FCNCAnalysis::AnalyzeEvent(int ievt, Pythia8::Pythia* pythia8, Pythia8::Pythia* pythia8PU, int nPU){
    if(fDebug) cout << "FCNCAnalysis::AnalyzeEvent Begin " << endl;

    // -------------------------
    if(fDebug) cout << "FCNCAnalysis::AnalyzeEvent Event Number " << ievt << endl;
    
    // reset branches 
    ResetBranches();
    

    // new event-----------------------
    fTEventNumber = ievt;
    fTNPV         = nPU;
    particlesForJets.clear();
    particlesForTruthJets.clear();
    Bosons.clear();
    PUtracks.clear();
    HStracks.clear();
    


    // get new event and return if failed
    bool ok(true);
    ok = GetPythiaParticles(pythia8,   hardscatter);
    ok = GetPythiaParticles(pythia8PU, pileup,     nPU);
    if (!ok) return;
    
    /// ------------------ Caluculate Rho ------------------------------
    fastjet::JetMedianBackgroundEstimator bge(fastjet::SelectorAbsRapMax(1.5), fastjet::JetDefinition(fastjet::kt_algorithm,0.4), fastjet::AreaDefinition(fastjet::active_area_explicit_ghosts));
    bge.set_particles(particlesForJets);
    fastjet::Subtractor theSubtractor(&bge);

    /*
    // truth jets  -----------------------------------------------
    fastjet::ClusterSequence csTruthSmall(particlesForTruthJets, fastjet::JetDefinition(fastjet::antikt_algorithm, 0.4));
    vector<fastjet::PseudoJet> myTruthJetsSmall = fastjet::sorted_by_pt(csTruthSmall.inclusive_jets(5.0)); 
    for(unsigned int iJ=0; iJ<myTruthJetsSmall.size(); ++iJ){
        if(fTNJetsSmallTruthFilled == MaxNJetSmallR) continue;
        fTJsmallTruthPt  [fTNJetsSmallTruthFilled] = myTruthJetsSmall[iJ].pt();
        fTJsmallTruthEta [fTNJetsSmallTruthFilled] = myTruthJetsSmall[iJ].eta();
        fTJsmallTruthPhi [fTNJetsSmallTruthFilled] = myTruthJetsSmall[iJ].phi();
        fTJsmallTruthM   [fTNJetsSmallTruthFilled] = myTruthJetsSmall[iJ].m();
        fTNJetsSmallTruthFilled++;
    }


    // small R jets: ATLAS Style ------------------------------------------
    fastjet::AreaDefinition        area_def(fastjet::active_area_explicit_ghosts, fastjet::GhostedAreaSpec(5.));
    fastjet::ClusterSequenceArea   cs(particlesForJets, fastjet::JetDefinition(fastjet::antikt_algorithm, 0.4), area_def);

    fastjet::Selector sel_jets_small    = fastjet::SelectorPtMin(20) * fastjet::SelectorAbsRapMax(2.5);

    vector<fastjet::PseudoJet> myJets_beforeSubtraction = fastjet::sorted_by_pt(cs.inclusive_jets()); //was 10
    vector<fastjet::PseudoJet> myJets                   = sel_jets_small(theSubtractor(myJets_beforeSubtraction));

    // CMS style jvertexjet
    fastjet::contrib::VertexJets    vertexjetCMS(SelectorPileupTrack(), SelectorHardScatterTrack());
    vertexjetCMS.set_tot_n_pu_tracks      (PUtracks.size());
    vertexjetCMS.set_ghost_scale_factor   (1);
    vertexjetCMS.set_corrJVF_scale_factor (0.01);
    vertexjetCMS.set_corrJVF_cut          (-1);


    if(fDebug){ cout << ">>>>> Small R jets <<<<< " << endl;}
    for (unsigned int ij = 0; ij < myJets.size(); ij++){
        if(fTNJetsSmallRFilled == MaxNJetSmallR) {if(fDebug) cout << "Warning: More than " << MaxNJetSmallR << " small R jets" << endl; continue;}
        float truthmatchpt = -1;

        fastjet::PseudoJet  resultjet       = vertexjetCMS(myJets[ij]);

        if(fDebug){
        cout << "myJets[ij] pt " << myJets[ij].pt() 
             << " corrJVF " << resultjet.structure_of<fastjet::contrib::VertexJets>().corrJVF() << endl;
        }

        fTJsmallPt          [fTNJetsSmallRFilled] = myJets[ij].pt();
        fTJsmallEta         [fTNJetsSmallRFilled] = myJets[ij].eta();
        fTJsmallPhi         [fTNJetsSmallRFilled] = myJets[ij].phi();
        fTJsmallM           [fTNJetsSmallRFilled] = myJets[ij].m();
        fTJsmallCharge      [fTNJetsSmallRFilled] = tool->JetCharge(myJets[ij],0.5);
        fTJsmallIsHS        [fTNJetsSmallRFilled] = tool->TruthMatchDR(myJets[ij], myTruthJetsSmall, 0.4, truthmatchpt)? 1:0;
        fTJsmallTMpt        [fTNJetsSmallRFilled] = truthmatchpt;
        fTJsmalldRtruth5    [fTNJetsSmallRFilled] = tool->dRTruth(myJets[ij], myTruthJetsSmall, 5.);
        fTJsmalldRtruth10   [fTNJetsSmallRFilled] = tool->dRTruth(myJets[ij], myTruthJetsSmall, 10.);
        fTJsmallcorrJVF     [fTNJetsSmallRFilled] = resultjet.structure_of<fastjet::contrib::VertexJets>().corrJVF();
//        cout << "small R pt " << fTJsmallPt[fTNJetsSmallRFilled] << " fTJsmallTMpt " << fTJsmallTMpt[fTNJetsSmallRFilled] << " fTJsmallM " << fTJsmallM[fTNJetsSmallRFilled] << endl;

        
        fTNJetsSmallRFilled++;
    }
    */



    // ------------------------------------------------------------------------
    // large-R jets: Trimmed, ATLAS style ---------------------------------
    fastjet::AreaDefinition        area_def(fastjet::active_area_explicit_ghosts, fastjet::GhostedAreaSpec(5.));
    fastjet::ClusterSequenceArea csLargeR(particlesForJets, fastjet::JetDefinition(fastjet::antikt_algorithm, 1.0), area_def);
    vector<fastjet::PseudoJet> myJetsLargeR_beforeSubtraction = fastjet::sorted_by_pt(csLargeR.inclusive_jets());

    

    // VertexJet for largeR jets: corrJVFcut 0.6, trimming 0.5 
    fastjet::contrib::VertexJets    vertexjetLarge(SelectorPileupTrack(), SelectorHardScatterTrack(),fastjet::JetDefinition(fastjet::kt_algorithm, 0.3));
    vertexjetLarge.set_tot_n_pu_tracks      (PUtracks.size());
    vertexjetLarge.set_ghost_scale_factor   (1);
    vertexjetLarge.set_corrJVF_scale_factor (0.01);
    vertexjetLarge.set_corrJVF_cut          (0.6);
    vertexjetLarge.set_trimming_fcut        (0.05);
    vertexjetLarge.set_subtractor           (&theSubtractor);


    // trimmer
    fastjet::Filter trimmer (fastjet::JetDefinition(fastjet::kt_algorithm, 0.3), fastjet::SelectorPtFractionMin(0.05));
    trimmer.set_subtractor(&theSubtractor);

    // Jet loop
    int ij = 0;

    if(myJetsLargeR_beforeSubtraction.size() ==0 ) return;

    fastjet::PseudoJet subtracted = theSubtractor(myJetsLargeR_beforeSubtraction[0]);
    if(! (subtracted.pt() > 300 && fabs(subtracted.eta())<1.5)) return;


    fastjet::PseudoJet tj                        = trimmer(myJetsLargeR_beforeSubtraction[ij]);

    // different vertexjets configs:
    vertexjetLarge.set_corrJVF_cut(0.0); vertexjetLarge.set_trimming_fcut(0.00); fastjet::PseudoJet jet_corrJVFgroomed_00_000 = vertexjetLarge(myJetsLargeR_beforeSubtraction[ij]);
    vertexjetLarge.set_corrJVF_cut(0.0); vertexjetLarge.set_trimming_fcut(0.01); fastjet::PseudoJet jet_corrJVFgroomed_00_001 = vertexjetLarge(myJetsLargeR_beforeSubtraction[ij]);
    vertexjetLarge.set_corrJVF_cut(0.0); vertexjetLarge.set_trimming_fcut(0.02); fastjet::PseudoJet jet_corrJVFgroomed_00_002 = vertexjetLarge(myJetsLargeR_beforeSubtraction[ij]);
    vertexjetLarge.set_corrJVF_cut(0.0); vertexjetLarge.set_trimming_fcut(0.03); fastjet::PseudoJet jet_corrJVFgroomed_00_003 = vertexjetLarge(myJetsLargeR_beforeSubtraction[ij]);
    vertexjetLarge.set_corrJVF_cut(0.0); vertexjetLarge.set_trimming_fcut(0.04); fastjet::PseudoJet jet_corrJVFgroomed_00_004 = vertexjetLarge(myJetsLargeR_beforeSubtraction[ij]);
    vertexjetLarge.set_corrJVF_cut(0.0); vertexjetLarge.set_trimming_fcut(0.05); fastjet::PseudoJet jet_corrJVFgroomed_00_005 = vertexjetLarge(myJetsLargeR_beforeSubtraction[ij]);
    vertexjetLarge.set_corrJVF_cut(0.9); vertexjetLarge.set_trimming_fcut(0.00); fastjet::PseudoJet jet_corrJVFgroomed_06_000 = vertexjetLarge(myJetsLargeR_beforeSubtraction[ij]);
    vertexjetLarge.set_corrJVF_cut(0.9); vertexjetLarge.set_trimming_fcut(0.01); fastjet::PseudoJet jet_corrJVFgroomed_06_001 = vertexjetLarge(myJetsLargeR_beforeSubtraction[ij]);
    vertexjetLarge.set_corrJVF_cut(0.9); vertexjetLarge.set_trimming_fcut(0.02); fastjet::PseudoJet jet_corrJVFgroomed_06_002 = vertexjetLarge(myJetsLargeR_beforeSubtraction[ij]);
    vertexjetLarge.set_corrJVF_cut(0.9); vertexjetLarge.set_trimming_fcut(0.03); fastjet::PseudoJet jet_corrJVFgroomed_06_003 = vertexjetLarge(myJetsLargeR_beforeSubtraction[ij]);
    vertexjetLarge.set_corrJVF_cut(0.9); vertexjetLarge.set_trimming_fcut(0.04); fastjet::PseudoJet jet_corrJVFgroomed_06_004 = vertexjetLarge(myJetsLargeR_beforeSubtraction[ij]);
    vertexjetLarge.set_corrJVF_cut(0.9); vertexjetLarge.set_trimming_fcut(0.05); fastjet::PseudoJet jet_corrJVFgroomed_06_005 = vertexjetLarge(myJetsLargeR_beforeSubtraction[ij]);


    fTJlargeRPt                  = subtracted.pt();
    fTJlargeREta                 = subtracted.eta();
    fTJlargeRPhi                 = subtracted.phi();
    fTJlargeRM                   = subtracted.m();
    fTJlargeRMungroomed          = myJetsLargeR_beforeSubtraction[ij].m();
    fTJlargeRPtungroomed         = myJetsLargeR_beforeSubtraction[ij].pt();
    fTJlargeRCharge              = tool->JetCharge(subtracted,0.5);
    fTJlargeRMtrimmed            = tj.m();
    fTJlargeWplusMatch           = (tool->BosonMatch(subtracted, Bosons, 0.5, 24)      ? 1:0);
    fTJlargeWminusMatch          = (tool->BosonMatch(subtracted, Bosons, 0.5, -24)     ? 1:0);
    fTJlargeZMatch               = (tool->BosonMatch(subtracted, Bosons, 0.5, 23)      ? 1:0);
    fTJlargeRcorrJVFGr_00_000_M  = jet_corrJVFgroomed_00_000.m();
    fTJlargeRcorrJVFGr_00_001_M  = jet_corrJVFgroomed_00_001.m();
    fTJlargeRcorrJVFGr_00_002_M  = jet_corrJVFgroomed_00_002.m();
    fTJlargeRcorrJVFGr_00_003_M  = jet_corrJVFgroomed_00_003.m();
    fTJlargeRcorrJVFGr_00_004_M  = jet_corrJVFgroomed_00_004.m();
    fTJlargeRcorrJVFGr_00_005_M  = jet_corrJVFgroomed_00_005.m();

    fTJlargeRcorrJVFGr_06_000_M  = jet_corrJVFgroomed_06_000.m();
    fTJlargeRcorrJVFGr_06_001_M  = jet_corrJVFgroomed_06_001.m();
    fTJlargeRcorrJVFGr_06_002_M  = jet_corrJVFgroomed_06_002.m();
    fTJlargeRcorrJVFGr_06_003_M  = jet_corrJVFgroomed_06_003.m();
    fTJlargeRcorrJVFGr_06_004_M  = jet_corrJVFgroomed_06_004.m();
    fTJlargeRcorrJVFGr_06_005_M  = jet_corrJVFgroomed_06_005.m();

    vector<fastjet::PseudoJet> subjets = jet_corrJVFgroomed_00_000.pieces();
    for(int iSub=0; iSub<subjets.size(); ++iSub){
        if(fTNSubjetsFilled == MaxNSubjets) { continue;}
        fTSubjetCorrJVF [fTNSubjetsFilled] = subjets[iSub].structure_of<fastjet::contrib::VertexJets>().corrJVF();
        fTSubjetPt      [fTNSubjetsFilled] = subjets[iSub].pt();
        fTNSubjetsFilled++;
    }




    // Fill
    tT->Fill();

    if(fDebug) cout << "FCNCAnalysis::AnalyzeEvent End " << endl;
    return;
}



// declate branches
void FCNCAnalysis::DeclareBranches(){
   
   // Event Properties 
   tT->Branch("EventNumber",               &fTEventNumber,            "EventNumber/I");
   tT->Branch("NPV",                       &fTNPV,                    "NPV/I");

   // smallR jets truth
   tT->Branch("NJetsSmallTruthFilled",     &fTNJetsSmallTruthFilled, "NJetsSmallTruthFilled/I");
   tT->Branch("JsmallTruthPt",             &fTJsmallTruthPt,         "JsmallTruthPt[NJetsSmallTruthFilled]/F");
   tT->Branch("JsmallTruthEta",            &fTJsmallTruthEta,        "JsmallTruthEta[NJetsSmallTruthFilled]/F");
   tT->Branch("JsmallTruthPhi",            &fTJsmallTruthPhi,        "JsmallTruthPhi[NJetsSmallTruthFilled]/F");
   tT->Branch("JsmallTruthM",              &fTJsmallTruthM,          "JsmallTruthM[NJetsSmallTruthFilled]/F");

   // smallR jets
   tT->Branch("NJetsFilledSmallR",         &fTNJetsSmallRFilled,       "NJetsFilledSmallR/I");
   tT->Branch("JsmallPt",                  &fTJsmallPt,                "JsmallPt[NJetsFilledSmallR]/F");
   tT->Branch("JsmallEta",                 &fTJsmallEta,               "JsmallEta[NJetsFilledSmallR]/F");
   tT->Branch("JsmallPhi",                 &fTJsmallPhi,               "JsmallPhi[NJetsFilledSmallR]/F");
   tT->Branch("JsmallM",                   &fTJsmallM,                 "JsmallM[NJetsFilledSmallR]/F");
   tT->Branch("JsmallCharge",              &fTJsmallCharge,            "JsmallCharge[NJetsFilledSmallR]/F");
   tT->Branch("JsmallBtag",                &fTJsmallBtag,              "JsmallBtag[NJetsFilledSmallR]/I");
   tT->Branch("JsmallcorrJVF",             &fTJsmallcorrJVF,           "JsmallcorrJVF[NJetsFilledSmallR]/F");
   tT->Branch("JsmallJVF",                 &fTJsmallJVF,               "JsmallJVF[NJetsFilledSmallR]/F");
   tT->Branch("JsmallIsHS",                &fTJsmallIsHS,              "JsmallIsHS[NJetsFilledSmallR]/I");
   tT->Branch("JsmallTMpt",                &fTJsmallTMpt,              "JsmallTMpt[NJetsFilledSmallR]/F");
   tT->Branch("JsmalldRtruth5",            &fTJsmalldRtruth5,          "JsmalldRtruth5[NJetsFilledSmallR]/F");
   tT->Branch("JsmalldRtruth10",           &fTJsmalldRtruth10,         "JsmalldRtruth10[NJetsFilledSmallR]/F");
   
   // largeR jets
   tT->Branch("NJetsFilledLargeR",         &fTNJetsLargeRFilled,         "NJetsFilledLargeR/I");
   tT->Branch("JlargeRPt",                 &fTJlargeRPt,                 "JlargeRPt/F");
   tT->Branch("JlargeREta",                &fTJlargeREta,                "JlargeREta/F");
   tT->Branch("JlargeRPhi",                &fTJlargeRPhi,                "JlargeRPhi/F");
   tT->Branch("JlargeRM",                  &fTJlargeRM,                  "JlargeRM/F");
   tT->Branch("JlargeRMungroomed",         &fTJlargeRMungroomed,         "JlargeRMungroomed/F");
   tT->Branch("JlargeRPtungroomed",         &fTJlargeRPtungroomed,         "JlargeRPtungroomed/F");
   tT->Branch("JlargeRMtrimmed",           &fTJlargeRMtrimmed,           "JlargeRMtrimmed/F");
   tT->Branch("JlargeRCharge",             &fTJlargeRCharge,             "JlargeRCharge/F");
   tT->Branch("JlargeWplusMatch",          &fTJlargeWplusMatch,          "JlargeWplusMatch/I");
   tT->Branch("JlargeWminusMatch",         &fTJlargeWminusMatch,         "JlargeWminusMatch/I");
   tT->Branch("JlargeZMatch",              &fTJlargeZMatch,              "JlargeZMatch/I");

   tT->Branch("JlargeRcorrJVFGr_00_000_M", &fTJlargeRcorrJVFGr_00_000_M  ,"JlargeRcorrJVFGr_00_000_M/F"); 
   tT->Branch("JlargeRcorrJVFGr_00_001_M", &fTJlargeRcorrJVFGr_00_001_M  ,"JlargeRcorrJVFGr_00_001_M/F");
   tT->Branch("JlargeRcorrJVFGr_00_002_M", &fTJlargeRcorrJVFGr_00_002_M  ,"JlargeRcorrJVFGr_00_002_M/F");
   tT->Branch("JlargeRcorrJVFGr_00_003_M", &fTJlargeRcorrJVFGr_00_003_M  ,"JlargeRcorrJVFGr_00_003_M/F");
   tT->Branch("JlargeRcorrJVFGr_00_004_M", &fTJlargeRcorrJVFGr_00_004_M  ,"JlargeRcorrJVFGr_00_004_M/F");
   tT->Branch("JlargeRcorrJVFGr_00_005_M", &fTJlargeRcorrJVFGr_00_005_M  ,"JlargeRcorrJVFGr_00_005_M/F");
   tT->Branch("JlargeRcorrJVFGr_06_000_M", &fTJlargeRcorrJVFGr_06_000_M  ,"JlargeRcorrJVFGr_06_000_M/F");
   tT->Branch("JlargeRcorrJVFGr_06_001_M", &fTJlargeRcorrJVFGr_06_001_M  ,"JlargeRcorrJVFGr_06_001_M/F");
   tT->Branch("JlargeRcorrJVFGr_06_002_M", &fTJlargeRcorrJVFGr_06_002_M  ,"JlargeRcorrJVFGr_06_002_M/F");
   tT->Branch("JlargeRcorrJVFGr_06_003_M", &fTJlargeRcorrJVFGr_06_003_M  ,"JlargeRcorrJVFGr_06_003_M/F");
   tT->Branch("JlargeRcorrJVFGr_06_004_M", &fTJlargeRcorrJVFGr_06_004_M  ,"JlargeRcorrJVFGr_06_004_M/F");
   tT->Branch("JlargeRcorrJVFGr_06_005_M", &fTJlargeRcorrJVFGr_06_005_M  ,"JlargeRcorrJVFGr_06_005_M/F");
   
   tT->Branch("NSubjetsFilled",&fTNSubjetsFilled,"NSubjetsFilled/I");
   tT->Branch("SubJetPt",      &fTSubjetPt,     "SubJetPt[NSubjetsFilled]/F");
   tT->Branch("SubJetCorrJVF", &fTSubjetCorrJVF,"SubJetCorrJVF[NSubjetsFilled]/F");

   // Bosons
   tT->Branch("NBosonsFilled",             &fTNBosonsFilled,          "NBosonsFilled/I");
   tT->Branch("BosonID",                   &fTBosonID,                "BosonID[NBosonsFilled]/I");
   tT->Branch("BosonPt",                   &fTBosonPt,                "BosonPt[NBosonsFilled]/F");
   tT->Branch("BosonEta",                  &fTBosonEta,               "BosonEta[NBosonsFilled]/F");
   tT->Branch("BosonPhi",                  &fTBosonPhi,               "BosonPhi[NBosonsFilled]/F");
   tT->Branch("BosonM",                    &fTBosonM,                 "BosonM[NBosonsFilled]/F");


   tT->GetListOfBranches()->ls();
    
   return;
}


// resets vars
void FCNCAnalysis::ResetBranches(){
      // reset branches 
      fTEventNumber                 = -999;
      fTNPV                         = -999;
      fTNJetsSmallRFilled           = 0;
      fTNSmallVJFilled              = 0;
      fTNJetsLargeRFilled           = 0;
      fTNBosonsFilled               = 0;
      fTNJetsSmallTruthFilled       = 0;
      fTNSubjetsFilled              = 0;

      for(int ij=0; ij<MaxNSubjets; ++ij){
        fTSubjetPt      [ij] = -99;
        fTSubjetCorrJVF [ij] = -99;
      }

      for(int ij=0; ij<MaxNJetLargeR; ++ij){
        fTJsmallTruthPt    [ij] =-999;
        fTJsmallTruthEta   [ij] =-999;
        fTJsmallTruthPhi   [ij] =-999;
        fTJsmallTruthM     [ij] =-999;
      }

      for (int iP=0; iP < MaxNJetSmallR; ++iP){
          fTJsmallPt       [iP]= -999;
          fTJsmallPhi      [iP]= -999;
          fTJsmallEta      [iP]= -999;
          fTJsmallM        [iP]= -999;
          fTJsmallCharge   [iP]= -999;
          fTJsmallBtag     [iP]= -999;
          fTJsmallcorrJVF  [iP]= -999;
          fTJsmallJVF      [iP]= -999;
          fTJsmallIsHS     [iP]= -999;
          fTJsmallTMpt     [iP]= -999;
          fTJsmalldRtruth5 [iP]= -999;
          fTJsmalldRtruth10[iP]= -999;
      }

      fTJlargeRPt                 = -999;
      fTJlargeRPhi                = -999;
      fTJlargeREta                = -999;
      fTJlargeRM                  = -999;
      fTJlargeRMtrimmed           = -999;
      fTJlargeRMungroomed         = -999;
      fTJlargeRPtungroomed         = -999;
      fTJlargeRCharge             = -999;
      fTJlargeWplusMatch          = -999;
      fTJlargeWminusMatch         = -999;
      fTJlargeZMatch              = -999;

    fTJlargeRcorrJVFGr_00_000_M = -999;  
    fTJlargeRcorrJVFGr_00_001_M = -999;
    fTJlargeRcorrJVFGr_00_002_M = -999;
    fTJlargeRcorrJVFGr_00_003_M = -999;
    fTJlargeRcorrJVFGr_00_004_M = -999;
    fTJlargeRcorrJVFGr_00_005_M = -999;
    fTJlargeRcorrJVFGr_06_000_M = -999;
    fTJlargeRcorrJVFGr_06_001_M = -999;
    fTJlargeRcorrJVFGr_06_002_M = -999;
    fTJlargeRcorrJVFGr_06_003_M = -999;
    fTJlargeRcorrJVFGr_06_004_M = -999;
    fTJlargeRcorrJVFGr_06_005_M = -999;




      for (int iP=0; iP < MaxNBosons; ++iP){
          fTBosonPt  [iP]= -999;
          fTBosonEta [iP]= -999;
          fTBosonPhi [iP]= -999;
          fTBosonM   [iP]= -999;
          fTBosonID  [iP]= -999;
      }
}

// ------------------------------- 
// Ghost Matching
vector<fastjet::PseudoJet> FCNCAnalysis::AddGhosts(const vector<fastjet::PseudoJet> &jet_constits, const vector<fastjet::PseudoJet> & hs_tracks, const vector<fastjet::PseudoJet> & pu_tracks){

    std::vector<fastjet::PseudoJet> clustersandtracks(jet_constits);

    // first HS tracks
    for (int iT = 0; iT < hs_tracks.size(); ++iT){
        fastjet::PseudoJet ghost = kGhostScaleFact*hs_tracks[iT];
        ghost.set_user_info(new GhostUserInfo(true, false, true));
        clustersandtracks.push_back(ghost);
    }
    // then pileup tracks
    for (int iT = 0; iT < pu_tracks.size(); ++iT){
        fastjet::PseudoJet ghost = kGhostScaleFact*pu_tracks[iT];
        ghost.set_user_info(new GhostUserInfo(true, true, true));
        clustersandtracks.push_back(ghost);
    }
    
    return clustersandtracks;
}


