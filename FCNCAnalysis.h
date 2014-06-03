#ifndef  FCNCAnalysis_H
#define  FCNCAnalysis_H

#include <vector>
#include <math.h>
#include <string>

#include "fastjet/ClusterSequence.hh"
#include "fastjet/PseudoJet.hh"  
#include "fastjet/tools/Filter.hh"

#include "TFile.h"
#include "TTree.h"
#include "TClonesArray.h"
#include "TParticle.h"

#include "FCNCTools.h"
#include "myFastJetBase.h"
#include "Pythia8/Pythia.h"

using namespace std;

class FCNCAnalysis{
    public:
        FCNCAnalysis ();
        ~FCNCAnalysis ();
        
        void Begin();
        void AnalyzeEvent(int iEvt, Pythia8::Pythia *pythia8, Pythia8::Pythia* pythia8PU, int nPU);
        void End();
        void Debug(int debug){
            fDebug = debug;
        }
        void SetOutName(string outname){
            fOutName = outname;
        }

        enum event_type {
            hardscatter,
            pileup
        } event_type_;

    private:
        int  ftest;
        int  fDebug;
        int  nPU_;
        string fOutName;

        TFile *tF;
        TTree *tT;
        FCNCTools *tool;


        static const float kGhostScaleFact = 1E-20;

        // Tree Vars ---------------------------------------
        int              fTEventNumber;
        int              fTNPV;

        

        static const int MaxNJetSmallR = 20;
    
        int              fTNJetsSmallTruthFilled;
        float            fTJsmallTruthPt   [MaxNJetSmallR];
        float            fTJsmallTruthEta  [MaxNJetSmallR];
        float            fTJsmallTruthPhi  [MaxNJetSmallR];
        float            fTJsmallTruthM    [MaxNJetSmallR];


        int              fTNJetsSmallRFilled;
        float            fTJsmallPt        [MaxNJetSmallR];
        float            fTJsmallEta       [MaxNJetSmallR];
        float            fTJsmallPhi       [MaxNJetSmallR];
        float            fTJsmallM         [MaxNJetSmallR];
        float            fTJsmallCharge    [MaxNJetSmallR];
        float            fTJsmallcorrJVF   [MaxNJetSmallR];
        float            fTJsmallJVF       [MaxNJetSmallR];
        int              fTJsmallBtag      [MaxNJetSmallR];
        int              fTJsmallIsHS      [MaxNJetSmallR];
        float            fTJsmallTMpt      [MaxNJetSmallR];
        float            fTJsmalldRtruth5  [MaxNJetSmallR];
        float            fTJsmalldRtruth10 [MaxNJetSmallR];
        
        int              fTNSmallVJFilled;
        float            fTVJsmallPt        [MaxNJetSmallR];
        float            fTVJsmallEta       [MaxNJetSmallR];
        float            fTVJsmallPhi       [MaxNJetSmallR];
        float            fTVJsmallM         [MaxNJetSmallR];
        float            fTVJsmallcorrJVF   [MaxNJetSmallR];
        float            fTVJsmallJVF       [MaxNJetSmallR];
        int              fTVJsmallIsHS      [MaxNJetSmallR];
        float            fTVJsmallTMpt      [MaxNJetSmallR];
        float            fTVJsmalldRtruth5  [MaxNJetSmallR];
        float            fTVJsmalldRtruth10 [MaxNJetSmallR];

        static const int MaxNJetLargeR = 20;
        int              fTNJetsLargeRFilled;
        float            fTJlargeRPt                ;
        float            fTJlargeREta               ;
        float            fTJlargeRPhi               ;
        float            fTJlargeRM                 ;
        float            fTJlargeRMungroomed        ;
        float            fTJlargeRPtungroomed        ;
        float            fTJlargeRMtrimmed          ;
        float            fTJlargeRCharge            ;
        int              fTJlargeRBtag              ;
        int              fTJlargeWplusMatch         ;
        int              fTJlargeWminusMatch        ;
        int              fTJlargeZMatch             ;

        static const int MaxNSubjets = 20;
        int              fTNSubjetsFilled;
        float            fTSubjetCorrJVF[MaxNSubjets];
        float            fTSubjetPt     [MaxNSubjets];

        float fTJlargeRcorrJVFGr_00_000_M; 
        float fTJlargeRcorrJVFGr_00_001_M;
        float fTJlargeRcorrJVFGr_00_002_M;
        float fTJlargeRcorrJVFGr_00_003_M;
        float fTJlargeRcorrJVFGr_00_004_M;
        float fTJlargeRcorrJVFGr_00_005_M;
        float fTJlargeRcorrJVFGr_06_000_M;
        float fTJlargeRcorrJVFGr_06_001_M;
        float fTJlargeRcorrJVFGr_06_002_M;
        float fTJlargeRcorrJVFGr_06_003_M;
        float fTJlargeRcorrJVFGr_06_004_M;
        float fTJlargeRcorrJVFGr_06_005_M;
        // Bosons
        static const int MaxNBosons    = 5;
        int              fTNBosonsFilled;
        float            fTBosonPt [MaxNBosons];
        float            fTBosonEta[MaxNBosons];
        float            fTBosonPhi[MaxNBosons];
        float            fTBosonM  [MaxNBosons];
        int              fTBosonID [MaxNBosons];

        fastjet::JetDefinition     *m_jet_def;
        fastjet::JetDefinition     *m_jet_def_largeR_ALTAS;

    
        std::vector <fastjet::PseudoJet>           particlesForJets;
        std::vector <fastjet::PseudoJet>           particlesForTruthJets;
        std::vector <fastjet::PseudoJet>           Bosons;
        std::vector <fastjet::PseudoJet>           PUtracks;
        std::vector <fastjet::PseudoJet>           HStracks;


        std::vector<fastjet::PseudoJet> AddGhosts(const std::vector<fastjet::PseudoJet> &jet_constits, const std::vector<fastjet::PseudoJet> & hs_tracks, const std::vector<fastjet::PseudoJet> & pu_tracks);
        void                   DeclareBranches();
        void                   ResetBranches();
        bool                   GetPythiaParticles(Pythia8::Pythia* pythia8, event_type etype, int nPU=0);

        
};



#endif

