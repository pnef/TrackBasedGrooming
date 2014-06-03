#include <iostream>
#include <sstream>
#include <fstream>
#include <iomanip>
#include <math.h>
#include <vector>
#include <stdlib.h>
#include <stdio.h>

#include "TString.h"
#include "TSystem.h"
#include "TError.h"
//#include "TPythia8.h"
#include "TClonesArray.h"
#include "TParticle.h"
#include "TDatabasePDG.h"

#include "fastjet/PseudoJet.hh"  
#include "fastjet/ClusterSequence.hh"
#include "fastjet/Selector.hh"

#include "Pythia8/Pythia.h"

#include "FCNCTools.h"
#include "FCNCAnalysis.h"

#include "boost/program_options.hpp"

#include "TRandom3.h"

using std::cout;
using std::endl;
using std::string;
using std::map;
using namespace std;
namespace po = boost::program_options;

int getSeed(int optSeed){
    if (optSeed > -1) return optSeed;
    int timeSeed = time(NULL);
    return abs(((timeSeed*181)*((getpid()-83)*359))%104729);
}

int main(int argc, char* argv[]){
    // argument parsing  ------------------------
    cout << "Called as: ";
    for(int ii=0; ii<argc; ++ii){
        cout << argv[ii] << " ";
    }
    cout << endl;

    // agruments 
    int nEvents        = 0;
    int fDebug         = 0;
    int nPU            = 0;
    int nFlatPUMax     = 0;
    int seed           = -1;
    string infile = "/u/at/pnef/Work/Code/MadGraph/MG5_aMC_v2_0_1/HP_tt_20140314/Events/run_01/unweighted_events.lhe";
    string outName = "TopToHc.root";

    po::options_description desc("Allowed options");
    desc.add_options()
      ("help", "produce help message")
      ("NEvents",        po::value<int>(&nEvents)->default_value(10) ,    "Number of Events ")
      ("Debug",          po::value<int>(&fDebug) ->default_value(0) ,     "Debug flag")
      ("InFlie",         po::value<string>(&infile) ->default_value("/u/at/pnef/Work/Code/MadGraph/MG5_aMC_v2_0_1/20140421_SM_ttbar_allhadronic/Events/run_01/unweighted_events.lhe") ,     "input file")
      ("OutFile",        po::value<string>(&outName)->default_value("test.root"), "output file name")
      ("nPU",            po::value<int>(&nPU) ->default_value(0) ,     "n PU interactions")
      ("nFlatPUMax",     po::value<int>(&nFlatPUMax) ->default_value(0) ,     "generate flat PU distribution up to NPV = nFlatPUMax")
      ("seed",           po::value<int>(&seed) ->default_value(-1) ,     "random seed")
      ;
    po::variables_map vm;
    po::store(po::parse_command_line(argc, argv, desc), vm);
    po::notify(vm);

    if (vm.count("help")>0){
        cout << desc << "\n";
        return 1;
    }
    //------


   seed = getSeed(seed);
   // Configure and initialize pythia: LHE mode
   /*
   Pythia8::Pythia* pythia8 = new Pythia8::Pythia();
   pythia8->init(infile);
   pythia8->readString("24:onMode = off");
   pythia8->readString("24:onIfAny = 1 2 3 4");
   pythia8->readString("32:onMode = off");
   pythia8->readString("32:onIfAny = 1 2 3 4");
   */

   // Configure pythia for W'->WZ
   Pythia8::Pythia* pythia8 = new Pythia8::Pythia();
   pythia8->readString("Random:setSeed = on");
   std::stringstream ss; ss << "Random:seed = " << seed;
   cout << ss.str() << endl;
   pythia8->readString(ss.str());
   pythia8->readString("32:m0=1000");
   pythia8->readString("NewGaugeBoson:ffbar2gmZZprime= on");
   pythia8->readString("Zprime:gmZmode=3");
   pythia8->readString("32:onMode = off");
   pythia8->readString("32:onIfAny = 6");
   pythia8->readString("24:onMode = off");
   pythia8->readString("24:onIfAny = 1 2 3 4");
   pythia8->init(2212 /* p */, 2212 /* p */, 14000. /* TeV */); //this has to be the last line!



   // configure pythia for PU generation
   Pythia8::Pythia* pythia_MB = new Pythia8::Pythia();
   pythia_MB->readString("Random:setSeed = on");
   ss.clear(); ss.str(""); ss << "Random:seed = " << seed+1;
   cout << ss.str() << endl;
   pythia_MB->readString(ss.str());
   pythia_MB->readString("SoftQCD:nonDiffractive = on");
   pythia_MB->readString("HardQCD:all = off");
   pythia_MB->readString("PhaseSpace:pTHatMin  = .1");
   pythia_MB->readString("PhaseSpace:pTHatMax  = 20000");
   pythia_MB->init(2212,2212,14000);
   
   // FCNCAnalysis
   FCNCAnalysis * analysis = new FCNCAnalysis();
   analysis->SetOutName(outName);
   analysis->Begin();
   analysis->Debug(fDebug);


   // Event loop
   TRandom3 *rand = new TRandom3();
   for (Int_t iev = 0; iev < nEvents; iev++) {
      if (iev%10==0) std::cout << iev << std::endl;
      int nPUevents =0;
      if(nFlatPUMax>1) {
       nPUevents = (int) rand->Uniform(nFlatPUMax);
      }else if(nPU>0){
       nPUevents = nPU;
      }
      analysis->AnalyzeEvent(iev, pythia8, pythia_MB, nPUevents);
   }

   analysis->End();

   // that was it
   delete pythia8;
   delete analysis;

   return 0;
}
