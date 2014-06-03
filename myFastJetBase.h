#ifndef MYFASTJETBASE_H
#define MYFASTJETBASE_H

#include "fastjet/PseudoJet.hh"
#include "fastjet/Selector.hh"


class MyUserInfo : public fastjet::PseudoJet::UserInfoBase{
 public:
 MyUserInfo(const int & pdg_id_in,const int & pythia_id_in,  const double & charge_in, const int & origin_in) :
  _pdg_id(pdg_id_in),_pythia_id(pythia_id_in), _charge(charge_in), _origin(origin_in){}
  int pdg_id() const { return _pdg_id;}
  int pythia_id() const {return _pythia_id;}
  double charge() const { return _charge;}
  int origin() const { return _origin;}
 protected:
  int _pdg_id;         // the associated pdg id
  int _pythia_id;  // index in pythia.event
  double _charge;  // the particle charge
  int _origin; // origin of particle: -1 = nowhere, 0 = HS, >=1 PU interaction
};

class GhostUserInfo : public fastjet::PseudoJet::UserInfoBase{
 public:
    GhostUserInfo(const bool & is_ghost, const bool &is_pileup, const bool &is_track): _is_ghost(is_ghost), _is_pileup(is_pileup), _is_track(is_track){}
    bool is_ghost()   const { return _is_ghost;}
    bool is_pileup()  const { return _is_pileup;}
    bool is_track()   const { return _is_track;}
 protected:
    bool _is_ghost;  // true is ghost particle
    bool _is_pileup; // true is ghost particle
    bool _is_track;  // true is track 
};

// Selector
class SelectorWorkerPileupGhostTrack : public fastjet::SelectorWorker {
public:

  virtual bool pass(const fastjet::PseudoJet & particle) const {
  // we check that the user_info_ptr is non-zero so as to make
  // sure that explicit ghosts don't cause the selector to fail
    return (particle.has_user_info<GhostUserInfo>()  
            &&  particle.user_info<GhostUserInfo>().is_pileup() == true 
            &&  particle.user_info<GhostUserInfo>().is_ghost()  == true );
  }
  
  virtual string description() const {return "comes from pileup interaction";}
};

// Selector
class SelectorWorkerHardScatterGhostTrack : public fastjet::SelectorWorker {
public:

  virtual bool pass(const fastjet::PseudoJet & particle) const {
  // we check that the user_info_ptr is non-zero so as to make
  // sure that explicit ghosts don't cause the selector to fail
    return (particle.has_user_info<GhostUserInfo>()  
            &&  particle.user_info<GhostUserInfo>().is_pileup() == false
            &&  particle.user_info<GhostUserInfo>().is_ghost()  == true );
  }
  
  virtual string description() const {return "comes from pileup interaction";}
};

// Selector
class SelectorWorkerPileupTrack : public fastjet::SelectorWorker {
public:

  virtual bool pass(const fastjet::PseudoJet & particle) const {
  // we check that the user_info_ptr is non-zero so as to make
  // sure that explicit ghosts don't cause the selector to fail
    return (particle.has_user_info<GhostUserInfo>()  
            &&  particle.user_info<GhostUserInfo>().is_track()  == true 
            &&  particle.user_info<GhostUserInfo>().is_pileup() == true 
            &&  particle.user_info<GhostUserInfo>().is_ghost()  == false );
  }
  
  virtual string description() const {return "comes from pileup interaction";}
};

// Selector
class SelectorWorkerHardScatterTrack : public fastjet::SelectorWorker {
public:

  virtual bool pass(const fastjet::PseudoJet & particle) const {
  // we check that the user_info_ptr is non-zero so as to make
  // sure that explicit ghosts don't cause the selector to fail
    return (particle.has_user_info<GhostUserInfo>()  
            &&  particle.user_info<GhostUserInfo>().is_track()  == true
            &&  particle.user_info<GhostUserInfo>().is_pileup() == false
            &&  particle.user_info<GhostUserInfo>().is_ghost()  == false );
  }
  
  virtual string description() const {return "comes from pileup interaction";}
};

#endif
