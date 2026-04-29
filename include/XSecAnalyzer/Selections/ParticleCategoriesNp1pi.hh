#pragma once

// Standard library includes
#include <map>
#include <string>

// Enum used to label particles of interest, to make particle dependent plots (aiming to make BDT plotting)
enum ParticleCategoryNp1pi {

  // Unable to categorize (e.g., because the event is real data and thus
  // has no MC truth information)
  kParticleUnknown = -1,
  kMuon = 0, //
  kPion = 1, //
  kProton = 2, //
  kOther = 3 //

};

extern std::map< int, std::pair< std::string, int > > CC1muNp1pi_PARTICLE_MAP;
