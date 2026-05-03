#pragma once

// Standard library includes
#include <map>
#include <string>

// Enum used to label particles of interest, to make particle dependent plots (aiming to make BDT plotting)
enum EventCategoryNp1pi_part {

  // Unable to categorize (e.g., because the event is real data and thus
  // has no MC truth information)
  kParticleUnknown = 0,
  kMuon = 1, //
  kPion = 2, //
  kProton = 3, //
  kOther = 4 //

};

extern std::map< int, std::pair< std::string, int > > CC1muNp1pi_MAPart;
