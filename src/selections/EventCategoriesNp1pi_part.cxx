#include "TH1.h"
#include "TColor.h"
#include "XSecAnalyzer/Selections/EventCategoriesNp1pi_part.hh"

std::map<int, std::pair<std::string, int>> CC1muNp1pi_MAPart = {

  { kMuon,        { "Muon", TColor::GetColor("#8B0000") } }, // darkred
  { kPion,     { "Pion", TColor::GetColor("#DAA520") } }, // goldenrod
  { kProton,    { "Proton", TColor::GetColor("#B8860B") } }, // darkgoldenrod
  { kOther,    { "Other", TColor::GetColor("#F0E68C") } }, // khaki
  
};
