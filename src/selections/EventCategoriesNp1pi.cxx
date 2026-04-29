#include "TH1.h"
#include "TColor.h"
#include "XSecAnalyzer/Selections/EventCategoriesNp1pi.hh"

std::map<int, std::pair<std::string, int>> CC1muNp1pi_MAP = {

  { kUnknown,        { "Unknown", TColor::GetColor("#8B0000") } }, // darkred
  { kSignalCCQE,     { "Signal (CCQE)", TColor::GetColor("#DAA520") } }, // goldenrod
  { kSignalCCMEC,    { "Signal (CCMEC)", TColor::GetColor("#B8860B") } }, // darkgoldenrod
  { kSignalCCRES,    { "Signal (CCRES)", TColor::GetColor("#F0E68C") } }, // khaki
  { kSignalCCDIS,    { "Signal (CCDIS)", TColor::GetColor("#FFD700") } }, // gold
  { kSignalCCCOH,    { "Signal (CCCOH)", TColor::GetColor("#FF0000") } }, // red
  { kSignalOther,    { "Signal (Other)", TColor::GetColor("#FFA500") } }, // orange

  { kNuMuCCpi0,      { "#nu_{#mu} CCN#pi^{0}", TColor::GetColor("#483D8B") } }, // darkslateblue
  { kNuMuCC0piXp,    { "#nu_{#mu} CC0#piXp", TColor::GetColor("#B0C4DE") } },   // lightsteelblue
  { kNuMuCCOther,    { "Other #nu_{#mu} CC", TColor::GetColor("#7B68EE") } },   // mediumslateblue

  { kNuECC,          { "#nu_{e} CC", TColor::GetColor("#00008B") } }, // darkblue
  { kNC,             { "NC", TColor::GetColor("#FFA07A") } },         // lightsalmon
  { kOOFV,           { "OOFV", TColor::GetColor("#FF6347") } },        // tomato
  { kOther,          { "Other", TColor::GetColor("#8B0000") } }        // darkred
};
