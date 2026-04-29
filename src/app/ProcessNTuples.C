// Post-processing program for the MicroBooNE xsec_analyzer framework. This is
// currently designed for use with the PeLEE group's "searchingfornues" ntuples
//
// Updated 24 September 2024
// Steven Gardiner <gardiner@fnal.gov>
// Daniel Barrow <daniel.barrow@physics.ox.ac.uk>

// Updated 19/02/2026: Modified for adding weighting for FSI effects (Temporary solution) - M.L. Velazquez Fernandez

// Standard library includes
#include <cmath>
#include <iostream>
#include <map>
#include <memory>
#include <sstream>
#include <string>
#include <vector>

// ROOT includes
#include "TChain.h"
#include "TFile.h"
#include "TBranch.h"
#include "TParameter.h"
#include "TTree.h"
#include "TVector3.h"
#include "TH2D.h"

// XSecAnalyzer includes
#include "XSecAnalyzer/AnalysisEvent.hh"
#include "XSecAnalyzer/Branches.hh"
#include "XSecAnalyzer/Constants.hh"
#include "XSecAnalyzer/Functions.hh"
#include "XSecAnalyzer/Selections/SelectionBase.hh"
#include "XSecAnalyzer/Selections/SelectionFactory.hh"

bool ends_with(const std::string& str, const std::string& suffix) {
        if (str.length() < suffix.length()) return false;
        return str.compare(str.length() - suffix.length(), suffix.length(), suffix) == 0;
      }
      
void analyze( const std::vector< std::string >& in_file_names,
  const std::vector< std::string >& selection_names,
  const std::string& output_filename )
{
  std::cout << "\nRunning ProcessNTuples with options:\n";
  std::cout << "\toutput_filename: " << output_filename << '\n';
  std::cout << "\tinput_file_names:\n";
  for ( size_t i = 0u; i < in_file_names.size(); ++i ) {
    std::cout << "\t\t- " << in_file_names[i] << '\n';
  }
  std::cout << "\n\nselection names:\n";
  for ( const auto& sel_name : selection_names ) {
    std::cout << "\t\t- " << sel_name << '\n';
  }

  // Toggle FSI reweighting
  bool use_fsi_weighter = true;
  //bool force_reco_pion_truth = true;   // TOY STUDY

  TChain events_ch( "NeutrinoSelectionFilter" );
  TChain subruns_ch( "SubRun" );

  for ( const auto& f_name : in_file_names ) {
    events_ch.Add( f_name.c_str() );
    subruns_ch.Add( f_name.c_str() );
  }

  // OUTPUT TTREE
  // Make an output TTree for plotting (one entry per event)
  TFile* out_file = new TFile( output_filename.c_str(), "recreate" );
  out_file->cd();
  TTree* out_tree = new TTree( "stv_tree", "STV analysis tree" );

  // Get the total POT from the subruns TTree. Save it in the output
  // TFile as a TParameter<float>. Real data doesn't have this TTree,
  // so check that it exists first.
  float pot;
  float summed_pot = 0.;
  bool has_pot_branch = ( subruns_ch.GetBranch("pot") != nullptr );
  if ( has_pot_branch ) {
    subruns_ch.SetBranchAddress( "pot", &pot );
    for ( int se = 0; se < subruns_ch.GetEntries(); ++se ) {
      subruns_ch.GetEntry( se );
      summed_pot += pot;
    }
  }

  TParameter<float>* summed_pot_param =
    new TParameter<float>( "summed_pot", summed_pot );
  summed_pot_param->Write();

  // Load FSI ratio histograms
  TH2D* h_ratio_211  = nullptr;
  TH2D* h_ratio_m211 = nullptr;
  TH2D* h_ratio_111  = nullptr;

  if (use_fsi_weighter) {

    TFile f211("ratio_weights_211.root","READ");
    TFile fm211("ratio_weights_neg211.root","READ");
    TFile f111("ratio_weights_111.root","READ");

    if (f211.IsZombie() || fm211.IsZombie() || f111.IsZombie()) {
      std::cout << "Error opening ratio files\n";
      exit(1);
    }

    h_ratio_211  = (TH2D*)f211.Get("h_ratio_211");
    h_ratio_m211 = (TH2D*)fm211.Get("h_ratio_neg211");
    h_ratio_111  = (TH2D*)f111.Get("h_ratio_111");

    if (!h_ratio_211 || !h_ratio_m211 || !h_ratio_111) {
      std::cout << "Error retrieving ratio histograms\n";
      exit(1);
    }

    h_ratio_211->SetDirectory(0);
    h_ratio_m211->SetDirectory(0);
    h_ratio_111->SetDirectory(0);

    f211.Close();
    fm211.Close();
    f111.Close();
  }

  std::vector< std::unique_ptr<SelectionBase> > selections;

  SelectionFactory sf;
  for ( const auto& sel_name : selection_names ) {
    selections.emplace_back().reset( sf.CreateSelection(sel_name) );
  }

  out_file->cd();
  for ( auto& sel : selections ) {
    sel->setup( out_tree );
  }

  // EVENT LOOP
  // TChains can potentially be really big (and spread out over multiple
  // files). When that's the case, calling TChain::GetEntries() can be very
  // slow. I get around this by using a while loop instead of a for loop.
  bool created_output_branches = false;
  long events_entry = 0;

  while ( true ) {

    //if ( events_entry > 1000) break;

    if ( events_entry % 100000 == 0 ) {
      std::cout << "Processing event #" << events_entry << '\n';
    }

    //if (events_entry > 10000) break;

    // Create a new AnalysisEvent object. This will reset all analysis
    // variables for the current event.
    AnalysisEvent cur_event;

    // Set branch addresses for the member variables that will be read
    // directly from the Event TTree.
    set_event_branch_addresses( events_ch, cur_event );

    // TChain::LoadTree() returns the entry number that should be used with
    // the current TTree object, which (together with the TBranch objects
    // that it owns) doesn't know about the other TTrees in the TChain.
    // If the return value is negative, there was an I/O error, or we've
    // attempted to read past the end of the TChain.
    int local_entry = events_ch.LoadTree( events_entry );

    // If we've reached the end of the TChain (or encountered an I/O error),
    // then terminate the event loop
    if ( local_entry < 0 ) break;

    // Load all of the branches for which we've called
    // TChain::SetBranchAddress() above
    events_ch.GetEntry( events_entry );

    // ------------------------------------------------------------
    // APPLY FSI REWEIGHTING BEFORE OUTPUT BRANCHES ARE SET
    // ------------------------------------------------------------
    // Apply same mask as used for obtaining the weights.
    struct TopologyMask {
      bool require_cc;
      int min_protons;
      int n_pi_plus;
      int n_pi_minus;
      int n_pi0;
      int pion_pdg;     
      TH2D* hist;        
      std::string name;  
    };
    std::vector<TopologyMask> fsi_masks;

    if (use_fsi_weighter) {

      // CC1pi0 + >=1p
      fsi_masks.push_back({
        true,     
        1,        
        0,        
        0,       
        1,        
        111,     
        h_ratio_111,
        "CC1pi0_1p"
      });

      // CC1pi+ + >=1p
      fsi_masks.push_back({
        true,
        1,
        1,
        0,
        0,
        211,
        h_ratio_211,
        "CC1piPlus_1p"
      });

      // CC1pi- + >=1p
      fsi_masks.push_back({
        true,
        1,
        0,
        1,
        0,
        -211,
        h_ratio_m211,
        "CC1piMinus_1p"
      });
    }

    if (use_fsi_weighter && cur_event.mc_nu_daughter_pdg_) {

      int n_proton = 0;
      int n_pi_plus = 0;
      int n_pi_minus = 0;
      int n_pi0 = 0;

      int pion_index_for_mask = -1;

      bool is_cc = (cur_event.mc_nu_ccnc_ == 0);

      size_t n_daughters = cur_event.mc_nu_daughter_pdg_->size();

      // Count final-state particles
      for (size_t d = 0; d < n_daughters; ++d) {

        int pdg = cur_event.mc_nu_daughter_pdg_->at(d);

        if (pdg == 2212) ++n_proton;
        if (pdg == 211)  ++n_pi_plus;
        if (pdg == -211) ++n_pi_minus;
        if (pdg == 111)  ++n_pi0;
      }

      double fsi_weight = 1.0;
      int n_masks_matched = 0;

      for (auto& mask : fsi_masks) {

        bool match = true;

        if (mask.require_cc && !is_cc) match = false;
        if (n_proton < mask.min_protons) match = false;
        if (n_pi_plus  != mask.n_pi_plus)  match = false;
        if (n_pi_minus != mask.n_pi_minus) match = false;
        if (n_pi0      != mask.n_pi0)      match = false;

        if (!match) continue;

        ++n_masks_matched;

        // Find pion for kinematics
        for (size_t d = 0; d < n_daughters; ++d) {
          if (cur_event.mc_nu_daughter_pdg_->at(d) == mask.pion_pdg) {
            pion_index_for_mask = d;
            break;
          }
        }

        if (pion_index_for_mask < 0) continue;

        double px = cur_event.mc_nu_daughter_px_->at(pion_index_for_mask);
        double py = cur_event.mc_nu_daughter_py_->at(pion_index_for_mask);
        double pz = cur_event.mc_nu_daughter_pz_->at(pion_index_for_mask);
        double p = std::sqrt(px*px + py*py + pz*pz);
        if (p <= 0) continue;

        const double m_pi_ch = 0.13957;
        const double m_pi0   = 0.13498;
        double mass = (mask.pion_pdg == 111) ? m_pi0 : m_pi_ch;
        double E = std::sqrt(p*p + mass*mass);
        double KE = (E - mass) * 1000.0;
        double cos_theta = pz / p;

        TH2D* h = mask.hist;
        if (!h) continue;

        int binx = h->GetXaxis()->FindBin(KE);
        binx = std::max(1, std::min(binx, h->GetXaxis()->GetNbins()));
        int biny = h->GetYaxis()->FindBin(cos_theta);
        biny = std::max(1, std::min(biny, h->GetYaxis()->GetNbins()));

        double w = h->GetBinContent(binx, biny);

        if (std::isfinite(w) && w > 0.0)
          fsi_weight = w;
      }

      // Apply weight
      // if (fsi_weight != 1.0 && cur_event.mc_weights_map_) {
      //   // std::cout << "Applying FSI weight " << fsi_weight
      //             // << " to event " << events_entry << "\n";
      //   auto it = cur_event.mc_weights_map_->find("TunedCentralValue_UBGenie");
      //   //print if found
      //   std::cout << "DEBUG: Looking for 'TunedCentralValue_UBGenie' in mc_weights_map_\n";
      //   std::cout << "DEBUG: mc_weights_map_ keys:\n";
      //   std::vector<std::string> keys;
      //   for (const auto& kv : *cur_event.mc_weights_map_) {
      //     keys.push_back(kv.first);
      //   }
      //   for (const auto& key : keys) {
      //     std::cout << "  - " << key << '\n';
      //   }

      //   if (it != cur_event.mc_weights_map_->end()) {
      //     it->second.at(0) *= fsi_weight;
      //   }
      // }
      
      if (fsi_weight != 1.0 && cur_event.mc_weights_map_) {
        for (auto& [wgt_name, wgt_vec] : *cur_event.mc_weights_map_) {
          if (ends_with(wgt_name, "UBGenie")) {
            for (double& w : wgt_vec) {
              w *= fsi_weight;
            }
          }
        }
      }
    }
    // ------------------------------------------------------------

    bool create_them = false;
    if ( !created_output_branches ) {
      create_them = true;
      created_output_branches = true;
    }
    set_event_output_branch_addresses(*out_tree, cur_event, create_them );

    for ( auto& sel : selections ) {
      sel->apply_selection( &cur_event );
    }

    // We're done. Save the results and move on to the next event.
    out_tree->Fill();
    ++events_entry;
  }

  for ( auto& sel : selections ) {
    sel->summary();
  }
  std::cout << "Wrote output to:" << output_filename << std::endl;

  for ( auto& sel : selections ) {
    sel->final_tasks();
  }

  std::cout << "DEBUG 0" << std::endl;
  out_tree->Write();
  std::cout << "DEBUG 1" << std::endl;
  out_file->Close();
  std::cout << "DEBUG 2" << std::endl;
  delete out_file;
  std::cout << "DEBUG 3" << std::endl;
}

void analyzer( const std::string& in_file_name,
 const std::vector< std::string > selection_names,
 const std::string& output_filename)
{
  std::vector< std::string > in_files = { in_file_name };
  analyze( in_files, selection_names, output_filename );
}

int main( int argc, char* argv[] ) {

  if ( argc != 4 ) {
    std::cout << "Usage: " << argv[0]
      << " INPUT_PELEE_NTUPLE_FILE SELECTION_NAMES OUTPUT_FILE\n";
    return 1;
  }

  std::string input_file_name( argv[1] );
  std::string output_file_name( argv[3] );

  std::vector< std::string > selection_names;

  std::stringstream sel_ss( argv[2] );
  std::string sel_name;
  while ( std::getline(sel_ss, sel_name, ',') ) {
    selection_names.push_back( sel_name );
  }

  analyzer( input_file_name, selection_names, output_file_name );

  return 0;
}
