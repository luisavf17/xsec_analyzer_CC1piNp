// Standard library includes
#include <algorithm>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <map>
#include <sstream>
#include <string>
#include <vector>

// ROOT includes
#include "TAxis.h"
#include "TCanvas.h"
#include "TFile.h"
#include "THStack.h"
#include "TLegend.h"
#include "TLatex.h"
#include "TLine.h"
#include "TColor.h"
#include "TParameter.h"
#include "TList.h"
#include "TCollection.h"
#include "TPad.h"

// XSecAnalyzer includes
#include "XSecAnalyzer/FilePropertiesManager.hh"
#include "XSecAnalyzer/MCC9SystematicsCalculator.hh"
#include "XSecAnalyzer/PlotUtils.hh"
#include "XSecAnalyzer/SliceBinning.hh"
#include "XSecAnalyzer/SliceHistogram.hh"

using NFT = NtupleFileType;
const int FontStyle = 132;

//#define USE_FAKE_DATA ""

namespace {

  void set_mc_histogram_style( int event_category, TH1* mc_hist, int color ) {
    mc_hist->SetFillColor( color );
    mc_hist->SetLineColor( color );
    mc_hist->SetStats( false );
  }

  void set_ext_histogram_style( TH1* ext_hist ) {
    ext_hist->SetFillColor( 28 );
    ext_hist->SetLineColor( 28 );
    ext_hist->SetLineWidth( 2 );
    ext_hist->SetFillStyle( 3005 );
    ext_hist->SetStats( false );
  }

  void set_bnb_data_histogram_style( TH1* bnb_hist ) {

    bnb_hist->SetLineColor( kBlack );
    bnb_hist->SetLineWidth( 2 );
    bnb_hist->SetMarkerStyle( kFullCircle );
    bnb_hist->SetMarkerSize( 0.8 );
    bnb_hist->SetStats( false );

    bnb_hist->GetXaxis()->SetTitleOffset( 0.0 );
    bnb_hist->GetXaxis()->SetTitleSize( 0.0 );
    bnb_hist->GetYaxis()->SetTitleSize( 0.05 );
    bnb_hist->GetYaxis()->CenterTitle( true );
    bnb_hist->GetXaxis()->SetLabelSize( 0.0 );

    // This prevents the first y-axis label label (0) to be clipped by the
    // ratio plot
    bnb_hist->SetMinimum( 1e-3 );
  }

  void set_stat_err_histogram_style( TH1* stat_err_hist ) {
    stat_err_hist->SetFillColor( kBlack );
    stat_err_hist->SetLineColor( kBlack );
    stat_err_hist->SetLineWidth( 2 );
    stat_err_hist->SetFillStyle( 3004 );
  }

} // anonymous namespace

void tutorial_slice_plots(std::string FPM_Config,
                          std::string SYST_Config,
                          std::string SLICE_Config,
                          std::string Univ_Output,
                          std::string Plot_OutputDir) {

  // Counter to ensure plots aren't overwritten
  uint FileNameCounter = 0;
  uint DetSysFileNameCounter = 0;

  std::string Plot_Prefix = "SlicePlots";
  std::string Plot_Suffix = ".png";

  std::string DetPlot_Prefix = "PlotDetSys";
  std::string DetPlot_Suffix = ".png";

  std::string PlotFileName = Plot_OutputDir + "/" + Plot_Prefix + Form("_%i",FileNameCounter) + Plot_Suffix;

  std::cout << "\nRunning Slice_Plots with options:" << std::endl;
  std::cout << "\tFPM_Config: " << FPM_Config << std::endl;
  std::cout << "\tSYST_Config: " << SYST_Config << std::endl;
  std::cout << "\tSLICE_Config: " <<  SLICE_Config << std::endl;
  std::cout << "\tUniv_Output: " << Univ_Output << std::endl;
  std::cout << "\tPlot_OutputDir: " << Plot_OutputDir << std::endl;
  std::cout << "\t\tWith filename: " << PlotFileName << std::endl;
  std::cout << "\n" << std::endl;

  std::cout << "DEBUG0" << std::endl;

// #ifdef USE_FAKE_DATA
//   // Initialize the FilePropertiesManager and tell it to treat the NuWro
//   // MC ntuples as if they were data
//   auto& fpm = FilePropertiesManager::Instance();
//   fpm.load_file_properties( FPM_Config );
// #endif

  auto& fpm = FilePropertiesManager::Instance();
  fpm.load_file_properties( FPM_Config );

  std::cout << "DEBUG1" << std::endl;

  // ============================================================
  // Compute beam-on POT using FilePropertiesManager
  // ============================================================
  double total_data_POT = 0.0;

  try {

    const auto& data_norm_map = fpm.data_norm_map();

    for (const auto& pair : data_norm_map) {

      const std::string& file_name = pair.first;
      const auto& trig_pot = pair.second;

      if (fpm.get_ntuple_file_type(file_name) == NtupleFileType::kOnBNB) {
        total_data_POT += trig_pot.pot_;
      }
    }

  }
  catch (const std::exception& e) {
    std::cerr << "Exception while computing total_data_POT: "
              << e.what() << std::endl;
  }

  std::cout << "Total beam-on POT = "
            << total_data_POT << std::endl;

  std::cout << "DEBUG2" << std::endl;

  // Check that we can read the universe output file
  TFile* temp_file = new TFile(Univ_Output.c_str(), "read");
  if (!temp_file || temp_file->IsZombie()) {
    std::cerr << "Could not read file: " << Univ_Output << std::endl;
    throw;
  }
  delete temp_file;

  std::cout << "DEBUG3" << std::endl;

  auto* syst_ptr = new MCC9SystematicsCalculator(Univ_Output, SYST_Config);
  // syst_ptr->set_syst_mode(syst_ptr->SystMode::ForXSec);
  syst_ptr->set_syst_mode(syst_ptr->SystMode::VaryBackgroundAndSignalDirectly);
  auto& syst = *syst_ptr;

  std::cout << "DEBUG4" << std::endl;

  // Get access to the relevant histograms owned by the SystematicsCalculator
  // object. These contain the reco bin counts that we need to populate the
  // slices below.
  TH1D* reco_bnb_hist = syst.data_hists_.at( NFT::kOnBNB ).get();
  TH1D* reco_ext_hist = syst.data_hists_.at( NFT::kExtBNB ).get();

#ifdef USE_FAKE_DATA
  // Add the EXT to the "data" when working with fake data
  reco_bnb_hist->Add( reco_ext_hist );
#endif

  TH2D* category_hist = syst.cv_universe().hist_categ_.get();

  // Total MC+EXT prediction in reco bin space. Start by getting EXT.
  TH1D* reco_mc_plus_ext_hist = dynamic_cast< TH1D* >(
    reco_ext_hist->Clone("reco_mc_plus_ext_hist") );
  reco_mc_plus_ext_hist->SetDirectory( nullptr );

  // Add in the CV MC prediction
  reco_mc_plus_ext_hist->Add( syst.cv_universe().hist_reco_.get() );

  // Keys are covariance matrix types, values are CovMatrix objects that
  // represent the corresponding matrices
  auto* matrix_map_ptr = syst.get_covariances().release();
  auto& matrix_map = *matrix_map_ptr;
  std::cout << "\n================ XSEC-RELATED COVARIANCE KEYS ================\n";

  for (const auto& [key, cov] : matrix_map) {

    // Print anything that looks like it belongs to xsec
    if (key.find("xsec") != std::string::npos ||
        key.find("XSec") != std::string::npos ||
        key.find("genie") != std::string::npos ||
        key.find("GENIE") != std::string::npos ||
        key.find("multisim") != std::string::npos ||
        key.find("unisim") != std::string::npos) {

      std::cout << key << std::endl;
    }
  }

  std::cout << "==============================================================\n\n";
  auto* sb_ptr = new SliceBinning( SLICE_Config );
  auto& sb = *sb_ptr;

  for ( size_t sl_idx = 0u; sl_idx < sb.slices_.size(); ++sl_idx ) {

    const auto& slice = sb.slices_.at( sl_idx );

    // We now have all of the reco bin space histograms that we need as input.
    // Use them to make new histograms in slice space.
    SliceHistogram* slice_bnb = SliceHistogram::make_slice_histogram(
      *reco_bnb_hist, slice, &matrix_map.at("BNBstats") );

    SliceHistogram* slice_ext = SliceHistogram::make_slice_histogram(
      *reco_ext_hist, slice, &matrix_map.at("EXTstats") );
    slice_ext->hist_->SetTitle("EXT");

    SliceHistogram* slice_mc_plus_ext = SliceHistogram::make_slice_histogram(
      *reco_mc_plus_ext_hist, slice, &matrix_map.at("total") );

    std::string data_mc_agreement;
    SliceHistogram::Chi2Result chi2_result;

    try {
      chi2_result = slice_bnb->get_chi2(*slice_mc_plus_ext);

      std::cout << "Slice " << sl_idx << ": χ² = "
                << chi2_result.chi2_ << "/"
                << chi2_result.num_bins_
                << " bins, p-value = "
                << chi2_result.p_value_ << '\n';

      data_mc_agreement = Form("#chi^{2} = %.2f/%d bins, p = %.3f",
                               chi2_result.chi2_,
                               chi2_result.num_bins_,
                               chi2_result.p_value_);
    }
    catch (const std::exception& e) {
      data_mc_agreement = "";
    }

    // Build a stack of categorized central-value MC predictions plus the
    // extBNB contribution in slice space
    set_ext_histogram_style( slice_ext->hist_.get() );

    THStack* slice_pred_stack = new THStack( "mc+ext", "" );
    slice_pred_stack->Add( slice_ext->hist_.get() ); // extBNB

    const auto& sel_for_cat = syst.get_selection_for_categories();
    const auto& cat_map = sel_for_cat.category_map();

    // ================= DEBUG BLOCK =================
    std::cout << "\n--- DEBUG CATEGORY INFO ---\n";
    std::cout << "Slice index: " << sl_idx << std::endl;
    std::cout << "Category map size = "
              << cat_map.size() << std::endl;
    std::cout << "Category hist total integral = "
              << category_hist->Integral() << std::endl;
    std::cout << "Category hist X bins = "
              << category_hist->GetNbinsX() << std::endl;
    std::cout << "Category hist Y bins = "
              << category_hist->GetNbinsY() << std::endl;
    std::cout << "-----------------------------\n";

    // Go in reverse so that, if the signal is defined first in the map, it
    // ends up on top. Note that this index is one-based to match the ROOT
    // histograms
    int cat_bin_index = static_cast<int>(cat_map.size());
    int unique_proj_counter = 0;
    for ( auto iter = cat_map.crbegin(); iter != cat_map.crend(); ++iter )
    {
      int cat = iter->first;
      int color = iter->second.second;
      std::string label = iter->second.first;

      // make projection name unique to avoid ROOT name collisions
      std::string proj_name = Form("temp_mc_hist_sl%zu_cat%d_u%d", sl_idx, cat, unique_proj_counter);
      ++unique_proj_counter;

      TH1D* temp_mc_hist = category_hist->ProjectionY(proj_name.c_str(),
                                                      cat_bin_index,
                                                      cat_bin_index);
      temp_mc_hist->SetDirectory(nullptr);

      SliceHistogram* temp_slice_mc =
          SliceHistogram::make_slice_histogram(*temp_mc_hist, slice);

      set_mc_histogram_style(cat, temp_slice_mc->hist_.get(), color);
      temp_slice_mc->hist_->SetTitle(label.c_str());

      slice_pred_stack->Add(temp_slice_mc->hist_.get());

      --cat_bin_index;
    }

    std::cout << "DEBUG8" << std::endl;

    // ========================================
    // CANVAS
    // ========================================
    auto* c1 = new TCanvas("c1", "", 700, 620);

    // ========================================
    // TOP PAD (MAIN DISTRIBUTION)
    // ========================================
    TPad* pad1 = new TPad("pad1", "", 0.0, 0.3, 0.7, 1.0);
    pad1->SetBottomMargin(0.02);
    pad1->SetRightMargin(0.03);
    pad1->SetLeftMargin(0.13);
    pad1->SetGridx();
    pad1->Draw();
    pad1->cd();

    // Hide x axis on top pad
    slice_bnb->hist_->GetXaxis()->SetLabelSize(0.0);
    slice_bnb->hist_->GetXaxis()->SetTitleSize(0.0);

    slice_bnb->hist_->SetLineColor(kBlack);
    slice_bnb->hist_->SetLineWidth(2);
    slice_bnb->hist_->SetStats(false);

    double ymax = std::max(slice_bnb->hist_->GetMaximum(),
                           slice_mc_plus_ext->hist_->GetMaximum()) * 1.3;

    slice_bnb->hist_->GetYaxis()->SetRangeUser(0., ymax);

    // Format POT as A × 10^{B}
    int exponent = 0;
    double mantissa = 0.0;
    if (total_data_POT > 0.0) {
      exponent = static_cast<int>(std::floor(std::log10(total_data_POT)));
      mantissa = total_data_POT / std::pow(10., exponent);
      mantissa = std::round(mantissa * 10.0) / 10.0;
    }

    std::string y_axis_title =
        Form("Events / %.1f #times 10^{%d}", mantissa, exponent);

    slice_bnb->hist_->GetYaxis()->SetTitle(y_axis_title.c_str());
    slice_bnb->hist_->GetYaxis()->CenterTitle(true);
    slice_bnb->hist_->GetYaxis()->SetTitleFont(FontStyle);
    slice_bnb->hist_->GetYaxis()->SetLabelFont(FontStyle);
    slice_bnb->hist_->GetYaxis()->SetTitleSize(0.045);
    slice_bnb->hist_->GetYaxis()->SetLabelSize(0.045);
    slice_bnb->hist_->GetYaxis()->SetTitleOffset(1.2);

    slice_bnb->hist_->Draw("E1");

    // Stack already created earlier
    slice_pred_stack->Draw("hist same");

    // ========================================
    // FULL UNCERTAINTY BAND
    // ========================================
    TH1D* total_unc_band =
      dynamic_cast<TH1D*>(slice_mc_plus_ext->hist_->Clone("total_unc_band"));

    total_unc_band->SetFillColor(kGray+1);
    total_unc_band->SetFillStyle(3004);
    total_unc_band->SetLineWidth(0);
    total_unc_band->Draw("E2 same");

    slice_mc_plus_ext->hist_->SetLineWidth(2);
    slice_mc_plus_ext->hist_->SetLineColor(kGray+1);
    slice_mc_plus_ext->hist_->Draw("hist same");

    slice_bnb->hist_->Draw("E1 same");

// ========================================
// χ² TEXT
// ========================================
TLatex text;
text.SetNDC();
text.SetTextFont(FontStyle);
text.SetTextSize(0.045);

// Create a default/empty chi2 result
SliceHistogram::Chi2Result current_chi2_result;
bool chi2_valid = false;

try {
  current_chi2_result = slice_bnb->get_chi2(*slice_mc_plus_ext);
  chi2_valid = true;

  std::cout << "Slice " << sl_idx << ": χ² = "
            << current_chi2_result.chi2_ << "/"
            << current_chi2_result.num_bins_
            << " bins, p-value = "
            << current_chi2_result.p_value_ << '\n';

  data_mc_agreement = Form("#chi^{2} = %.2f/%d bins, p = %.3f",
                           current_chi2_result.chi2_,
                           current_chi2_result.num_bins_,
                           current_chi2_result.p_value_);
}
catch (const std::exception& e) {
  data_mc_agreement = "";
  // Optionally log the error
  std::cerr << "Could not compute chi2 for slice " << sl_idx 
            << ": " << e.what() << std::endl;
}

// Only draw chi2 if it was successfully computed
if (chi2_valid) {
  std::ostringstream chi2_ss;
  chi2_ss << "#chi^{2} = "
          << std::fixed << std::setprecision(2)
          << current_chi2_result.chi2_
          << " / " << current_chi2_result.num_bins_
          << ",  p = "
          << std::setprecision(3)
          << current_chi2_result.p_value_;
  
  text.DrawLatex(0.18, 0.85, chi2_ss.str().c_str());
} else {
  // Optionally draw something like "Chi² N/A"
  text.DrawLatex(0.18, 0.85, "#chi^{2} not available");
}

    // Repaint the frame on top
    gPad->RedrawAxis();

    c1->cd();

    // ========================================
    // LEGEND PAD (RIGHT)
    // ========================================
    TPad* pad3 = new TPad("pad3", "", 0.7, 0.5, 1.0, 0.94);
    pad3->SetRightMargin(0.06);
    pad3->SetLeftMargin(0.);
    pad3->Draw();
    pad3->cd();

    TLegend* lg = new TLegend(0.,0.,1.,1.);
    lg->SetTextFont(FontStyle);
    lg->SetTextSize(0.07);
    lg->SetBorderSize(0);

    lg->AddEntry(slice_bnb->hist_.get(), "Data (beam on)", "lp");
    lg->AddEntry(slice_mc_plus_ext->hist_.get(), "MC + EXT", "l");
    lg->AddEntry(total_unc_band, "Total Uncertainty", "f");

    // Use actual stack histograms for proper colors
    TList* stack_hists = slice_pred_stack->GetHists();
    if (stack_hists) {

      double total_pred_integral =
          slice_mc_plus_ext->hist_->Integral();

      TIter next(stack_hists);
      TObject* obj;

      while ((obj = next())) {
        TH1* h = dynamic_cast<TH1*>(obj);
        if (!h) continue;

        double cat_integral = h->Integral();
        if (cat_integral <= 0) continue;

        double frac = 0.0;
        if (total_pred_integral > 0)
          frac = cat_integral / total_pred_integral;

        double percent = frac * 100.0;

        std::string legend_label =
            Form("%s (%.1f%%)", h->GetTitle(), percent);

        lg->AddEntry(h, legend_label.c_str(), "f");
      }
    }

    lg->Draw();

    c1->cd();

  // ========================================
  // RATIO PAD
  // ========================================
  TPad* pad2 = new TPad("pad2", "", 0.0, 0.01, 0.7, 0.3);
  pad2->SetTopMargin(0.05);
  pad2->SetBottomMargin(0.38);
  pad2->SetRightMargin(0.03);
  pad2->SetLeftMargin(0.13);
  pad2->SetGridx();
  pad2->Draw();
  pad2->cd();

  // Create ratio histogram with zero errors
  TH1D* h_ratio = dynamic_cast<TH1D*>(slice_bnb->hist_->Clone("h_ratio"));
  TH1D* h_ratio_values = dynamic_cast<TH1D*>(slice_mc_plus_ext->hist_.get()->Clone("h_ratio_values"));

  // Set errors to zero for the denominator
  for(int i = 0; i <= h_ratio_values->GetNbinsX(); i++) {
      h_ratio_values->SetBinError(i, 0);
  }

  // Calculate ratio
  h_ratio->Divide(h_ratio_values);

  // Create error band histogram
  TH1D* h_ratio_error = dynamic_cast<TH1D*>(slice_mc_plus_ext->hist_.get()->Clone("h_ratio_error"));
  h_ratio_error->Divide(h_ratio_values);

  // Set styling for ratio plot
  h_ratio->SetStats(false);
  h_ratio->SetMarkerStyle(kFullCircle);
  h_ratio->SetMarkerSize(0.5);
  h_ratio->SetLineColor(kBlack);

  // Set axis labels and styling
  h_ratio->GetYaxis()->SetTitle("Ratio");
  h_ratio->GetYaxis()->SetTitleSize(0.11);
  h_ratio->GetYaxis()->SetTitleOffset(0.45);
  h_ratio->GetYaxis()->SetLabelSize(0.11);
  h_ratio->GetYaxis()->SetNdivisions(505);
  h_ratio->GetYaxis()->SetRangeUser(0, 2);

  h_ratio->SetBit(TH1::kNoTitle);

  h_ratio->GetXaxis()->SetLabelSize(0.11);
  h_ratio->GetXaxis()->SetTitleSize(0.11);
  h_ratio->GetXaxis()->SetTitleOffset(1);
  h_ratio->GetXaxis()->SetTickLength(0.03);
  h_ratio->GetXaxis()->SetTitle(slice_bnb->hist_->GetXaxis()->GetTitle());
  h_ratio->GetXaxis()->CenterTitle(true);
  h_ratio->GetXaxis()->SetTitleFont(FontStyle);
  h_ratio->GetXaxis()->SetLabelFont(FontStyle);

  // Set error band styling
  h_ratio_error->SetFillColor(kGray+2);
  h_ratio_error->SetFillStyle(3004);
  h_ratio_error->SetLineWidth(0);

  // Draw ratio and error band
  h_ratio->Draw();
  h_ratio_error->Draw("e2 same");

    // Unity line
    TLine* line = new TLine(h_ratio->GetXaxis()->GetXmin(),1.0,
                            h_ratio->GetXaxis()->GetXmax(),1.0);
    line->SetLineStyle(2);
    line->Draw();

    c1->Update();

    // ========================================
    // SAVE MAIN PLOT
    // ========================================
    PlotFileName = Plot_OutputDir + "/" + Plot_Prefix +
                   Form("_%i",FileNameCounter) + Plot_Suffix;
    c1->SaveAs(PlotFileName.c_str());
    FileNameCounter++;

    std::cout << "DEBUG9" << std::endl;

    // Get the binning and axis labels for the current slice by cloning the
    // (empty) histogram owned by the Slice object
    TH1* slice_hist = dynamic_cast< TH1* >(
      slice.hist_->Clone("slice_hist") );

    slice_hist->SetDirectory( nullptr );

    // Keys are labels, values are fractional uncertainty histograms
    auto* fr_unc_hists = new std::map< std::string, TH1* >();
    auto& frac_uncertainty_hists = *fr_unc_hists;

    // Separate container for detector systematics only
    std::map< std::string, TH1* > det_frac_uncertainty_hists;

    // XSEC-only breakdown
    std::map< std::string, TH1* > xsec_frac_uncertainty_hists;
        
    const std::vector< std::string > cov_mat_keys = {
      "total",
      "detVar_total",
      "flux",
      "reint",
      "xsec_total",
      "POT",
      "numTargets",
      "MCstats",
      "EXTstats",
      "BNBstats"
    };

    // Detector-only plot order
    const std::vector< std::string > det_mat_keys = {
      "detVarLYatten",
      "detVarLYdown",
      "detVarLYrayl",
      "detVarRecomb2",
      "detVarSCE",
      "detVarWMAngleXZ",
      "detVarWMAngleYZ",
      "detVarWMX",
      "detVarWMYZ",
      "detVar_total"
    };

    const std::vector< int > cov_mat_colors = {
      TColor::GetColor("#D4B308"),
      TColor::GetColor("#FF6347"),
      TColor::GetColor("#483D8B"),
      TColor::GetColor("#B0C4DE"),
      TColor::GetColor("#D4B308"),
      TColor::GetColor("#FF6347"),
      TColor::GetColor("#483D8B"),
      TColor::GetColor("#B0C4DE"),
      kOrange,
      12,
      kYellow
    };

    const std::vector<int> cov_mat_linestyle = {1,1,1,1,2,2,2,2,2,2};

    // Loop over the various systematic uncertainties
    int color = 1;
    int color_counter = 0;
    for ( const auto& pair : matrix_map ) {

      const auto& key = pair.first;
      const auto& cov_matrix = pair.second;

      SliceHistogram* slice_for_syst = SliceHistogram::make_slice_histogram(
        *reco_mc_plus_ext_hist, slice, &cov_matrix );

      // The SliceHistogram object already set the bin errors appropriately
      // based on the slice covariance matrix. Just change the bin contents
      // for the current histogram to be fractional uncertainties. Also set
      // the "uncertainties on the uncertainties" to zero.
      // TODO: revisit this last bit, possibly assign bin errors here
      for ( const auto& bin_pair : slice.bin_map_ ) {
        int global_bin_idx = bin_pair.first;
        double y = slice_for_syst->hist_->GetBinContent( global_bin_idx );
        double err = slice_for_syst->hist_->GetBinError( global_bin_idx );
        double frac = 0.;
        if ( y > 0. ) frac = err / y;
        slice_for_syst->hist_->SetBinContent( global_bin_idx, frac );
        slice_for_syst->hist_->SetBinError( global_bin_idx, 0. );
      }

      // --- PRINT FULL BREAKDOWN (ALL SYSTEMATICS) ---
      std::cout << "---- Slice " << sl_idx
                << " : fractional uncertainty in bin #1 ----\n";

      double frac_bin1 = slice_for_syst->hist_->GetBinContent(1) * 100.0;
      std::cout << key << " = " << frac_bin1 << " %\n";

      std::cout << "DEBUG10" << std::endl;

      // Keep ALL detector variations, including detVar_total,
      // for the detector-only plot
      if (key.rfind("detVar", 0) == 0) {
        det_frac_uncertainty_hists[key] = slice_for_syst->hist_.get();
      }

      // do same for xsec variations for xsec-only plot
      // Keep all xsec-related terms for the dedicated xsec breakdown plot
      if (key.rfind("xsec_", 0) == 0) {
        xsec_frac_uncertainty_hists[key] = slice_for_syst->hist_.get();
      }

      // For the main plot, only keep the selected list
      auto cbegin = cov_mat_keys.cbegin();
      auto cend = cov_mat_keys.cend();
      auto iter = std::find( cbegin, cend, key );
      if ( iter == cend ) continue;

      frac_uncertainty_hists[ key ] = slice_for_syst->hist_.get();

      if ( color <= 9 ) ++color;
      if ( color == 5 ) ++color;
      if ( color >= 10 ) color += 11;

      slice_for_syst->hist_->SetLineColor( cov_mat_colors[color_counter] );
      slice_for_syst->hist_->SetLineStyle( cov_mat_linestyle[color_counter] );
      slice_for_syst->hist_->SetLineWidth( 2 );
      color_counter += 1;
    }

    std::cout << "DEBUG11" << std::endl;

    // ========================================
    // SYSTEMATICS BREAKDOWN PLOT
    // ========================================
    TCanvas* c2 = new TCanvas("c2", "", 700, 580);

    // Create legend at top with 5 columns
    TLegend* lg2 = new TLegend(0.1, 0.92, 0.9, 0.99);
    lg2->SetNColumns(5);
    lg2->SetColumnSeparation(0.05);
    lg2->SetBorderSize(0);
    lg2->SetTextSize(0.035);
    lg2->SetTextFont(FontStyle);

    auto* total_frac_err_hist = frac_uncertainty_hists.at( "total" );
    // add margin on the left
    total_frac_err_hist->GetYaxis()->SetTitleOffset(1.1);
    total_frac_err_hist->SetStats( false );
    total_frac_err_hist->GetYaxis()->SetRangeUser( 0.,
      total_frac_err_hist->GetMaximum() * 1.05 );
    total_frac_err_hist->SetLineColor( kBlack );
    total_frac_err_hist->SetLineWidth( 2 );
    total_frac_err_hist->SetLineStyle( 1 );
    total_frac_err_hist->GetXaxis()->SetLabelSize(0.035);
    total_frac_err_hist->GetYaxis()->SetLabelSize(0.035);
    total_frac_err_hist->GetXaxis()->SetTitleSize(0.04);
    total_frac_err_hist->GetYaxis()->SetTitleSize(0.04);
    total_frac_err_hist->GetXaxis()->SetLabelFont(FontStyle);
    total_frac_err_hist->GetYaxis()->SetLabelFont(FontStyle);
    total_frac_err_hist->GetXaxis()->SetTitleFont(FontStyle);
    total_frac_err_hist->GetYaxis()->SetTitleFont(FontStyle);

    // Use the same total_data_POT computed earlier
    double pot = total_data_POT;
    std::cout << "On-beam POT = " << pot << std::endl;

    // total_frac_err_hist is fractional uncertainty.
    total_frac_err_hist->GetYaxis()->SetTitle("Fractional Uncertainty");
    total_frac_err_hist->GetYaxis()->CenterTitle(true);
    total_frac_err_hist->Draw( "hist" );

    lg2->AddEntry( total_frac_err_hist, "total", "l" );

    for ( auto& pair : frac_uncertainty_hists ) {
      const auto& name = pair.first;
      TH1* hist = pair.second;
      // We already plotted the "total" one above
      if ( name == "total" ) continue;

      std::string legend_name;

      if (name == "detVar_total"){
        legend_name = "detvar";
      }
      else if (name == "xsec_total"){
        legend_name = "xsec";
      }
      else {
        legend_name = name;
      }

      lg2->AddEntry( hist, legend_name.c_str(), "l" );
      hist->Draw( "same hist" );

      std::cout << name << " frac err in bin #1 = "
                << hist->GetBinContent( 1 )*100. << "%\n";
    }

    std::cout << "DEBUG12" << std::endl;

    lg2->Draw( "same" );

    std::cout << "Total frac error in bin #1 = "
              << total_frac_err_hist->GetBinContent( 1 )*100. << "%\n";

    PlotFileName = Plot_OutputDir + "/" + Plot_Prefix +
                   Form("_%i",FileNameCounter) + Plot_Suffix;
    c2->SaveAs(PlotFileName.c_str());
    FileNameCounter += 1;

    // ========================================
    // DETECTOR BREAKDOWN PLOT
    // ========================================
    if (!det_frac_uncertainty_hists.empty()) {

      TCanvas* cDet = new TCanvas(
        Form("cDet_%zu", sl_idx), "", 700, 620 );

      TLegend* lgDet = new TLegend(0.10, 0.85, 0.90, 0.99);
      lgDet->SetNColumns(3);
      lgDet->SetColumnSeparation(0.045);
      lgDet->SetBorderSize(0);
      lgDet->SetTextSize(0.035);
      lgDet->SetTextFont(FontStyle);

      // Same palette as your main plot
      const int col_gold   = TColor::GetColor("#D4B308");
      const int col_tomato = TColor::GetColor("#FF6347");
      const int col_blue   = TColor::GetColor("#483D8B");
      const int col_steel  = TColor::GetColor("#B0C4DE");

      std::map<std::string, int> det_color_map = {
        {"detVarLYatten",   col_gold},
        {"detVarSCE",       col_gold},

        {"detVarLYdown",    col_tomato},
        {"detVarWMAngleXZ", col_tomato},

        {"detVarLYrayl",    col_blue},
        {"detVarWMAngleYZ", col_blue},

        {"detVarRecomb2",   col_steel},
        {"detVarWMX",       col_steel},

        {"detVarWMYZ",      kRed+1},

        {"detVar_total",    kBlack}
      };

      // Track how many times each color is used → assign style
      std::map<int, int> color_usage_counter;

      // Find max
      double det_ymax = 0.0;
      for (const auto& name : det_mat_keys) {
        auto it = det_frac_uncertainty_hists.find(name);
        if (it == det_frac_uncertainty_hists.end()) continue;
        det_ymax = std::max(det_ymax, it->second->GetMaximum());
      }
      if (det_ymax <= 0.0) det_ymax = 1.0;

      bool first_det = true;

      for (const auto& name : det_mat_keys) {
        auto it = det_frac_uncertainty_hists.find(name);
        if (it == det_frac_uncertainty_hists.end()) continue;

        TH1* hist = it->second;
        hist->SetStats(false);
        hist->SetTitle("");

        int color = kBlack;
        if (det_color_map.count(name)) {
          color = det_color_map[name];
        }

        int linestyle = 1;

        if (name == "detVar_total") {
          // keep total as solid black
          linestyle = 1;
        } else {
          int count = color_usage_counter[color];

          // alternate styles within same color
          if (count == 0) linestyle = 1;      // solid
          else if (count == 1) linestyle = 2; // dashed
          else linestyle = 3;                 // dotted (fallback)

          color_usage_counter[color]++;
        }

        hist->SetLineColor(color);
        hist->SetLineStyle(linestyle);
        hist->SetLineWidth(name == "detVar_total" ? 3 : 2);

        hist->GetYaxis()->SetTitle("Fractional Uncertainty");
        hist->GetYaxis()->CenterTitle(true);
        hist->GetXaxis()->SetLabelSize(0.035);
        hist->GetYaxis()->SetLabelSize(0.035);
        hist->GetXaxis()->SetTitleSize(0.04);
        hist->GetYaxis()->SetTitleSize(0.04);
        hist->GetXaxis()->SetLabelFont(FontStyle);
        hist->GetYaxis()->SetLabelFont(FontStyle);
        hist->GetXaxis()->SetTitleFont(FontStyle);
        hist->GetYaxis()->SetTitleFont(FontStyle);
        hist->GetYaxis()->SetRangeUser(0., det_ymax * 1.15);

        if (first_det) {
          //add top margin so that the legend goers above the plot
          cDet->SetTopMargin(0.18);
          cDet->SetLeftMargin(0.11);
          cDet->SetRightMargin(0.07);
          hist->Draw("hist");
          first_det = false;
        } else {
          hist->Draw("same hist");
        }

        lgDet->AddEntry(hist, name.c_str(), "l");
      }

      lgDet->Draw("same");
      gPad->RedrawAxis();
      cDet->Update();

      std::string detPlotName =
        Plot_OutputDir + "/" + DetPlot_Prefix +
        Form("_%u", DetSysFileNameCounter) + DetPlot_Suffix;

      cDet->SaveAs(detPlotName.c_str());
      ++DetSysFileNameCounter;
    }

    // ========================================
    // XSEC BREAKDONW PLOT
    // ========================================
    if (!xsec_frac_uncertainty_hists.empty()) {

      TCanvas* cXsec = new TCanvas("cXsec", "", 700, 620);

      TLegend* lgXsec = new TLegend(0.10, 0.84, 0.90, 0.99);
      lgXsec->SetNColumns(3);
      lgXsec->SetColumnSeparation(0.045);
      lgXsec->SetBorderSize(0);
      lgXsec->SetTextSize(0.033);
      lgXsec->SetTextFont(FontStyle);

      const int cGold  = TColor::GetColor("#D4B308");
      const int cTom   = TColor::GetColor("#FF6347");
      const int cBlue  = TColor::GetColor("#483D8B");
      const int cSteel = TColor::GetColor("#B0C4DE");

      // Same style family you were already using: gold / tomato / blue / steel
      const std::vector<std::string> xsec_order = {
        "xsec_total",
        "xsec_multi",
        "xsec_unisim",
        "xsec_AxFFCCQEshape",
        "xsec_DecayAngMEC",
        "xsec_NormCCCOH",
        "xsec_NormNCCOH",
        "xsec_RPA_CCQE",
        "xsec_ThetaDelta2NRad",
        "xsec_Theta_Delta2Npi",
        "xsec_VecFFCCQEshape",
        "xsec_XSecShape_CCMEC",
        "xsec_xsr_scc_Fa3_SCC",
        "xsec_xsr_scc_Fv3_SCC"
      };

      double xsec_ymax = 0.0;
      for (const auto& name : xsec_order) {
        auto it = xsec_frac_uncertainty_hists.find(name);
        if (it == xsec_frac_uncertainty_hists.end()) continue;
        xsec_ymax = std::max(xsec_ymax, it->second->GetMaximum());
      }
      if (xsec_ymax <= 0.0) xsec_ymax = 1.0;

      int xsec_idx = 0;
      bool first_xsec = true;

      for (const auto& name : xsec_order) {
        auto it = xsec_frac_uncertainty_hists.find(name);
        if (it == xsec_frac_uncertainty_hists.end()) continue;

        TH1* hist = it->second;
        hist->SetStats(false);
        hist->SetTitle("");

        int color = kBlack;
        int linestyle = 1;
        int linewidth = 2;

        if (name == "xsec_total") {
          color = kBlack;
          linestyle = 1;
          linewidth = 3;
        } else if (name == "xsec_multi") {
          color = cGold;
          linestyle = 1;
        } else if (name == "xsec_unisim") {
          color = cTom;
          linestyle = 2;
        } else {
          const std::vector<int> palette = {cGold, cTom, cBlue, cSteel};
          const std::vector<int> styles   = {1, 2, 3};
          color = palette[xsec_idx % palette.size()];
          linestyle = styles[(xsec_idx / palette.size()) % styles.size()];
          ++xsec_idx;
        }

        hist->SetLineColor(color);
        hist->SetLineStyle(linestyle);
        hist->SetLineWidth(linewidth);

        hist->GetYaxis()->SetTitle("Fractional Uncertainty");
        hist->GetYaxis()->CenterTitle(true);
        hist->GetXaxis()->SetLabelSize(0.035);
        hist->GetYaxis()->SetLabelSize(0.035);
        hist->GetXaxis()->SetTitleSize(0.04);
        hist->GetYaxis()->SetTitleSize(0.04);
        hist->GetXaxis()->SetLabelFont(FontStyle);
        hist->GetYaxis()->SetLabelFont(FontStyle);
        hist->GetXaxis()->SetTitleFont(FontStyle);
        hist->GetYaxis()->SetTitleFont(FontStyle);
        hist->GetYaxis()->SetRangeUser(0., xsec_ymax * 1.15);

        if (first_xsec) {
          hist->Draw("hist");
          first_xsec = false;
        } else {
          hist->Draw("same hist");
        }

        lgXsec->AddEntry(hist, name.c_str(), "l");
      }

      lgXsec->Draw("same");
      gPad->RedrawAxis();
      cXsec->Update();

      std::string xsecPlotName =
        Plot_OutputDir + "/xsecVars_" + Form("%zu", sl_idx) + ".png";

      cXsec->SaveAs(xsecPlotName.c_str());
    }
    
  }

}

int main(int argc, char* argv[]) {
  if ( argc != 6 ) {
    std::cout << "Usage: Slice_Plots FPM_Config"
              << " SYST_Config SLICE_Config Univ_Output Plot_OutputDir\n";
    return 1;
  }

  std::string list_file_name( argv[1] );
  std::string univmake_config_file_name( argv[2] );
  std::string output_file_name( argv[3] );

  std::string FPM_Config( argv[1] );
  std::string SYST_Config( argv[2] );
  std::string SLICE_Config( argv[3] );
  std::string Univ_Output( argv[4] );
  std::string Plot_OutputDir( argv[5] );

  tutorial_slice_plots(FPM_Config, SYST_Config, SLICE_Config, Univ_Output, Plot_OutputDir);
  return 0;
}