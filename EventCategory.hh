#pragma once

#include "TH1.h"

// Enum used to label event categories of interest for analysis plots
enum EventCategory {

  // Unable to categorize (e.g., because the event is real data and thus
  // has no MC truth information)
  kUnknown = 0,

  // Signal
  kNumuCC1PiChargedNonGolden = 1,
  kNumuCC1PiChargedGolden = 2,

  // MC Background
  kExternal = 3,
  kDirt = 4,
  kNonFiducial = 5,
  kNue = 6,
  kNC = 7,
  kNumuCC0Pi = 8,
  kNumuCC1PiZero = 9,
  kNumuCCOther = 10,
  kNumuCC1PiNonSignal = 11,
  kNumuCC0PiSignal = 12
};

// Singleton class that helps manipulate EventCategory enum values
class EventCategoryInterpreter {

  public:

    // This is a singleton class, so we'll delete the copy constructor
    // the move constructor, and the assignment operators
    EventCategoryInterpreter( const EventCategoryInterpreter& ) = delete;
    EventCategoryInterpreter( EventCategoryInterpreter&& ) = delete;
    EventCategoryInterpreter& operator=( const EventCategoryInterpreter& )
      = delete;
    EventCategoryInterpreter& operator=( EventCategoryInterpreter&& )
      = delete;

    // Get a const reference to the singleton instance of the
    // EventCategoryInterpreter
    inline static const EventCategoryInterpreter& Instance() {

      // Create the EventCategoryInterpreter object using a static variable.
      // This ensures that the singleton instance is only created once.
      static std::unique_ptr<EventCategoryInterpreter>
        the_instance( new EventCategoryInterpreter() );

      // Return a reference to the singleton instance
      return *the_instance;
    }

    inline const std::map< EventCategory, std::string >& label_map() const
      { return event_category_to_label_map_; }

    inline std::string label( EventCategory ec ) const
      { return event_category_to_label_map_.at( ec ); }

    inline int color_code( EventCategory ec ) const
      { return event_category_to_color_map_.at( ec ); }

    inline void set_mc_histogram_style( EventCategory ec, TH1* mc_hist ) const
    {
      int color = color_code( ec );
      mc_hist->SetFillColor( color );
      mc_hist->SetLineColor( color );
      mc_hist->SetStats( false );
    }

    inline void set_ext_histogram_style( TH1* ext_hist ) const {
      ext_hist->SetFillColor( 28 );
      ext_hist->SetLineColor( 28 );
      ext_hist->SetLineWidth( 2 );
      ext_hist->SetFillStyle( 3005 );
      ext_hist->SetStats( false );
    }

    inline void set_bnb_data_histogram_style( TH1* bnb_hist ) const {

      bnb_hist->SetLineColor( kBlack );
      bnb_hist->SetLineWidth( 3 );
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

    inline void set_stat_err_histogram_style( TH1* stat_err_hist ) const {
      stat_err_hist->SetFillColor( kBlack );
      stat_err_hist->SetLineColor( kBlack );
      stat_err_hist->SetLineWidth( 2 );
      stat_err_hist->SetFillStyle( 3004 );
    }

  private:

    EventCategoryInterpreter() {}

    std::map< EventCategory, std::string > event_category_to_label_map_ = {
      { kUnknown, "Unknown/Data" },
      { kNumuCC1PiChargedGolden, "#nu_{#mu} CC1#pi^{#pm} for p_{#pi}" },
      { kNumuCC1PiChargedNonGolden, "#nu_{#mu} CC1#pi^{#pm} Other" },
      { kNumuCC1PiNonSignal, "#nu_{#mu} CC1#pi^{#pm} Non-Signal" },
      { kNumuCC0PiSignal, "#nu_{#mu} CC0#pi Constraint Signal" },
      { kNumuCC0Pi, "#nu_{#mu} CC0#pi Other" },
      { kNumuCC1PiZero, "#nu_{#mu} CC1#pi^{0}" },
      { kNumuCCOther, "Other #nu_{#mu} CC" },
      { kNue, "#nu_{e}" },
      { kNC, "NC" },
      { kDirt, "Dirt" },
      { kNonFiducial, "Non-Fiducial" },
      { kExternal, "External" }
    };

    std::map< EventCategory, Color_t > event_category_to_color_map_ = {
      { kUnknown, kGray },
      { kNumuCC1PiChargedGolden,  kGreen },
      { kNumuCC1PiChargedNonGolden, kViolet},
      { kNumuCC1PiNonSignal , kOrange + 5 },
      { kNumuCC0PiSignal, kYellow - 6},
      { kNumuCC0Pi, kYellow },
      { kNumuCC1PiZero, kOrange - 4 },
      { kNumuCCOther, kOrange - 3 },
      { kNue, kBlue },
      { kNC, kBlue+10 },
      { kDirt, kOrange+3 },
      { kNonFiducial, kBlack },
      { kExternal, kBlack }
    };
};
