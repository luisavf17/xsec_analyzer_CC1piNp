#pragma once

// Standard library includes
#include <string>

// XSecAnalyzer includes
#include "./SelectionBase.hh"

class SelectionFactory {
  public:
    SelectionFactory();
    SelectionBase* CreateSelection( const std::string& selection_name );
};
