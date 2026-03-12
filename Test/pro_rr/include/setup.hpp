#include "globals.h"

RoadRunner* rrg = nullptr;

void rr_setup(std::string antimony_file) {
    if (rrg != nullptr) { delete rrg; rrg = nullptr; }
    
    grn_results.clear();
    gene_names.clear();
    
    loadAntimonyFile(antimony_file.c_str());
  
    string sbml = getSBMLString(NULL);
    freeAll();
    
    rrg = new RoadRunner(sbml);
    rrg->setIntegrator("cvode");
    
    std::vector<string> selections = {"G0", "G1", "G2", "G3", "G4", "G5", "G6", "G7", "G8", "G9"};
    rrg->setSelections(selections);
}
