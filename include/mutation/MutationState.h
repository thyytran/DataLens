#pragma once
#include <string>
#include <map>
#include <mutex>
#include <atomic>
#include <vector>
#include "../alphamissense/AlphaMissenseSummary.h"

// Sidechain atom classification — populated once when mutant loads
std::vector<std::string> lostAtoms;    // in WT only
std::vector<std::string> gainedAtoms;  // in mutant only  
std::vector<std::string> sharedAtoms;  // in both


struct MutationPanel {
    std::map<char, ChainAmSummary> chains;  // chain -> per-residue AM data
    std::atomic<bool> loading{ false };
    std::atomic<int>  loadedCount{ 0 };
    int               totalCount = 0;
    std::mutex        dataMutex;
    bool              showWindow = false;
    std::string       statusMsg;
};

struct MutationFocusData {
    bool        show = false;
    std::string variantId = "";
    std::string wtAA = "";
    std::string mutAA = "";
    std::string chain = "";
    int         position = -1;
    float       amScore = 0.0f;
    std::string amClass = "";
    float       ddg = 0.0f;
    std::string ddgInterp = "";
    std::string wtProps = "";
    std::string mutProps = "";
    bool        analyzing = false;
    std::string secondaryStructure;   
    std::string uniprotId;        // "P38398"
    int         uniprotPosition = 0;

    // Sidechain atom name sets — computed once when mutant loads, read everywhere
    std::vector<std::string> lostAtoms;    // in WT only  (orange in 3D)
    std::vector<std::string> gainedAtoms;  // in mutant only (magenta in 3D)
    std::vector<std::string> sharedAtoms;  // in both (blue in 3D)
};

struct ManualMutState {
    int  selectedChainIdx = 0;
    int  position = 0;
    char mutAA = 'A';
    bool runManual = false;
    int selectedResidueIdx = 0;
    // UniProt variant fetch
    bool fetchUniProtVariants = false;
    std::string uniprotIdToFetch = "";
    int positionToFetch = 0;

};

extern MutationPanel    mutationPanel;
extern MutationFocusData mutFocus;