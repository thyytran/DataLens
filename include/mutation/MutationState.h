#pragma once
#include <string>
#include <map>
#include <mutex>
#include <atomic>
#include <vector>
#include "../alphamissense/AlphaMissenseSummary.h"



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
};

struct ManualMutState {
    int  selectedChainIdx = 0;
    int  position = 0;
    char mutAA = 'A';
    bool runManual = false;
    int selectedResidueIdx = 0;
};

extern MutationPanel    mutationPanel;
extern MutationFocusData mutFocus;