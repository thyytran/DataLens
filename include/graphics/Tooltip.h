#pragma once
#include <string>
#include <vector>
struct HoverDisplayInfo {
    std::string residueName;  // "GLY", "ALA" - for tooltip
    int residueNum;           // 45 - for Selection
    char chain;               // 'A' - for Selection
    std::string atomName;     // "CA" - for tooltip
    std::string element;      // "C" - for tooltip
};

std::vector<HoverDisplayInfo> atomIndexToDisplayInfo;