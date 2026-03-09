#pragma once
// ── StructureComparison.h ────────────────────────────────────────────────────

#pragma once
#include <string>
#include <vector>
#include <unordered_map>
#include <optional>
#include "bio/MoleculeData.h"
#include "imgui/imgui.h"

struct ResidueComparison {
    int         residueNum;
    std::string chain;
    std::string residueName;
    float       rmsd;          // per-residue CA RMSD
    bool        isMutationSite;
};

struct StructureComparison {
    bool                           active = false;
    std::string                    wtLabel;
    std::string                    mutLabel;
    float                          globalRmsd = 0.0f;
    float                          maxRmsd = 0.001f;  
    int                            mutResNum = 0;
    std::string                    mutChain;
    std::vector<ResidueComparison> residues;
    std::vector<int>               sortedByRmsd;
    std::vector<int>               sortedByResNum;

    enum class SortMode { BY_RESIDUE, BY_RMSD };
    SortMode sortMode = SortMode::BY_RESIDUE;
    float    rmsdFilter = 0.0f;
    int      hoveredIdx = -1;
};

void buildComparison(
    StructureComparison& cmp,
    const MoleculeData* wtMol,
    const MoleculeData* mutMol,
    const std::string& wtLabel,
    const std::string& mutLabel,
    int                  mutResNum,
    const std::string& mutChain
);

void renderComparisonPanel(StructureComparison& cmp);
