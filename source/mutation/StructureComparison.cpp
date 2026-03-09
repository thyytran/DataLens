// ── StructureComparison.cpp ──────────────────────────────────────────────────

#include "mutation/StructureComparison.h"
#include <cmath>
#include <algorithm>
#include <sstream>
#include <numeric>

// Build a {chain+residueNum -> CA position} map from a MoleculeData
static std::unordered_map<std::string, Vec3> buildCAMap(const MoleculeData* mol) {
    std::unordered_map<std::string, Vec3> map;
    if (!mol) return map;
    for (const auto& atom : mol->atoms) {
        if (atom.name == "CA") {
            std::string key = std::string(1, atom.chain) + "_" + std::to_string(atom.residueNum);
            map[key] = Vec3(atom.coords.getX(), atom.coords.getY(), atom.coords.getZ());
        }
    }
    return map;
}

static float dist3(const Vec3& a, const Vec3& b) {
    float dx = a.getX() - b.getX();
    float dy = a.getY() - b.getY();
    float dz = a.getZ() - b.getZ();
    return std::sqrt(dx * dx + dy * dy + dz * dz);
}

// Call this after mutant PDB is loaded into mutantMoleculeData
void buildComparison(
    StructureComparison& cmp,
    const MoleculeData* wtMol,
    const MoleculeData* mutMol,
    const std::string& wtLabel,
    const std::string& mutLabel,
    int                   mutResNum,
    const std::string& mutChain)
{
    cmp = {};   // reset
    if (!wtMol || !mutMol) return;

    cmp.active = true;
    cmp.wtLabel = wtLabel;
    cmp.mutLabel = mutLabel;
    cmp.mutResNum = mutResNum;
    cmp.mutChain = mutChain;

    auto wtMap = buildCAMap(wtMol);
    auto mutMap = buildCAMap(mutMol);

    float sumSq = 0.0f;
    int   count = 0;

    for (const auto& [key, wtPos] : wtMap) {
        auto mutIt = mutMap.find(key);
        if (mutIt == mutMap.end()) continue;

        // Parse key back to chain + resNum
        size_t sep = key.find('_');
        std::string chain = key.substr(0, sep);
        int         resNum = std::stoi(key.substr(sep + 1));

        float d = dist3(wtPos, mutIt->second);
        sumSq += d * d;
        count++;

        // Get residue name from WT
        std::string resName = "???";
        for (const auto& atom : wtMol->atoms) {
            if (atom.residueNum == resNum &&
                std::string(1, atom.chain) == chain &&
                atom.name == "CA")
            {
                resName = atom.residueName;
                break;
            }
        }

        cmp.residues.push_back({
            resNum,
            chain,
            resName,
            d,
            (resNum == mutResNum && chain == mutChain)
            });
    }

    cmp.globalRmsd = (count > 0) ? std::sqrt(sumSq / count) : 0.0f;

    // Build sorted index lists
    cmp.sortedByResNum.resize(cmp.residues.size());
    std::iota(cmp.sortedByResNum.begin(), cmp.sortedByResNum.end(), 0);
    std::sort(cmp.sortedByResNum.begin(), cmp.sortedByResNum.end(),
        [&](int a, int b) {
            return cmp.residues[a].residueNum < cmp.residues[b].residueNum;
        });

    cmp.sortedByRmsd = cmp.sortedByResNum;
    std::sort(cmp.sortedByRmsd.begin(), cmp.sortedByRmsd.end(),
        [&](int a, int b) {
            return cmp.residues[a].rmsd > cmp.residues[b].rmsd;
        });

    cmp.maxRmsd = 0.001f;
    for (const auto& r : cmp.residues)
        if (r.rmsd > cmp.maxRmsd) cmp.maxRmsd = r.rmsd;
}