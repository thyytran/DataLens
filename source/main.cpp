#define NOMINMAX

// ── Standard library ───────────────────────────────────────────────────────
#include <algorithm>
#include <deque>
#include <iostream>
#include <optional>
#include <set>
#include <stdexcept>
#include <string>
#include <thread>
#include <utility>
#include <vector>

// ── Third-party ────────────────────────────────────────────────────────────
#include "imgui/imgui.h"
#include "imgui/imgui_impl_glfw.h"
#include "imgui/imgui_impl_opengl3.h"
#include "nlohmann/json.hpp"

// ── Math ───────────────────────────────────────────────────────────────────
#include "math/Vec.h"
#include "math/Mat.h"
#include "math/MathUtils.h"

// ── Graphics ───────────────────────────────────────────────────────────────
#include "graphics/Color.h"
#include "graphics/ConnectorTemplate.h"
#include "graphics/Model.h"
#include "graphics/RayCaster.h"
#include "graphics/Shader.h"
#include "graphics/SphereTemplate.h"
#include "graphics/Tooltip.h"
#include "graphics/Window.h"
#include "graphics/RepresentationRenderer.h"
#include "graphics/RepresentationType.h"

// ── Biology / data ─────────────────────────────────────────────────────────
#include "bio/Protein.h"
#include "bio/PDBFile.h"
#include "bio/MoleculeData.h"
#include "bio/ConnectorType.h"

// ── Mutation pipeline ──────────────────────────────────────────────────────
#include "alphamissense/AlphaMissenseDB.h"
#include "mapper/PDBToUniProtMapper.h"
#include "mutation/MutationAnalyzer.h"
#include "mutation/MutationSuggester.h"
#include "mutation/StructureComparison.h"

// ── Application ────────────────────────────────────────────────────────────
#include "assistant/ChatWindow.h"
#include "assistant/ClaudeHandler.h"
#include "ColorCommand.h"
#include "FetchPDB.h"
#include "Input.h"
#include "Parser.h"
#include "ResourceManager.h"
#include "Selection.h"
#include "../include/mutation/MutationState.h"

// ── Config ────────────────────────────────────────────────────────────
#include "Config.h"
#include "filesystem"

using json = nlohmann::json;

// ── Configuration ──────────────────────────────────────────────────────────

#ifdef _DEBUG
AppConfig config = loadConfig("config.json");  // relative to x64/Debug/
#else
AppConfig config = loadConfig("config.json");  // relative to x64/Release/
#endif

//static const std::string BACKEND_URL = "http://localhost:8000";
static const std::string BACKEND_URL = config.backendUrl;

std::string apiKey = config.aiApiKey;
// std::string apiKey = "config.aiApiKey";

// ── Application control ────────────────────────────────────────────────────
bool shouldExit = false;

// ── Window & ImGui shared state ────────────────────────────────────────────
GLFWwindow* mainWindowPtr = nullptr;
std::atomic<bool> mainWindowReady(false);
ImFontAtlas* g_sharedFontAtlas = nullptr;
ImGuiContext* g_mainImGuiContext = nullptr;

// ── Chat ───────────────────────────────────────────────────────────────────
std::shared_ptr<ChatWindow> chatWindowPtr;
std::vector<std::string> chatHistory;
char chatInputBuffer[512] = "";
bool chatProcessing = false;
std::mutex chatMutex;

// ── Molecule & structure ───────────────────────────────────────────────────
MoleculeData* currentMolData = nullptr;
PDBStructureData currentStructureData;
PDBUniProtMapper pdbMapper;
std::string currentPdbId = "";
std::string lastChainId;
json allChains;
json* pAllChains = &allChains;
std::vector<std::pair<int, std::string>> chainResidues;  // {residueNum, residueName}
std::vector<HoverDisplayInfo> atomToResidueMap;
std::vector<std::string> commands;

// ── Shaders for representations ─────────────────────────────────────────────
Shader* ribbonShader = nullptr;
Shader* surfaceShader = nullptr;
RepresentationRenderer representationRenderer;

// ── Mutation pipeline ──────────────────────────────────────────────────────
AlphaMissenseDB      alphaMissenseDB;
MutationSuggester* mutationSuggester = nullptr;
MutationAnalyzer* mutationAnalyzer = nullptr;
std::atomic<bool>    mutationAnalysisRunning(false);
ManualMutState manualMut;
MutationPanel mutationPanel;

// ── Comparison odule ──────────────────────────────────────────────────────
StructureComparison structureComparison;
MoleculeData* mutantMoleculeData = nullptr;
std::string mutantPdbContent = "";
std::string mutantPdbPath = "";

// -- Shader programs ---------

const char* ribbonVS = R"(#version 330 core
layout (location = 0) in vec3 aPos;
layout (location = 1) in vec3 aNormal;
layout (location = 2) in vec2 aTexCoord;
layout (location = 3) in vec4 aColor;
out vec3 FragPos;
out vec3 Normal;
out vec2 TexCoord;
out vec4 Color;
uniform mat4 u_model;
uniform mat4 u_view;
uniform mat4 u_projection;
void main() {
    vec4 worldPos = u_model * vec4(aPos, 1.0);
    FragPos = worldPos.xyz;
    Normal = mat3(transpose(inverse(u_model))) * aNormal;
    TexCoord = aTexCoord;
    Color = aColor;
    gl_Position = u_projection * u_view * worldPos;
})";

const char* ribbonFS = R"(#version 330 core
in vec3 FragPos;
in vec3 Normal;
in vec2 TexCoord;
in vec4 Color;
out vec4 FragColor;
uniform vec3 u_lightDir;
uniform vec3 u_cameraPos;
void main() {
    vec3 n = normalize(Normal);
    vec3 l = normalize(u_lightDir);
    vec3 v = normalize(u_cameraPos - FragPos);
    if (dot(n, v) < 0.0) n = -n;
    float diff = max(dot(n, l), 0.0);
    vec3 h = normalize(l + v);
    float spec = pow(max(dot(n, h), 0.0), 32.0);
    vec3 result = (0.15 + 0.7 * diff) * Color.rgb + 0.4 * spec * vec3(1.0);
    FragColor = vec4(result, Color.a);
})";

const char* surfaceVS = R"(#version 330 core
layout (location = 0) in vec3 aPos;
layout (location = 1) in vec3 aNormal;
layout (location = 2) in vec4 aColor;
out vec3 FragPos;
out vec3 Normal;
out vec4 Color;
uniform mat4 u_model;
uniform mat4 u_view;
uniform mat4 u_projection;
void main() {
    vec4 worldPos = u_model * vec4(aPos, 1.0);
    FragPos = worldPos.xyz;
    Normal = mat3(transpose(inverse(u_model))) * aNormal;
    Color = aColor;
    gl_Position = u_projection * u_view * worldPos;
})";

const char* surfaceFS = R"(#version 330 core
in vec3 FragPos;
in vec3 Normal;
in vec4 Color;
out vec4 FragColor;
uniform vec3 u_lightDir;
uniform vec3 u_cameraPos;
void main() {
    vec3 n = normalize(Normal);
    vec3 l = normalize(u_lightDir);
    vec3 v = normalize(u_cameraPos - FragPos);
    if (dot(n, v) < 0.0) n = -n;
    float diff = max(dot(n, l), 0.0);
    vec3 h = normalize(l + v);
    float spec = pow(max(dot(n, h), 0.0), 16.0);
    vec3 result = (0.2 + 0.6 * diff) * Color.rgb + 0.3 * spec * vec3(1.0);
    FragColor = vec4(result, 0.8);
})";

// -- Shader program ends ---

/*
* The chat window shares OpenGL resources (font atlas, ImGui context) with the main
* window. It must not initialize until the main window has completed its first frame,
* otherwise glfwMakeContextCurrent and ImGui setup will race or crash.
*/
void chatWindowThread() {
	// Wait for main window to be created
	while (!mainWindowPtr || !mainWindowReady) {
		std::this_thread::sleep_for(std::chrono::milliseconds(100));
	}
	std::cout << "Starting chat window after main window is ready..." << std::endl;

	// Reuse the font atlas already uploaded to the GPU by the main ImGui context.
	// Creating a second atlas would duplicate texture memory and cause glyph mismatches
	// between the two windows.
	ImFontAtlas* const sharedFonts = ImGui::GetIO().Fonts;  // Get from main context

	ChatWindow chatWindow(1000, 800, mainWindowPtr);
	chatWindowPtr = std::shared_ptr<ChatWindow>(&chatWindow, [](ChatWindow*) {});

	chatWindow.run();

	chatWindowPtr.reset();
}

void addChatMessage(const std::string& sender, const std::string& message) {
	std::lock_guard<std::mutex> lock(chatMutex);
	chatHistory.push_back(sender + ": " + message);
}

std::string getPrimaryUniProtID(const std::string& pdbID, char chain) {
	HTTPConnection conn;
	std::string url = "https://data.rcsb.org/rest/v1/core/uniprot/" +
		pdbID + "/" + std::string(1, chain);

	if (conn.get(url)) {
		json response = json::parse(conn.response);
		return response["rcsb_uniprot_container_identifiers"][0]
			["uniprot_id"].get<std::string>();
	}
	return "";
}

// Count how many residues are within a distance threshold
int countNearbyResidues(const Atom& targetAtom, MoleculeData* moleculeData, float threshold) {
	int count = 0;
	int targetResNum = targetAtom.residueNum;
	char targetChain = targetAtom.chain;

	/* A residue's own atoms are always within threshold distance of each other, so
	* counting them would inflate the result. We compare both residueNum and chain
	* because the same sequence number can appear on different chains.
	*/
	for (const auto& atom : moleculeData->atoms) {

		// Skip same residue
		if (atom.residueNum == targetResNum && atom.chain == targetChain) {
			continue;
		}

		// Euclidean distance in Angstroms between the two atom coordinate positions.
		// PDB coordinates are stored in Angstroms, so no unit conversion is needed.
		// std::sqrt is intentional here — squared distance would change threshold semantics.		
		float dx = atom.coords.getX() - targetAtom.coords.getX();
		float dy = atom.coords.getY() - targetAtom.coords.getY();
		float dz = atom.coords.getZ() - targetAtom.coords.getZ();
		float distance = std::sqrt(dx * dx + dy * dy + dz * dz);

		if (distance < threshold) {
			count++;
		}
	}

	return count;
}

// Get chemical properties of amino acid
std::string getResidueProperties(const std::string& resName) {
	// Nonpolar side chains. These residues avoid water and tend to cluster in the
	// protein core. PRO is included despite its rigidity because its side chain
	// is also nonpolar.
	if (resName == "ALA" || resName == "VAL" || resName == "ILE" ||
		resName == "LEU" || resName == "MET" || resName == "PHE" ||
		resName == "TRP" || resName == "PRO") {
		return "Hydrophobic";
	}

	// Basic residues carrying a net positive charge at physiological pH (~7.4).
	// HIS is borderline (pKa ~6) and may be neutral in some environments,
	// but is conventionally grouped here.
	else if (resName == "LYS" || resName == "ARG" || resName == "HIS") {
		return "Positively charged";
	}

	// Acidic residues carrying a net negative charge at physiological pH.
	// Both have carboxylate side chains that are fully deprotonated under normal conditions.
	else if (resName == "ASP" || resName == "GLU") {
		return "Negatively charged";
	}

	// Uncharged but hydrophilic residues. They form hydrogen bonds via -OH, -SH,
	// -NH2, or -CONH2 groups. CYS is included here but can behave as nonpolar
	// depending on its redox state (free thiol vs. disulfide bond).
	else if (resName == "SER" || resName == "THR" || resName == "ASN" ||
		resName == "GLN" || resName == "TYR" || resName == "CYS") {
		return "Polar";
	}

	// GLY has no side chain — just a hydrogen — making it uniquely flexible.
	// It introduces conformational freedom and is often found in tight turns or
	// where backbone flexibility is structurally required.
	else if (resName == "GLY") {
		return "Small/flexible (glycine)";
	}
	else {
		return "Unknown";
	}
}

void renderComparisonPanel(StructureComparison& cmp) {
	if (!cmp.active) return;

	ImGui::SetNextWindowSize({ 480, 560 }, ImGuiCond_FirstUseEver);
	ImGui::SetNextWindowPos({ 20, 300 }, ImGuiCond_FirstUseEver);

	if (!ImGui::Begin("Structure Comparison", &cmp.active)) {
		ImGui::End();
		return;
	}

	// Header
	ImGui::TextColored({ 0.4f, 0.8f, 1.0f, 1.0f }, "%s", cmp.wtLabel.c_str());
	ImGui::SameLine();
	ImGui::TextColored({ 0.6f, 0.6f, 0.6f, 1.0f }, "vs");
	ImGui::SameLine();
	ImGui::TextColored({ 1.0f, 0.55f, 0.3f, 1.0f }, "%s", cmp.mutLabel.c_str());

	ImGui::Spacing();

	// Global RMSD badge
	ImVec4 rmsdCol = cmp.globalRmsd < 0.5f ? ImVec4(0.3f, 0.9f, 0.4f, 1.0f)
		: cmp.globalRmsd < 1.5f ? ImVec4(0.9f, 0.8f, 0.2f, 1.0f)
		: ImVec4(1.0f, 0.35f, 0.3f, 1.0f);
	ImGui::TextColored({ 0.7f, 0.7f, 0.7f, 1.0f }, "Global CA RMSD:");
	ImGui::SameLine();
	ImGui::TextColored(rmsdCol, "%.3f A", cmp.globalRmsd);
	ImGui::SameLine();
	ImGui::TextColored({ 0.5f, 0.5f, 0.5f, 1.0f }, "(%zu residues)", cmp.residues.size());

	ImGui::Spacing();
	ImGui::Separator();
	ImGui::Spacing();

	// Controls   
	ImGui::Text("Sort:");
	ImGui::SameLine();
	if (ImGui::RadioButton("By Residue",
		cmp.sortMode == StructureComparison::SortMode::BY_RESIDUE))
		cmp.sortMode = StructureComparison::SortMode::BY_RESIDUE;
	ImGui::SameLine();
	if (ImGui::RadioButton("By RMSD",
		cmp.sortMode == StructureComparison::SortMode::BY_RMSD))
		cmp.sortMode = StructureComparison::SortMode::BY_RMSD;

	ImGui::Text("Min RMSD:");
	ImGui::SameLine();
	ImGui::SetNextItemWidth(120);
	ImGui::SliderFloat("##rmsdfilter", &cmp.rmsdFilter, 0.0f, 3.0f, "%.2f A");

	ImGui::Spacing();

	// RMSD bar: inline proportional bar scaled to cmp.maxRmsd; color encodes severity
	if (ImGui::CollapsingHeader("RMSD Overview", ImGuiTreeNodeFlags_DefaultOpen)) {
		ImVec2 chartOrigin = ImGui::GetCursorScreenPos();
		float chartW = ImGui::GetContentRegionAvail().x;
		const int BAR_COUNT = std::min((int)cmp.residues.size(), 40);
		float barW = chartW / std::max(BAR_COUNT, 1);
		float chartH = 60.0f;
		float maxRmsd = 0.001f;

		for (const auto& r : cmp.residues)
			maxRmsd = std::max(maxRmsd, r.rmsd);

		ImDrawList* dl = ImGui::GetWindowDrawList();
		dl->AddRectFilled(
			chartOrigin,
			{ chartOrigin.x + chartW, chartOrigin.y + chartH },
			IM_COL32(20, 20, 30, 200), 4.0f);

		const auto& order = cmp.sortedByRmsd; // always show highest first in chart
		for (int i = 0; i < BAR_COUNT && i < (int)order.size(); i++) {
			const auto& r = cmp.residues[order[i]];
			float   h = (r.rmsd / maxRmsd) * (chartH - 4.0f);
			float   x0 = chartOrigin.x + i * barW + 1.0f;
			float   y1 = chartOrigin.y + chartH - 2.0f;
			float   y0 = y1 - h;

			ImU32 col = r.isMutationSite
				? IM_COL32(255, 120, 40, 255)
				: (r.rmsd > 1.0f ? IM_COL32(255, 80, 80, 200)
					: r.rmsd > 0.3f ? IM_COL32(255, 200, 60, 200)
					: IM_COL32(60, 180, 100, 200));

			dl->AddRectFilled({ x0, y0 }, { x0 + barW - 2.0f, y1 }, col, 2.0f);

			// Per-bar hover detection uses screen-space rect testing rather than ImGui item
			// hover, because the bars are drawn directly via ImDrawList and are not ImGui items.
			if (ImGui::IsMouseHoveringRect(
				{ x0, chartOrigin.y }, { x0 + barW, chartOrigin.y + chartH }))
			{
				ImGui::SetTooltip("%s %d%s\nRMSD: %.3f A",
					r.residueName.c_str(), r.residueNum, r.chain.c_str(), r.rmsd);
			}
		}

		// Axis label
		dl->AddText(
			{ chartOrigin.x + 2, chartOrigin.y + 2 },
			IM_COL32(150, 150, 150, 200),
			"Top affected residues (CA RMSD)");

		ImGui::Dummy({ chartW, chartH });
		ImGui::Spacing();
	}

	// Residue diff table    
	if (ImGui::CollapsingHeader("Per-Residue Detail", ImGuiTreeNodeFlags_DefaultOpen)) {
		const auto& order = (cmp.sortMode == StructureComparison::SortMode::BY_RMSD)
			? cmp.sortedByRmsd : cmp.sortedByResNum;

		ImGuiTableFlags tableFlags =
			ImGuiTableFlags_Borders |
			ImGuiTableFlags_RowBg |
			ImGuiTableFlags_ScrollY |
			ImGuiTableFlags_SizingStretchProp;

		ImGui::SetNextWindowContentSize({ 0, 0 });
		if (ImGui::BeginTable("##cmptable", 4, tableFlags, { 0, 280 })) {
			ImGui::TableSetupScrollFreeze(0, 1);
			ImGui::TableSetupColumn("Chain", ImGuiTableColumnFlags_WidthFixed, 45);
			ImGui::TableSetupColumn("Residue", ImGuiTableColumnFlags_WidthFixed, 70);
			ImGui::TableSetupColumn("Name", ImGuiTableColumnFlags_WidthFixed, 50);
			ImGui::TableSetupColumn("RMSD (A)", ImGuiTableColumnFlags_WidthStretch);
			ImGui::TableHeadersRow();

			for (int idx : order) {
				const auto& r = cmp.residues[idx];
				if (r.rmsd < cmp.rmsdFilter) continue;

				ImGui::TableNextRow();

				// Row highlight for mutation site
				if (r.isMutationSite) {
					ImGui::TableSetBgColor(ImGuiTableBgTarget_RowBg0,
						IM_COL32(80, 35, 10, 180));
				}

				// Chain
				ImGui::TableSetColumnIndex(0);
				ImGui::TextColored({ 0.7f,0.7f,0.7f,1.0f }, "%s", r.chain.c_str());

				// Residue number
				ImGui::TableSetColumnIndex(1);
				if (r.isMutationSite)
					ImGui::TextColored({ 1.0f,0.55f,0.3f,1.0f }, "%d *", r.residueNum);
				else
					ImGui::Text("%d", r.residueNum);

				// Residue name
				ImGui::TableSetColumnIndex(2);
				ImGui::TextColored({ 0.6f,0.85f,1.0f,1.0f }, "%s", r.residueName.c_str());

				// RMSD bar
				ImGui::TableSetColumnIndex(3);
				float  cellW = ImGui::GetContentRegionAvail().x;
				float barFrac = r.rmsd / cmp.maxRmsd;
				if (barFrac > 1.0f) barFrac = 1.0f;

				ImVec2 barMin = ImGui::GetCursorScreenPos();
				ImVec2 barMax = { barMin.x + cellW * barFrac, barMin.y + 13.0f };

				ImU32 barCol = r.rmsd > 1.0f ? IM_COL32(255, 80, 80, 180)
					: r.rmsd > 0.3f ? IM_COL32(255, 200, 60, 180)
					: IM_COL32(60, 180, 100, 180);
				ImGui::GetWindowDrawList()->AddRectFilled(barMin, barMax, barCol, 2.0f);
				ImGui::SetCursorPosX(ImGui::GetCursorPosX() + 4);
				ImGui::TextColored({ 1,1,1,0.9f }, "%.3f", r.rmsd);
			}
			ImGui::EndTable();
		}
	}

	ImGui::End();
}

// Convert 3-letter residue name to 1-letter amino acid code
// Returns '?' if unknown
char threeLetterToOne(const std::string& three) {
	static const std::unordered_map<std::string, char> lookup = {
		{"ALA",'A'}, {"ARG",'R'}, {"ASN",'N'}, {"ASP",'D'}, {"CYS",'C'},
		{"GLN",'Q'}, {"GLU",'E'}, {"GLY",'G'}, {"HIS",'H'}, {"ILE",'I'},
		{"LEU",'L'}, {"LYS",'K'}, {"MET",'M'}, {"PHE",'F'}, {"PRO",'P'},
		{"SER",'S'}, {"THR",'T'}, {"TRP",'W'}, {"TYR",'Y'}, {"VAL",'V'},

		// Common variants
		{"HSD",'H'}, {"HSE",'H'}, {"HSP",'H'},  // Histidine protonation states
		{"MSE",'M'},  // Selenomethionine
		{"SEC",'U'}   // Selenocysteine
	};

	std::string upper = three;
	for (auto& c : upper) c = std::toupper(c);

	auto it = lookup.find(upper);
	return (it != lookup.end()) ? it->second : '?';
}

// Format AlphaMissense class string with color/indicator for console
std::string formatAmClass(const std::string& amClass, float score) {
	if (amClass == "pathogenic")  return "[PATH  " + std::to_string(score).substr(0, 5) + "]";
	if (amClass == "benign") return "[BEN   " + std::to_string(score).substr(0, 5) + "]";
	return   "[AMB   " + std::to_string(score).substr(0, 5) + "]";
}

// Run mutation analysis in background thread — called via std::thread().detach()
void runMutationAnalysis(
	std::string pdbId,
	std::string chain,
	int pdbResidueNum,
	std::string residueName3Letter
) {
	mutationAnalysisRunning = true;

	char wtAA = threeLetterToOne(residueName3Letter);
	if (wtAA == '?') {
		std::cerr << "\n[!] Cannot analyze residue \"" << residueName3Letter << "\".\n"
			<< " This residue type is non-standard or not yet supported.\n"
			<< " Try selecting a standard amino acid residue (e.g. ALA, GLY, LYS).\n\n";
		mutationAnalysisRunning = false;
		return;
	}

	std::cout << "\n╔══════════════════════════════════════════════╗\n";
	std::cout << "║  MUTATION ANALYSIS: " << residueName3Letter << " " << pdbResidueNum
		<< " (Chain " << chain << ")\n";
	std::cout << "╚══════════════════════════════════════════════╝\n";

	HTTPConnection http;

	// STEP 1: PDB residue → UniProt position  
	std::string mappingUrl = BACKEND_URL + "/api/pdb/" + pdbId + "/chain/" + chain
		+ "/residue/" + std::to_string(pdbResidueNum);

	std::cout << "[1/3] Mapping PDB residue to UniProt coordinates...\n";

	if (!http.get(mappingUrl)) {
		std::cerr << "\n[!] Could not reach the mapping service.\n"
			<< " Make sure the DataLens backend is running:\n"
			<< " cd backend && uvicorn main:app --reload\n"
			<< " Expected at: " << BACKEND_URL << "\n\n";
		mutationAnalysisRunning = false;
		return;
	}

	std::string uniprotId;
	int uniprotPosition = -1;

	try {
		auto mappingJson = json::parse(http.getResponse());

		if (mappingJson.contains("detail")) {
			// 404 / 400 from API
			std::cout << "\n[!] No PDB-UniProt mapping found for residue " << pdbResidueNum
				<< " (chain " << chain << ") in " << pdbId << ".\n"
				<< " Detail: " << mappingJson["detail"].get<std::string>() << "\n"
				<< " The residue may be a ligand, solvent, or outside the SIFTS mapping.\n\n";
			mutationAnalysisRunning = false;
			return;
		}

		uniprotId = mappingJson.value("uniprot_id", "");
		uniprotPosition = mappingJson.value("uniprot_residue", -1);

		std::cout << " PDB " << pdbResidueNum << " → UniProt " << uniprotId
			<< " position " << uniprotPosition << "\n";
	}
	catch (const json::exception& e) {
		std::cerr << "\n[!] Failed to parse the mapping response from the backend.\n"
			<< " This may indicate a server-side error or an unexpected API change.\n"
			<< " Details: " << e.what() << "\n\n";
		mutationAnalysisRunning = false;
		return;
	}

	if (uniprotId.empty() || uniprotPosition < 0) {
		std::cerr << "\n[!] The backend returned incomplete mapping data for residue "
			<< pdbResidueNum << " (chain " << chain << ").\n"
			<< " UniProt ID or position was missing. Try a different residue.\n\n";
		mutationAnalysisRunning = false;
		return;
	}

	// STEP 2: All 19 AlphaMissense predictions at this position 
	std::string amUrl = BACKEND_URL + "/api/alphamissense/position/"
		+ uniprotId + "/" + std::to_string(uniprotPosition);

	std::cout << "[2/3] Fetching AlphaMissense predictions...\n";

	if (!http.get(amUrl)) {
		std::cerr << "\n[!] Could not reach the AlphaMissense service.\n"
			<< " Make sure the backend is running at " << BACKEND_URL << "\n"
			<< " and the AlphaMissense database is loaded (216M+ variants).\n\n";
		mutationAnalysisRunning = false;
		return;
	}

	try {
		auto amJson = json::parse(http.getResponse());

		if (amJson.contains("detail")) {
			std::cout << "\n[!] No AlphaMissense predictions found for "
				<< uniprotId << " position " << uniprotPosition << ".\n"
				<< " Detail: " << amJson["detail"].get<std::string>() << "\n"
				<< " This protein or position may not be covered by AlphaMissense v1.0.\n\n";
			mutationAnalysisRunning = false;
			return;
		}

		auto& predictions = amJson["predictions"];
		int count = amJson.value("prediction_count", 0);
		std::string refAA = amJson.value("reference_aa", std::string(1, wtAA));

		std::cout << "\n┌ AlphaMissense: " << uniprotId << " position "
			<< uniprotPosition << " (" << refAA << ")    ┐\n";
		std::cout << "│  " << count << " substitutions ranked by pathogenicity:\n";
		std::cout << "├    ┤\n";

		// Track most pathogenic for STEP 3
		std::string topMutAA = "";
		float topScore = -1.0f;

		int shown = 0;
		for (auto& pred : predictions) {
			std::string variant = pred.value("protein_variant", "?");  // e.g. K28H
			float score = pred.value("am_pathogenicity", 0.0f);
			std::string amClass = pred.value("am_class", "ambiguous");
			std::string altAA = pred.value("alternate_aa", "?");

			// Color code (console): PATH=!, AMB=~, BEN=.
			char indicator = (amClass == "pathogenic") ? '!' :
				(amClass == "benign") ? '.' : '~';

			std::cout << "│  " << indicator << " " << variant
				<< "  " << formatAmClass(amClass, score) << "\n";

			if (score > topScore) {
				topScore = score;
				topMutAA = altAA;
			}

			// Only show top 10 to avoid wall of text
			if (++shown >= 10) {
				std::cout << "│  ... and " << (count - 10) << " more\n";
				break;
			}
		}
		std::cout << "└    ┘\n";

		// Legend
		std::cout << "  Legend: ! pathogenic  ~ ambiguous  . benign\n";

		// STEP 3: Full mutation report for most pathogenic substitution  
		if (!topMutAA.empty() && topScore > 0.5f) {
			std::cout << "\n[3/3] Running full analysis for most pathogenic: "
				<< wtAA << chain << pdbResidueNum << topMutAA
				<< " (AM score: " << topScore << ")...\n";

			// /analyze_mutation requires PDB coordinates (not UniProt) for FoldX, which operates
			// directly on the structure file. pdb_id + chain + position uniquely identifies the
			// residue in 3D space; wt_aa is included as a sanity check against the backend's
			// own residue lookup.
			json requestBody = {
				{"pdb_id",   pdbId},
				{"chain", chain},
				{"position", pdbResidueNum},
				{"wt_aa", std::string(1, wtAA)},
				{"mut_aa",   topMutAA},
				{"alphamissense_score", topScore}
			};

			std::string analyzeUrl = BACKEND_URL + "/api/analyze_mutation";

			if (!http.post(analyzeUrl, requestBody.dump())) {
				std::cout << "\n[i] Full structural analysis is unavailable.\n"
					<< " FoldX may not be installed or configured on the backend.\n"
					<< " AlphaMissense predictions above are still valid.\n";
			}
			else {
				try {
					auto reportJson = json::parse(http.getResponse());

					if (reportJson.contains("detail")) {
						std::cout << " Analysis error: "
							<< reportJson["detail"].get<std::string>() << "\n";
					}
					else {
						// Print summary  
						std::cout << "\n┌ Full Mutation Report    ┐\n";

						if (reportJson.contains("summary")) {
							std::cout << "│  SUMMARY:\n";
							std::cout << "│  " << reportJson["summary"].get<std::string>() << "\n";
						}

						// Location context
						if (reportJson.contains("location_context")) {
							auto& loc = reportJson["location_context"];
							std::cout << "│\n│  LOCATION:  "
								<< loc.value("description", "unknown") << "\n";
							if (loc.contains("sasa")) {
								std::cout << "│  SASA: " << loc["sasa"].get<float>() << " Å²\n";
							}
						}

						// Secondary structure
						if (reportJson.contains("structural_impact")) {
							auto& si = reportJson["structural_impact"];
							std::cout << "│  SEC.STRUCT: "
								<< si.value("secondary_structure_description", "unknown") << "\n";
						}

						// Predictions section
						if (reportJson.contains("predictions")) {
							auto& preds = reportJson["predictions"];
							if (preds.contains("alphamissense")) {
								std::cout << "│\n│  ALPHAMISSENSE: score="
									<< preds["alphamissense"].value("score", 0.0f)
									<< "  class="
									<< preds["alphamissense"].value("classification", "?") << "\n";
							}
							if (preds.contains("foldx") && !preds["foldx"].is_null()) {
								std::cout << "│  FOLDX ΔΔG:  "
									<< preds["foldx"].value("ddg", 0.0f) << " kcal/mol\n";
							}
						}

						// Top interpretations
						if (reportJson.contains("interpretations") &&
							reportJson["interpretations"].is_array()) {
							std::cout << "│\n│  KEY FINDINGS:\n";
							int iCount = 0;
							for (auto& interp : reportJson["interpretations"]) {
								std::cout << "│   "
									<< interp.value("icon", "•") << " "
									<< "[" << interp.value("severity", "?") << "] "
									<< interp.value("text", "") << "\n";
								if (++iCount >= 5) break;  // Limit output
							}
						}

						// Load mutant structure for comparison 
						if (reportJson.contains("mutant_pdb") &&
							reportJson["mutant_pdb"].is_string() &&
							!reportJson["mutant_pdb"].get<std::string>().empty())
						{
							mutantPdbContent = reportJson["mutant_pdb"].get<std::string>();
							commands.push_back("loadmutant");
							std::cout << "\n[COMPARISON] Mutant PDB received — queued for comparison.\n";
						}
						else {
							std::cout << "\n[i] No mutant PDB in response — structure comparison unavailable.\n";
						}

						std::cout << "└    ┘\n";

						// Clinical assessment (if available)
						if (reportJson.contains("clinical_assessment")) {
							auto& ca = reportJson["clinical_assessment"];
							std::cout << "\n  CLINICAL:  "
								<< ca.value("overall_assessment", "N/A") << "\n";
							std::cout << "  CONFIDENCE: "
								<< ca.value("confidence", "N/A") << "\n";
						}
					}
				}
				catch (const json::exception& e) {
					std::cerr << "[MUTATION] Failed to parse analysis report: " << e.what() << "\n";
				}
			}
		}
		else {
			std::cout << "\n[i] Skipped full FoldX analysis — top AlphaMissense score is "
				<< topScore << " (threshold: 0.5).\n"
				<< " All substitutions at this position appear benign or ambiguous.\n\n";
		}
	}
	catch (const json::exception& e) {
		std::cerr << "\n[!] Failed to parse the AlphaMissense response from the backend.\n"
			<< " The server may have returned an unexpected format.\n"
			<< " Details: " << e.what() << "\n\n";
	}

	std::cout << "\n[DONE] Press M to show menu again.\n\n";
	mutationAnalysisRunning = false;
}

void fetchAllMutationsForPDB(
	const std::string& pdbId,
	const std::map<char, std::vector<int>>& chainResidues  // chain -> pdb residue list
) {
	mutationPanel.loading = true;
	mutationPanel.loadedCount = 0;
	mutationPanel.totalCount = 0;
	{
		std::lock_guard<std::mutex> lock(mutationPanel.dataMutex);
		mutationPanel.chains.clear();
		mutationPanel.statusMsg = "Fetching chain mappings...";
	}

	HTTPConnection http;

	// Step 1: which chains have AM data  
	std::string chainsUrl = BACKEND_URL + "/api/alphamissense/pdb/"
		+ pdbId + "/analyzable-chains";

	if (!http.get(chainsUrl)) {
		std::lock_guard<std::mutex> lock(mutationPanel.dataMutex);
		mutationPanel.statusMsg = "Could not reach backend. Is it running at " + BACKEND_URL + "?";
		mutationPanel.loading = false;
		return;
	}

	// Map: chain char -> uniprot_id, pdb_start, uniprot_start
	struct ChainInfo {
		std::string uniprotId;
		int pdbStart;
		int uniprotStart;
	};
	std::map<char, ChainInfo> analyzableChains;

	try {
		auto j = json::parse(http.getResponse());
		if (j.contains("detail") || !j.contains("chains")) {
			std::lock_guard<std::mutex> lock(mutationPanel.dataMutex);
			mutationPanel.statusMsg = "No AlphaMissense data available for this structure.";
			mutationPanel.loading = false;
			return;
		}
		for (auto& ch : j["chains"]) {
			std::string chainStr = ch.value("chain_id", "");
			if (chainStr.empty()) continue;
			char c = chainStr[0];
			analyzableChains[c] = {
				ch.value("uniprot_id", ""),
				ch.value("pdb_start", 0),
				ch.value("uniprot_start", 0)
			};
		}
	}
	catch (const json::exception& e) {
		std::lock_guard<std::mutex> lock(mutationPanel.dataMutex);
		mutationPanel.statusMsg = std::string("Parse error: ") + e.what();
		mutationPanel.loading = false;
		return;
	}

	// Count total residues to fetch
	int total = 0;
	for (auto& [chain, info] : analyzableChains) {
		auto it = chainResidues.find(chain);
		if (it != chainResidues.end()) total += (int)it->second.size();
	}
	mutationPanel.totalCount = total;

	{
		std::lock_guard<std::mutex> lock(mutationPanel.dataMutex);
		mutationPanel.statusMsg = "Loading 0 / " + std::to_string(total) + " residues...";
	}

	// Step 2: per-residue worst score    
	for (auto& [chain, chainInfo] : analyzableChains) {
		auto it = chainResidues.find(chain);
		if (it == chainResidues.end()) continue;

		ChainAmSummary chainSummary;
		chainSummary.uniprotId = chainInfo.uniprotId;

		for (int pdbRes : it->second) {

			// SIFTS maps PDB residue numbers to UniProt positions using a per-chain linear
			// offset. The offset is not always zero — PDB files often start numbering at a
			// non-1 value (e.g. after signal peptide cleavage or engineered constructs).
			// AlphaMissense is indexed by UniProt position, so this conversion is mandatory.
			int uniprotPos = chainInfo.uniprotStart + (pdbRes - chainInfo.pdbStart);

			std::string posUrl = BACKEND_URL + "/api/alphamissense/position/"
				+ chainInfo.uniprotId + "/" + std::to_string(uniprotPos);

			ResidueAmSummary summary;
			summary.pdbResidueNum = pdbRes;
			summary.uniprotPos = uniprotPos;

			if (http.get(posUrl)) {
				try {
					auto j = json::parse(http.getResponse());
					if (!j.contains("detail") && j.contains("predictions")) {
						for (auto& pred : j["predictions"]) {
							float score = pred.value("am_pathogenicity", 0.0f);
							if (score > summary.worstScore) {
								summary.worstScore = score;
								summary.worstClass = pred.value("am_class", "");
								summary.worstVariant = pred.value("protein_variant", "");
							}
						}
					}
				}
				catch (...) {}
			}

			{
				std::lock_guard<std::mutex> lock(mutationPanel.dataMutex);
				chainSummary.residues.push_back(summary);
				mutationPanel.loadedCount++;
				mutationPanel.statusMsg = "Loading " + std::to_string(mutationPanel.loadedCount.load())
					+ " / " + std::to_string(total) + " residues...";
			}
		}

		chainSummary.loaded = true;
		{
			std::lock_guard<std::mutex> lock(mutationPanel.dataMutex);
			mutationPanel.chains[chain] = std::move(chainSummary);
		}
	}

	{
		std::lock_guard<std::mutex> lock(mutationPanel.dataMutex);
		mutationPanel.statusMsg = "Loaded " + std::to_string(total) + " residues.";
		mutationPanel.showWindow = true;
	}
	mutationPanel.loading = false;
}


// In your ImGui panel render function
void renderVariantSelector(
	const json& variantData, // AM data
	const std::set<int>& positions,   // AM positions
	int& selectedPosition,
	std::string& selectedVariantId,
	bool& runAnalysis,
	const json& allChains,   // from /api/pdb/{id}/all-residues
	ManualMutState& manualMut
) {
	ImGui::Begin("Variant Selector");

	float totalHeight = ImGui::GetContentRegionAvail().y;
	float sectionHeight = (totalHeight / 2.0f) - 20.0f;

	// ═══════════════════════════════════════════════════════
	// SECTION 1 — AlphaMissense Predictions
	// ═══════════════════════════════════════════════════════
	if (ImGui::CollapsingHeader(
		"AlphaMissense Predictions", ImGuiTreeNodeFlags_DefaultOpen))
	{
		static char searchBuf[32] = "";
		ImGui::SetNextItemWidth(-1);
		ImGui::InputText("Search pos##am", searchBuf, sizeof(searchBuf));
		ImGui::Separator();

		if (variantData.empty()) {
			ImGui::TextDisabled("No AlphaMissense data for this structure.");
		}
		else {
			// Left: position list
			ImGui::BeginChild("##am_positions", ImVec2(140, sectionHeight), true);
			for (int pos : positions) {
				std::string label = std::to_string(pos);
				if (strlen(searchBuf) > 0 &&
					label.find(searchBuf) == std::string::npos) continue;

				bool isSel = (selectedPosition == pos);
				if (ImGui::Selectable(label.c_str(), isSel)) {
					selectedPosition = pos;
					selectedVariantId = "";
				}
			}
			ImGui::EndChild();

			ImGui::SameLine();

			// Renders the per-position AlphaMissense substitution table and the manual mutation
			// input panel. Split into two collapsing sections: Handler 1 (ray-cast click from
			// the 3D view populates variantData) and Handler 2 (manual residue + AA input).
			// allChains provides the chain/residue dropdown data fetched on structure load.
			ImGui::BeginChild("##am_subs", ImVec2(0, sectionHeight), true);
			if (selectedPosition < 0) {
				ImGui::TextDisabled("<- Select a position");
			}
			else {
				if (ImGui::BeginTable("##am_vars", 5,
					ImGuiTableFlags_Borders |
					ImGuiTableFlags_RowBg |
					ImGuiTableFlags_ScrollY |
					ImGuiTableFlags_SizingFixedFit,
					ImVec2(0, sectionHeight - 8.0f)))
				{
					ImGui::TableSetupScrollFreeze(0, 1);
					ImGui::TableSetupColumn("Variant", ImGuiTableColumnFlags_WidthFixed, 75);
					ImGui::TableSetupColumn("Sub", ImGuiTableColumnFlags_WidthFixed, 50);
					ImGui::TableSetupColumn("Score", ImGuiTableColumnFlags_WidthStretch);
					ImGui::TableSetupColumn("Class", ImGuiTableColumnFlags_WidthFixed, 75);
					ImGui::TableSetupColumn("FoldX", ImGuiTableColumnFlags_WidthFixed, 42);
					ImGui::TableHeadersRow();

					for (auto& v : variantData) {
						if (v.value("position", 0) != selectedPosition) continue;

						std::string varId = v.value("variant_id", "");
						std::string ref = v.value("ref_aa", "");
						std::string alt = v.value("alt_aa", "");
						float  score = v.value("pathogenicity_score", 0.0f);
						std::string cls = v.value("classification", "ambiguous");
						bool   foldxReady = v.value("foldx_ready", false);

						ImVec4 color = (cls == "pathogenic")
							? ImVec4(0.9f, 0.2f, 0.2f, 1.f)
							: (cls == "benign")
							? ImVec4(0.2f, 0.8f, 0.2f, 1.f)
							: ImVec4(0.9f, 0.7f, 0.1f, 1.f);

						ImGui::TableNextRow();

						ImGui::TableSetColumnIndex(0);
						bool rowSel = (selectedVariantId == varId);
						if (ImGui::Selectable(varId.c_str(), rowSel,
							ImGuiSelectableFlags_SpanAllColumns))
							selectedVariantId = varId;

						ImGui::TableSetColumnIndex(1);
						ImGui::Text("%s>%s", ref.c_str(), alt.c_str());

						ImGui::TableSetColumnIndex(2);
						ImGui::PushStyleColor(ImGuiCol_PlotHistogram, color);
						ImGui::ProgressBar(score, ImVec2(-1, 10),
							std::to_string(score).substr(0, 5).c_str());
						ImGui::PopStyleColor();

						ImGui::TableSetColumnIndex(3);
						ImGui::TextColored(color, "%s", cls.c_str());

						ImGui::TableSetColumnIndex(4);
						if (foldxReady)
							ImGui::TextColored(ImVec4(0.2f, 0.9f, 0.4f, 1.f), "YES");
						else
							ImGui::TextColored(ImVec4(0.4f, 0.4f, 0.4f, 1.f), "--");
					}
					ImGui::EndTable();
				}
			}
			ImGui::EndChild();

			// Analyze button for Section 1
			if (!selectedVariantId.empty()) {
				ImGui::Separator();
				ImGui::Text("Selected: %s", selectedVariantId.c_str());
				ImGui::SameLine();

				bool canAnalyze = false;
				for (auto& v : variantData) {
					if (v.value("variant_id", "") == selectedVariantId) {
						canAnalyze = v.value("foldx_ready", false);
						break;
					}
				}

				if (canAnalyze) {
					if (ImGui::Button("Analyze -> ##am"))
						runAnalysis = true;
				}
				else {
					ImGui::BeginDisabled();
					ImGui::Button("Analyze ->##am");
					ImGui::EndDisabled();
					ImGui::SameLine();
					ImGui::TextColored(
						ImVec4(0.5f, 0.5f, 0.5f, 1.f), "(no PDB mapping found)");
				}
			}
		}
	}

	ImGui::Spacing();

	// ═══════════════════════════════════════════════════════
	// SECTION 2 — All Chains (FoldX Manual Mutation)
	// ═══════════════════════════════════════════════════════
	if (ImGui::CollapsingHeader(
		"All Chains - FoldX Manual", ImGuiTreeNodeFlags_DefaultOpen))
	{
		if (allChains.empty()) {
			ImGui::TextDisabled("No chain data loaded.");
		}
		else {
			// Build chain labels for combo
			std::vector<std::string> chainLabels;
			for (auto& ch : allChains) {
				std::string label =
					"Chain " + ch.value("chain_id", "?") +
					" (" + ch.value("uniprot_id", "?") + ")" +
					"  pos " + std::to_string(ch.value("pdb_start", 0)) +
					"-" + std::to_string(ch.value("pdb_end", 0));
				chainLabels.push_back(label);
			}

			// Chain selector
			ImGui::Text("Chain:");
			ImGui::SetNextItemWidth(-1);
			if (ImGui::BeginCombo("##chain_sel",
				chainLabels[manualMut.selectedChainIdx].c_str()))
			{
				for (int i = 0; i < (int)chainLabels.size(); ++i) {
					bool isSel = (manualMut.selectedChainIdx == i);
					if (ImGui::Selectable(chainLabels[i].c_str(), isSel))
						manualMut.selectedChainIdx = i;
					if (isSel) ImGui::SetItemDefaultFocus();
				}
				ImGui::EndCombo();
			}

			// Show residue range for selected chain
			auto& selChain = allChains[manualMut.selectedChainIdx];
			int pdbStart = selChain.value("pdb_start", 0);
			int pdbEnd = selChain.value("pdb_end", 0);
			ImGui::TextDisabled("Residue range: %d - %d", pdbStart, pdbEnd);

			ImGui::Spacing();

			// Position input
			ImGui::Text("Position (PDB):");
			ImGui::SetNextItemWidth(100);
			ImGui::InputInt("##manual_pos", &manualMut.position, 1);
			if (manualMut.position < pdbStart) manualMut.position = pdbStart;
			if (manualMut.position > pdbEnd)   manualMut.position = pdbEnd;

			// Mutant AA selector
			ImGui::Text("Mutate to:");
			static const char* AA_LIST =
				"A\0R\0N\0D\0C\0Q\0E\0G\0H\0I\0"
				"L\0K\0M\0F\0P\0S\0T\0W\0Y\0V\0";
			static int mutAAIdx = 0;
			ImGui::SetNextItemWidth(80);
			if (ImGui::Combo("##mut_aa", &mutAAIdx, AA_LIST)) {
				const char* aas = "ARNDCQEGHILKMFPSTWYV";
				manualMut.mutAA = aas[mutAAIdx];
			}

			ImGui::Spacing();

			// Mutation preview
			ImGui::TextColored(ImVec4(0.8f, 0.8f, 0.2f, 1.f),
				"Mutation: ?%s%d%c",
				selChain.value("chain_id", "?").c_str(),
				manualMut.position,
				manualMut.mutAA);
			ImGui::TextDisabled("(wt_aa resolved from PDB at runtime)");

			ImGui::Spacing();

			if (ImGui::Button("Run FoldX ->##manual", ImVec2(-1, 0)))
				manualMut.runManual = true;
		}
	}

	ImGui::End();
}

void applyMutationHighlight(
	MoleculeData* wtMol,
	MoleculeData* mutMol,
	const std::string& chainId,
	int position,
	bool showingMutant
) {
	if (!wtMol || !mutMol) return;
	char chain = chainId.empty() ? 'A' : chainId[0];

	// Count atoms to determine which residue is bigger
	int wtAtomCount = 0, mutAtomCount = 0;
	std::string wtResName, mutResName;

	for (const auto& atom : wtMol->atoms) {
		if (atom.residueNum == position && atom.chain == chain) {
			wtAtomCount++;
			wtResName = atom.residueName;
		}
	}
	for (const auto& atom : mutMol->atoms) {
		if (atom.residueNum == position && atom.chain == chain) {
			mutAtomCount++;
			mutResName = atom.residueName;
		}
	}

	Selection mutSel;
	mutSel.residue = position;
	mutSel.chain = chain;

	if (showingMutant) {
		// Viewing mutant — color mutation site magenta
		// If mutant is SMALLER (like ALA vs LYS), make it larger so visible
		Color magenta = Color::fromByte(255, 0, 255);
		Model::setAtomColor(&magenta, &mutSel);
		float radius = (mutAtomCount < wtAtomCount) ? 0.45f : 0.30f;
		Model::setAtomRadius(radius, &mutSel);
	}
	else {
		// Viewing wildtype — color mutation site orange
		// If WT is LARGER (like LYS vs ALA), make sidechain prominent
		Color orange = Color::fromByte(255, 140, 0);
		Model::setAtomColor(&orange, &mutSel);
		float radius = (wtAtomCount > mutAtomCount) ? 0.40f : 0.30f;
		Model::setAtomRadius(radius, &mutSel);
	}

	std::cout << "[COMPARE] " << wtResName << "(" << wtAtomCount << " atoms)"
		<< " -> " << mutResName << "(" << mutAtomCount << " atoms)"
		<< " at pos " << position << " chain " << chainId << "\n";
}

void renderMutationComparison(
	MoleculeData* wtMol,
	MoleculeData* mutMol,
	const std::string& chainId,
	int position
) {
	if (!wtMol || !mutMol) return;

	char chain = chainId.empty() ? 'A' : chainId[0];

	// Collect atom names at this position from both structures
	std::set<std::string> wtAtomNames, mutAtomNames;

	for (const auto& atom : wtMol->atoms) {
		if (atom.residueNum == position && atom.chain == chain)
			wtAtomNames.insert(atom.name);
	}
	for (const auto& atom : mutMol->atoms) {
		if (atom.residueNum == position && atom.chain == chain)
			mutAtomNames.insert(atom.name);
	}

	// Backbone atoms — present in both
	std::set<std::string> backbone = { "N", "CA", "C", "O" };

	// Lost atoms  = in WT but not in mutant (e.g. CG,CD,CE,NZ from LYS)
	std::set<std::string> lostAtoms, gainedAtoms, sharedSidechain;

	for (auto& name : wtAtomNames) {
		if (backbone.count(name)) continue;
		if (!mutAtomNames.count(name))
			lostAtoms.insert(name);   // in WT only
		else
			sharedSidechain.insert(name);  // in both (e.g. CB)
	}
	for (auto& name : mutAtomNames) {
		if (backbone.count(name)) continue;
		if (!wtAtomNames.count(name))
			gainedAtoms.insert(name); // in mutant only (new atoms)
	}

	// Color the WT structure at mutation site    
	// Backbone → cyan
	Color cyan = Color::fromByte(0, 255, 255);
	// Shared sidechain (e.g. CB) → white
	Color white = Color::fromByte(220, 220, 220);
	// Lost sidechain → red + enlarged (K's long chain)
	Color red = Color::fromByte(255, 50, 50);

	for (const auto& atom : wtMol->atoms) {
		if (atom.residueNum != position || atom.chain != chain) continue;

		Selection atomSel;
		atomSel.residue = position;
		atomSel.chain = chain;
		// We need per-atom coloring — use element/name selection

		if (backbone.count(atom.name)) {
			// backbone — cyan
		}
		else if (lostAtoms.count(atom.name)) {
			// lost sidechain — red enlarged
		}
		else {
			// shared sidechain — white
		}
	}

	// Apply colors using residue selection (your existing API):
	Selection resSel;
	resSel.residue = position;
	resSel.chain = chain;

	// Step 1: Color entire residue cyan first (backbone color)
	Model::setAtomColor(&cyan, &resSel);
	Model::setAtomRadius(0.22f, &resSel);

	// Step 2: Now color sidechain atoms
	// Since your Selection supports element filter, use it:
	// For lost atoms (K sidechain) — red + bigger
	// Note: your system colors by element, so color by iterating
	// The cleanest approach with your existing API:

	std::cout << "[COMPARE] WT residue: " << wtAtomNames.size() << " atoms\n";
	std::cout << "[COMPARE] Mutant residue: " << mutAtomNames.size() << " atoms\n";
	std::cout << "[COMPARE] Lost atoms: ";
	for (auto& a : lostAtoms) std::cout << a << " ";
	std::cout << "\n";
	std::cout << "[COMPARE] Shared sidechain: ";
	for (auto& a : sharedSidechain) std::cout << a << " ";
	std::cout << "\n";
}

void renderMutationFocusPanel(MutationFocusData& f) {
	if (!f.show) return;

	ImGui::SetNextWindowPos(ImVec2(10.0f, 230.0f), ImGuiCond_FirstUseEver);
	ImGui::SetNextWindowSize(ImVec2(280.0f, 320.0f), ImGuiCond_FirstUseEver);
	ImGui::Begin("Mutation Focus", &f.show);

	// Header   
	ImGui::TextColored(ImVec4(1.0f, 0.4f, 0.8f, 1.0f),
		"%s", f.variantId.c_str());
	ImGui::SameLine();
	ImGui::TextDisabled("Chain %s · Pos %d", f.chain.c_str(), f.position);
	ImGui::Separator();

	// Amino acid change   
	ImGui::Spacing();
	ImGui::Text("Substitution");
	ImGui::Spacing();

	// WT box
	ImGui::PushStyleColor(ImGuiCol_ChildBg, ImVec4(0.1f, 0.2f, 0.1f, 1.0f));
	ImGui::BeginChild("##wt", ImVec2(115, 60), true);
	ImGui::TextColored(ImVec4(0.4f, 1.0f, 0.4f, 1.0f), "WT: %s", f.wtAA.c_str());
	ImGui::TextWrapped("%s", f.wtProps.c_str());
	ImGui::EndChild();
	ImGui::PopStyleColor();

	ImGui::SameLine();
	ImGui::SetCursorPosY(ImGui::GetCursorPosY() + 20);
	ImGui::TextColored(ImVec4(1.0f, 0.8f, 0.0f, 1.0f), " ->");
	ImGui::SameLine();
	ImGui::SetCursorPosY(ImGui::GetCursorPosY() - 20);

	// Mutant box
	ImGui::PushStyleColor(ImGuiCol_ChildBg, ImVec4(0.2f, 0.1f, 0.1f, 1.0f));
	ImGui::BeginChild("##mut", ImVec2(115, 60), true);
	ImGui::TextColored(ImVec4(1.0f, 0.4f, 0.4f, 1.0f), "MUT: %s", f.mutAA.c_str());
	ImGui::TextWrapped("%s", f.mutProps.c_str());
	ImGui::EndChild();
	ImGui::PopStyleColor();

	ImGui::Spacing();
	ImGui::Separator();

	// AlphaMissense  
	ImGui::Spacing();
	ImGui::Text("AlphaMissense");
	ImVec4 amColor = (f.amClass == "pathogenic")
		? ImVec4(0.9f, 0.2f, 0.2f, 1.f)
		: (f.amClass == "benign")
		? ImVec4(0.2f, 0.9f, 0.2f, 1.f)
		: ImVec4(0.9f, 0.7f, 0.1f, 1.f);
	ImGui::PushStyleColor(ImGuiCol_PlotHistogram, amColor);
	ImGui::ProgressBar(f.amScore, ImVec2(-1, 12),
		std::to_string(f.amScore).substr(0, 5).c_str());
	ImGui::PopStyleColor();
	ImGui::TextColored(amColor, "%s", f.amClass.c_str());

	ImGui::Spacing();
	ImGui::Separator();

	// FoldX ΔΔG   
	ImGui::Spacing();
	ImGui::Text("FoldX Stability");
	if (f.analyzing) {
		ImGui::TextColored(ImVec4(1.0f, 1.0f, 0.4f, 1.0f), "Calculating...");
	}
	else {
		ImVec4 ddgColor = (f.ddg > 2.0f) ? ImVec4(0.9f, 0.2f, 0.2f, 1.f) :
			(f.ddg > 0.5f) ? ImVec4(0.9f, 0.6f, 0.2f, 1.f) :
			(f.ddg > -0.5f) ? ImVec4(0.8f, 0.8f, 0.8f, 1.f) :
			ImVec4(0.2f, 0.9f, 0.4f, 1.f);
		ImGui::TextColored(ddgColor, "DDG: %+.3f kcal/mol", f.ddg);
		ImGui::TextColored(ddgColor, "%s", f.ddgInterp.c_str());
	}

	ImGui::Spacing();
	ImGui::Separator();

	// Key difference    
	ImGui::Spacing();
	ImGui::TextDisabled("Key structural change:");

	// Highlight what was lost/gained
	if (f.wtProps != f.mutProps) {
		ImGui::TextColored(ImVec4(1.0f, 0.6f, 0.2f, 1.0f),
			"%s -> %s", f.wtProps.c_str(), f.mutProps.c_str());
	}
	// Specific well-known changes
	if (f.wtAA == "S" || f.wtAA == "SER") {
		ImGui::TextWrapped("Loss of hydroxyl (-OH)\nremoves H-bond donor");
	}
	else if (f.mutAA == "P" || f.mutAA == "PRO") {
		ImGui::TextWrapped("Proline introduced:\nrestricts backbone flexibility");
	}
	else if (f.mutAA == "G" || f.mutAA == "GLY") {
		ImGui::TextWrapped("Glycine introduced:\nincreases local flexibility");
	}

	ImGui::End();
}

/*
 * Runs the graphics rendering loop for 3D visualization
 */
void renderingThread() {
	// In renderingThread(), add these variables:
	std::optional<size_t> hoveredAtomIndex;
	std::optional<size_t> selectedAtomIndex;
	bool showContextMenu = false;
	ImVec2 contextMenuPos;

	static ImVec2 colorLegendPos;
	static ImVec2 colorLegendSize;

	// MutationFocusData mutFocus;
	MoleculeData* mutantMoleculeData = nullptr;   // ← add
	bool showingMutant = false; // ← add

	// Initialize GLFW
	if (!ResourceManager::initGLFW(3, 3)) {
		std::cin.get();
		return;
	}

	// Create main visualization window
	Window window("Datalens - Protein Visualization", 2000, 1330, false);
	mainWindowPtr = window.get(); // Store reference for chat window

	// Initialize GLEW
	if (!ResourceManager::initGLEW()) {
		std::cin.get();
		return;
	}

	// Initialize OpenGL
	ResourceManager::initOpenGL();
	//RepresentationRenderer representationRenderer;

	// Load default shaders
	Shader::loadDefaultShaders();

	// Load representation shaders
	representationRenderer.initialize();
	
	ribbonShader = new Shader(ribbonVS, ribbonFS, true);
	surfaceShader = new Shader(surfaceVS, surfaceFS, true);

	if (ribbonShader) {
		std::cout << "Ribbon shader created, ID: " << ribbonShader->getID() << std::endl;
	}
	else {
		std::cout << "Ribbon shader FAILED" << std::endl;
	}

	if (surfaceShader) {
		std::cout << "Surface shader created, ID: " << surfaceShader->getID() << std::endl;
	}
	else {
		std::cout << "Surface shader FAILED" << std::endl;
	}

	// Prepare Model
	SphereTemplate sphereTemplate;
	Model::setSphereTemplate(&sphereTemplate);

	ConnectorTemplate connectorTemplate;
	Model::setConnectorTemplate(&connectorTemplate);

	Model::genBuffers();
	Model::fillSphereTemplateBuffer();

	// Initialize ImGui - CORRECT ORDER
	IMGUI_CHECKVERSION();
	ImGui::CreateContext();

	ImGuiContext* mainImGuiContext = ImGui::GetCurrentContext();  // Store main context

	std::cout << "ImGui Version: " << IMGUI_VERSION << std::endl;

	ImGuiIO& io = ImGui::GetIO(); (void)io;

	io.DisplaySize = ImVec2((float)window.getWidth(), (float)window.getHeight());

	ImGui::StyleColorsDark();

	// Initialize platform/renderer backends IN CORRECT ORDER
	ImGui_ImplGlfw_InitForOpenGL(window.get(), true);
	ImGui_ImplOpenGL3_Init("#version 330");
	Input::initCallbacks(&window);   // ← add this, AFTER ImGui init


	std::cout << "ImGui backends initialized" << std::endl;

	// Store context globally for chat window to use
	g_mainImGuiContext = mainImGuiContext;
	addChatMessage("Assistant", "Hello! I'm your AI assistant. Ask me about protein structures!");

	MoleculeData* moleculeData = nullptr;
	Camera camera(Vec3(0.0f, 0.0f, 15.0f));
	Selection selection;

	// Claude handler variables
	ClaudeHandler claudeHandler;
	claudeHandler.loadConfig("config.json");

	if (!claudeHandler.loadConfig("config.json")) {
		claudeHandler.setApiKey(config.aiApiKey);  // fallback to already-loaded key
	}

	bool claudePending = false;

	json variantData = json::array();
	std::set<int> variantPositions;
	int selectedPosition = -1;
	std::string selectedVariantId = "";
	bool runAnalysis = false;
	bool showVariantPanel = false;
	MutationFocusData mutFocus;

	double prevMouseX = 0.0;
	double prevMouseY = 0.0;
	bool prevMousePosSet = false;

	// Variables for click detection
	bool prevLeftMouseButtonPressed = false;
	double clickStartX = 0.0;
	double clickStartY = 0.0;
	const double DRAG_THRESHOLD = 5.0; // pixels - distinguish click from drag

	// Build font atlas for main context by rendering one dummy frame
	// This must happen AFTER all OpenGL initialization is complete
	//std::cout << "Building font atlas for main window..." << std::endl;

	// Ensure we're using main context
	ImGui::SetCurrentContext(mainImGuiContext);

	// Render a minimal frame to trigger font atlas build
	glfwMakeContextCurrent(mainWindowPtr);

	ImGuiIO& initIO = ImGui::GetIO();

	// IMPORTANT: Reset frame state before entering main loop
	// This ensures any residual state from font atlas building is cleared
	// std::cout << "Resetting frame state before main loop..." << std::endl;
	glfwMakeContextCurrent(mainWindowPtr);
	ImGui::SetCurrentContext(mainImGuiContext);


	// Main rendering loop
	while (!window.shouldClose() && !shouldExit) {

		ImGui::SetCurrentContext(mainImGuiContext);

		// Handle window input
		Input::pollInput(&window);

		if (Input::keyPressed(&window, Key::ESCAPE) && selectedAtomIndex && moleculeData) {
			const Atom& atom = moleculeData->atoms[*selectedAtomIndex];
			Selection clearSel;
			clearSel.residue = atom.residueNum;
			clearSel.chain = atom.chain;
			Model::setAtomRadius(SphereTemplate::DEFAULT_RADIUS, &clearSel);
			Model::colorAtomsDefault(&clearSel);

			selectedAtomIndex = std::nullopt;
			selection.reset();
			showContextMenu = false;

			std::cout << "Selection cleared (ESC)" << std::endl;
		}

		camera.zoom((float)Input::getMouseScrollOffsetY());

		double mouseX = Input::getMouseX();
		double mouseY = Input::getMouseY();

		// Check current mouse button state
		bool leftMouseButtonPressed = Input::mouseButtonPressed(&window, MouseButton::LEFT);

		static bool prevRightMouseButtonPressed = false;
		bool rightMouseButtonPressed = Input::mouseButtonPressed(&window, MouseButton::RIGHT);

		//      
		// M-CLICK DETECTION (for context menu)
		//      
		static bool prevMKeyPressed = false;
		bool mKeyPressed = Input::keyPressed(&window, Key::M);

		// Detect right-click release
		if (mKeyPressed && !prevMKeyPressed) {
			if (selectedAtomIndex && moleculeData) {
				showContextMenu = true;
				// std::cout << "ColorLegendPos X: " << colorLegendPos.x << std::endl;
				// std::cout << "ColorLegendPos Y: " << colorLegendPos.y << std::endl;

				contextMenuPos = ImVec2(
					colorLegendPos.x,
					colorLegendPos.y + colorLegendSize.y + 5.0f  // 5px gap
				);
			}
		}

		// Handle keyboard shortcuts when context menu is shown
		if (showContextMenu) {
			static bool prevAKeyPressed = false;
			static bool prevSKeyPressed = false;
			static bool prevDKeyPressed = false;

			bool aKeyPressed = Input::keyPressed(&window, Key::A);
			bool sKeyPressed = Input::keyPressed(&window, Key::S);
			bool dKeyPressed = Input::keyPressed(&window, Key::D);

			// A - Mutation Analysis
			if (aKeyPressed && !prevAKeyPressed) {
				if (selectedAtomIndex && moleculeData) {
					const Atom& atom = moleculeData->atoms[*selectedAtomIndex];

					char refAA = threeLetterToOne(atom.residueName);
					if (refAA == '?') {
						std::cerr << "\n[!] Cannot analyze residue \"" << atom.residueName << "\" — non-standard type.\n"
							<< " Select a standard amino acid residue to run mutation analysis.\n\n";
					}
					else if (currentPdbId.empty()) {
						std::cerr << "\n[!] No PDB ID is associated with the loaded structure.\n"
							<< " Please reload the structure with: fetch pdb <pdb_id>\n\n";
					}
					else if (mutationAnalysisRunning) {
						std::cout << "\n[i] Mutation analysis is already in progress.\n"
							<< " Please wait for it to complete before selecting another residue.\n\n";
					}
					else {
						std::string capPdbId = currentPdbId;
						std::string capChain = std::string(1, atom.chain);
						int capResNum = atom.residueNum;
						char   capRefAA = refAA;

						std::cout << "[MUTATION] Fetching variants for " << atom.residueName
							<< atom.residueNum << " Chain " << atom.chain << "...\n";

						bool* showVariantPanelPtr = &showVariantPanel;


						std::thread([=, &variantData, &variantPositions, &selectedPosition, &selectedVariantId, &showVariantPanel]() {
							mutationAnalysisRunning = true;
							HTTPConnection http;

							// Step 1: PDB residue → UniProt ID + position  
							std::string mappingUrl = BACKEND_URL + "/api/pdb/" + capPdbId
								+ "/chain/" + capChain
								+ "/residue/" + std::to_string(capResNum);

							if (!http.get(mappingUrl)) {
								std::cerr << "[MUTATION] Failed to reach mapping API.\n";
								mutationAnalysisRunning = false;
								return;
							}

							std::string uniprotId;
							int uniprotPos = -1;

							try {
								auto mappingJson = json::parse(http.getResponse());

								if (mappingJson.contains("detail")) {
									std::cout << "[MUTATION] No mapping: "
										<< mappingJson["detail"].get<std::string>() << "\n";
									mutationAnalysisRunning = false;
									return;
								}

								uniprotId = mappingJson.value("uniprot_id", "");
								uniprotPos = mappingJson.value("uniprot_residue", -1);

								std::cout << "[1/2] PDB " << capResNum << " → "
									<< uniprotId << " pos " << uniprotPos << "\n";
							}
							catch (const json::exception& e) {
								std::cerr << "[MUTATION] Mapping parse error: " << e.what() << "\n";
								mutationAnalysisRunning = false;
								return;
							}

							if (uniprotId.empty() || uniprotPos < 0) {
								std::cerr << "[MUTATION] Invalid mapping data.\n";
								mutationAnalysisRunning = false;
								return;
							}

							// Step 2: All 19 substitutions from PostgreSQL   
							std::string amUrl = BACKEND_URL + "/api/alphamissense/position/"
								+ uniprotId + "/" + std::to_string(uniprotPos);

							if (!http.get(amUrl)) {
								std::cerr << "[MUTATION] Failed to reach AlphaMissense API.\n";
								mutationAnalysisRunning = false;
								return;
							}

							try {
								auto amJson = json::parse(http.getResponse());

								if (amJson.contains("detail")) {
									std::cout << "[MUTATION] No AM data: "
										<< amJson["detail"].get<std::string>() << "\n";
									mutationAnalysisRunning = false;
									return;
								}

								int count = amJson.value("prediction_count", 0);

								std::cout << "\n╔══════════════════════════════════════════════════╗\n";
								std::printf("║  %c%d (Chain %s)  |  %s pos %d\n",
									capRefAA, capResNum, capChain.c_str(),
									uniprotId.c_str(), uniprotPos);
								std::cout << "╠══════════════════════════════════════════════════╣\n";

								int pathCount = 0, ambCount = 0, benCount = 0;

								// Reset panel state
								variantData = json::array();
								variantPositions.clear();
								selectedPosition = -1;
								selectedVariantId = "";

								for (auto& pred : amJson["predictions"]) {
									std::string variant = pred.value("protein_variant", "?");
									float  score = pred.value("am_pathogenicity", 0.0f);
									std::string amClass = pred.value("am_class", "ambiguous");
									std::string altAA = pred.value("alternate_aa", "?");

									char indicator = (amClass == "pathogenic") ? '!' :
										(amClass == "benign") ? '.' : '~';

									std::printf("║  %c  %-8s  %.4f  %s\n",
										indicator, variant.c_str(), score, amClass.c_str());

									if (indicator == '!') pathCount++;
									else if (indicator == '~') ambCount++;
									else  benCount++;

									variantData.push_back({
										{"variant_id",  variant},
										{"position",  uniprotPos},
										{"ref_aa", std::string(1, capRefAA)},
										{"alt_aa", altAA},
										{"pathogenicity_score", score},
										{"classification", amClass},
										{ "foldx_ready", true }   // A-key only fires for mapped chains

										});
									variantPositions.insert(uniprotPos);
								}

								std::cout << "╠══════════════════════════════════════════════════╣\n";
								std::printf("║  ! pathogenic: %d   ~ ambiguous: %d   . benign: %d\n",
									pathCount, ambCount, benCount);
								std::cout << "╚══════════════════════════════════════════════════╝\n";
								std::cout << "  → Select a variant in the panel, then click Analyze\n\n";

								// Open the variant selector panel
								showVariantPanel = true;
							}
							catch (const json::exception& e) {
								std::cerr << "[MUTATION] AM parse error: " << e.what() << "\n";
								mutationAnalysisRunning = false;
								return;
							}

							// Step 3 removed — FoldX now fires from "Analyze ->" in the variant panel

							mutationAnalysisRunning = false;
							}).detach();
					}
				}
				else {
					std::cout << "[MUTATION] Select a residue first, then press A.\n";
				}
				// showContextMenu = false;
			}

			// S - Show Structure
			if (sKeyPressed && !prevSKeyPressed) {
				if (selectedAtomIndex && moleculeData) {
					const Atom& atom = moleculeData->atoms[*selectedAtomIndex];

					std::cout << "\n========== STRUCTURAL ANALYSIS ==========" << std::endl;
					std::cout << "Residue:  " << atom.residueName << " " << atom.residueNum
						<< " (Chain " << atom.chain << ")" << std::endl;
					std::cout << "Atom:  " << atom.name << " (" << atom.element << ")" << std::endl;
					std::cout << "==========================================" << std::endl;

					// 1. Count nearby residues (simple contact analysis)
					int nearbyResidues = countNearbyResidues(atom, moleculeData, 8.0f);
					std::cout << "Nearby Residues (8Å): " << nearbyResidues << std::endl;

					// 2. Chemical properties from residue name
					std::string properties = getResidueProperties(atom.residueName);
					std::cout << "Chemical Properties:  " << properties << std::endl;

					// 3. Estimate burial based on contact count
					std::string burial;
					if (nearbyResidues > 15) {
						burial = "Buried (core)";
					}
					else if (nearbyResidues > 8) {
						burial = "Partially buried";
					}
					else {
						burial = "Surface exposed";
					}
					std::cout << "Estimated Location:   " << burial << std::endl;

					std::cout << "\nMutation Impact:" << std::endl;
					if (nearbyResidues > 15) {
						std::cout << "  • High constraint - mutations likely disruptive" << std::endl;
					}
					else {
						std::cout << "  • Lower constraint - more mutation tolerant" << std::endl;
					}

					std::cout << "=========================================\n" << std::endl;
				}
			}

			// D - Deselect
			if (dKeyPressed && !prevDKeyPressed) {
				if (selectedAtomIndex && moleculeData) {
					const Atom& atom = moleculeData->atoms[*selectedAtomIndex];
					Selection clearSel;
					clearSel.residue = atom.residueNum;
					clearSel.chain = atom.chain;
					Model::colorAtomsDefault(&clearSel);
					Model::setAtomRadius(SphereTemplate::DEFAULT_RADIUS, &clearSel);

					selectedAtomIndex = std::nullopt;
					std::cout << "\n[DESELECTED]" << std::endl;
				}
				showContextMenu = false;
			}

			prevAKeyPressed = aKeyPressed;
			prevSKeyPressed = sKeyPressed;
			prevDKeyPressed = dKeyPressed;
		}

		prevMKeyPressed = mKeyPressed;

		// Read commands sent from console
		for (size_t i = 0; i < commands.size(); ++i) {
			std::vector<std::string> commandWords =
				Parser::split(Parser::lowercase(commands[0]), ' ');

			std::cout << "DEBUG: Command has " << commandWords.size() << " words" << std::endl;
			if (!commandWords.empty()) {
				std::cout << "DEBUG: First word = '" << commandWords[0] << "'" << std::endl;
			}

			if (commandWords.size() >= 2 && commandWords[0] == "load") {
				std::string pdb_file = commandWords[1];

				std::cout << "Loading local PDB: " << pdb_file << "\n";

				try {
					// Check if file exists
					std::ifstream check(pdb_file);
					if (!check.good()) {
						std::cerr << "\n[!] File not found: \"" << pdb_file << "\"\n"
							<< " Check that the path is correct and the file exists.\n\n";
						commands.erase(commands.begin());
						continue;
					}
					check.close();

				}
				catch (const std::exception& e) {
					std::cerr << "\n[!] Error loading PDB file: " << e.what() << "\n"
						<< " Make sure the file is a valid PDB format.\n\n";
					delete moleculeData;
					moleculeData = nullptr;
				}
			}

			if (commandWords.size() == 3 && commandWords[0] == "fetch") {
				if (commandWords[1] == "pdb") {
					showContextMenu = false;  // ADD THIS

					// CLEAR HOVER AND SELECTION STATES
					if (hoveredAtomIndex && moleculeData) {
						const Atom& atom = moleculeData->atoms[*hoveredAtomIndex];
						Selection clearSel;
						clearSel.residue = atom.residueNum;
						clearSel.chain = atom.chain;
						Model::setAtomRadius(SphereTemplate::DEFAULT_RADIUS, &clearSel);
						Model::colorAtomsDefault(&clearSel);

					}
					if (selectedAtomIndex && moleculeData) {
						const Atom& atom = moleculeData->atoms[*selectedAtomIndex];
						Selection clearSel;
						clearSel.residue = atom.residueNum;
						clearSel.chain = atom.chain;
						Model::setAtomRadius(SphereTemplate::DEFAULT_RADIUS, &clearSel);
						Model::colorAtomsDefault(&clearSel);
					}
					hoveredAtomIndex = std::nullopt;
					selectedAtomIndex = std::nullopt;

					std::string pdbId = commandWords[2];
					currentPdbId = pdbId;

					// Validate PDB ID format (must be exactly 4 alphanumeric characters)
					if (commandWords[2].size() == 4) {
						bool validId = true;
						for (char c : commandWords[2]) {
							if (!std::isalnum(c)) { validId = false; break; }
						}
						if (!validId) {
							std::cerr << "\n[!] Invalid PDB ID \"" << commandWords[2] << "\".\n"
								<< " PDB IDs must be exactly 4 alphanumeric characters.\n"
								<< " Example: fetch pdb 1TUP\n\n";
							commands.erase(commands.begin());
							continue;
						}
					}

					delete moleculeData;
					moleculeData = nullptr;
					currentStructureData = PDBStructureData();

					std::string url;
					if (commandWords[2].size() == 4) {
						url = "http://files.rcsb.org/view/" + commandWords[2] + ".pdb";
						std::cout << "\nFetching " << commandWords[2] << " from RCSB PDB...\n";
					}
					else {
						url = commandWords[2];
						std::cout << "\nLoading from URL: " << url << "\n";
					}

					moleculeData = new PDBFile(url);
					if (moleculeData) {
						std::cout << "PDBFile created successfully\n";
					}

					// Fetch both PDB structure and mutation data
					bool success = FetchPDB::fetchComplete(
						pdbId,
						moleculeData,
						currentStructureData
					);

					// Build chain -> residue list from loaded structure
					std::map<char, std::vector<int>> chainResidues;
					std::map<char, std::set<int>> seen;
					for (const auto& atom : moleculeData->atoms) {
						if (seen[atom.chain].insert(atom.residueNum).second) {
							chainResidues[atom.chain].push_back(atom.residueNum);
						}
					}
					// Sort each chain's residue list
					for (auto& [ch, vec] : chainResidues)
						std::sort(vec.begin(), vec.end());

					// Launch background fetch
					mutationPanel.showWindow = false;
					std::string capPdbId = currentPdbId;
					std::thread([capPdbId, chainResidues]() {
						fetchAllMutationsForPDB(capPdbId, chainResidues);
						}).detach();

					// ADD THIS — pre-load all variants for the variant selector panel
					variantData = json::array();
					variantPositions.clear();

					bool* pShowVariantPanel = &showVariantPanel;
					bool* pVariantLoading = nullptr; // optional

					std::thread([capPdbId, &variantData, &variantPositions, pShowVariantPanel]() {
						HTTPConnection http;
						std::string url = BACKEND_URL + "/api/alphamissense/variants/" + capPdbId;
						std::cout << "[VARIANTS] Fetching all variants for " << capPdbId << "...\n";

						if (!http.get(url)) {
							std::cerr << "[VARIANTS] Failed to reach variants endpoint.\n";
							return;
						}

						try {
							auto j = json::parse(http.getResponse());
							if (j.contains("detail")) {
								std::cerr << "[VARIANTS] " << j["detail"].get<std::string>() << "\n";
								return;
							}

							variantData = j["variants"];
							;

							variantPositions.clear();
							for (auto& v : variantData) {
								variantPositions.insert(v["position"].get<int>());
								std::string chainId = v.value("chain_id", "A");
								bool foldxReady = false;
								{
									std::lock_guard<std::mutex> lock(mutationPanel.dataMutex);
									foldxReady = !chainId.empty() &&
										mutationPanel.chains.count(chainId[0]) > 0;
								}
								v["foldx_ready"] = foldxReady;
							}

							std::cout << "[VARIANTS] Loaded " << variantData.size()
								<< " variants across " << variantPositions.size()
								<< " positions.\n";
							*pShowVariantPanel = true;
						}
						catch (const json::exception& e) {
							std::cerr << "[VARIANTS] Parse error: " << e.what() << "\n";
						}
						}).detach();

					std::thread([capPdbId]() {
						HTTPConnection http;
						std::string url = BACKEND_URL + "/api/pdb/" + capPdbId + "/all-residues";
						std::cout << "[CHAINS] Fetching all chains for " << capPdbId << "...\n";

						if (!http.get(url)) {
							std::cerr << "[CHAINS] Failed to reach endpoint.\n";
							return;
						}

						try {
							auto j = json::parse(http.getResponse());
							if (j.contains("detail")) {
								std::cerr << "[CHAINS] " << j["detail"].get<std::string>() << "\n";
								return;
							}
							allChains = j.value("chains", json::array());   // ← just use it directly
							std::cout << "[CHAINS] Loaded " << allChains.size() << " chains.\n";
						}
						catch (const json::exception& e) {
							std::cerr << "[CHAINS] Parse error: " << e.what() << "\n";
						}
						}).detach();

					// Print first atom position for reference:
					if (success && moleculeData && !moleculeData->atoms.empty()) {
						auto& atom = moleculeData->atoms[0];
						std::cout << "First atom at: " << atom.coords.getX() << ", "
							<< atom.coords.getY() << ", " << atom.coords.getZ() << std::endl;

						// Calculate center
						Vec3 center(0, 0, 0);
						for (const auto& atom : moleculeData->atoms) {
							center = center + atom.coords;
						}
						center = center / (float)moleculeData->atoms.size();

						std::cout << "Original protein center: " << center.getX() << ", "
							<< center.getY() << ", " << center.getZ() << std::endl;

						// Shift all atoms to origin
						for (auto& atom : moleculeData->atoms) {
							const_cast<Vec3&>(atom.coords) = atom.coords - center;
						}

						// Camera should look at origin (where atoms now are)
						camera.setTarget(Vec3(0, 0, 0));  // ← Look at origin, not old center!

						Model::loadMoleculeData(moleculeData);
						std::cout << "Model::loadMoleculeData called\n";

						representationRenderer.loadMoleculeData(moleculeData);
						representationRenderer.setRepresentationEnabled(RepresentationType::RIBBON, true);

						selection.reset();
						std::cout << "Loaded protein: " << pdbId << std::endl;

					}
					else if (!success) {
						std::cerr << "\n[!] Failed to load structure \"" << pdbId << "\".\n"
							<< " Possible causes:\n"
							<< " - The PDB ID does not exist (check https://rcsb.org)\n"
							<< " - No internet connection\n"
							<< " - The RCSB server is temporarily unavailable\n\n";
					}

					// NEW: Suggest mutations after loading
					if (mutationSuggester && alphaMissenseDB.isReady()) {
						// Assume chain A for now (you can make this configurable)
						auto suggestions = mutationSuggester->suggestMutations(
							pdbId,
							moleculeData,
							20  // Show top 20 mutations
						);
					}

					// Notify chat window
					if (chatWindowPtr) {
						chatWindowPtr->addSystemMessage("Loaded protein structure: " + pdbId);
					}

				}
			}
			if (commandWords.size() == 1 && commandWords[0] == "listchains") {
				if (!moleculeData) {
					std::cerr << "\n[!] No protein structure is loaded.\n"
						<< " Load one first with: fetch pdb <pdb_id>\n"
						<< " Example: fetch pdb 1TUP\n\n";
				}
				else {
					std::set<char> uniqueChains;
					for (const auto& atom : moleculeData->atoms) {
						uniqueChains.insert(atom.chain);
					}

					std::cout << "\n========== CHAINS IN STRUCTURE ==========" << std::endl;
					std::cout << "Available chains: ";
					for (char chain : uniqueChains) {
						std::cout << chain << " ";
					}
					std::cout << "\n";
					std::cout << "Total chains: " << uniqueChains.size() << std::endl;
					std::cout << "========================================\n" << std::endl;
				}
			}

			else if (commandWords.size() >= 2 && commandWords[0] == "select") {
				if (!moleculeData) {
					std::cerr << "\n[!] No protein structure is loaded.\n"
						<< " Load one first with: fetch pdb <pdb_id>\n\n";
				}
				else {
					// Parse selection using existing selection syntax
					std::vector<std::string> selectionQuery(
						commandWords.begin() + 1,  // Skip "select"
						commandWords.end()
					);

					selection.parseQuery(selectionQuery);

					// Highlight the residue
					Color highlightColor = Color::fromName("white");
					Model::setAtomRadius(0.4f, &selection);
					Model::setAtomColor(&highlightColor, &selection);

					std::cout << "Selection applied and highlighted" << std::endl;
					selection.print();  // Show what was selected
				}
			}
			else if (commandWords.size() >= 4 && commandWords[0] == "mutate") {
				if (!moleculeData) {
					std::cerr << "\n[!] No protein structure is loaded.\n"
						<< " Load one first with: fetch pdb <pdb_id>\n"
						<< " Example: fetch pdb 1TUP\n\n";
				}
				else if (!mutationAnalyzer) {
					std::cerr << "\n[!] The mutation analyzer failed to initialize.\n"
						<< " Check that the backend is running and try restarting DataLens.\n\n";
				}
				else if (currentPdbId.empty()) {
					std::cerr << "\n[!] No PDB ID is associated with the current structure.\n"
						<< " Please reload the structure with: fetch pdb <pdb_id>\n\n";
				}
				else {
					try {
						char chain = commandWords[1][0];
						int position = std::stoi(commandWords[2]);
						char mutantAA = std::toupper(commandWords[3][0]);

						auto result = mutationAnalyzer->analyzeMutation(
							currentPdbId,
							moleculeData,
							chain,
							position,
							mutantAA
						);

						if (result) {
							result->print();
						}
					}
					catch (const std::exception& e) {
						std::cerr << "\n[!] Invalid mutate arguments: " << e.what() << "\n"
							<< " Usage:   mutate <chain> <position> <amino_acid>\n"
							<< " Example: mutate A 175 H\n"
							<< " Chain must be a single letter (e.g. A, B).\n"
							<< " Position must be an integer.\n"
							<< " Amino acid must be a single letter (e.g. H, G, P).\n\n";
					}
				}
			}

			else if (commandWords.size() == 2 && commandWords[0] == "rotate") {
				const float TOTAL_ROTATION_DEGREES = 20.0f;
				const int NUMBER_OF_STEPS = 10;

				if (commandWords[1] == "x") {
					for (int i = 0; i < NUMBER_OF_STEPS; i++) {
						Model::rotate(Vec3(TOTAL_ROTATION_DEGREES / NUMBER_OF_STEPS, 0.0f, 0.0f));
					}
				}
				else if (commandWords[1] == "y") {
					for (int i = 0; i < NUMBER_OF_STEPS; i++) {
						Model::rotate(Vec3(0.0f, TOTAL_ROTATION_DEGREES / NUMBER_OF_STEPS, 0.0f));
					}
				}
				else if (commandWords[1] == "z") {
					for (int i = 0; i < NUMBER_OF_STEPS; i++) {
						Model::rotate(Vec3(0.0f, 0.0f, TOTAL_ROTATION_DEGREES / NUMBER_OF_STEPS));
					}
				}
				else {
					std::cerr << "\n[!] Unknown rotation axis \"" << commandWords[1] << "\".\n"
						<< " Valid axes are: x, y, z\n"
						<< " Example: rotate x\n\n";
				}
			}

			else if (commandWords.size() == 1 && commandWords[0] == "chains") {
				if (currentPdbId.empty()) {
					std::cerr << "\n[!] No protein structure is loaded.\n"
						<< " Load one first with: fetch pdb <pdb_id>\n"
						<< " Example: fetch pdb 1TUP\n\n";
				}
				else if (pdbMapper.isReady()) {
					pdbMapper.displayPDBInfo(currentPdbId);
				}
				else {
					std::cerr << "\n[!] The PDB-UniProt mapper is not ready.\n"
						<< " Try reloading the structure or restarting DataLens.\n\n";
				}
			}

			// NEW COMMAND: Clear selection and reset highlighting
			else if (commandWords.size() == 1 && commandWords[0] == "clearselection") {
				if (selectedAtomIndex && moleculeData) {
					const Atom& atom = moleculeData->atoms[*selectedAtomIndex];
					Selection clearSel;
					clearSel.residue = atom.residueNum;
					clearSel.chain = atom.chain;

					// RESET SIZE AND COLOR
					Model::setAtomRadius(SphereTemplate::DEFAULT_RADIUS, &clearSel);
					Model::colorAtomsDefault(&clearSel);
				}
				selectedAtomIndex = std::nullopt;
				selection.reset();

				std::cout << "Selection cleared" << std::endl;
			}

			// In the command parsing loop, add:
			else if (commandWords.size() == 2 && commandWords[0] == "repr") {
				if (commandWords[1] == "ribbon") {
					representationRenderer.toggleRepresentation(RepresentationType::RIBBON);
				}
				else if (commandWords[1] == "surface") {
					representationRenderer.toggleRepresentation(RepresentationType::SURFACE);
				}
				else if (commandWords[1] == "ballandstick") {
					representationRenderer.toggleRepresentation(RepresentationType::BALL_AND_STICK);
				}
				else if (commandWords[1] == "vdw") {
					representationRenderer.setSurfaceType(SurfaceTemplate::SurfaceType::VAN_DER_WAALS);
				}
				else if (commandWords[1] == "sas") {
					representationRenderer.setSurfaceType(SurfaceTemplate::SurfaceType::SOLVENT_ACCESSIBLE);
				}
				else if (commandWords[1] == "ses") {
					representationRenderer.setSurfaceType(SurfaceTemplate::SurfaceType::SOLVENT_EXCLUDED);
				}
			}

			// NEW COMMAND: Highlight currently selected residue
			else if (commandWords.size() >= 2 && commandWords[0] == "highlight") {
				if (commandWords.size() == 2) {
					try {
						Color highlightColor = Color::fromName(commandWords[1]);
						Model::setAtomRadius(0.35f, &selection);
						Model::setAtomColor(&highlightColor, &selection);
						std::cout << "Highlighted selection with color: " << commandWords[1] << std::endl;
					}
					catch (...) {
						std::cerr << "\n[!] Unknown color name \"" << commandWords[1] << "\".\n"
							<< " Try a basic color: red, green, blue, white, yellow, cyan, magenta, orange.\n\n";
					}
				}
				else {
					// Default yellow highlight
					Color highlightColor = Color::fromName("yellow");
					Model::setAtomRadius(0.35f, &selection);
					Model::setAtomColor(&highlightColor, &selection);
					std::cout << "Highlighted selection" << std::endl;
				}
			}

			// line 2178 — was: commandWords.size() == 1 (correct count, wrong path source)
			else if (commandWords.size() == 1 && commandWords[0] == "loadmutant") {
				if (!moleculeData) {
					std::cerr << "[COMPARISON] No WT structure loaded.\n";
				}
				else if (mutantPdbContent.empty()) {
					std::cerr << "[COMPARISON] No mutant PDB available. Run an analysis first.\n";
				}
				else {
					delete mutantMoleculeData;
					mutantMoleculeData = new PDBFile(mutantPdbContent, true);

					if (mutantMoleculeData && !mutantMoleculeData->atoms.empty()) {
						buildComparison(
							structureComparison,
							moleculeData,
							mutantMoleculeData,
							currentPdbId + " (WT)",
							mutFocus.variantId,
							mutFocus.position,
							mutFocus.chain
						);
						std::cout << "[COMPARISON] Built. Global RMSD: "
							<< structureComparison.globalRmsd << " Å\n";
					}
					else {
						std::cerr << "[COMPARISON] Failed to parse mutant PDB from: "
							<< mutantPdbPath << "\n";
					}
				}
			}

			else if (commandWords.size() == 1 && commandWords[0] == "leave") {
				// Clear hover and selection WITH SIZE RESET
				if (hoveredAtomIndex && moleculeData) {
					const Atom& atom = moleculeData->atoms[*hoveredAtomIndex];
					Selection clearSel;
					clearSel.residue = atom.residueNum;
					clearSel.chain = atom.chain;
					Model::setAtomRadius(SphereTemplate::DEFAULT_RADIUS, &clearSel);
					Model::colorAtomsDefault(&clearSel);
				}
				if (selectedAtomIndex && moleculeData) {
					const Atom& atom = moleculeData->atoms[*selectedAtomIndex];
					Selection clearSel;
					clearSel.residue = atom.residueNum;
					clearSel.chain = atom.chain;
					Model::setAtomRadius(SphereTemplate::DEFAULT_RADIUS, &clearSel);
					Model::colorAtomsDefault(&clearSel);
				}
				hoveredAtomIndex = std::nullopt;
				selectedAtomIndex = std::nullopt;

				delete moleculeData;
				moleculeData = nullptr;
				currentMolData = nullptr;
				currentPdbId = "";
				currentStructureData = PDBStructureData();

				// Mutant structure
				delete mutantMoleculeData;
				mutantMoleculeData = nullptr;
				showingMutant = false;

				Model::reset();
				camera.reset();
				selection.reset();
				std::cout << "Structure unloaded" << std::endl;
			}

			else if (commandWords.size() == 1 && commandWords[0] == "setdefault") {
				selection.reset();
				Model::setAtomRadius(SphereTemplate::DEFAULT_RADIUS, &selection);
				Model::setConnectorRadius(
					ConnectorTemplate::DEFAULT_RADIUS,
					&selection,
					ConnectorType::BACKBONE
				);
				Model::setConnectorRadius(
					0.0f,
					&selection,
					ConnectorType::DISULFIDE_BOND
				);
				Model::colorAtomsDefault(&selection);
				Model::colorConnectorsDefault(&selection, ConnectorType::BACKBONE);
				Model::colorConnectorsDefault(&selection, ConnectorType::DISULFIDE_BOND);

				std::cout << "Reset to default view" << std::endl;
			}
			else {
				std::cerr << "\n[!] Unknown command: \"" << commandWords[0] << "\"\n\n"
					<< "  Available commands:\n"
					<< " fetch pdb <id> Load a protein (e.g. fetch pdb 1TUP)\n"
					<< " rotate <x|y|z> Rotate model along an axis\n"
					<< " select r=<range> c=<chain>  Select residues (e.g. select r=1:50 c=A)\n"
					<< " highlight [color]   Highlight selection (e.g. highlight yellow)\n"
					<< " color atom <r> <g> <b>   Set atom color\n"
					<< " mutate <chain> <pos> <aa>   Run mutation analysis (e.g. mutate A 175 H)\n"
					<< " chains    List chains in current structure\n"
					<< " listchains  List all chain identifiers\n"
					<< " clearselection Clear current selection\n"
					<< " setdefault  Reset all visual settings\n"
					<< " leave  Unload current structure\n"
					<< " exit   Quit DataLens\n\n";
			}

			commands.erase(commands.begin());
		}

		// Check window size is valid before starting ImGui frame
		float windowWidth = (float)window.getWidth();
		float windowHeight = (float)window.getHeight();

		if (windowWidth <= 0.0f || windowHeight <= 0.0f) {
			std::cerr << "Warning: Invalid window size, skipping frame" << std::endl;
			continue;  // Skip this frame entirely
		}

		// Poll Claude response — fires once when worker thread finishes
		if (claudePending && claudeHandler.isReady()) {
			std::string reply = claudeHandler.getResponse();
			if (!reply.empty()) {
				addChatMessage("Assistant", reply);
			}
			claudePending = false;
		}

		// Start ImGui frame
		ImGui_ImplOpenGL3_NewFrame();
		ImGui_ImplGlfw_NewFrame();

		// Get IO from current context and set display size
		ImGuiIO& io = ImGui::GetIO();
		io.DisplaySize = ImVec2(windowWidth, windowHeight);

		ImGui::NewFrame();
		//imguiWantsMouse = ImGui::GetIO().WantCaptureMouse;  // ← refresh after NewFrame

		// Check ImGui vs Menu
		// == Mouse handling starts ==
		// Check if ImGui wants to capture mouse input (hovering/clicking UI)
		bool imguiWantsMouse = io.WantCaptureMouse;

		// Detect mouse button press (start of potential click)
		if (leftMouseButtonPressed && !prevLeftMouseButtonPressed && !imguiWantsMouse) {
			clickStartX = mouseX;
			clickStartY = mouseY;
		}

		if (prevMousePosSet && !imguiWantsMouse) {
			if (leftMouseButtonPressed) {
				if (Input::keyPressed(&window, Key::LEFT_CTRL) ||
					Input::keyPressed(&window, Key::RIGHT_CTRL)) {
					camera.move(
						Vec3(
							float(prevMouseX - mouseX) / 100,
							float(prevMouseY - mouseY) / -100,
							0.0f
						)
					);
				}
				else {
					Model::rotate(
						Vec3(
							float(prevMouseY - mouseY) / 100,
							float(prevMouseX - mouseX) / 100,
							0.0f
						)
					);
				}
			}
			else if (Input::mouseButtonPressed(&window, MouseButton::RIGHT)) {
				Model::rotate(Vec3(0.0f, 0.0f, float(prevMouseX - mouseX) / 100));
			}
		}

		/*
		 * HOVER WITH SIZE RESET
		 *
		 * Behavior:
		 * - HOVER: Gold color + LARGER size
		 * - CLICK: Orange color + LARGER size (persists)
		 * - MOUSE AWAY: Return to NORMAL size + color (unless selected)
		 */

		if (moleculeData != nullptr && !showContextMenu && !imguiWantsMouse) {
			bool isMouseDragging = (leftMouseButtonPressed &&
				(Input::keyPressed(&window, Key::LEFT_CTRL) ||
					Input::keyPressed(&window, Key::RIGHT_CTRL) ||
					prevMousePosSet));

			// Only hover when not dragging
			if (!isMouseDragging && !leftMouseButtonPressed) {

				// Unproject (mouseX, mouseY) from screen space through the inverse view-projection
				// matrix to produce a normalized world-space ray originating at the camera position.
				Ray ray = RayCaster::screenToRay(mouseX, mouseY, &window, &camera);

				// Test the ray against each atom's bounding sphere in world space, returning the
				// instance buffer index of the nearest hit, or std::nullopt if nothing was struck.
				auto hit = RayCaster::castRay(ray, Model::getAtomSpheres(), Model::getModelMatrix());

				if (hit && hit->atomIndex < moleculeData->atoms.size()) {
					size_t atomIndex = hit->atomIndex;
					const Atom& currentAtom = moleculeData->atoms[atomIndex];

					// Check if this atom is part of the selected RESIDUE (not just same atom index)
					bool isInSelectedResidue = false;
					if (selectedAtomIndex) {
						const Atom& selectedAtom = moleculeData->atoms[*selectedAtomIndex];
						isInSelectedResidue = (currentAtom.residueNum == selectedAtom.residueNum &&
							currentAtom.chain == selectedAtom.chain);
					}

					// If hovering over the SELECTED RESIDUE
					if (isInSelectedResidue) {
						hoveredAtomIndex = atomIndex;
					}

					// If hovering over an UNSELECTED residue, apply cyan hover
					else {
						// Clear previous hover only if it's different and NOT selected
						if (hoveredAtomIndex && *hoveredAtomIndex != atomIndex) {
							const Atom& prevAtom = moleculeData->atoms[*hoveredAtomIndex];

							// Check if previous hover was in selected residue
							bool prevWasInSelectedResidue = false;
							if (selectedAtomIndex) {
								const Atom& selectedAtom = moleculeData->atoms[*selectedAtomIndex];
								prevWasInSelectedResidue = (prevAtom.residueNum == selectedAtom.residueNum &&
									prevAtom.chain == selectedAtom.chain);
							}

							if (!prevWasInSelectedResidue) {
								Selection clearSel;
								clearSel.residue = prevAtom.residueNum;
								clearSel.chain = prevAtom.chain;
								Model::setAtomRadius(SphereTemplate::DEFAULT_RADIUS, &clearSel);
								Model::colorAtomsDefault(&clearSel);
							}
						}

						// Apply cyan hover to the new unselected residue
						if (!hoveredAtomIndex || *hoveredAtomIndex != atomIndex) {
							hoveredAtomIndex = atomIndex;

							Selection hoverSel;
							hoverSel.residue = currentAtom.residueNum;
							hoverSel.chain = currentAtom.chain;
							Color cyan = Color::fromByte(0, 255, 255);
							Model::setAtomRadius(0.30f, &hoverSel);
							Model::setAtomColor(&cyan, &hoverSel);
						}
					}
				}
				else {
					// Mouse not over any atom - clear hover (but NEVER the selection)
					if (hoveredAtomIndex) {
						const Atom& prevAtom = moleculeData->atoms[*hoveredAtomIndex];

						// Check if hover was in selected residue
						bool wasInSelectedResidue = false;
						if (selectedAtomIndex) {
							const Atom& selectedAtom = moleculeData->atoms[*selectedAtomIndex];
							wasInSelectedResidue = (prevAtom.residueNum == selectedAtom.residueNum &&
								prevAtom.chain == selectedAtom.chain);
						}

						if (!wasInSelectedResidue) {
							Selection clearSel;
							clearSel.residue = prevAtom.residueNum;
							clearSel.chain = prevAtom.chain;
							Model::setAtomRadius(SphereTemplate::DEFAULT_RADIUS, &clearSel);
							Model::colorAtomsDefault(&clearSel);
						}
					}
					hoveredAtomIndex = std::nullopt;
				}
			}
		}

		// Detect mouse button release (end of click)
		if (prevLeftMouseButtonPressed && !leftMouseButtonPressed && !imguiWantsMouse) {
			// Calculate distance moved during click
			double deltaX = mouseX - clickStartX;
			double deltaY = mouseY - clickStartY;
			double dragDistance = std::sqrt(deltaX * deltaX + deltaY * deltaY);

			if (dragDistance < DRAG_THRESHOLD &&
				!Input::keyPressed(&window, Key::LEFT_CTRL) &&
				!Input::keyPressed(&window, Key::RIGHT_CTRL) &&
				moleculeData != nullptr) {

				Ray ray = RayCaster::screenToRay(mouseX, mouseY, &window, &camera);
				auto hit = RayCaster::castRay(ray, Model::getAtomSpheres(), Model::getModelMatrix());

				if (hit) {
					const Atom& clickedAtom = moleculeData->atoms[hit->atomIndex];

					// Check if clicking on a different atom than currently selected
					bool isDifferentAtom = (!selectedAtomIndex || *selectedAtomIndex != hit->atomIndex);

					// Only clear previous selection if clicking a DIFFERENT atom
					if (isDifferentAtom && selectedAtomIndex) {
						const Atom& prevAtom = moleculeData->atoms[*selectedAtomIndex];
						Selection clearSel;
						clearSel.residue = prevAtom.residueNum;
						clearSel.chain = prevAtom.chain;
						Model::colorAtomsDefault(&clearSel);
						Model::setAtomRadius(SphereTemplate::DEFAULT_RADIUS, &clearSel);
					}

					// Update selected index
					selectedAtomIndex = hit->atomIndex;

					std::cout << "\n========== RESIDUE SELECTED ==========" << std::endl;
					std::cout << "Residue:  " << clickedAtom.residueName
						<< " " << clickedAtom.residueNum << std::endl;
					std::cout << "Chain: " << clickedAtom.chain << std::endl;
					std::cout << "Atom:  " << clickedAtom.name
						<< " (" << clickedAtom.element << ")" << std::endl;
					std::cout << "Position: ("
						<< clickedAtom.coords.getX() << ", "
						<< clickedAtom.coords.getY() << ", "
						<< clickedAtom.coords.getZ() << ")" << std::endl;
					std::cout << "=========================================\n" << std::endl;

					// Update selection to this residue
					std::vector<std::string> selectionQuery = {
						"r=" + std::to_string(clickedAtom.residueNum),
						"c=" + std::string(1, clickedAtom.chain)
					};
					selection.parseQuery(selectionQuery);

					// Apply ORANGE highlight for selection
					Color selectColor = Color::fromByte(255, 100, 0);
					Model::setAtomRadius(0.2f, &selection);
					Model::setAtomColor(&selectColor, &selection);
				}
				else {
					std::cout << "No atom selected (clicked on empty space)" << std::endl;
				}

				std::cout << "Ray origin: " << ray.origin.getX() << ", "
					<< ray.origin.getY() << ", " << ray.origin.getZ() << std::endl;
				std::cout << "Ray direction: " << ray.direction.getX() << ", "
					<< ray.direction.getY() << ", " << ray.direction.getZ() << std::endl;
			}
		}
		// == Mouse handling ends ==

		prevMouseX = mouseX;
		prevMouseY = mouseY;
		prevMousePosSet = true;
		prevLeftMouseButtonPressed = leftMouseButtonPressed;

		// ImGui Tooltip, before Legend window
		// Hover tooltip
		if (hoveredAtomIndex && moleculeData &&
			*hoveredAtomIndex < moleculeData->atoms.size()) {

			const Atom& hoveredAtom = moleculeData->atoms[*hoveredAtomIndex];

			ImGui::BeginTooltip();

			ImGui::TextColored(ImVec4(0.0f, 1.0f, 1.0f, 1.0f),
				"%s %d", hoveredAtom.residueName.c_str(),
				hoveredAtom.residueNum);
			ImGui::Text("Chain: %c", hoveredAtom.chain);
			ImGui::Text("Atom: %s (%s)", hoveredAtom.name.c_str(),
				hoveredAtom.element.c_str());

			if (selectedAtomIndex && *selectedAtomIndex == *hoveredAtomIndex) {
				ImGui::Separator();
				ImGui::TextColored(ImVec4(1.0f, 0.4f, 0.0f, 1.0f), "SELECTED!");
				ImGui::TextColored(ImVec4(0.7f, 0.7f, 0.7f, 1.0f), "Press M for options");
			}
			else {
				ImGui::Separator();
				ImGui::TextColored(ImVec4(0.7f, 0.7f, 0.7f, 1.0f), "Click to select");
			}

			ImGui::EndTooltip();
		}

		/*
		* CONTEXT MENU for selected residue
		*/
		if (showContextMenu && selectedAtomIndex && moleculeData) {
			const Atom& selectedAtom = moleculeData->atoms[*selectedAtomIndex];

			ImGui::SetNextWindowPos(contextMenuPos, ImGuiCond_Always);

			ImGui::Begin("Residue Menu", &showContextMenu);

			ImGui::Text("%s %d (Chain %c)",
				selectedAtom.residueName.c_str(),
				selectedAtom.residueNum,
				selectedAtom.chain);
			ImGui::Separator();

			if (ImGui::MenuItem("Mutation Analysis", "A")) {
				std::cout << "Mutate clicked!" << std::endl;
				showContextMenu = false;
			}

			if (ImGui::MenuItem("Show Structure", "S")) {
				std::cout << "Show info clicked!" << std::endl;
				showContextMenu = false;
			}

			if (ImGui::MenuItem("Deselect", "D")) {
				std::cout << "Deselect clicked!" << std::endl;
				showContextMenu = false;
			}

			ImGui::End();
		}

		// Close context menu if clicked elsewhere
		if (showContextMenu && leftMouseButtonPressed) {
			showContextMenu = false;
		}

		// Create Legend Window
		ImGui::SetNextWindowPos(ImVec2(10, 10), ImGuiCond_FirstUseEver);
		ImGui::SetNextWindowSize(ImVec2(280, 200), ImGuiCond_FirstUseEver);
		ImGui::Begin("Color Legend", nullptr, ImGuiWindowFlags_AlwaysAutoResize);

		ImGui::Text("Atom/Structure Colors");
		ImGui::Separator();

		// Create color legend table
		if (ImGui::BeginTable("LegendTable", 2, ImGuiTableFlags_Borders | ImGuiTableFlags_RowBg)) {
			ImGui::TableSetupColumn("Color", ImGuiTableColumnFlags_WidthFixed, 60.0f);
			ImGui::TableSetupColumn("Meaning", ImGuiTableColumnFlags_WidthStretch);
			ImGui::TableHeadersRow();

			// Red
			ImGui::TableNextRow();
			ImGui::TableSetColumnIndex(0);
			ImGui::PushStyleColor(ImGuiCol_Text, ImVec4(1.0f, 0.0f, 0.0f, 1.0f));
			ImGui::Text("Red");
			ImGui::PopStyleColor();
			ImGui::TableSetColumnIndex(1);
			ImGui::Text("Oxygen Atoms");

			// Blue
			ImGui::TableNextRow();
			ImGui::TableSetColumnIndex(0);
			ImGui::PushStyleColor(ImGuiCol_Text, ImVec4(0.0f, 0.5f, 1.0f, 1.0f));
			ImGui::Text("Blue");
			ImGui::PopStyleColor();
			ImGui::TableSetColumnIndex(1);
			ImGui::Text("Nitrogen Atoms");

			// Gray
			ImGui::TableNextRow();
			ImGui::TableSetColumnIndex(0);
			ImGui::PushStyleColor(ImGuiCol_Text, ImVec4(0.6f, 0.6f, 0.6f, 1.0f));
			ImGui::Text("Gray");
			ImGui::PopStyleColor();
			ImGui::TableSetColumnIndex(1);
			ImGui::Text("Carbon/Backbone");

			// Yellow
			ImGui::TableNextRow();
			ImGui::TableSetColumnIndex(0);
			ImGui::PushStyleColor(ImGuiCol_Text, ImVec4(1.0f, 1.0f, 0.0f, 1.0f));
			ImGui::Text("Yellow");
			ImGui::PopStyleColor();
			ImGui::TableSetColumnIndex(1);
			ImGui::Text("Sulfur/SS Bonds");

			ImGui::EndTable();
		}

		ImGui::Separator();
		ImGui::TextWrapped("Rotate: Left Mouse Drag");
		ImGui::TextWrapped("Pan: Ctrl + Left Mouse");
		ImGui::TextWrapped("Zoom: Mouse Wheel");

		// Store position and size for context menu positioning
		colorLegendPos = ImGui::GetWindowPos();
		colorLegendSize = ImGui::GetWindowSize();

		// Amino Acid Legend table starts
		ImGui::Spacing();
		ImGui::Separator();
		ImGui::Spacing();

		ImGui::TextColored(ImVec4(0.9f, 0.85f, 0.6f, 1.0f), "Amino Acid Codes");
		ImGui::Spacing();

		// 20 standard amino acids: {3-letter, 1-letter, full name}
		static const std::tuple<const char*, char, const char*> AA_TABLE[] = {
			{"ALA", 'A', "Alanine"},
			{"ARG", 'R', "Arginine"},
			{"ASN", 'N', "Asparagine"},
			{"ASP", 'D', "Aspartate"},
			{"CYS", 'C', "Cysteine"},
			{"GLN", 'Q', "Glutamine"},
			{"GLU", 'E', "Glutamate"},
			{"GLY", 'G', "Glycine"},
			{"HIS", 'H', "Histidine"},
			{"ILE", 'I', "Isoleucine"},
			{"LEU", 'L', "Leucine"},
			{"LYS", 'K', "Lysine"},
			{"MET", 'M', "Methionine"},
			{"PHE", 'F', "Phenylalanine"},
			{"PRO", 'P', "Proline"},
			{"SER", 'S', "Serine"},
			{"THR", 'T', "Threonine"},
			{"TRP", 'W', "Tryptophan"},
			{"TYR", 'Y', "Tyrosine"},
			{"VAL", 'V', "Valine"},
		};

		// Render in 2 columns
		if (ImGui::BeginTable("##aa_legend", 2, ImGuiTableFlags_SizingStretchSame)) {
			for (int i = 0; i < 20; i++) {
				ImGui::TableNextColumn();
				ImGui::TextColored(ImVec4(0.5f, 0.85f, 1.0f, 1.0f),
					"%s", std::get<0>(AA_TABLE[i])); // ALA
				ImGui::SameLine();
				ImGui::TextColored(ImVec4(1.0f, 1.0f, 1.0f, 0.5f), "=");
				ImGui::SameLine();
				ImGui::TextColored(ImVec4(1.0f, 0.75f, 0.4f, 1.0f),
					"%c", std::get<1>(AA_TABLE[i])); // A
				ImGui::SameLine();
				ImGui::TextColored(ImVec4(0.75f, 0.75f, 0.75f, 1.0f),
					"(%s)", std::get<2>(AA_TABLE[i]));  // (Alanine)
			}
			ImGui::EndTable();
		}
		ImGui::Spacing();
		// Amino legend table ends

		ImGui::End(); // Legend window ends


		ImGui::SetNextWindowPos(
			ImVec2((float)window.getWidth() - 420.0f, (float)window.getHeight() - 380.0f),
			ImGuiCond_FirstUseEver
		);
		ImGui::SetNextWindowSizeConstraints(ImVec2(300.0f, 200.0f), ImVec2(600.0f, 500.0f));
		ImGui::SetNextWindowSize(ImVec2(400.0f, 360.0f), ImGuiCond_FirstUseEver);
		ImGui::Begin("AI Assistant");

		ImGui::TextColored(ImVec4(0.4f, 0.8f, 1.0f, 1.0f), "AI Chat Assistant");

		// Show current mutation context as a subtitle if available
		if (mutFocus.show && !mutFocus.variantId.empty()) {
			ImGui::SameLine();
			ImGui::TextDisabled(" | %s  DDG %+.2f  AM %.2f",
				mutFocus.variantId.c_str(),
				mutFocus.ddg,
				mutFocus.amScore);
		}

		ImGui::Separator();
		ImGui::Spacing();

		// ── Chat history scroll region ────────────────────────────────
		ImGui::BeginChild("ChatHistory", ImVec2(0, -100), true);
		{
			std::lock_guard<std::mutex> lock(chatMutex);
			for (const auto& msg : chatHistory) {
				if (msg.find("You:") == 0) {
					ImGui::TextColored(ImVec4(0.4f, 0.8f, 1.0f, 1.0f), "%s", msg.c_str());
				}
				else if (msg.find("Assistant:") == 0) {
					// Wrap long responses — TextWrapped respects window width
					ImGui::PushStyleColor(ImGuiCol_Text, ImVec4(0.4f, 1.0f, 0.4f, 1.0f));
					ImGui::TextWrapped("%s", msg.c_str());
					ImGui::PopStyleColor();
				}
				else {
					// System messages (e.g. protein loaded notification)
					ImGui::PushStyleColor(ImGuiCol_Text, ImVec4(0.6f, 0.6f, 0.6f, 1.0f));
					ImGui::TextWrapped("%s", msg.c_str());
					ImGui::PopStyleColor();
				}
				ImGui::Spacing();
			}
		}

		// Animated thinking indicator
		if (claudePending) {
			static float dotTimer = 0.0f;
			dotTimer += ImGui::GetIO().DeltaTime;
			int dots = (int)(dotTimer * 2.0f) % 4;   // cycles 0-3
			std::string thinking = "Thinking" + std::string(dots, '.');
			ImGui::TextColored(ImVec4(1.0f, 1.0f, 0.4f, 1.0f), "%s", thinking.c_str());
		}

		ImGui::SetScrollHereY(1.0f);  // auto-scroll to latest message
		ImGui::EndChild();

		// ── Input area ────────────────────────────────────────────────
		ImGui::Separator();
		ImGui::Spacing();

		// Submit on Enter (but allow Shift+Enter for newline)
		bool submitOnEnter = ImGui::IsKeyPressed(ImGuiKey_Enter)
			&& !ImGui::GetIO().KeyShift
			&& ImGui::IsItemFocused();

		ImGui::SetNextItemWidth(-1);
		ImGui::InputTextMultiline(
			"##chat_input",
			chatInputBuffer, sizeof(chatInputBuffer),
			ImVec2(-1.0f, 40.0f),
			ImGuiInputTextFlags_EnterReturnsTrue
		);

		bool sendPressed = ImGui::Button("Send", ImVec2(80, 0));

		// "Ask about selected mutation" shortcut button
		if (mutFocus.show && !mutFocus.variantId.empty()) {
			ImGui::SameLine();
			if (ImGui::Button("Ask about mutation", ImVec2(0, 0))) {
				std::snprintf(
					chatInputBuffer, sizeof(chatInputBuffer),
					"Explain the clinical significance of %s (DDG %+.2f kcal/mol, AlphaMissense %s %.2f)",
					mutFocus.variantId.c_str(),
					mutFocus.ddg,
					mutFocus.amClass.c_str(),
					mutFocus.amScore
				);
			}
		}

		// ── Send logic ────────────────────────────────────────────────
		if ((sendPressed || submitOnEnter) && strlen(chatInputBuffer) > 0 && !claudePending) {
			std::string userMessage(chatInputBuffer);
			addChatMessage("You", userMessage);
			chatInputBuffer[0] = '\0';

			// ── Build MutationContext from current mutFocus state ─────
			MutationContext ctx;
			ctx.proteinName = currentPdbId.empty() ? std::nullopt
				: std::optional<std::string>(currentPdbId);

			if (mutFocus.show && !mutFocus.variantId.empty()) {
				if (!mutFocus.wtAA.empty())  ctx.wildTypeAA = mutFocus.wtAA;
				if (!mutFocus.mutAA.empty()) ctx.mutantAA = mutFocus.mutAA;
				if (mutFocus.position > 0)   ctx.residueNum = mutFocus.position;
				if (!mutFocus.chain.empty()) ctx.chain = mutFocus.chain;

				if (!mutFocus.analyzing) {
					ctx.amScore = mutFocus.amScore;
					ctx.amClass = mutFocus.amClass;
					ctx.ddg = mutFocus.ddg;
				}
			}

			// Replay last N turns for multi-turn context (cap at 6 to stay within tokens)
			{
				std::lock_guard<std::mutex> lock(chatMutex);
				const int MAX_HISTORY_TURNS = 6;
				int start = std::max(0, (int)chatHistory.size() - MAX_HISTORY_TURNS * 2);
				for (int i = start; i < (int)chatHistory.size(); ++i) {
					const auto& entry = chatHistory[i];
					if (entry.find("You: ") == 0)
						ctx.history.push_back({ "user", entry.substr(5) });
					else if (entry.find("Assistant: ") == 0)
						ctx.history.push_back({ "assistant", entry.substr(11) });
				}
			}

			// ── Fire async query ──────────────────────────────────────
			bool sent = claudeHandler.query(ctx, userMessage, nullptr); // poll-based, no callback
			if (sent) {
				claudePending = true;
			}
			else {
				addChatMessage("System", "A query is already in flight. Please wait.");
			}
		}

		ImGui::End();  // AI Assistant end

		// Mutation Panel Starts 
		if (mutationPanel.showWindow || mutationPanel.loading) {
			ImGui::SetNextWindowSize(ImVec2(420, 520), ImGuiCond_FirstUseEver);
			ImGui::SetNextWindowPos(ImVec2(20, 20), ImGuiCond_FirstUseEver);
			ImGui::Begin("Mutation Landscape", &mutationPanel.showWindow);

			// Status / progress bar
			{
				std::lock_guard<std::mutex> lock(mutationPanel.dataMutex);
				ImGui::TextUnformatted(mutationPanel.statusMsg.c_str());
			}

			if (mutationPanel.loading && mutationPanel.totalCount > 0) {
				float progress = (float)mutationPanel.loadedCount.load()
					/ (float)mutationPanel.totalCount;
				ImGui::ProgressBar(progress, ImVec2(-1, 0));
			}

			ImGui::Separator();

			// Legend
			ImGui::TextColored(ImVec4(1, 0, 0, 1), "■"); ImGui::SameLine();
			ImGui::Text("Pathogenic (>=0.564)"); ImGui::SameLine();
			ImGui::Spacing(); ImGui::SameLine();
			ImGui::TextColored(ImVec4(1, 1, 0, 1), "■"); ImGui::SameLine();
			ImGui::Text("Ambiguous"); ImGui::SameLine();
			ImGui::Spacing(); ImGui::SameLine();
			ImGui::TextColored(ImVec4(0, 1, 0, 1), "■"); ImGui::SameLine();
			ImGui::Text("Benign");
			ImGui::Separator();

			// Per-chain tables
			std::lock_guard<std::mutex> lock(mutationPanel.dataMutex);
			for (auto& [chain, chainData] : mutationPanel.chains) {
				std::string header = "Chain " + std::string(1, chain)
					+ "  (" + chainData.uniprotId + ")";
				if (!ImGui::CollapsingHeader(header.c_str(), ImGuiTreeNodeFlags_DefaultOpen))
					continue;

				if (ImGui::BeginTable("##mut", 4,
					ImGuiTableFlags_Borders |
					ImGuiTableFlags_RowBg |
					ImGuiTableFlags_ScrollY |
					ImGuiTableFlags_SizingStretchSame,
					ImVec2(0, 300)))
				{
					ImGui::TableSetupColumn("Residue", ImGuiTableColumnFlags_WidthFixed, 70);
					ImGui::TableSetupColumn("UniProt", ImGuiTableColumnFlags_WidthFixed, 70);
					ImGui::TableSetupColumn("Score", ImGuiTableColumnFlags_WidthFixed, 60);
					ImGui::TableSetupColumn("Worst", ImGuiTableColumnFlags_WidthStretch);
					ImGui::TableHeadersRow();

					for (auto& res : chainData.residues) {
						if (res.worstScore < 0) continue; // no data

						ImGui::TableNextRow();

						// Color row by classification
						ImVec4 color;
						if (res.worstScore >= 0.564f)
							color = ImVec4(1.0f, 0.3f, 0.3f, 1.0f);   // red
						else if (res.worstScore >= 0.340f)
							color = ImVec4(1.0f, 1.0f, 0.3f, 1.0f);   // yellow
						else
							color = ImVec4(0.3f, 1.0f, 0.3f, 1.0f);   // green

						// Residue number — clickable to select in 3D view
						ImGui::TableSetColumnIndex(0);
						std::string label = std::to_string(res.pdbResidueNum)
							+ "##" + std::string(1, chain)
							+ std::to_string(res.pdbResidueNum);
						ImGui::PushStyleColor(ImGuiCol_Text, color);
						if (ImGui::Selectable(label.c_str(), false,
							ImGuiSelectableFlags_SpanAllColumns)) {
							// Highlight this residue in 3D view
							Selection sel;
							sel.residue = res.pdbResidueNum;
							sel.chain = chain;
							Color highlight = Color::fromByte(255, 165, 0); // orange
							Model::setAtomColor(&highlight, &sel);
							Model::setAtomRadius(0.04f, &sel);
							selectedAtomIndex = std::nullopt; // clear so A-key works fresh
						}
						ImGui::PopStyleColor();

						ImGui::TableSetColumnIndex(1);
						ImGui::TextColored(color, "%d", res.uniprotPos);

						ImGui::TableSetColumnIndex(2);
						ImGui::TextColored(color, "%.3f", res.worstScore);

						ImGui::TableSetColumnIndex(3);
						ImGui::TextColored(color, "%s", res.worstVariant.c_str());
					}
					ImGui::EndTable();
				}
			}
			ImGui::End();
		}
		// Mutation Panel Ends

		renderMutationFocusPanel(mutFocus);
		renderComparisonPanel(structureComparison);


		// Mutant Viewer   
		if (mutantMoleculeData != nullptr) {
			ImGui::SetNextWindowPos(ImVec2(10.0f, 560.0f), ImGuiCond_FirstUseEver);
			ImGui::SetNextWindowSize(ImVec2(280.0f, 220.0f), ImGuiCond_FirstUseEver);
			ImGui::Begin("Mutant Viewer");

			ImGui::TextColored(ImVec4(1.0f, 0.4f, 0.8f, 1.0f),
				"Mutation: %s", mutFocus.variantId.c_str());
			ImGui::Separator();

			// DDG
			ImVec4 ddgColor = (mutFocus.ddg > 2.0f) ? ImVec4(0.9f, 0.2f, 0.2f, 1.f) :
				(mutFocus.ddg > 0.5f) ? ImVec4(0.9f, 0.6f, 0.2f, 1.f) :
				(mutFocus.ddg > -0.5f) ? ImVec4(0.8f, 0.8f, 0.8f, 1.f) :
				ImVec4(0.2f, 0.9f, 0.4f, 1.f);
			ImGui::TextColored(ddgColor, "DDG: %+.3f kcal/mol", mutFocus.ddg);
			ImGui::TextColored(ddgColor, "%s", mutFocus.ddgInterp.c_str());

			ImGui::Spacing();
			ImGui::Separator();
			ImGui::Spacing();

			// Toggle buttons
			ImGui::Text("Viewing:");
			ImGui::SameLine();

			if (!showingMutant) {
				ImGui::PushStyleColor(ImGuiCol_Button, ImVec4(0.1f, 0.4f, 0.1f, 1.f));
				ImGui::Button("Wildtype##active", ImVec2(90, 0));
				ImGui::PopStyleColor();
				ImGui::SameLine();
				if (ImGui::Button("Mutant##switch", ImVec2(90, 0))) {
					Model::loadMoleculeData(mutantMoleculeData);
					applyMutationHighlight(
						moleculeData, mutantMoleculeData,
						mutFocus.chain, mutFocus.position, true
					);
					showingMutant = true;
				}
			}
			else {
				if (ImGui::Button("Wildtype##switch", ImVec2(90, 0))) {
					Model::loadMoleculeData(moleculeData);
					applyMutationHighlight(
						moleculeData, mutantMoleculeData,
						mutFocus.chain, mutFocus.position, false
					);
					showingMutant = false;
				}
				ImGui::SameLine();
				ImGui::PushStyleColor(ImGuiCol_Button, ImVec4(0.4f, 0.1f, 0.4f, 1.f));
				ImGui::Button("Mutant##active", ImVec2(90, 0));
				ImGui::PopStyleColor();
			}

			ImGui::Spacing();
			ImGui::Separator();
			ImGui::Spacing();

			// Residue info
			ImGui::Text("Site: Chain %s · Pos %d",
				mutFocus.chain.c_str(), mutFocus.position);
			ImGui::Text("%s -> %s  (%s -> %s)",
				mutFocus.wtAA.c_str(), mutFocus.mutAA.c_str(),
				mutFocus.wtProps.c_str(), mutFocus.mutProps.c_str());

			ImGui::Spacing();

			if (ImGui::Button("Close##mutviewer", ImVec2(-1, 0))) {
				if (showingMutant && moleculeData) {
					Model::loadMoleculeData(moleculeData);
					showingMutant = false;
				}
				delete mutantMoleculeData;
				mutantMoleculeData = nullptr;
			}

			ImGui::End();
		}

		// Render Variant Selector
		if (showVariantPanel && !variantData.empty()) {

			renderVariantSelector(
				variantData, variantPositions,
				selectedPosition, selectedVariantId, runAnalysis,
				allChains, manualMut
			);
		}

		// Handler 1 - AM guided (Section 1)
		if (runAnalysis && !selectedVariantId.empty()) {
			runAnalysis = false;

			for (auto& v : variantData) {
				if (v.value("variant_id", "") == selectedVariantId) {

					std::string capPdbId = currentPdbId;
					std::string capRef = v.value("ref_aa", "");
					std::string capAlt = v.value("alt_aa", "");
					float  capScore = v.value("pathogenicity_score", 0.0f);
					std::string capClassification = v.value("classification", "ambiguous");
					int capUniPos = v.value("position", 0);
					std::string capVariantId = selectedVariantId;

					// Get chain + PDB residue from selected atom
					std::string capChain = "A";
					int capResNum = capUniPos;

					if (selectedAtomIndex && moleculeData) {
						const Atom& atom = moleculeData->atoms[*selectedAtomIndex];
						std::lock_guard<std::mutex> lock(mutationPanel.dataMutex);
						if (mutationPanel.chains.count(atom.chain)) {
							capChain = std::string(1, atom.chain);
							capResNum = atom.residueNum;
						}
					}

					if (capAlt.empty() || capRef.empty()) {
						std::cerr << "[ANALYZE] Missing ref/alt for " << capVariantId << "\n";
						break;
					}

					// AA property lookup
					static const std::unordered_map<char, std::string> oneToThree = {
						{'A',"ALA"},{'R',"ARG"},{'N',"ASN"},{'D',"ASP"},{'C',"CYS"},
						{'Q',"GLN"},{'E',"GLU"},{'G',"GLY"},{'H',"HIS"},{'I',"ILE"},
						{'L',"LEU"},{'K',"LYS"},{'M',"MET"},{'F',"PHE"},{'P',"PRO"},
						{'S',"SER"},{'T',"THR"},{'W',"TRP"},{'Y',"TYR"},{'V',"VAL"}
					};
					auto wtIt = oneToThree.find(capRef.empty() ? '?' : capRef[0]);
					auto mutIt = oneToThree.find(capAlt.empty() ? '?' : capAlt[0]);

					// Populate mutFocus immediately with AM data
					mutFocus.show = true;
					mutFocus.analyzing = true;
					mutFocus.variantId = capVariantId;
					mutFocus.wtAA = capRef;
					mutFocus.mutAA = capAlt;
					mutFocus.chain = capChain;
					mutFocus.position = capResNum;
					mutFocus.amScore = capScore;
					mutFocus.amClass = capClassification;
					mutFocus.ddg = 0.0f;
					mutFocus.ddgInterp = "Calculating...";
					mutFocus.wtProps = (wtIt != oneToThree.end())
						? getResidueProperties(wtIt->second) : "Unknown";
					mutFocus.mutProps = (mutIt != oneToThree.end())
						? getResidueProperties(mutIt->second) : "Unknown";

					std::cout << "[ANALYZE] Sending: "
						<< capRef << capChain << capResNum << capAlt
						<< "  AM=" << capScore << "\n";

					std::thread([=, &mutFocus]() {
						HTTPConnection http;
						json body = {
							{"pdb_id", capPdbId},
							{"chain",  capChain},
							{"position",  capResNum},
							{"wt_aa",  capRef},
							{"mut_aa", capAlt},
							{"alphamissense_score",  capScore}
						};

						if (http.post(BACKEND_URL + "/api/analyze_mutation", body.dump())) {
							try {
								auto report = json::parse(http.getResponse());

								// Update mutFocus with results
								if (report.contains("predictions")) {
									auto& p = report["predictions"];
									if (p.contains("foldx") && !p["foldx"].is_null()) {
										mutFocus.ddg = p["foldx"].value("ddg", 0.0f);
										mutFocus.ddgInterp = p["foldx"].value("interpretation", "");
									}
									if (p.contains("alphamissense"))
										mutFocus.amClass = p["alphamissense"].value("interpretation", "");
								}

								std::cout << "\n   Report: " << capVariantId << "   \n";
								if (report.contains("summary"))
									std::cout << "  " << report["summary"].get<std::string>() << "\n";
								std::cout << "  FoldX DDG: " << mutFocus.ddg
									<< " kcal/mol (" << mutFocus.ddgInterp << ")\n";
								std::cout << "  \n\n";

								// Load mutant PDB into viewer
								if (report.contains("mutant_pdb") && !report["mutant_pdb"].is_null()) {
									std::string mutantContent = report["mutant_pdb"].get<std::string>();

									// Write to temp file
									std::string tempPath = "./temp/" + capPdbId + "_mutant.pdb";

									// Create temp directory if needed
									std::filesystem::create_directories("./temp");

									std::ofstream f(tempPath);
									if (f.is_open()) {
										f << mutantContent;
										f.close();

										// Load into viewer — runs on render thread via command queue
										commands.push_back("loadmutant " + tempPath);
										std::cout << "[MUTANT] Queued load: " << tempPath << "\n";
									}
								}
							}
							catch (const std::exception& e) {
								std::cerr << "[ANALYZE] Parse error: " << e.what() << "\n";
								mutFocus.ddgInterp = "Parse error";
							}
						}
						mutFocus.analyzing = false;
						}).detach();

					break;
				}
			}
		}

		// Handler 2 — Manual FoldX (Section 2)
		if (manualMut.runManual && !allChains.empty()) {
			manualMut.runManual = false;

			// Guard: don't fire if no mutation AA entered
			if (manualMut.mutAA == '\0') {
				std::cerr << "[MANUAL] No mutation AA entered, aborting.\n";
				return;
			}

			auto& selChain = allChains[manualMut.selectedChainIdx];
			std::string capPdbId = currentPdbId;
			std::string capChainId = selChain.value("chain_id", "A");
			char   capMutAA = manualMut.mutAA;
			int capPdbPos = 0;
			std::string capWtAA = "?";

			// Primary: resolve from residue dropdown selection
			if (!chainResidues.empty() && manualMut.selectedResidueIdx < (int)chainResidues.size()) {
				capPdbPos = chainResidues[manualMut.selectedResidueIdx].first;
				std::string resName = chainResidues[manualMut.selectedResidueIdx].second;
				char one = threeLetterToOne(resName);
				if (one != '?') capWtAA = std::string(1, one);

				std::cout << "[MANUAL] Resolved from residue dropdown: chain=" << capChainId
					<< " residueNum=" << capPdbPos
					<< " residueName=" << resName
					<< " wt_aa=" << capWtAA << "\n";
			}
			// Fallback: resolve from selected atom in viewer
			else if (selectedAtomIndex && moleculeData) {
				const Atom& atom = moleculeData->atoms[*selectedAtomIndex];
				capChainId = std::string(1, atom.chain);
				capPdbPos = atom.residueNum;
				char one = threeLetterToOne(atom.residueName);
				if (one != '?') capWtAA = std::string(1, one);

				std::cout << "[MANUAL] Resolved from selected atom: chain=" << capChainId
					<< " residueNum=" << capPdbPos
					<< " residueName=" << atom.residueName
					<< " wt_aa=" << capWtAA << "\n";
			}
			// Last resort: offset math
			else {
				int capPdbStart = selChain.value("uniprot_start", 0);
				int capUniStart = selChain.value("pdb_start", 0);
				int capUniPos = manualMut.position;
				capPdbPos = capPdbStart + (capUniPos - capUniStart);

				std::cout << "[MANUAL] No residue/atom selected, using offset math: "
					<< "pdbStart=" << capPdbStart
					<< " uniStart=" << capUniStart
					<< " uniPos=" << capUniPos
					<< " -> capPdbPos=" << capPdbPos << "\n";

				if (moleculeData) {
					// Pass 1: CA atom
					for (const auto& atom : moleculeData->atoms) {
						if (atom.residueNum == capPdbPos &&
							std::string(1, atom.chain) == capChainId &&
							atom.name == "CA")
						{
							char one = threeLetterToOne(atom.residueName);
							if (one != '?') { capWtAA = std::string(1, one); break; }
						}
					}
					// Pass 2: any atom
					if (capWtAA == "?") {
						for (const auto& atom : moleculeData->atoms) {
							if (atom.residueNum == capPdbPos &&
								std::string(1, atom.chain) == capChainId)
							{
								char one = threeLetterToOne(atom.residueName);
								if (one != '?') { capWtAA = std::string(1, one); break; }
							}
						}
					}
				}
			}

			// Abort if wt_aa still unresolved
			if (capWtAA == "?") {
				std::cerr << "[MANUAL] Could not resolve wt_aa at chain="
					<< capChainId << " pdbPos=" << capPdbPos << " — aborting.\n";
				mutFocus.analyzing = false;
				return;
			}

			std::string capVariantId = capWtAA + capChainId +
				std::to_string(capPdbPos) +
				std::string(1, capMutAA);

			// AA properties
			static const std::unordered_map<char, std::string> oneToThree = {
				{'A',"ALA"},{'R',"ARG"},{'N',"ASN"},{'D',"ASP"},{'C',"CYS"},
				{'Q',"GLN"},{'E',"GLU"},{'G',"GLY"},{'H',"HIS"},{'I',"ILE"},
				{'L',"LEU"},{'K',"LYS"},{'M',"MET"},{'F',"PHE"},{'P',"PRO"},
				{'S',"SER"},{'T',"THR"},{'W',"TRP"},{'Y',"TYR"},{'V',"VAL"}
			};
			auto wtIt = oneToThree.find(capWtAA[0]);
			auto mutIt = oneToThree.find(capMutAA);

			// Populate mutFocus immediately
			mutFocus.show = true;
			mutFocus.analyzing = true;
			mutFocus.variantId = capVariantId;
			mutFocus.wtAA = capWtAA;
			mutFocus.mutAA = std::string(1, capMutAA);
			mutFocus.chain = capChainId;
			mutFocus.position = capPdbPos;
			mutFocus.amScore = 0.0f;
			mutFocus.amClass = "N/A";
			mutFocus.ddg = 0.0f;
			mutFocus.ddgInterp = "Calculating...";
			mutFocus.wtProps = (wtIt != oneToThree.end()) ? getResidueProperties(wtIt->second) : "Unknown";
			mutFocus.mutProps = (mutIt != oneToThree.end()) ? getResidueProperties(mutIt->second) : "Unknown";

			std::cout << "[MANUAL] Sending: " << capVariantId
				<< "  PDB pos=" << capPdbPos
				<< "  wt_aa=" << capWtAA << "\n";

			std::thread([=, &mutFocus]() {
				HTTPConnection http;
				json body = {
					{"pdb_id",   capPdbId},
					{"chain", capChainId},
					{"position", capPdbPos},
					{"wt_aa", capWtAA},
					{"mut_aa",   std::string(1, capMutAA)}
				};

				if (http.post(BACKEND_URL + "/api/analyze_mutation", body.dump())) {
					try {
						auto report = json::parse(http.getResponse());

						// FoldX predictions
						if (report.contains("predictions")) {
							auto& p = report["predictions"];

							if (p.contains("foldx") && !p["foldx"].is_null()) {
								mutFocus.ddg = p["foldx"].value("ddg", 0.0f);
								mutFocus.ddgInterp = p["foldx"].value("interpretation", "");
							}

							// AlphaMissense — backend returns this, we were ignoring it
							if (p.contains("alphamissense") && !p["alphamissense"].is_null()) {
								mutFocus.amClass = p["alphamissense"].value("interpretation", "N/A");
							}
						}

						// Summary
						if (report.contains("summary"))
							std::cout << "  " << report["summary"].get<std::string>() << "\n";

						std::cout << "\n-- Manual FoldX: " << capVariantId << " --\n";
						std::cout << "  DDG: " << mutFocus.ddg << " (" << mutFocus.ddgInterp << ")\n";
						//std::cout << "  AM class: " << mutFocus.amClass << "\n";
						std::cout << "-------------------------\n\n";

						// Load mutant structure for comparison 
						if (report.contains("mutant_pdb") &&
							report["mutant_pdb"].is_string() &&
							!report["mutant_pdb"].get<std::string>().empty())
						{
							mutantPdbContent = report["mutant_pdb"].get<std::string>();
							commands.push_back("loadmutant");
							std::cout << "\n[COMPARISON] Mutant PDB received — queued for comparison.\n";
						}
						else {
							std::cout << "\n[i] No mutant PDB in response — structure comparison unavailable.\n";
						}
					}
					catch (const std::exception& e) {
						std::cerr << "[MANUAL] Parse error: " << e.what() << "\n";
						mutFocus.ddgInterp = "Parse error";
					}
				}
				else {
					std::cerr << "[MANUAL] POST failed\n";
					mutFocus.ddgInterp = "Request failed";
				}
				mutFocus.analyzing = false;
				}).detach();
		}

		ImGui::Render();

		// ImGui::EndFrame();

		Window::clear(0, 0, 0);
		Model::render(Shader::getSphereDefault(), Shader::getConnectorDefault(), &window, &camera);

		float identityMatrix[16] = {
			1.0f, 0.0f, 0.0f, 0.0f,
			0.0f, 1.0f, 0.0f, 0.0f,
			0.0f, 0.0f, 1.0f, 0.0f,
			0.0f, 0.0f, 0.0f, 1.0f
			};

		// std::cout << "Ribbon vertices: " << representationRenderer.getRibbonVertexCount() << std::endl;
		// std::cout << "Surface vertices: " << representationRenderer.getSurfaceVertexCount() << std::endl;

		
		if (representationRenderer.getRibbonVertexCount() > 0 ||
			representationRenderer.getSurfaceVertexCount() > 0) {
			representationRenderer.render(ribbonShader, surfaceShader, &window, &camera, identityMatrix);
			//continue;
		}
		

		// Render ImGui draw data
		//ImGui_ImplOpenGL3_RenderDrawData(ImGui::GetDrawData());
		ImDrawData* draw_data = ImGui::GetDrawData();
		if (draw_data != nullptr) {
			ImGui_ImplOpenGL3_RenderDrawData(draw_data);
		}

		window.swapBuffers();

		// chatWindowThread polls this flag before initializing its own ImGui context.
		// Setting it here (end of first rendered frame) guarantees the main window's
		// font atlas and OpenGL state are fully committed before the chat window
		// attempts to share them.
		mainWindowReady = true;
	}
	ImGui::SetCurrentContext(mainImGuiContext);

	representationRenderer.cleanup();

	
	if (ribbonShader) {
		delete ribbonShader;
		ribbonShader = nullptr;
	}
	if (surfaceShader) {
		delete surfaceShader;
		surfaceShader = nullptr;
	}
	

	// Order matters: OpenGL backend must be shut down before GLFW backend, and both
	// before DestroyContext. Reversing this order causes a use-after-free on the
	// ImDrawData that GLFW holds a reference to.
	ImGui_ImplOpenGL3_Shutdown();
	ImGui_ImplGlfw_Shutdown();
	ImGui::DestroyContext();

	Shader::freeResources();
	ResourceManager::freeResources();

}

void buildAtomToResidueMapping(MoleculeData* moleculeData)
{
	atomIndexToDisplayInfo.clear();

	// Builds a flat index from atom array position → display info (residue name,
	// number, chain, atom name). This mirrors the GPU instance buffer layout, where
	// each instance index corresponds to the same atom index in moleculeData->atoms.
	// Required so the ray-cast hit index can be resolved to hover tooltip data
	// without traversing the full atom list on every frame.
	for (const Atom& atom : moleculeData->atoms) {
		HoverDisplayInfo info;
		info.residueName = atom.residueName;
		info.residueNum = atom.residueNum;
		info.chain = atom.chain;
		info.atomName = atom.name;
		info.element = atom.element;

		atomIndexToDisplayInfo.push_back(info);
	}

	std::cout << "Built hover mapping for " << atomIndexToDisplayInfo.size()
		<< " atoms" << std::endl;
}

/*
 *  Main function of the program.
 *  This function initializes the OpenAI API, creates a graphics thread,
 *  reads commands from the console, and waits for the graphics thread to finish.
 *  @return 0 if the program executed successfully.
 */
int main() {
	openai::start(apiKey); // Do not include real API key in production code

	std::cout << "=== Datalens Protein Visualization ===" << std::endl;
	std::cout << "Starting visualization and chat windows..." << std::endl;
	std::cout << "\nCommands:" << std::endl;
	std::cout << "  fetch pdb <pdb_id>  - Load a protein structure" << std::endl;
	std::cout << "  rotate <x|y|z> - Rotate model" << std::endl;
	std::cout << "  setdefault  - Reset to default view" << std::endl;
	std::cout << "  leave  - Clear current model" << std::endl;
	std::cout << "  mutate <original_AA> <position> <mutated_AA>  - Mutate current PDB with args" << std::endl;
	std::cout << "  exit   - Exit application" << std::endl;
	std::cout << "\nUse the separate chat window for AI assistance!\n" << std::endl;

	std::thread graphicsThread(renderingThread);

	// Handle console commands
	std::string command;
	while (std::getline(std::cin, command)) {
		if (command == "exit") {
			shouldExit = true;
			std::cout << "====> Exiting....\n" << std::endl;
			break;
		}
		if (command.size() > 0) {
			commands.push_back(command);
		}
	}

	graphicsThread.join();

	return 0;
}