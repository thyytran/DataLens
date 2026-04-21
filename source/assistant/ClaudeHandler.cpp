#include "assistant/ClaudeHandler.h"
#include "assistant/VariantDatabaseFetcher.h"

#include <iostream>
#include <fstream>
#include <sstream>
#include <iomanip>
#include <curl/curl.h>

// ─── Constants ───────────────────────────────────────────────────────────────

const std::string ClaudeHandler::MODEL = "claude-opus-4-6";
const int         ClaudeHandler::MAX_TOKENS = 1024;
const std::string ClaudeHandler::CLAUDE_API_URL =
"https://api.anthropic.com/v1/messages";

// ─── MutationContext helpers ──────────────────────────────────────────────────

std::string MutationContext::buildSystemPrompt() const {
    std::ostringstream ss;
    ss << "You are an expert protein bioinformatics assistant embedded in DataLens, "
        << "a real-time GPU-accelerated protein mutation analysis platform. "
        << "You interpret computational predictions and curated database evidence "
        << "to help researchers understand the biological significance of missense mutations.\n\n"

        << "Evidence hierarchy (highest to lowest reliability):\n"
        << "1. ClinVar   — expert-reviewed clinical classification (strongest)\n"
        << "2. COSMIC    — somatic tumor observations (strong functional evidence)\n"
        << "3. UniProt   — curated functional site annotation at the residue\n"
        << "4. gnomAD    — population allele frequency (low AF supports pathogenicity)\n"
        << "5. AlphaMissense — computational pathogenicity score (0=benign, 1=pathogenic)\n"
        << "6. FoldX ΔΔG — computational stability change (kcal/mol)\n\n"

        << "Thresholds:\n"
        << "AlphaMissense: <0.34 benign · 0.34-0.564 ambiguous · >0.564 pathogenic\n"
        << "FoldX: >2 kcal/mol destabilizing · <-1 stabilizing · -1 to 2 roughly neutral\n"
        << "gnomAD: AF >0.001 suggests benign · AF <0.0001 supports pathogenicity\n\n"

        << "REQUIRED RESPONSE FORMAT — always use this exact structure:\n\n"
        << "**Summary**\n"
        << "One or two sentences with the key clinical/functional takeaway.\n\n"
        << "**Evidence**\n"
        << "List each available evidence source on its own line in hierarchy order:\n"
        << "  [***] ClinVar: <classification> — <condition if known>\n"
        << "  [***] COSMIC: <N tumor samples> (<cancer types>)\n"
        << "  [**]  UniProt: <site annotation>\n"
        << "  [**]  gnomAD: AF=<value> — <interpretation>\n"
        << "  [*]   AlphaMissense: <score> (<class>)\n"
        << "  [*]   FoldX ΔΔG: <value> kcal/mol (<interpretation>)\n"
        << "Only include lines for evidence that was actually provided.\n"
        << "Star ratings: [***] = curated fact  [**] = curated annotation  [*] = prediction\n\n"
        << "**Mechanism**\n"
        << "1-3 sentences explaining the structural or functional reason this mutation matters. "
        << "Reference the secondary structure context if provided.\n\n"
        << "**Confidence**\n"
        << "One sentence: state whether evidence sources agree or conflict, "
        << "and flag if only computational predictions are available (no curated evidence).\n\n"
        << "**Disclaimer**\n"
        << "Always end with exactly this line: "
        << "'For research use only. Do not use for clinical or diagnostic decisions.'\n\n"
        << "Rules:\n"
        << "- Never omit the Disclaimer section\n"
        << "- Never recommend a clinical action\n"
        << "- Use 'predicted' not 'is' when describing computational results\n"
        << "- If no curated evidence (ClinVar/COSMIC/UniProt) is available, "
        << "say so explicitly in the Confidence section\n"
        << "- Add a disclaimer for user - do not base this for clinical decision.\n";

    return ss.str();
}

std::string MutationContext::buildUserMessage(const std::string& userQuery) const {
    std::ostringstream ss;

    // ── Variant identity ──────────────────────────────────────────────────────
    bool hasMutationData = residueNum || wildTypeAA || mutantAA || amScore || ddg;

    if (hasMutationData) {
        ss << "=== Current Variant ===\n";
        if (proteinName) ss << "Protein:  " << *proteinName << "\n";
        if (uniprotId)   ss << "UniProt:  " << *uniprotId << "\n";
        if (pdbId)       ss << "PDB:      " << *pdbId << "\n";
        if (chain)       ss << "Chain:    " << *chain << "\n";

        if (wildTypeAA && mutantAA && residueNum)
            ss << "Mutation: " << *wildTypeAA << *residueNum << *mutantAA << "\n";
        else if (residueNum)
            ss << "Residue:  " << *residueNum << "\n";

        // Secondary structure with structural reasoning hints
        if (secondaryStructure) {
            ss << "Location: residue sits in a " << *secondaryStructure << "\n";
            if (*secondaryStructure == "helix")
                ss << "  (consider: backbone H-bond disruption, helix-dipole effects, "
                << "packing against adjacent helices)\n";
            else if (*secondaryStructure == "sheet")
                ss << "  (consider: inter-strand H-bond loss, hydrophobic core packing, "
                << "edge vs. buried strand)\n";
            else if (*secondaryStructure == "loop" || *secondaryStructure == "coil")
                ss << "  (consider: surface exposure, active-site or binding-interface proximity, "
                << "conformational flexibility)\n";
        }
        ss << "\n";
    }

    // ── Curated database evidence ─────────────────────────────────────────────
    bool hasDatabaseEvidence = clinvarClass || uniprotSiteAnnotation ||
        gnomadAF || cosmicCount;

    if (hasDatabaseEvidence) {
        ss << "=== Curated Database Evidence ===\n";

        // ClinVar
        if (clinvarClass) {
            ss << "ClinVar classification: " << *clinvarClass;
            if (clinvarStars)    ss << " (" << *clinvarStars << "-star review)";
            if (clinvarId)       ss << "  [" << *clinvarId << "]";
            ss << "\n";
            if (clinvarCondition)
                ss << "  Associated condition: " << *clinvarCondition << "\n";
        }

        // UniProt site annotation
        if (uniprotSiteAnnotation) {
            ss << "UniProt site annotation: " << *uniprotSiteAnnotation << "\n";
            ss << "  (this residue position has curated functional importance)\n";
        }

        // gnomAD population frequency
        if (gnomadAF) {
            ss << "gnomAD allele frequency: ";
            // Format as scientific notation for very small values
            if (*gnomadAF < 0.0001f && *gnomadAF > 0.0f) {
                ss << std::scientific << std::setprecision(2) << *gnomadAF;
                ss << std::defaultfloat;
            }
            else {
                ss << std::fixed << std::setprecision(6) << *gnomadAF;
                ss << std::defaultfloat;
            }
            if (gnomadAC)      ss << " (AC=" << *gnomadAC << ")";
            if (gnomadPopmax)  ss << ", highest in " << *gnomadPopmax << " population";
            ss << "\n";

            // Interpret frequency for Claude
            if (*gnomadAF == 0.0f)
                ss << "  [not observed in gnomAD — supports pathogenicity]\n";
            else if (*gnomadAF < 0.0001f)
                ss << "  [extremely rare — consistent with pathogenic variant]\n";
            else if (*gnomadAF < 0.001f)
                ss << "  [rare in general population]\n";
            else
                ss << "  [present at appreciable frequency — may suggest benign/VUS]\n";
        }

        // COSMIC somatic observations
        if (cosmicCount) {
            ss << "COSMIC somatic observations: " << *cosmicCount << " tumor samples";
            if (cosmicId)          ss << "  [" << *cosmicId << "]";
            ss << "\n";
            if (cosmicCancerTypes) ss << "  Cancer types: " << *cosmicCancerTypes << "\n";
            if (*cosmicCount > 100)
                ss << "  [high recurrence — strong evidence of driver/functional consequence]\n";
            else if (*cosmicCount > 10)
                ss << "  [moderate recurrence — likely functionally relevant in cancer]\n";
        }

        ss << "\n";
    }

    // ── Computational predictions ─────────────────────────────────────────────
    bool hasPredictions = amScore || ddg || uniprotPolyphenPred || uniprotSiftPred;

    if (hasPredictions) {
        ss << "=== Computational Predictions ===\n";

        if (amScore) {
            ss << "AlphaMissense: " << std::fixed << std::setprecision(3) << *amScore;
            if (amClass) ss << " (" << *amClass << ")";
            ss << "\n";
        }

        if (ddg) {
            ss << "FoldX ΔΔG: " << std::fixed << std::setprecision(2) << *ddg << " kcal/mol";
            if (*ddg > 2.0f)       ss << " [destabilizing]";
            else if (*ddg < -1.0f) ss << " [stabilizing]";
            else                   ss << " [roughly neutral]";
            ss << "\n";
        }


        // UniProt variant predictions (PolyPhen, SIFT)
        if (uniprotPolyphenPred || uniprotSiftPred) {
            ss << "\n-- UniProt Variant Predictions --\n";

            if (uniprotPolyphenPred) {
                ss << "PolyPhen: " << *uniprotPolyphenPred;
                if (uniprotPolyphenScore)
                    ss << " (score: " << std::fixed << std::setprecision(3) << *uniprotPolyphenScore << ")";
                ss << "\n";
            }

            if (uniprotSiftPred) {
                ss << "SIFT: " << *uniprotSiftPred;
                if (uniprotSiftScore)
                    ss << " (score: " << std::fixed << std::setprecision(3) << *uniprotSiftScore << ")";
                ss << "\n";
            }

            if (uniprotConsequence)
                ss << "Consequence type: " << *uniprotConsequence << "\n";

            if (uniprotCodon)
                ss << "Codon change: " << *uniprotCodon << "\n";

            if (uniprotSomatic)
                ss << "Somatic: " << (*uniprotSomatic ? "Yes" : "No") << "\n";
        }

        // UniProt clinical/disease data
        if (uniprotClinicalSignificance) {
            ss << "UniProt clinical significance: " << *uniprotClinicalSignificance << "\n";
        }

        if (uniprotDiseaseAssociation) {
            ss << "UniProt disease associations: " << *uniprotDiseaseAssociation << "\n";
        }

        // After the other database evidence sections (around line 169):
        // UniProt variant summary (pre-formatted)
        if (uniprotVariantSummary && !uniprotVariantSummary->empty()) {
            ss << "\n-- UniProt Variant Data --\n";
            ss << *uniprotVariantSummary << "\n";
        }

        ss << "\n";
    }

    // ── Memory context from previous sessions ─────────────────────────────────
    if (!memoryContext.empty()) {
        ss << "=== Session Memory ===\n" << memoryContext << "\n\n";
    }

    // ── User question ─────────────────────────────────────────────────────────
    ss << "=== User Question ===\n" << userQuery;
    return ss.str();
}

// ─── ClaudeHandler ────────────────────────────────────────────────────────────

ClaudeHandler::ClaudeHandler() {}

ClaudeHandler::~ClaudeHandler() {
    if (workerThread.joinable()) {
        workerThread.join();
    }
}

bool ClaudeHandler::loadConfig(const std::string& configPath) {
    std::ifstream file(configPath);
    if (!file.is_open()) {
        std::cerr << "ClaudeHandler: cannot open " << configPath << "\n";
        return false;
    }

    try {
        json cfg;
        file >> cfg;
        apiKey = cfg.at("api_key").get<std::string>();
    }
    catch (const std::exception& e) {
        std::cerr << "ClaudeHandler: config parse error: " << e.what() << "\n";
        return false;
    }

    if (apiKey.empty() || apiKey == "YOUR_CLAUDE_API_KEY_HERE") {
        std::cerr << "ClaudeHandler: claude api_key not set in " << configPath << "\n";
        return false;
    }

    return true;
}

void ClaudeHandler::setApiKey(const std::string& key) {
    apiKey = key;
}

void ClaudeHandler::setCosmicToken(const std::string& bearerToken) {
    cosmicToken = bearerToken;
}

bool ClaudeHandler::query(
    const MutationContext& context,
    const std::string& userMessage,
    std::function<void(std::string)> onComplete
) {
    if (requestActive.load()) return false;

    if (workerThread.joinable()) {
        workerThread.join();
    }

    requestActive.store(true);
    resultReady.store(false);

    workerThread = std::thread([this, context, userMessage, onComplete,
        cosmicTok = cosmicToken]() mutable {

            // fetchAll writes database results into ctx — needs a non-const copy.
            // The lambda already captures context by value; bind a non-const ref to it.
            MutationContext ctx = context;

            VariantDatabaseFetcher fetcher;
            fetcher.setTimeoutSeconds(8);
            if (!cosmicTok.empty()) fetcher.setCosmicToken(cosmicTok);
            fetcher.fetchAll(ctx);   // enriches ctx in-place

            // ── Build messages array ──────────────────────────────────────────────
            json messages = json::array();

            for (const auto& [role, content] : ctx.history) {
                messages.push_back({ {"role", role}, {"content", content} });
            }

            std::string fullUserContent = ctx.buildUserMessage(userMessage);
            messages.push_back({ {"role", "user"}, {"content", fullUserContent} });

            // ── Assemble request body ─────────────────────────────────────────────
            json body = {
                { "model",      MODEL      },
                { "max_tokens", MAX_TOKENS },
                { "system",     ctx.buildSystemPrompt() },
                { "messages",   messages   }
            };

            // ── POST and extract text ─────────────────────────────────────────────
            std::string rawResponse = postToClaudeAPI(body);
            std::string reply;

            if (!rawResponse.empty()) {
                try {
                    json parsed = json::parse(rawResponse);

                    if (parsed.value("type", "") == "error" && parsed.contains("error")) {
                        std::string errType = parsed["error"].value("type", "unknown");
                        std::string errMsg = parsed["error"].value("message", rawResponse);
                        reply = "[API error – " + errType + ": " + errMsg + "]";
                        std::cerr << "ClaudeHandler API error: " << errType << ": " << errMsg << "\n";
                    }
                    else if (parsed.contains("content")
                        && parsed["content"].is_array()
                        && !parsed["content"].empty()
                        && parsed["content"][0].contains("text")
                        && parsed["content"][0]["text"].is_string())
                    {
                        reply = parsed["content"][0]["text"].get<std::string>();
                    }
                    else {
                        std::cerr << "ClaudeHandler: unexpected response shape:\n"
                            << parsed.dump(2) << "\n";
                        reply = "[Unexpected response format]";
                    }
                }
                catch (const std::exception& e) {
                    std::cerr << "ClaudeHandler parse error: " << e.what()
                        << "\nRaw: " << rawResponse.substr(0, 300) << "\n";
                    reply = "[Parse error: " + std::string(e.what()) + "]";
                }
            }
            else {
                reply = "[No response from Claude API]";
            }

            // ── Store + notify ────────────────────────────────────────────────────
            {
                std::lock_guard<std::mutex> lock(resultMutex);
                pendingResult = reply;
            }
            resultReady.store(true);
            requestActive.store(false);

            if (onComplete) onComplete(reply);
        });

    return true;
}

bool ClaudeHandler::isReady() const {
    return resultReady.load();
}

bool ClaudeHandler::isInflight() const {
    return requestActive.load();
}

std::string ClaudeHandler::getResponse() {
    if (!resultReady.load()) return "";
    resultReady.store(false);
    std::lock_guard<std::mutex> lock(resultMutex);
    return pendingResult;
}

// ─── Private ──────────────────────────────────────────────────────────────────

size_t ClaudeHandler::curlWriteCallback(
    char* ptr, size_t size, size_t nmemb, void* userdata
) {
    auto* buf = static_cast<std::string*>(userdata);
    buf->append(ptr, size * nmemb);
    return size * nmemb;
}

std::string ClaudeHandler::postToClaudeAPI(const json& requestBody) {
    CURL* curl = curl_easy_init();
    if (!curl) {
        std::cerr << "ClaudeHandler: curl_easy_init failed\n";
        return "";
    }

    std::string responseBuffer;
    std::string bodyStr = requestBody.dump(
        -1, ' ', false, nlohmann::json::error_handler_t::replace
    );

    struct curl_slist* headers = nullptr;
    std::string authHeader = "x-api-key: " + apiKey;
    headers = curl_slist_append(headers, "Content-Type: application/json");
    headers = curl_slist_append(headers, authHeader.c_str());
    headers = curl_slist_append(headers, "anthropic-version: 2023-06-01");

    curl_easy_setopt(curl, CURLOPT_URL, CLAUDE_API_URL.c_str());
    curl_easy_setopt(curl, CURLOPT_HTTPHEADER, headers);
    curl_easy_setopt(curl, CURLOPT_POSTFIELDS, bodyStr.c_str());
    curl_easy_setopt(curl, CURLOPT_POSTFIELDSIZE, (long)bodyStr.size());
    curl_easy_setopt(curl, CURLOPT_WRITEFUNCTION, curlWriteCallback);
    curl_easy_setopt(curl, CURLOPT_WRITEDATA, &responseBuffer);
    curl_easy_setopt(curl, CURLOPT_SSL_VERIFYPEER, 1L);
    curl_easy_setopt(curl, CURLOPT_SSL_VERIFYHOST, 2L);
    curl_easy_setopt(curl, CURLOPT_CAINFO, "cacert.pem");

    CURLcode result = curl_easy_perform(curl);
    if (result != CURLE_OK) {
        std::cerr << "ClaudeHandler: request failed: "
            << curl_easy_strerror(result) << "\n";
        responseBuffer.clear();
    }

    curl_slist_free_all(headers);
    curl_easy_cleanup(curl);
    return responseBuffer;
}