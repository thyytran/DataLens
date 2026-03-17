#include "assistant/ClaudeHandler.h"

#include <iostream>
#include <fstream>
#include <sstream>
#include <curl/curl.h>

// ─── Constants ───────────────────────────────────────────────────────────────

const std::string ClaudeHandler::MODEL = "claude-opus-4-6";
const int         ClaudeHandler::MAX_TOKENS = 1024;
const std::string ClaudeHandler::CLAUDE_API_URL =
"https://api.anthropic.com/v1/messages";

// ─── MutationContext helpers ─────────────────────────────────────────────────

std::string MutationContext::buildSystemPrompt() const {
    std::ostringstream ss;
    ss << "You are an expert protein bioinformatics assistant embedded in DataLens, "
        << "a real-time GPU-accelerated protein mutation analysis platform. "
        << "You interpret computational predictions to help researchers understand "
        << "the biological significance of missense mutations.\n\n"
        << "Prediction tools in use:\n"
        << "- AlphaMissense: pathogenicity scoring (0=benign, 1=pathogenic; "
        << "threshold: <0.34 benign, >0.564 pathogenic, ambiguous in between)\n"
        << "- FoldX: structural stability change (ΔΔG kcal/mol; "
        << ">2 = significantly destabilizing, <-1 = stabilizing, "
        << "-1 to 2 = roughly neutral)\n\n"
        << "Response style:\n"
        << "- Lead with the key clinical/functional takeaway in 1-2 sentences\n"
        << "- Then explain mechanism (structural, functional, or evolutionary reasoning)\n"
        << "- Flag uncertainty honestly; distinguish computational prediction from "
        << "experimental evidence\n"
        << "- Be concise — this appears in a GUI panel next to a 3D protein viewer\n";

    return ss.str();
}

std::string MutationContext::buildUserMessage(const std::string& userQuery) const {
    std::ostringstream ss;

    // Build mutation identifier
    bool hasMutationData = residueNum || wildTypeAA || mutantAA ||
        amScore || ddg;

    if (hasMutationData) {
        ss << "=== Current Variant Context ===\n";

        if (proteinName) ss << "Protein:   " << *proteinName << "\n";
        if (uniprotId)   ss << "UniProt:   " << *uniprotId << "\n";
        if (pdbId)       ss << "PDB:       " << *pdbId << "\n";
        if (chain)       ss << "Chain:     " << *chain << "\n";

        if (wildTypeAA && mutantAA && residueNum) {
            ss << "Mutation:  " << *wildTypeAA << *residueNum << *mutantAA << "\n";
        }
        else if (residueNum) {
            ss << "Residue:   " << *residueNum << "\n";
        }

        if (amScore) {
            ss << "AlphaMissense score: " << *amScore;
            if (amClass) ss << " (" << *amClass << ")";
            ss << "\n";
        }

        if (ddg) {
            ss << "FoldX ΔΔG: " << *ddg << " kcal/mol";
            if (*ddg > 2.0f)       ss << " [destabilizing]";
            else if (*ddg < -1.0f) ss << " [stabilizing]";
            else                   ss << " [roughly neutral]";
            ss << "\n";
        }

        ss << "\n";
    }

    ss << "=== User Question ===\n" << userQuery;
    return ss.str();
}

// ─── ClaudeHandler ───────────────────────────────────────────────────────────

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

bool ClaudeHandler::query(
    const MutationContext& context,
    const std::string& userMessage,
    std::function<void(std::string)> onComplete
) {
    if (requestActive.load()) {
        return false;  // one request at a time
    }

    // Join any finished previous thread
    if (workerThread.joinable()) {
        workerThread.join();
    }

    requestActive.store(true);
    resultReady.store(false);

    workerThread = std::thread([this, context, userMessage, onComplete]() {

        // ── Build messages array ──────────────────────────────────────────────
        json messages = json::array();

        // Replay conversation history
        for (const auto& [role, content] : context.history) {
            messages.push_back({ {"role", role}, {"content", content} });
        }

        // Current user turn with mutation context injected
        std::string fullUserContent = context.buildUserMessage(userMessage);
        messages.push_back({ {"role", "user"}, {"content", fullUserContent} });

        // ── Assemble request body ─────────────────────────────────────────────
        json body = {
            { "model",      MODEL      },
            { "max_tokens", MAX_TOKENS },
            { "system",     context.buildSystemPrompt() },
            { "messages",   messages   }
        };

        // ── POST and extract text ─────────────────────────────────────────────
        std::string rawResponse = postToClaudeAPI(body);
        std::string reply;

        if (!rawResponse.empty()) {
            try {
                json parsed = json::parse(rawResponse);

                // Claude error response: { "type": "error", "error": { "type": "...", "message": "..." } }
                if (parsed.value("type", "") == "error" && parsed.contains("error")) {
                    std::string errType = parsed["error"].value("type", "unknown");
                    std::string errMsg = parsed["error"].value("message", rawResponse);
                    reply = "[API error – " + errType + ": " + errMsg + "]";
                    std::cerr << "ClaudeHandler API error: " << errType << ": " << errMsg << "\n";
                }
                // Normal response: { "content": [ { "type": "text", "text": "..." } ] }
                else if (parsed.contains("content")
                    && parsed["content"].is_array()
                    && !parsed["content"].empty()
                    && parsed["content"][0].contains("text")
                    && parsed["content"][0]["text"].is_string())
                {
                    reply = parsed["content"][0]["text"].get<std::string>();
                }
                else {
                    // Unexpected shape — log raw for debugging
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

        if (onComplete) {
            onComplete(reply);
        }
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

void ClaudeHandler::setApiKey(const std::string& key) {
    apiKey = key;
}

// ─── Private ─────────────────────────────────────────────────────────────────

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
    // std::string bodyStr = requestBody.dump();
    std::string bodyStr = requestBody.dump(-1, ' ', false, nlohmann::json::error_handler_t::replace);

    // ── Headers ───────────────────────────────────────────────────────────────
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