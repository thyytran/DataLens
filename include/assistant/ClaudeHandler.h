#pragma once

#include <string>
#include <thread>
#include <mutex>
#include <atomic>
#include <functional>
#include <vector>
#include <optional>

#include "nlohmann/json.hpp"

using json = nlohmann::json;

/*
 *  Holds all protein/mutation context sent to Claude as grounding data.
 *  Fill whichever fields are available; unset fields are omitted from the prompt.
 */
struct MutationContext {
    std::optional<std::string> proteinName;    // e.g. "TP53"
    std::optional<std::string> uniprotId;      // e.g. "P04637"
    std::optional<std::string> pdbId;          // e.g. "7LMK"
    std::optional<std::string> chain;
    std::optional<int>         residueNum;
    std::optional<std::string> wildTypeAA;     // single-letter, e.g. "R"
    std::optional<std::string> mutantAA;       // single-letter, e.g. "H"

    // AlphaMissense
    std::optional<float>       amScore;        // 0.0 – 1.0
    std::optional<std::string> amClass;        // "pathogenic" | "benign" | "ambiguous"

    // FoldX
    std::optional<float>       ddg;            // ΔΔG in kcal/mol (positive = destabilizing)

    // Injected by DataLensMemory — previous sessions, mutations, prior summaries
    std::string memoryContext;

    // Conversation history for multi-turn chat
    std::vector<std::pair<std::string, std::string>> history; // {role, content}

    std::string buildSystemPrompt() const;
    std::string buildUserMessage(const std::string& userQuery) const;
};

/*
 *  Async, non-blocking Claude API handler.
 *
 *  Usage:
 *      handler.query(context, "What does this mutation mean clinically?",
 *          [](const std::string& reply) {
 *              // called on worker thread — guard shared state with mutex
 *          });
 *
 *  Alternatively, poll with isReady() / getResponse().
 */
class ClaudeHandler {
public:
    static const std::string MODEL;
    static const int         MAX_TOKENS;
    static const std::string CLAUDE_API_URL;

    ClaudeHandler();
    ~ClaudeHandler();

    /* Load API key from config.json — call once at startup */
    bool loadConfig(const std::string& configPath = "config.json");

    /* Direct setter — fallback if config file path is unreliable */
    void setApiKey(const std::string& key);

    /*
     *  Fire-and-forget async query.
     *  onComplete is called on the worker thread when the response arrives.
     *  Returns false immediately if a query is already in flight.
     */
    bool query(
        const MutationContext& context,
        const std::string& userMessage,
        std::function<void(std::string)>    onComplete
    );

    /* Polling API — alternative to callback */
    bool        isReady()     const;   // true when a result is waiting
    bool        isInflight()  const;   // true while request is pending
    std::string getResponse();         // clears ready flag; returns "" if not ready

private:
    std::string  apiKey;

    mutable std::mutex  resultMutex;
    std::string         pendingResult;
    std::atomic<bool>   resultReady{ false };
    std::atomic<bool>   requestActive{ false };

    std::thread workerThread;

    /* Synchronous HTTP POST — runs on worker thread */
    std::string postToClaudeAPI(const json& requestBody);

    static size_t curlWriteCallback(char* ptr, size_t size, size_t nmemb, void* userdata);
};