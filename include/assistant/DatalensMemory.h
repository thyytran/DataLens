#pragma once

#include <string>
#include <vector>
#include <mutex>
#include <chrono>
#include <optional>
#include <fstream>
#include <iostream>
#include <filesystem>

#include "nlohmann/json.hpp"

using json = nlohmann::json;

// ─── Individual record types ──────────────────────────────────────────────────

struct ChatEntry {
    std::string role;       // "user" | "assistant"
    std::string content;
    std::string timestamp;  // ISO-8601
    std::string pdbContext; // which protein was loaded when this was said
};

struct MutationRecord {
    std::string timestamp;
    std::string pdbId;
    std::string uniprotId;
    std::string chain;
    int         pdbPosition = 0;
    int         uniprotPosition = 0;
    std::string wtAA;
    std::string mutAA;
    std::string variantId;   // e.g. "R175H"

    // Predictions
    float       amScore = 0.0f;
    std::string amClass;     // "pathogenic" | "benign" | "ambiguous"
    float       ddg = 0.0f;
    std::string ddgInterp;   // e.g. "Highly destabilizing"

    // Free-text notes (Claude summaries, user annotations)
    std::string claudeSummary;
};

struct SessionRecord {
    std::string timestamp;
    std::string pdbId;
    std::string uniprotId;
    int         mutationsAnalyzed = 0;
    std::string notes;
};

// ─── DataLensMemory ───────────────────────────────────────────────────────────

/*
 *  Persistent memory store for DataLens.
 *  Saves to memory.json alongside config.json on every write.
 *
 *  Thread-safe: all public methods lock memoryMutex.
 *
 *  Claude integration:
 *    Call buildMemoryContext() to get a formatted string injected
 *    into MutationContext::buildSystemPrompt().
 */
class DataLensMemory {
public:
    static const int MAX_CHAT_HISTORY = 40;  // turns kept in memory.json
    static const int MAX_MUTATION_HISTORY = 50; // mutation records kept
    static const int CLAUDE_CHAT_WINDOW = 10;  // turns sent to Claude per query
    static const int CLAUDE_MUT_WINDOW = 5;   // recent mutations sent to Claude

    explicit DataLensMemory(const std::string& path = "memory.json");

    // ── Load / save ──────────────────────────────────────────────────────────
    bool load();
    bool save() const;

    // ── Chat history ─────────────────────────────────────────────────────────
    void addChat(const std::string& role,
        const std::string& content,
        const std::string& pdbContext = "");

    // Returns last N turns as {role, content} pairs for Claude messages array
    std::vector<std::pair<std::string, std::string>>
        getRecentChatPairs(int n = CLAUDE_CHAT_WINDOW) const;

    // ── Mutation history ─────────────────────────────────────────────────────
    void addMutation(const MutationRecord& record);

    // Look up a previous analysis for the same variant (exact match)
    std::optional<MutationRecord>
        findMutation(const std::string& pdbId,
            const std::string& variantId) const;

    // Returns the N most recent records
    std::vector<MutationRecord>
        getRecentMutations(int n = CLAUDE_MUT_WINDOW) const;

    // ── Session history ───────────────────────────────────────────────────────
    void startSession(const std::string& pdbId, const std::string& uniprotId = "");
    void endSession(const std::string& notes = "");

    // ── Claude context builder ────────────────────────────────────────────────
    /*
     *  Returns a compact text block injected into the system prompt.
     *  Covers: proteins seen this session, recent mutations analyzed,
     *  and any prior Claude summaries for the current variant.
     */
    std::string buildMemoryContext(
        const std::string& currentPdbId = "",
        const std::string& currentVariant = "") const;

    // ── Attach Claude summary to most recent matching mutation ────────────────
    void attachClaudeSummary(const std::string& variantId,
        const std::string& summary);

    // ── Accessors ─────────────────────────────────────────────────────────────
    const std::vector<MutationRecord>& getMutations() const { return mutations; }
    const std::vector<ChatEntry>& getChats()     const { return chats; }

private:
    std::string filePath;
    mutable std::mutex memoryMutex;

    std::vector<ChatEntry>     chats;
    std::vector<MutationRecord> mutations;
    std::vector<SessionRecord>  sessions;

    std::string currentPdbId;
    SessionRecord activeSession;
    bool sessionOpen = false;

    static std::string nowISO();
    void trimToLimits();
};