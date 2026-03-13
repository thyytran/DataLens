#include "assistant/DataLensMemory.h"

#include <algorithm>
#include <sstream>
#include <iomanip>
#include <ctime>

// ─── Helpers ──────────────────────────────────────────────────────────────────

std::string DataLensMemory::nowISO() {
    auto now = std::chrono::system_clock::now();
    std::time_t t = std::chrono::system_clock::to_time_t(now);
    std::ostringstream ss;
    ss << std::put_time(std::gmtime(&t), "%Y-%m-%dT%H:%M:%SZ");
    return ss.str();
}

// ─── Constructor ─────────────────────────────────────────────────────────────

DataLensMemory::DataLensMemory(const std::string& path) : filePath(path) {}

// ─── Load / Save ─────────────────────────────────────────────────────────────

bool DataLensMemory::load() {
    std::lock_guard<std::mutex> lock(memoryMutex);

    std::ifstream file(filePath);
    if (!file.is_open()) {
        std::cerr << "[Memory] No existing memory file at \"" << filePath
            << "\" — starting fresh.\n";
        return false;
    }

    try {
        json j;
        file >> j;

        // Chat history
        chats.clear();
        for (const auto& c : j.value("chat_history", json::array())) {
            ChatEntry e;
            e.role = c.value("role", "");
            e.content = c.value("content", "");
            e.timestamp = c.value("timestamp", "");
            e.pdbContext = c.value("pdb_context", "");
            chats.push_back(std::move(e));
        }

        // Mutation history
        mutations.clear();
        for (const auto& m : j.value("mutation_history", json::array())) {
            MutationRecord r;
            r.timestamp = m.value("timestamp", "");
            r.pdbId = m.value("pdb_id", "");
            r.uniprotId = m.value("uniprot_id", "");
            r.chain = m.value("chain", "");
            r.pdbPosition = m.value("pdb_position", 0);
            r.uniprotPosition = m.value("uniprot_position", 0);
            r.wtAA = m.value("wt_aa", "");
            r.mutAA = m.value("mut_aa", "");
            r.variantId = m.value("variant_id", "");
            r.amScore = m.value("am_score", 0.0f);
            r.amClass = m.value("am_class", "");
            r.ddg = m.value("ddg", 0.0f);
            r.ddgInterp = m.value("ddg_interp", "");
            r.claudeSummary = m.value("claude_summary", "");
            mutations.push_back(std::move(r));
        }

        // Session history
        sessions.clear();
        for (const auto& s : j.value("session_history", json::array())) {
            SessionRecord sr;
            sr.timestamp = s.value("timestamp", "");
            sr.pdbId = s.value("pdb_id", "");
            sr.uniprotId = s.value("uniprot_id", "");
            sr.mutationsAnalyzed = s.value("mutations_analyzed", 0);
            sr.notes = s.value("notes", "");
            sessions.push_back(std::move(sr));
        }

        std::cerr << "[Memory] Loaded: " << chats.size() << " chat entries, "
            << mutations.size() << " mutations, "
            << sessions.size() << " sessions.\n";
        return true;
    }
    catch (const std::exception& e) {
        std::cerr << "[Memory] Load error: " << e.what() << "\n";
        return false;
    }
}

bool DataLensMemory::save() const {
    // Caller must hold memoryMutex or call under lock
    json j;

    json chatArr = json::array();
    for (const auto& c : chats) {
        chatArr.push_back({
            {"role",        c.role},
            {"content",     c.content},
            {"timestamp",   c.timestamp},
            {"pdb_context", c.pdbContext}
            });
    }
    j["chat_history"] = chatArr;

    json mutArr = json::array();
    for (const auto& m : mutations) {
        mutArr.push_back({
            {"timestamp",        m.timestamp},
            {"pdb_id",           m.pdbId},
            {"uniprot_id",       m.uniprotId},
            {"chain",            m.chain},
            {"pdb_position",     m.pdbPosition},
            {"uniprot_position", m.uniprotPosition},
            {"wt_aa",            m.wtAA},
            {"mut_aa",           m.mutAA},
            {"variant_id",       m.variantId},
            {"am_score",         m.amScore},
            {"am_class",         m.amClass},
            {"ddg",              m.ddg},
            {"ddg_interp",       m.ddgInterp},
            {"claude_summary",   m.claudeSummary}
            });
    }
    j["mutation_history"] = mutArr;

    json sessArr = json::array();
    for (const auto& s : sessions) {
        sessArr.push_back({
            {"timestamp",          s.timestamp},
            {"pdb_id",             s.pdbId},
            {"uniprot_id",         s.uniprotId},
            {"mutations_analyzed", s.mutationsAnalyzed},
            {"notes",              s.notes}
            });
    }
    j["session_history"] = sessArr;

    try {
        std::ofstream file(filePath);
        if (!file.is_open()) {
            std::cerr << "[Memory] Cannot write to \"" << filePath << "\"\n";
            return false;
        }
        file << j.dump(2);
        return true;
    }
    catch (const std::exception& e) {
        std::cerr << "[Memory] Save error: " << e.what() << "\n";
        return false;
    }
}

// ─── Chat ─────────────────────────────────────────────────────────────────────

void DataLensMemory::addChat(
    const std::string& role,
    const std::string& content,
    const std::string& pdbContext
) {
    std::lock_guard<std::mutex> lock(memoryMutex);

    ChatEntry e;
    e.role = role;
    e.content = content;
    e.timestamp = nowISO();
    e.pdbContext = pdbContext.empty() ? currentPdbId : pdbContext;
    chats.push_back(std::move(e));

    trimToLimits();
    save();
}

std::vector<std::pair<std::string, std::string>>
DataLensMemory::getRecentChatPairs(int n) const {
    std::lock_guard<std::mutex> lock(memoryMutex);

    std::vector<std::pair<std::string, std::string>> result;
    int start = std::max(0, (int)chats.size() - n);
    for (int i = start; i < (int)chats.size(); ++i) {
        // Skip system-style entries (no role / empty content)
        if (chats[i].role.empty() || chats[i].content.empty()) continue;
        result.push_back({ chats[i].role, chats[i].content });
    }
    return result;
}

// ─── Mutations ────────────────────────────────────────────────────────────────

void DataLensMemory::addMutation(const MutationRecord& record) {
    std::lock_guard<std::mutex> lock(memoryMutex);

    MutationRecord r = record;
    if (r.timestamp.empty()) r.timestamp = nowISO();

    // Update existing record for same pdbId + variantId rather than duplicating
    for (auto& existing : mutations) {
        if (existing.pdbId == r.pdbId && existing.variantId == r.variantId) {
            existing = r;
            save();
            return;
        }
    }

    mutations.push_back(std::move(r));

    if (sessionOpen) {
        activeSession.mutationsAnalyzed++;
    }

    trimToLimits();
    save();
}

std::optional<MutationRecord>
DataLensMemory::findMutation(
    const std::string& pdbId,
    const std::string& variantId
) const {
    std::lock_guard<std::mutex> lock(memoryMutex);

    // Search from most recent backwards
    for (int i = (int)mutations.size() - 1; i >= 0; --i) {
        if (mutations[i].pdbId == pdbId && mutations[i].variantId == variantId)
            return mutations[i];
    }
    return std::nullopt;
}

std::vector<MutationRecord>
DataLensMemory::getRecentMutations(int n) const {
    std::lock_guard<std::mutex> lock(memoryMutex);

    int start = std::max(0, (int)mutations.size() - n);
    return std::vector<MutationRecord>(mutations.begin() + start, mutations.end());
}

// ─── Sessions ─────────────────────────────────────────────────────────────────

void DataLensMemory::startSession(
    const std::string& pdbId,
    const std::string& uniprotId
) {
    std::lock_guard<std::mutex> lock(memoryMutex);

    if (sessionOpen) {
        // Close previous session silently
        sessions.push_back(activeSession);
    }

    activeSession = SessionRecord();
    activeSession.timestamp = nowISO();
    activeSession.pdbId = pdbId;
    activeSession.uniprotId = uniprotId;
    currentPdbId = pdbId;
    sessionOpen = true;
    save();
}

void DataLensMemory::endSession(const std::string& notes) {
    std::lock_guard<std::mutex> lock(memoryMutex);

    if (!sessionOpen) return;
    activeSession.notes = notes;
    sessions.push_back(activeSession);
    sessionOpen = false;
    save();
}

// ─── Claude context builder ───────────────────────────────────────────────────

std::string DataLensMemory::buildMemoryContext(
    const std::string& currentPdbId_,
    const std::string& currentVariant
) const {
    std::lock_guard<std::mutex> lock(memoryMutex);

    std::ostringstream ss;
    ss << "=== DataLens Memory Context ===\n";

    // 1. Session history — proteins the researcher has worked with
    if (!sessions.empty()) {
        ss << "\nPrevious sessions (most recent first):\n";
        int shown = 0;
        for (int i = (int)sessions.size() - 1; i >= 0 && shown < 3; --i, ++shown) {
            const auto& s = sessions[i];
            ss << "  [" << s.timestamp.substr(0, 10) << "] "
                << s.pdbId;
            if (!s.uniprotId.empty()) ss << " (" << s.uniprotId << ")";
            ss << " — " << s.mutationsAnalyzed << " mutations analyzed";
            if (!s.notes.empty()) ss << "; " << s.notes;
            ss << "\n";
        }
    }

    // 2. Recent mutation analyses — what has been explored already
    int mutStart = std::max(0, (int)mutations.size() - CLAUDE_MUT_WINDOW);
    if (mutStart < (int)mutations.size()) {
        ss << "\nRecent mutation analyses:\n";
        for (int i = mutStart; i < (int)mutations.size(); ++i) {
            const auto& m = mutations[i];
            ss << "  " << m.variantId
                << " [" << m.pdbId << " chain " << m.chain << "]"
                << "  AM=" << m.amScore << " (" << m.amClass << ")"
                << "  DDG=" << m.ddg << " kcal/mol (" << m.ddgInterp << ")\n";
        }
    }

    // 3. If we have a prior Claude summary for THIS variant, include it
    if (!currentVariant.empty() && !currentPdbId_.empty()) {
        for (int i = (int)mutations.size() - 1; i >= 0; --i) {
            const auto& m = mutations[i];
            if (m.pdbId == currentPdbId_ && m.variantId == currentVariant
                && !m.claudeSummary.empty())
            {
                ss << "\nPrevious AI summary for " << currentVariant << ":\n"
                    << "  " << m.claudeSummary << "\n";
                break;
            }
        }
    }

    ss << "=== End Memory Context ===\n";
    return ss.str();
}

// ─── Attach summary ───────────────────────────────────────────────────────────

void DataLensMemory::attachClaudeSummary(
    const std::string& variantId,
    const std::string& summary
) {
    std::lock_guard<std::mutex> lock(memoryMutex);

    // Attach to most recent matching record
    for (int i = (int)mutations.size() - 1; i >= 0; --i) {
        if (mutations[i].variantId == variantId) {
            mutations[i].claudeSummary = summary;
            save();
            return;
        }
    }
    // No matching record — nothing to attach to
}

// ─── Trim ─────────────────────────────────────────────────────────────────────

void DataLensMemory::trimToLimits() {
    // Keep only the most recent N entries (called under lock)
    if ((int)chats.size() > MAX_CHAT_HISTORY) {
        chats.erase(chats.begin(),
            chats.begin() + (chats.size() - MAX_CHAT_HISTORY));
    }
    if ((int)mutations.size() > MAX_MUTATION_HISTORY) {
        mutations.erase(mutations.begin(),
            mutations.begin() + (mutations.size() - MAX_MUTATION_HISTORY));
    }
}