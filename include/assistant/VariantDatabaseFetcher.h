#pragma once

/*
 *  VariantDatabaseFetcher
 *  ──────────────────────
 *  Populates the four curated-evidence fields of MutationContext from public REST APIs:
 *    - ClinVar  (NCBI E-utilities)
 *    - UniProt  (REST API — functional site annotations)
 *    - gnomAD   (GraphQL API)
 *    - COSMIC   (REST API v3.4 — requires an account token)
 *
 *  All calls are synchronous.  Call fetchAll() inside a background thread
 *  (e.g. before invoking ClaudeHandler::query()) so the UI never blocks.
 *
 *  Minimal usage:
 *      VariantDatabaseFetcher fetcher;
 *      fetcher.setCosmicToken("Bearer <base64_email:password>");
 *      fetcher.fetchAll(context);   // fills context in-place
 *      handler.query(context, userMsg, onComplete);
 */

#include "assistant/ClaudeHandler.h"  // for MutationContext
#include <string>

class VariantDatabaseFetcher {
public:
    VariantDatabaseFetcher();
    ~VariantDatabaseFetcher() = default;

    /*
     *  COSMIC requires HTTP Basic auth encoded as base64("email:password").
     *  Pass the full "Bearer <token>" header value here.
     *  If not set, COSMIC queries are skipped.
     */
    void setCosmicToken(const std::string& bearerToken);

    /*
     *  Timeout in seconds for each individual HTTP request (default: 8 s).
     *  Keeps the UI responsive if a database is temporarily slow.
     */
    void setTimeoutSeconds(long seconds);

    /*
     *  Run all four fetches in sequence and populate ctx in-place.
     *  Fields already set in ctx are not overwritten.
     *  Returns true if at least one database returned data.
     */
    bool fetchAll(MutationContext& ctx);

    /* Individual fetchers — call selectively if needed */
    bool fetchClinVar(MutationContext& ctx);
    bool fetchUniProt(MutationContext& ctx);
    bool fetchGnomAD(MutationContext& ctx);
    bool fetchCOSMIC(MutationContext& ctx);
    bool fetchUniProtVariant(MutationContext& ctx);

private:
    std::string cosmicToken;
    long        timeoutSeconds = 8;

    /* Low-level HTTP helpers */
    std::string httpGet(const std::string& url,
        const std::string& extraHeader = "");
    std::string httpPost(const std::string& url,
        const std::string& body,
        const std::string& contentType = "application/json");

    static size_t curlWriteCallback(char* ptr, size_t size,
        size_t nmemb, void* userdata);

    /*
     *  Resolves a gene symbol for ClinVar / gnomAD / COSMIC queries.
     *  1. Returns ctx.proteinName immediately if already set (no HTTP).
     *  2. Otherwise calls UniProt REST with ctx.uniprotId, caches the result
     *     back into ctx.proteinName so subsequent fetchers reuse it for free.
     *  Returns empty string if neither source is available.
     */
    std::string resolveGeneName(MutationContext& ctx);

    /* Amino-acid helpers */
    static std::string oneToThree(const std::string& one); // "R" -> "Arg"
};