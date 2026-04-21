#include "assistant/VariantDatabaseFetcher.h"

#include <iostream>
#include <sstream>
#include <map>
#include <curl/curl.h>
#include "nlohmann/json.hpp"

using json = nlohmann::json;

// ─── Construction ─────────────────────────────────────────────────────────────

VariantDatabaseFetcher::VariantDatabaseFetcher() {}

void VariantDatabaseFetcher::setCosmicToken(const std::string& bearerToken) {
    cosmicToken = bearerToken;
}

void VariantDatabaseFetcher::setTimeoutSeconds(long seconds) {
    timeoutSeconds = seconds;
}

// ─── Public API ───────────────────────────────────────────────────────────────

bool VariantDatabaseFetcher::fetchAll(MutationContext& ctx) {
    bool any = false;
    any |= fetchClinVar(ctx);
    // any |= fetchUniProt(ctx);           // site annotations
    any |= fetchUniProtVariant(ctx);    // variant predictions (PolyPhen, SIFT, etc.)
    any |= fetchGnomAD(ctx);
    any |= fetchCOSMIC(ctx);
    return any;
}

// ─── Gene name resolver ───────────────────────────────────────────────────────
//
//  Shared by fetchClinVar, fetchGnomAD, and fetchCOSMIC.
//  If ctx.proteinName is already set, returns it immediately — no HTTP call.
//  Otherwise queries UniProt with ctx.uniprotId and caches the result back
//  into ctx.proteinName so the next fetcher gets it for free.

std::string VariantDatabaseFetcher::resolveGeneName(MutationContext& ctx) {

    // Helper: does a string look like a UniProt accession?
    // Format: [A-NR-Z][0-9][A-Z][A-Z0-9]{2}[0-9]  (6 chars)
    // or      [OPQ][0-9][A-Z0-9]{3}[0-9]           (6 chars)
    // Simple check: 6 chars, first is a letter, rest are alphanumeric
    auto looksLikeAccession = [](const std::string& s) -> bool {
        if (s.size() != 6) return false;
        if (!std::isupper((unsigned char)s[0])) return false;
        for (size_t i = 1; i < s.size(); ++i)
            if (!std::isalnum((unsigned char)s[i])) return false;
        return true;
        };

    // If proteinName is set and is NOT a UniProt accession, use it directly
    if (ctx.proteinName && !ctx.proteinName->empty()) {
        if (!looksLikeAccession(ctx.proteinName.value()))
            return ctx.proteinName.value();

        std::cerr << "VariantDatabaseFetcher::resolveGeneName: "
            << "proteinName looks like a UniProt accession (\""
            << ctx.proteinName.value() << "\") — will look up gene name\n";
        // Fall through to UniProt lookup
    }

    // Resolve UniProt accession to use for lookup — check multiple fields
    std::string accession;

    if (ctx.uniprotId && !ctx.uniprotId->empty()) {
        accession = ctx.uniprotId.value();
    }
    else if (ctx.chain && !ctx.chain->empty()
        && looksLikeAccession(*ctx.chain)) {
        // Chain field is being used to carry a UniProt accession
        accession = ctx.chain.value();
        std::cerr << "VariantDatabaseFetcher::resolveGeneName: "
            << "reading UniProt accession from ctx.chain (\""
            << accession << "\")\n";
        // Cache it into the right field
        ctx.uniprotId = accession;
    }
    else if (ctx.proteinName && looksLikeAccession(ctx.proteinName.value())) {
        accession = ctx.proteinName.value();
    }

    if (accession.empty()) return "";

    std::string url =
        "https://rest.uniprot.org/uniprotkb/" + accession +
        "?fields=gene_names&format=json";

    std::string raw = httpGet(url);
    if (raw.empty()) return "";

    try {
        json entry = json::parse(raw);
        if (entry.contains("genes")
            && entry["genes"].is_array()
            && !entry["genes"].empty())
        {
            auto& first = entry["genes"][0];
            if (first.contains("geneName")
                && first["geneName"].contains("value"))
            {
                std::string gene = first["geneName"]["value"].get<std::string>();
                ctx.proteinName = gene;   // cache gene symbol for subsequent fetchers
                ctx.uniprotId = accession;  // ensure uniprotId is set correctly
                return gene;
            }
        }
    }
    catch (...) {}

    return "";
}

// ─── ClinVar ──────────────────────────────────────────────────────────────────
//
//  Strategy: NCBI E-utilities esearch → esummary (JSON).
//
//  Query strategy (tried in order until one returns results):
//    1. "BRCA1 L1701P"          — gene + one-letter shorthand (space-separated,
//                                  matches ClinVar web search behavior exactly)
//    2. "BRCA1[gene] L1701P"    — with gene field tag for precision
//    3. "p.Leu1701Pro BRCA1"    — three-letter HGVS + gene
//    4. "p.Leu1701Pro"          — broadest fallback (no gene filter)
//
//  URL encoding: spaces → %20, brackets → %5B %5D, dots are safe.
//  NEVER use raw '+' as space — NCBI's CGI interprets '+' as literal plus,
//  not as space, in the term= parameter outside of HTML form POST context.
//
//  Docs: https://www.ncbi.nlm.nih.gov/clinvar/docs/maintenance_use/
//        https://www.ncbi.nlm.nih.gov/clinvar/docs/help/

// Minimal percent-encoder for ClinVar term= values.
// Encodes only characters that break NCBI's query parser.
// Dots are intentionally NOT encoded — p.Leu1701Pro must arrive as-is.
// Brackets ARE encoded — [gene] becomes %5Bgene%5D in the URL,
// but NCBI decodes them correctly on its end.
static std::string urlEncodeTerm(const std::string& s) {
    std::ostringstream out;
    for (unsigned char c : s) {
        switch (c) {
        case ' ': out << "%20"; break;
        case '[': out << "%5B"; break;
        case ']': out << "%5D"; break;
        case '(': out << "%28"; break;
        case ')': out << "%29"; break;
        case '+': out << "%2B"; break;
            // '.' is safe in URLs — do NOT encode, required for p.Leu1701Pro
        default:  out << c;    break;
        }
    }
    return out.str();
}

bool VariantDatabaseFetcher::fetchClinVar(MutationContext& ctx) {
    if (!ctx.wildTypeAA || !ctx.mutantAA || !ctx.residueNum) return false;
    if (ctx.clinvarClass) return true;  // already populated

    int         pos = ctx.residueNum.value();
    std::string wt = ctx.wildTypeAA.value();
    std::string mut = ctx.mutantAA.value();

    // Resolve gene symbol from proteinName or uniprotId (UniProt lookup if needed)
    std::string gene = resolveGeneName(ctx);

    // One-letter shorthand: "L1701P"
    std::string shorthand = wt + std::to_string(pos) + mut;
    // Three-letter HGVS:   "p.Leu1701Pro"
    std::string hgvsP = "p." + oneToThree(wt) + std::to_string(pos) + oneToThree(mut);

    // Build candidate queries in priority order.
    // Start with shorthand alone — always runs even if gene resolution fails.
    // Gene-prefixed queries narrow the results when available.
    std::vector<std::string> queries;
    queries.push_back(shorthand);                              // "L1701P" — always first

    if (!gene.empty() && gene.size() > 4) {
        // Only use gene prefix if it looks like a real gene symbol,
        // not a PDB ID (4 chars) or UniProt accession (6 chars, caught elsewhere)
        queries.push_back(gene + "[gene] AND " + hgvsP);      // "BRCA1[gene] AND p.Leu1701Pro"
        queries.push_back(gene + "[gene] AND " + shorthand);  // "BRCA1[gene] AND L1701P"
    }
    queries.push_back(hgvsP);                                 // "p.Leu1701Pro" — broadest fallback

    std::string uid;
    for (const auto& q : queries) {
        // Fetch up to 5 results — we need to validate each one matches our
        // exact protein change, since ClinVar's protein_change field is a
        // comma-separated list of all isoform aliases (L1701P may appear in
        // a record whose canonical title is L1729Pro).
        std::string esearchUrl =
            "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi"
            "?db=clinvar&retmode=json&retmax=5&term=" + urlEncodeTerm(q);

        std::cerr << "ClinVar esearch: " << esearchUrl << "\n";

        std::string raw = httpGet(esearchUrl);
        if (raw.empty()) continue;

        std::vector<std::string> candidates;
        try {
            json esearch = json::parse(raw);
            auto& idlist = esearch["esearchresult"]["idlist"];
            for (auto& id : idlist)
                candidates.push_back(id.get<std::string>());
        }
        catch (...) { continue; }

        if (candidates.empty()) continue;

        // Fetch all candidate summaries in one call (comma-separated IDs)
        std::string allIds;
        for (size_t i = 0; i < candidates.size(); ++i) {
            if (i > 0) allIds += ",";
            allIds += candidates[i];
        }

        std::string sumUrl =
            "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esummary.fcgi"
            "?db=clinvar&retmode=json&id=" + allIds;

        std::string sumRaw = httpGet(sumUrl);
        if (sumRaw.empty()) continue;

        try {
            json sumJson = json::parse(sumRaw);
            auto& resultMap = sumJson["result"];

            for (const auto& id : candidates) {
                if (!resultMap.contains(id)) continue;
                auto& rec = resultMap[id];

                // Validate: protein_change must contain our exact shorthand
                // e.g. "L1701P" must appear as a standalone token in the
                // comma-separated protein_change field, not just as a substring.
                if (rec.contains("protein_change") && rec["protein_change"].is_string()) {
                    std::string pc = rec["protein_change"].get<std::string>();
                    // Search for shorthand as a word boundary token
                    // comma-delimited list: ", L1701P," or starts with "L1701P,"
                    bool matched = false;
                    std::istringstream ss(pc);
                    std::string token;
                    while (std::getline(ss, token, ',')) {
                        // trim spaces
                        size_t s = token.find_first_not_of(' ');
                        size_t e = token.find_last_not_of(' ');
                        if (s != std::string::npos)
                            token = token.substr(s, e - s + 1);
                        if (token == shorthand) { matched = true; break; }
                    }
                    if (!matched) {
                        std::cerr << "ClinVar: skipping uid=" << id
                            << " (protein_change has no exact " << shorthand << ")\n";
                        continue;
                    }
                }

                uid = id;
                std::cerr << "ClinVar: matched uid=" << uid
                    << " query=" << q << "\n";
                break;
            }
        }
        catch (...) { continue; }

        if (!uid.empty()) break;
    }

    if (uid.empty()) {
        std::cerr << "VariantDatabaseFetcher::fetchClinVar: "
            << "no record found for " << (gene.empty() ? "" : gene + " ")
            << shorthand << "\n";
        return false;
    }

    // ── Fetch final esummary for the validated UID ────────────────────────────
    std::string esummaryUrl =
        "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esummary.fcgi"
        "?db=clinvar&retmode=json&id=" + uid;

    std::string esummaryRaw = httpGet(esummaryUrl);
    if (esummaryRaw.empty()) return false;

    try {
        json esummary = json::parse(esummaryRaw);
        auto& result = esummary["result"][uid];

        // ── Clinical classification ───────────────────────────────────────────
        auto extractClass = [&](const json& node) -> bool {
            if (!node.is_object()) return false;
            if (node.contains("description") && node["description"].is_string()) {
                std::string desc = node["description"].get<std::string>();
                if (!desc.empty()) ctx.clinvarClass = desc;
            }
            if (node.contains("review_status") && node["review_status"].is_string()) {
                std::string rs = node["review_status"].get<std::string>();
                if (rs.find("practice guideline") != std::string::npos) ctx.clinvarStars = 4;
                else if (rs.find("expert panel") != std::string::npos) ctx.clinvarStars = 3;
                else if (rs.find("multiple submitters") != std::string::npos) ctx.clinvarStars = 2;
                else if (rs.find("single submitter") != std::string::npos) ctx.clinvarStars = 1;
                else                                                           ctx.clinvarStars = 0;
            }
            // Trait is nested inside germline_classification.trait_set
            // NOT at result["trait_set"] — confirmed from live API response
            if (node.contains("trait_set") && node["trait_set"].is_array()) {
                for (auto& trait : node["trait_set"]) {
                    if (trait.contains("trait_name") && trait["trait_name"].is_string()) {
                        std::string name = trait["trait_name"].get<std::string>();
                        if (name != "not specified" && name != "not provided" && !name.empty()) {
                            ctx.clinvarCondition = name;
                            break;
                        }
                    }
                }
            }
            return ctx.clinvarClass.has_value();
            };

        bool found = false;
        if (result.contains("germline_classification"))
            found = extractClass(result["germline_classification"]);
        if (!found && result.contains("clinical_significance"))
            found = extractClass(result["clinical_significance"]);

        // Accession
        if (result.contains("accession") && result["accession"].is_string())
            ctx.clinvarId = result["accession"].get<std::string>();
    }
    catch (const std::exception& e) {
        std::cerr << "VariantDatabaseFetcher::fetchClinVar parse error: "
            << e.what() << "\n";
        return false;
    }

    return ctx.clinvarClass.has_value();
}

// ─── UniProt Site Annotations ────────────────────────────────────────────────
//
//  Strategy: UniProt REST API features endpoint.
//  Fetches site annotations (active sites, binding sites, domains) for a residue.

bool VariantDatabaseFetcher::fetchUniProt(MutationContext& ctx) {
    if (!ctx.uniprotResidueNum && !ctx.residueNum)
        return false;
    if (ctx.uniprotSiteAnnotation)
        return true;  // already populated

    resolveGeneName(ctx);
    if (!ctx.uniprotId || ctx.uniprotId->empty())
        return false;

    std::string accession = ctx.uniprotId.value();
    int targetPos = ctx.uniprotResidueNum
        ? ctx.uniprotResidueNum.value()
        : ctx.residueNum.value();

    std::string url =
        "https://rest.uniprot.org/uniprotkb/" + accession +
        "?fields=ft_act_site,ft_binding,ft_site,ft_region,ft_domain"
        "&format=json";

    std::string raw = httpGet(url);
    if (raw.empty()) return false;

    std::vector<std::string> hits;

    try {
        json entry = json::parse(raw);
        if (!entry.contains("features")) return false;

        for (auto& feature : entry["features"]) {
            if (!feature.contains("location")) continue;
            auto& loc = feature["location"];

            int start = -1, end = -1;
            if (loc.contains("start") && loc["start"].contains("value"))
                start = loc["start"]["value"].get<int>();
            if (loc.contains("end") && loc["end"].contains("value"))
                end = loc["end"]["value"].get<int>();

            if (start < 0) continue;
            if (end < 0) end = start;

            if (targetPos < start || targetPos > end) continue;

            std::string featureType = feature.value("type", "");
            std::string description = feature.value("description", "");

            std::string annotation = featureType;
            if (!description.empty())
                annotation += " (" + description + ")";
            hits.push_back(annotation);
        }
    }
    catch (const std::exception& e) {
        std::cerr << "fetchUniProt parse error: " << e.what() << "\n";
        return false;
    }

    if (hits.empty()) return false;

    std::ostringstream ss;
    for (size_t i = 0; i < hits.size(); ++i) {
        if (i > 0) ss << "; ";
        ss << hits[i];
    }
    ctx.uniprotSiteAnnotation = ss.str();
    return true;
}

// ─── UniProt Variant ─────────────────────────────────────────────────────────
//
//  Strategy: UniProt Proteins API /variation/{accession} endpoint.
//  Fetches all variants for a UniProt accession, then filters by position
//  and mutant amino acid to find the exact variant.
//
//  Docs: https://www.ebi.ac.uk/proteins/api/doc/#/variation

bool VariantDatabaseFetcher::fetchUniProtVariant(MutationContext& ctx) {
    // Need position and mutation info
    if (!ctx.mutantAA || (!ctx.uniprotResidueNum && !ctx.residueNum))
        return false;

    // Already populated?
    if (ctx.uniprotVariantId)
        return true;

    // Resolve UniProt accession
    resolveGeneName(ctx);
    if (!ctx.uniprotId || ctx.uniprotId->empty()) {
        std::cerr << "[fetchUniProtVariant] no uniprotId\n";
        return false;
    }

    std::string accession = ctx.uniprotId.value();
    int targetPos = ctx.uniprotResidueNum
        ? ctx.uniprotResidueNum.value()
        : ctx.residueNum.value();
    std::string targetMut = ctx.mutantAA.value();

    std::string url = "https://www.ebi.ac.uk/proteins/api/variation/" + accession;

    std::cerr << "[UniProtVariant] Fetching " << accession
        << " pos=" << targetPos << " mut=" << targetMut << "\n";

    std::string raw = httpGet(url);
    if (raw.empty()) return false;

    try {
        json response = json::parse(raw);

        if (!response.contains("features") || !response["features"].is_array())
            return false;

        // Helper for safe string extraction
        auto safeString = [](const json& obj, const std::string& key) -> std::string {
            if (!obj.contains(key)) return "";
            if (obj[key].is_string()) return obj[key].get<std::string>();
            return "";
            };

        for (auto& feature : response["features"]) {
            if (safeString(feature, "type") != "VARIANT") continue;

            // Parse position (API returns as string)
            std::string beginStr = safeString(feature, "begin");
            int pos = 0;
            try { pos = std::stoi(beginStr); }
            catch (...) { continue; }

            if (pos != targetPos) continue;

            // Check mutation matches
            std::string altSeq = safeString(feature, "alternativeSequence");
            if (altSeq.empty() || altSeq[0] != targetMut[0]) continue;

            // ══════════════════════════════════════════════════════════════
            // Found matching variant — extract all data
            // ══════════════════════════════════════════════════════════════

            std::cerr << "[UniProtVariant] MATCH: " << safeString(feature, "wildType")
                << " -> " << altSeq << " at " << pos << "\n";

            // Basic info
            ctx.uniprotVariantId = safeString(feature, "ftId");
            ctx.uniprotConsequence = safeString(feature, "consequenceType");
            ctx.uniprotCodon = safeString(feature, "codon");
            ctx.uniprotSourceType = safeString(feature, "sourceType");

            // Somatic flag
            if (feature.contains("somatic") && feature["somatic"].is_boolean())
                ctx.uniprotSomatic = feature["somatic"].get<bool>();

            // Predictions (PolyPhen, SIFT)
            if (feature.contains("predictions") && feature["predictions"].is_array()) {
                for (auto& pred : feature["predictions"]) {
                    std::string algName = safeString(pred, "predAlgorithmNameType");
                    std::string predVal = safeString(pred, "predictionValType");
                    float score = -1.0f;
                    if (pred.contains("score") && pred["score"].is_number())
                        score = pred["score"].get<float>();

                    if (algName.find("PolyPhen") != std::string::npos ||
                        algName.find("polyphen") != std::string::npos) {
                        ctx.uniprotPolyphenPred = predVal;
                        if (score >= 0) ctx.uniprotPolyphenScore = score;
                    }
                    else if (algName.find("SIFT") != std::string::npos ||
                        algName.find("sift") != std::string::npos) {
                        ctx.uniprotSiftPred = predVal;
                        if (score >= 0) ctx.uniprotSiftScore = score;
                    }
                }
            }

            // Clinical significances
            if (feature.contains("clinicalSignificances") &&
                feature["clinicalSignificances"].is_array() &&
                !feature["clinicalSignificances"].empty()) {
                std::ostringstream ss;
                bool first = true;
                for (auto& clin : feature["clinicalSignificances"]) {
                    std::string type = safeString(clin, "type");
                    if (!type.empty()) {
                        if (!first) ss << ", ";
                        ss << type;
                        first = false;
                    }
                }
                if (!first) ctx.uniprotClinicalSignificance = ss.str();
            }

            // Disease associations
            if (feature.contains("association") && feature["association"].is_array() &&
                !feature["association"].empty()) {
                std::ostringstream ss;
                bool first = true;
                for (auto& assoc : feature["association"]) {
                    std::string name = safeString(assoc, "name");
                    if (!name.empty()) {
                        if (!first) ss << "; ";
                        ss << name;
                        first = false;
                    }
                }
                if (!first) ctx.uniprotDiseaseAssociation = ss.str();
            }

            return true;  // Found and populated
        }

        std::cerr << "[UniProtVariant] No matching variant found at pos "
            << targetPos << " -> " << targetMut << "\n";
    }
    catch (const std::exception& e) {
        std::cerr << "[fetchUniProtVariant] parse error: " << e.what() << "\n";
        return false;
    }

    return false;
}

// ─── gnomAD ───────────────────────────────────────────────────────────────────
//
//  Strategy: gnomAD GraphQL API (v4, GRCh38).
//  We query by gene symbol + protein change using the variantSearch endpoint.
//
//  Docs: https://gnomad.broadinstitute.org/api
//
//  Note: gnomAD v4 GraphQL schema changed field names vs v2/v3.
//  We request exome + genome combined AF.

bool VariantDatabaseFetcher::fetchGnomAD(MutationContext& ctx) {
    if (!ctx.wildTypeAA || !ctx.mutantAA || !ctx.residueNum) return false;
    if (ctx.gnomadAF) return true;  // already populated

    // Resolve gene symbol from proteinName or uniprotId (UniProt lookup if needed)
    std::string gene = resolveGeneName(ctx);
    if (gene.empty()) {
        std::cerr << "VariantDatabaseFetcher::fetchGnomAD: "
            << "no gene symbol (set proteinName or uniprotId)\n";
        return false;
    }

    int         pos = ctx.residueNum.value();
    std::string wt = ctx.wildTypeAA.value();
    std::string mut = ctx.mutantAA.value();

    // Convert single-letter AA codes to three-letter for HGVS protein notation
    std::string wtThree = oneToThree(wt);
    std::string mutThree = oneToThree(mut);

    // gnomAD uses HGVS p. notation: p.Arg248Trp
    std::string hgvsP = "p." + wtThree + std::to_string(pos) + mutThree;

    // GraphQL query — search by gene + protein consequence
    std::string gql = R"({
  "query": "query VariantSearch($gene: String!, $query: String!) { variantSearch(dataset: gnomad_r4, query: $query) { variants { variant_id exome { af ac } genome { af ac } popmax_population } } }",
  "variables": { "gene": ")" + gene + R"(", "query": ")" + hgvsP + R"(" }
})";

    std::string raw = httpPost(
        "https://gnomad.broadinstitute.org/api",
        gql,
        "application/json"
    );

    if (raw.empty()) return false;

    try {
        json resp = json::parse(raw);
        auto& variants = resp["data"]["variantSearch"]["variants"];
        if (!variants.is_array() || variants.empty()) return false;

        // Take the first match
        auto& v = variants[0];

        // Prefer combined AF: exome AF if available, else genome
        float af = 0.0f;
        int   ac = 0;

        if (v.contains("exome") && !v["exome"].is_null()
            && v["exome"].contains("af") && !v["exome"]["af"].is_null()) {
            af = v["exome"]["af"].get<float>();
            ac = v["exome"].value("ac", 0);
        }
        else if (v.contains("genome") && !v["genome"].is_null()
            && v["genome"].contains("af") && !v["genome"]["af"].is_null()) {
            af = v["genome"]["af"].get<float>();
            ac = v["genome"].value("ac", 0);
        }

        ctx.gnomadAF = af;
        ctx.gnomadAC = ac;

        if (v.contains("popmax_population") && !v["popmax_population"].is_null())
            ctx.gnomadPopmax = v["popmax_population"].get<std::string>();
    }
    catch (const std::exception& e) {
        std::cerr << "VariantDatabaseFetcher::fetchGnomAD parse error: "
            << e.what() << "\n";
        return false;
    }

    return ctx.gnomadAF.has_value();
}

// ─── COSMIC ───────────────────────────────────────────────────────────────────
//
//  Strategy: COSMIC REST API v3.4 — /mutations endpoint filtered by gene + AA change.
//  Requires HTTP Basic Auth (base64-encoded "email:password") passed as Bearer token.
//
//  Gene symbol resolution (in priority order):
//    1. ctx.proteinName       — used directly if set (e.g. "TP53")
//    2. ctx.uniprotId         — if proteinName is absent, resolve gene symbol from
//                               UniProt REST API (/uniprotkb/{accession}?fields=gene_names)
//
//  Docs: https://cancer.sanger.ac.uk/cosmic/download/cosmic/v99/documentation

bool VariantDatabaseFetcher::fetchCOSMIC(MutationContext& ctx) {
    if (cosmicToken.empty()) return false;
    if (!ctx.wildTypeAA || !ctx.mutantAA || !ctx.residueNum) return false;
    if (ctx.cosmicCount) return true;

    // Resolve gene symbol from proteinName or uniprotId (UniProt lookup if needed)
    std::string gene = resolveGeneName(ctx);
    if (gene.empty()) {
        std::cerr << "VariantDatabaseFetcher::fetchCOSMIC: "
            << "no gene symbol (set proteinName or uniprotId)\n";
        return false;
    }

    // ── Build HGVS protein notation and query ─────────────────────────────────
    int         pos = ctx.residueNum.value();
    std::string wt = ctx.wildTypeAA.value();
    std::string mut = ctx.mutantAA.value();

    std::string hgvsP = "p." + oneToThree(wt) + std::to_string(pos) + oneToThree(mut);

    std::string url =
        "https://cancer.sanger.ac.uk/api/v3.4/mutations"
        "?gene_name=" + gene +
        "&mutation_aa=" + hgvsP +
        "&format=json&limit=20";  // fetch 20 so cancer-type aggregation is meaningful

    std::string authHeader = "Authorization: " + cosmicToken;
    std::string raw = httpGet(url, authHeader);
    if (raw.empty()) return false;

    try {
        json resp = json::parse(raw);

        if (!resp.contains("results") || !resp["results"].is_array()) return false;
        if (resp["results"].empty()) return false;

        int totalCount = resp.value("count", 0);
        ctx.cosmicCount = totalCount;

        auto& first = resp["results"][0];
        if (first.contains("id"))
            ctx.cosmicId = "COSV" + std::to_string(first["id"].get<int>());

        std::map<std::string, int> cancerCounts;
        for (auto& r : resp["results"]) {
            if (r.contains("primary_site") && r["primary_site"].is_string())
                cancerCounts[r["primary_site"].get<std::string>()]++;
        }

        if (!cancerCounts.empty()) {
            std::ostringstream ss;
            bool firstEntry = true;
            for (auto& [cancer, count] : cancerCounts) {
                if (!firstEntry) ss << ", ";
                ss << cancer << " (" << count << ")";
                firstEntry = false;
            }
            ctx.cosmicCancerTypes = ss.str();
        }
    }
    catch (const std::exception& e) {
        std::cerr << "VariantDatabaseFetcher::fetchCOSMIC parse error: "
            << e.what() << "\n";
        return false;
    }

    return ctx.cosmicCount.has_value();
}

// ─── HTTP helpers ─────────────────────────────────────────────────────────────

size_t VariantDatabaseFetcher::curlWriteCallback(
    char* ptr, size_t size, size_t nmemb, void* userdata
) {
    auto* buf = static_cast<std::string*>(userdata);
    buf->append(ptr, size * nmemb);
    return size * nmemb;
}

std::string VariantDatabaseFetcher::httpGet(
    const std::string& url,
    const std::string& extraHeader
) {
    CURL* curl = curl_easy_init();
    if (!curl) return "";

    std::string response;
    struct curl_slist* headers = nullptr;
    headers = curl_slist_append(headers, "Accept: application/json");
    if (!extraHeader.empty())
        headers = curl_slist_append(headers, extraHeader.c_str());

    curl_easy_setopt(curl, CURLOPT_URL, url.c_str());
    curl_easy_setopt(curl, CURLOPT_HTTPHEADER, headers);
    curl_easy_setopt(curl, CURLOPT_WRITEFUNCTION, curlWriteCallback);
    curl_easy_setopt(curl, CURLOPT_WRITEDATA, &response);
    curl_easy_setopt(curl, CURLOPT_TIMEOUT, timeoutSeconds);
    curl_easy_setopt(curl, CURLOPT_SSL_VERIFYPEER, 1L);
    curl_easy_setopt(curl, CURLOPT_SSL_VERIFYHOST, 2L);
    curl_easy_setopt(curl, CURLOPT_CAINFO, "cacert.pem");
    curl_easy_setopt(curl, CURLOPT_FOLLOWLOCATION, 1L);

    CURLcode res = curl_easy_perform(curl);
    if (res != CURLE_OK) {
        std::cerr << "VariantDatabaseFetcher GET failed (" << url << "): "
            << curl_easy_strerror(res) << "\n";
        response.clear();
    }

    curl_slist_free_all(headers);
    curl_easy_cleanup(curl);
    return response;
}

std::string VariantDatabaseFetcher::httpPost(
    const std::string& url,
    const std::string& body,
    const std::string& contentType
) {
    CURL* curl = curl_easy_init();
    if (!curl) return "";

    std::string response;
    struct curl_slist* headers = nullptr;
    std::string ctHeader = "Content-Type: " + contentType;
    headers = curl_slist_append(headers, ctHeader.c_str());
    headers = curl_slist_append(headers, "Accept: application/json");

    curl_easy_setopt(curl, CURLOPT_URL, url.c_str());
    curl_easy_setopt(curl, CURLOPT_HTTPHEADER, headers);
    curl_easy_setopt(curl, CURLOPT_POSTFIELDS, body.c_str());
    curl_easy_setopt(curl, CURLOPT_POSTFIELDSIZE, (long)body.size());
    curl_easy_setopt(curl, CURLOPT_WRITEFUNCTION, curlWriteCallback);
    curl_easy_setopt(curl, CURLOPT_WRITEDATA, &response);
    curl_easy_setopt(curl, CURLOPT_TIMEOUT, timeoutSeconds);
    curl_easy_setopt(curl, CURLOPT_SSL_VERIFYPEER, 1L);
    curl_easy_setopt(curl, CURLOPT_SSL_VERIFYHOST, 2L);
    curl_easy_setopt(curl, CURLOPT_CAINFO, "cacert.pem");

    CURLcode res = curl_easy_perform(curl);
    if (res != CURLE_OK) {
        std::cerr << "VariantDatabaseFetcher POST failed (" << url << "): "
            << curl_easy_strerror(res) << "\n";
        response.clear();
    }

    curl_slist_free_all(headers);
    curl_easy_cleanup(curl);
    return response;
}

// ─── Amino-acid helpers ───────────────────────────────────────────────────────

std::string VariantDatabaseFetcher::oneToThree(const std::string& one) {
    static const std::map<char, std::string> table = {
        {'A', "Ala"}, {'R', "Arg"}, {'N', "Asn"}, {'D', "Asp"},
        {'C', "Cys"}, {'Q', "Gln"}, {'E', "Glu"}, {'G', "Gly"},
        {'H', "His"}, {'I', "Ile"}, {'L', "Leu"}, {'K', "Lys"},
        {'M', "Met"}, {'F', "Phe"}, {'P', "Pro"}, {'S', "Ser"},
        {'T', "Thr"}, {'W', "Trp"}, {'Y', "Tyr"}, {'V', "Val"},
        {'*', "Ter"}
    };
    if (one.empty()) return "";
    auto it = table.find(one[0]);
    return (it != table.end()) ? it->second : one;
}