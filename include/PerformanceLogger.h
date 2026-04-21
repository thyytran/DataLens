#pragma once
// Minimal crash-proof PerformanceLogger

#include <fstream>
#include <chrono>
#include <string>
#include <vector>
#include <mutex>
#include <iostream>

class PerformanceLogger {
public:
    enum class MetricType {
        API_ALPHAMISSENSE_QUERY,
        API_FOLDX_ANALYSIS,
        API_LANDSCAPE_LOAD,
        LOAD_PDB_FETCH,
        LOAD_MUTANT_MODEL,
        RENDER_FRAME,
        UI_MUTATION_FOCUS,
        UI_VARIANT_SELECTOR
    };

private:
    struct Sample {
        MetricType type;
        double ms;
        std::string ctx;
        int atomCount;
    };

    std::vector<Sample> data;
    std::recursive_mutex mtx;  // Changed to recursive_mutex

    std::chrono::high_resolution_clock::time_point timerStart[8];
    bool timerActive[8] = { false };
    std::string timerCtx[8];
    int timerAtomCount[8] = { 0 };

    // Throttle: only sample RENDER_FRAME every 5 seconds
    std::chrono::steady_clock::time_point lastRenderSample = std::chrono::steady_clock::now();
    static constexpr double RENDER_SAMPLE_INTERVAL_SEC = 5.0;

    const char* name(MetricType t) {
        switch (t) {
        case MetricType::API_ALPHAMISSENSE_QUERY: return "API_AlphaMissense";
        case MetricType::API_FOLDX_ANALYSIS:      return "API_FoldX";
        case MetricType::API_LANDSCAPE_LOAD:      return "API_Landscape";
        case MetricType::LOAD_PDB_FETCH:          return "Load_PDB";
        case MetricType::LOAD_MUTANT_MODEL:       return "Load_Mutant";
        case MetricType::RENDER_FRAME:            return "Render_Frame";
        case MetricType::UI_MUTATION_FOCUS:       return "UI_MutFocus";
        case MetricType::UI_VARIANT_SELECTOR:     return "UI_Variant";
        default: return "Unknown";
        }
    }

public:
    void beginTimer(MetricType t, const std::string& ctx = "", int atomCount = 0, int = 0) {
        int i = (int)t;
        if (i >= 0 && i < 8) {
            timerStart[i] = std::chrono::high_resolution_clock::now();
            timerCtx[i] = ctx;
            timerAtomCount[i] = atomCount;
            timerActive[i] = true;
        }
    }

    void endTimer(MetricType t) {
        int i = (int)t;
        if (i >= 0 && i < 8 && timerActive[i]) {
            auto end = std::chrono::high_resolution_clock::now();
            double ms = std::chrono::duration<double, std::milli>(end - timerStart[i]).count();
            timerActive[i] = false;

            // Throttle RENDER_FRAME to every 5 seconds
            if (t == MetricType::RENDER_FRAME) {
                auto now = std::chrono::steady_clock::now();
                double elapsed = std::chrono::duration<double>(now - lastRenderSample).count();
                if (elapsed < RENDER_SAMPLE_INTERVAL_SEC) {
                    return;  // Skip this sample
                }
                lastRenderSample = now;
            }

            std::lock_guard<std::recursive_mutex> lock(mtx);
            data.push_back({ t, ms, timerCtx[i], timerAtomCount[i] });
        }
    }

    void recordSample(MetricType t, double ms, const std::string& ctx = "", int atomCount = 0, int = 0) {
        std::lock_guard<std::recursive_mutex> lock(mtx);
        data.push_back({ t, ms, ctx, atomCount });
    }

    size_t sampleCount() {
        std::lock_guard<std::recursive_mutex> lock(mtx);
        return data.size();
    }

    void exportAll(const std::string& prefix = "datalens_perf") {
        std::cout << "[PerfLog] Starting export...\n";

        std::vector<Sample> copy;
        {
            std::lock_guard<std::recursive_mutex> lock(mtx);
            copy = data;  // Copy under lock
        }

        std::cout << "[PerfLog] Copied " << copy.size() << " samples\n";

        if (copy.empty()) {
            std::cout << "[PerfLog] No data to export\n";
            return;
        }

        // Count per type
        int counts[8] = { 0 };
        double sums[8] = { 0 };
        double mins[8] = { 999999,999999,999999,999999,999999,999999,999999,999999 };
        double maxs[8] = { 0 };
        int lastAtomCount = 0;

        for (const auto& s : copy) {
            int i = (int)s.type;
            if (i >= 0 && i < 8) {
                counts[i]++;
                sums[i] += s.ms;
                if (s.ms < mins[i]) mins[i] = s.ms;
                if (s.ms > maxs[i]) maxs[i] = s.ms;
                if (s.type == MetricType::RENDER_FRAME && s.atomCount > 0) {
                    lastAtomCount = s.atomCount;
                }
            }
        }

        // Write summary to console
        std::cout << "\n=== PERFORMANCE SUMMARY ===\n";
        if (lastAtomCount > 0) {
            std::cout << "Atoms: " << lastAtomCount << "\n";
        }
        for (int i = 0; i < 8; i++) {
            if (counts[i] > 0) {
                double avg = sums[i] / counts[i];
                std::cout << name((MetricType)i) << ": "
                    << "n=" << counts[i]
                    << " avg=" << avg << "ms"
                        << " min=" << mins[i] << "ms"
                        << " max=" << maxs[i] << "ms";
                    if (i == (int)MetricType::RENDER_FRAME && avg > 0) {
                        std::cout << " (" << (1000.0 / avg) << " FPS)";
                    }
                    std::cout << "\n";
            }
        }
        std::cout << "===========================\n\n";

        // Try to write CSV
        std::ofstream f(prefix + "_summary.csv");
        if (f.is_open()) {
            f << "Metric,Count,Avg_ms,Min_ms,Max_ms,FPS\n";
            for (int i = 0; i < 8; i++) {
                if (counts[i] > 0) {
                    double avg = sums[i] / counts[i];
                    double fps = (i == (int)MetricType::RENDER_FRAME && avg > 0) ? 1000.0 / avg : 0;
                    f << name((MetricType)i) << "," << counts[i] << ","
                        << avg << "," << mins[i] << "," << maxs[i] << "," << fps << "\n";
                }
            }
            f.close();
            std::cout << "[PerfLog] Wrote " << prefix << "_summary.csv\n";
        }
        else {
            std::cout << "[PerfLog] Could not write CSV file\n";
        }
    }

    void clear() {
        std::lock_guard<std::recursive_mutex> lock(mtx);
        data.clear();
    }
};