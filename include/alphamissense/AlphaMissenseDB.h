#pragma once

#include <string>
#include <unordered_map>
#include <optional>
#include <fstream>
#include <iostream>
#include "Parser.h"  // Your existing Parser for split()

class AlphaMissenseDB {
    private:
        struct Score {
            float pathogenicity = 0.0f;
            std::string classification;
        };

        std::unordered_map<std::string, Score> scores;
        bool isLoaded = false;

    public:
        bool loadFromTSV(const std::string& filepath) {
            std::ifstream file(filepath);

            if (!file.is_open()) {
                std::cerr << "Err > Failed to open AlphaMissense database: "
                    << filepath << "\n\n";
                std::cerr << "      Make sure the file exists and path is correct.\n\n";
                return false;
            }

            std::cout << "Loading AlphaMissense database...\n";

            std::string line;

            // Read header and print it to see the format
            if (std::getline(file, line)) {
                std::cout << "Header: " << line << "\n\n";
            }

            int lineCount = 0;
            int errorCount = 0;

            while (std::getline(file, line)) {
                // Skip empty lines
                if (line.empty()) continue;

                std::vector<std::string> fields = Parser::split(line, '\t');

                // Debug: Print first few lines to see the format
                if (lineCount < 3) {
                    std::cout << "Line " << lineCount << " has " << fields.size() << " fields:\n";
                    for (size_t i = 0; i < fields.size(); ++i) {
                        std::cout << "  [" << i << "] = '" << fields[i] << "'\n";
                    }
                    std::cout << "\n";
                }

                if (fields.size() >= 4) {
                    try {
                        std::string uniprotID = fields[0];
                        std::string variant = fields[1];

                        // Validate and convert score
                        if (fields[2].empty()) {
                            errorCount++;
                            continue;
                        }

                        auto pathogenicity = std::stof(fields[2]);
                        std::string classification = fields[3];

                        std::string key = uniprotID + "_" + variant;
                        scores[key] = { pathogenicity, classification };

                        lineCount++;

                        if (lineCount % 1000000 == 0) {
                            std::cout << "  Loaded " << lineCount / 1000000 << "M entries...\n";
                        }
                    }
                    catch (const std::invalid_argument&) {
                        errorCount++;
                        if (errorCount <= 5) {
                            std::cerr << "Err > Invalid number format in line " << lineCount
                                << ": '" << fields[2] << "'\n";
                        }
                    }
                    catch (const std::out_of_range&) {
                        errorCount++;
                        if (errorCount <= 5) {
                            std::cerr << "Err > Number out of range in line " << lineCount
                                << ": '" << fields[2] << "'\n";
                        }
                    }
                }
            }

            file.close();
            isLoaded = true;
            std::cout << "Loaded " << lineCount
                << " AlphaMissense predictions\n\n";
            return true;
        }

        std::optional<Score> getScore(
            const std::string& uniprotID,
            char refAA,
            int position,
            char altAA
        ) {
            std::string variant = std::string(1, refAA) +
                std::to_string(position) +
                std::string(1, altAA);
            std::string key = uniprotID + "_" + variant;

            auto it = scores.find(key);
            if (it != scores.end()) {
                return it->second;
            }
            return std::nullopt;
        }
        bool isReady() const { return isLoaded; }

    };
