#pragma once
#include <string>

struct AppConfig {
    std::string aiApiKey;
    std::string backendUrl;  // move BACKEND_URL here too
};

AppConfig loadConfig(const std::string& path);