#include "Config.h"
#include "nlohmann/json.hpp"
#include <fstream>

AppConfig loadConfig(const std::string& path) {
    std::ifstream f(path);
    auto j = nlohmann::json::parse(f);
    return { j["api_key"], j["backend_url"] };
}