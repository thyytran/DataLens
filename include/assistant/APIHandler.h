#pragma once

#include <iostream>
#include <string>
#include <curl/curl.h>
#include <nlohmann/json.hpp>

using json = nlohmann::json;


size_t write_callback(void* contents, size_t size, size_t nmemb, std::string* output);
std::string get_openai_response(const std::string& prompt, const std::string& api_key);

