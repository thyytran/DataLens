#pragma once
#include <assistant/APIHandler.h>

size_t write_callback(void* contents, size_t size, size_t nmemb, std::string* output) {
    size_t total_size = size * nmemb;
    output->append((char*)contents, total_size);
    return total_size;
}

std::string get_openai_response(const std::string& prompt, const std::string& api_key) {
    CURL* curl;
    CURLcode res;
    std::string response_string;

    curl_global_init(CURL_GLOBAL_ALL);
    curl = curl_easy_init();
    if (curl) {

        // 1. Prepare the JSON payload
        json payload = {
            {"model", "gpt-3.5-turbo"},
            {"messages", {{ "role", "user"}, {"content", prompt}}}
        };
        std::string payload_str = payload.dump();

        // 2. Set the URL
        curl_easy_setopt(curl, CURLOPT_URL, "https://api.openai.com/v1/chat/completions");

        // 3. Set headers
        struct curl_slist* headers = NULL;
        headers = curl_slist_append(headers, "Content-Type: application/json");
        headers = curl_slist_append(headers, ("Authorization: Bearer " + api_key).c_str());
        curl_easy_setopt(curl, CURLOPT_HTTPHEADER, headers);

        // 4. Set the request method to POST and provide the payload
        curl_easy_setopt(curl, CURLOPT_POST, 1L);
        curl_easy_setopt(curl, CURLOPT_POSTFIELDS, payload_str.c_str());

        // 5. Set the callback function
        curl_easy_setopt(curl, CURLOPT_WRITEFUNCTION, write_callback);
        curl_easy_setopt(curl, CURLOPT_WRITEDATA, &response_string);

        // 6. Perform the request
        res = curl_easy_perform(curl);

        // 7. Handle errors
        if (res != CURLE_OK) {
            std::cerr << "curl_easy_perform() failed: " << curl_easy_strerror(res) << std::endl;
            curl_easy_cleanup(curl);
            curl_global_cleanup();
            return "Error: Could not connect to OpenAI";
        }

        curl_easy_cleanup(curl);

        // Parse the JSON response
        try {
            json response_json = json::parse(response_string);
            return response_json["choices"][0]["message"]["content"].get<std::string>();
        }
        catch (const json::parse_error& e) {
            std::cerr << "JSON parse error: " << e.what() << std::endl;
            curl_global_cleanup();
            return "Error: Could not parse OpenAI response";
        }
    }
    else {
        std::cerr << "curl_easy_init() failed" << std::endl;
        curl_global_cleanup();
        return "Error: Could not initialize cURL";
    }

    curl_global_cleanup();
}