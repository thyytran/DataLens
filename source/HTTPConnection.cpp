#include "HTTPConnection.h"

bool HTTPConnection::curlDidGlobalInit = false;

/*
 *  Callback function for libcurl to write received data to a string.
 *  This function is called by libcurl when it receives data from the server.
 *  It appends the received data to the provided string.
 *  @param contents A pointer to the received data.
 *  @param size The size of each data element.
 *  @param nmemb The number of data elements.
 *  @param userData A pointer to the string where the data should be written.
 *  @return The total number of bytes written.
 */
size_t HTTPConnection::writeCallback(char *contents, size_t size, size_t nmemb, void *userData) {
	std::string *pageContents = (std::string*)userData;
	pageContents->append(contents);
	return size * nmemb;
}

/*
 *  Initializes a libcurl easy handle.
 *  This function creates a new libcurl easy handle and checks if the initialization was successful.
 *  @return True if the initialization was successful, false otherwise.
 */
bool HTTPConnection::initCurl() {
	curl = curl_easy_init();
	if (!curl) {
		std::cerr << "Err > Failed to initialize libcurl.\n\n";
		return 0;
	}
	return 1;
}

/*
 *  Constructor for the HTTPConnection class.
 *  This constructor initializes the libcurl global state if it hasn't been initialized yet.
 */
HTTPConnection::HTTPConnection() {
	if (!curlDidGlobalInit) {
		curl_global_init(CURL_GLOBAL_ALL);
		curlDidGlobalInit = true;
	}
}

/*
 *  Destructor
 */
HTTPConnection::~HTTPConnection() {
	// curl_easy_cleanup is called in get/post methods
}

/*
 *  Performs an HTTP GET request to the specified URL.
 *  This function sends an HTTP GET request to the specified URL and stores the response in the response string.
 *  @param url The URL to send the GET request to.
 *  @return True if the request was successful, false otherwise.
 */
bool HTTPConnection::get(std::string url) {
	response.clear();

	if (!initCurl()) {
		return 0;
	}

	curl_easy_setopt(curl, CURLOPT_URL, url.c_str());
	curl_easy_setopt(curl, CURLOPT_HTTP_VERSION, (long)CURL_HTTP_VERSION_1_0);

	//Configure SSL
	curl_easy_setopt(curl, CURLOPT_SSL_VERIFYPEER, 1L);
	curl_easy_setopt(curl, CURLOPT_SSL_VERIFYHOST, 2L);
	curl_easy_setopt(curl, CURLOPT_CAINFO, "cacert.pem");

	//Write to response string
	curl_easy_setopt(curl, CURLOPT_WRITEFUNCTION, writeCallback);
	curl_easy_setopt(curl, CURLOPT_WRITEDATA, &response);


	CURLcode result = curl_easy_perform(curl);
	if (result != CURLE_OK) {
		std::cerr << "Err > Failed to make request: " << curl_easy_strerror(result) << "\n\n";
		curl_easy_cleanup(curl);
		return 0;
	}

	curl_easy_cleanup(curl);
	return 1;

}


/*
 *  Performs an HTTP POST request with JSON data.
 
bool HTTPConnection::post(std::string url, std::string jsonData) {
    response.clear();

    if (!initCurl()) {
        return false;
    }

    // Set URL
    curl_easy_setopt(curl, CURLOPT_URL, url.c_str());
    curl_easy_setopt(curl, CURLOPT_HTTP_VERSION, (long)CURL_HTTP_VERSION_1_1);

    // Set POST
    curl_easy_setopt(curl, CURLOPT_POST, 1L);
    curl_easy_setopt(curl, CURLOPT_POSTFIELDS, jsonData.c_str());
    curl_easy_setopt(curl, CURLOPT_POSTFIELDSIZE, jsonData.length());

    // Set headers
    struct curl_slist *headers = NULL;
    headers = curl_slist_append(headers, "Content-Type: application/json");
    headers = curl_slist_append(headers, "Accept: application/json");
    curl_easy_setopt(curl, CURLOPT_HTTPHEADER, headers);

    // Configure SSL
    curl_easy_setopt(curl, CURLOPT_SSL_VERIFYPEER, 1L);
    curl_easy_setopt(curl, CURLOPT_SSL_VERIFYHOST, 2L);
    curl_easy_setopt(curl, CURLOPT_CAINFO, "cacert.pem");

    // Write callback
    curl_easy_setopt(curl, CURLOPT_WRITEFUNCTION, writeCallback);
    curl_easy_setopt(curl, CURLOPT_WRITEDATA, &response);

    // Perform request
    CURLcode result = curl_easy_perform(curl);

    // Cleanup
    curl_slist_free_all(headers);
    curl_easy_cleanup(curl);

    if (result != CURLE_OK) {
        std::cerr << "Err > POST request failed: " << curl_easy_strerror(result) << "\n\n";
        return false;
    }

    return true;
}
*/

bool HTTPConnection::post(std::string url, const std::string& body, const std::string& contentType) {
    response.clear();
    if (!initCurl()) return false;

    curl_easy_setopt(curl, CURLOPT_URL, url.c_str());
    curl_easy_setopt(curl, CURLOPT_POST, 1L);
    curl_easy_setopt(curl, CURLOPT_POSTFIELDS, body.c_str());
    curl_easy_setopt(curl, CURLOPT_POSTFIELDSIZE, (long)body.size());

    struct curl_slist* headers = nullptr;
    std::string ct = "Content-Type: " + contentType;
    headers = curl_slist_append(headers, ct.c_str());
    curl_easy_setopt(curl, CURLOPT_HTTPHEADER, headers);

    curl_easy_setopt(curl, CURLOPT_SSL_VERIFYPEER, 1L);
    curl_easy_setopt(curl, CURLOPT_SSL_VERIFYHOST, 2L);
    curl_easy_setopt(curl, CURLOPT_CAINFO, "cacert.pem");
    curl_easy_setopt(curl, CURLOPT_WRITEFUNCTION, writeCallback);
    curl_easy_setopt(curl, CURLOPT_WRITEDATA, &response);

    CURLcode result = curl_easy_perform(curl);
    curl_slist_free_all(headers);
    curl_easy_cleanup(curl);

    if (result != CURLE_OK) {
        std::cerr << "Err > POST failed: " << curl_easy_strerror(result) << "\n\n";
        return false;
    }
    return true;
}