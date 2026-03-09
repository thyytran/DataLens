#pragma once

#include <string>
#include <iostream>

#include <curl/curl.h>

class HTTPConnection {
private:
	// Static member to track whether curl_global_init has been called
	static bool curlDidGlobalInit;

	// CURL object, pointer to the libcurl easy handle
	CURL *curl = nullptr;

	/*
	 *  Callback function used by libcurl to write the data received from the HTTP request
	 *  @param contents  Pointer to the received data
	 *  @param size      Size of each data block
	 *  @param nmemb     Number of data blocks
	 *  @param userData  Pointer to user-provided data (in this case, a string to store the response)
	 *  @return          The total number of bytes written
	 */
	static size_t writeCallback(char *contents, size_t size, size_t nmemb, void *userData);

	/*
	* Initializes the CURL object and sets basic options.
	* @return true if initialization was successful, false otherwise.
	*/
	bool initCurl();

public:
	// Public member to store the HTTP response
	std::string response;

	// Constructor for the HTTPConnection class
	HTTPConnection();
	~HTTPConnection();

	/*
	 * Performs an HTTP GET request to the specified URL.
	 * @param url The URL to send the GET request to.
	 * @return true if the request was successful, false otherwise.
	 */
	bool get(std::string url);

	// post(std::string url, std::string jsonData);  // NEW
	bool post(std::string url, const std::string& body, const std::string& contentType = "application/json");

	// Response access
	std::string getResponse() const { return response; }  // NEW
	void clearResponse() { response.clear(); }  // NEW
};
