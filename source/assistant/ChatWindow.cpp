#include "assistant/ChatWindow.h"
#include <iostream>

extern ImGuiContext* g_mainImGuiContext;

ChatWindow::ChatWindow(int width, int height, GLFWwindow* mainWindow)
    : window(nullptr), sharedContext(mainWindow), shouldClose(false),
    isProcessingMessage(false), imguiContext(nullptr) {

    inputBuffer[0] = '\0';

    // Create window with shared OpenGL context
    // glfwWindowHint(GLFW_CONTEXT_VERSION_MAJOR, 3);
    // glfwWindowHint(GLFW_CONTEXT_VERSION_MINOR, 3);
    // glfwWindowHint(GLFW_OPENGL_PROFILE, GLFW_OPENGL_CORE_PROFILE);

    // Create window WITHOUT shared OpenGL context
    glfwWindowHint(GLFW_CONTEXT_VERSION_MAJOR, 3);
    glfwWindowHint(GLFW_CONTEXT_VERSION_MINOR, 3);
    glfwWindowHint(GLFW_OPENGL_PROFILE, GLFW_OPENGL_CORE_PROFILE);
    glfwWindowHint(GLFW_VISIBLE, GLFW_TRUE);

    window = glfwCreateWindow(width, height, "AI Chat Assistant", nullptr, nullptr);
    if (!window) {
        std::cerr << "Failed to create chat window" << std::endl;
        return;
    }
    
    glfwShowWindow(window); // just added
    glfwMakeContextCurrent(window);
    glfwSwapInterval(1);

    // Initialize GLEW for THIS context (no static check)
    // Each independent OpenGL context needs its own GLEW init
    GLenum err = glewInit();
    if (err != GLEW_OK) {
        std::cerr << "Chat window GLEW init failed: " << glewGetErrorString(err) << std::endl;
        return;
    }
    std::cout << "Chat window GLEW initialized" << std::endl;

    // Initialize ImGui
    // initImGui();

    // Add welcome message
    ChatMessage welcome;
    welcome.sender = "Assistant";
    welcome.message = "Hello! I'm your AI assistant for protein structure analysis. "
        "Ask me anything about protein structures, mutations, or analysis techniques!";
    welcome.timestamp = std::chrono::duration_cast<std::chrono::milliseconds>(
        std::chrono::system_clock::now().time_since_epoch()
    ).count();
    chatHistory.push_back(welcome);
}

ChatWindow::~ChatWindow() {
    if (window) {
        glfwMakeContextCurrent(window);

        // Restore this window's ImGui context before cleanup
        if (imguiContext) {
            ImGui::SetCurrentContext(imguiContext);
            ImGui_ImplOpenGL3_Shutdown();
            ImGui_ImplGlfw_Shutdown();
            ImGui::DestroyContext(imguiContext);
        }

        glfwDestroyWindow(window);
    }
}

/*
void ChatWindow::initImGui() {
    glfwMakeContextCurrent(window);

    IMGUI_CHECKVERSION();
    imguiContext = ImGui::CreateContext(nullptr);

    if (!imguiContext) {
        std::cerr << "ERROR: Failed to create ImGui context!" << std::endl;
        isInitialized = false;
        return;
    }

    ImGui::SetCurrentContext(imguiContext);
    std::cout << "Chat ImGui context created" << std::endl;

    int width, height;
    glfwGetFramebufferSize(window, &width, &height);

    ImGuiIO& io = ImGui::GetIO();
    io.DisplaySize = ImVec2((float)width, (float)height);
    ImGui::StyleColorsDark();

    // DON'T manually build font atlas - let backends handle it

    // Initialize backends - they will build font atlas on first NewFrame()
    std::cout << "Initializing GLFW backend..." << std::endl;
    if (!ImGui_ImplGlfw_InitForOpenGL(window, true)) {
        std::cerr << "ERROR: GLFW backend init failed!" << std::endl;
        isInitialized = false;
        return;
    }
    std::cout << "GLFW backend OK" << std::endl;

    std::cout << "Initializing OpenGL3 backend..." << std::endl;
    if (!ImGui_ImplOpenGL3_Init("#version 330")) {
        std::cerr << "ERROR: OpenGL3 backend init failed!" << std::endl;
        ImGui_ImplGlfw_Shutdown();
        isInitialized = false;
        return;
    }
    std::cout << "OpenGL3 backend OK" << std::endl;

    isInitialized = true;
    std::cout << "Chat window ImGui fully initialized" << std::endl;
}
*/

void ChatWindow::initImGui() {
    std::cout << ">>> Initializing chat ImGui <<<" << std::endl;

    glfwMakeContextCurrent(window);

    IMGUI_CHECKVERSION();
    imguiContext = ImGui::CreateContext(nullptr);

    if (!imguiContext) {
        std::cerr << "ERROR: Failed to create context!" << std::endl;
        isInitialized = false;
        return;
    }

    ImGui::SetCurrentContext(imguiContext);
    std::cout << "Context created: " << imguiContext << std::endl;

    int width, height;
    glfwGetFramebufferSize(window, &width, &height);

    ImGuiIO& io = ImGui::GetIO();
    io.BackendPlatformUserData = nullptr;
    io.BackendRendererUserData = nullptr;
    io.BackendPlatformName = nullptr;
    io.BackendRendererName = nullptr;
    io.DisplaySize = ImVec2((float)width, (float)height);
    io.IniFilename = nullptr;

    ImGui::StyleColorsDark();

    std::cout << "Initializing GLFW backend..." << std::endl;
    if (!ImGui_ImplGlfw_InitForOpenGL(window, false)) {
        std::cerr << "GLFW backend failed!" << std::endl;
        isInitialized = false;
        return;
    }

    std::cout << "Initializing OpenGL backend..." << std::endl;
    if (!ImGui_ImplOpenGL3_Init("#version 330")) {
        std::cerr << "OpenGL backend failed!" << std::endl;
        isInitialized = false;
        return;
    }

    std::cout << "Chat ImGui initialized successfully!" << std::endl;
    isInitialized = true;
}

void ChatWindow::cleanupImGui() {
    if (imguiContext) {
        glfwMakeContextCurrent(window);
        ImGui::SetCurrentContext(imguiContext);

        ImGui_ImplOpenGL3_Shutdown();
        ImGui_ImplGlfw_Shutdown();
        ImGui::DestroyContext(imguiContext);
        imguiContext = nullptr;

        std::cout << "Chat window ImGui cleaned up" << std::endl;
    }
}

void ChatWindow::renderUI() {
    // DON'T set context here - it's already set in run()
    // ImGui::SetCurrentContext(imguiContext);  // REMOVE THIS LINE

    int display_w, display_h;
    glfwGetFramebufferSize(window, &display_w, &display_h);

    ImGui::SetNextWindowPos(ImVec2(0, 0));
    ImGui::SetNextWindowSize(ImVec2((float)display_w, (float)display_h));

    ImGui::Begin("Chat", nullptr,
        ImGuiWindowFlags_NoTitleBar |
        ImGuiWindowFlags_NoResize |
        ImGuiWindowFlags_NoMove |
        ImGuiWindowFlags_NoCollapse);

    // Header with larger text
    ImGui::PushStyleColor(ImGuiCol_Text, ImVec4(0.4f, 0.8f, 1.0f, 1.0f));
    ImGui::SetWindowFontScale(1.2f);
    ImGui::TextWrapped("AI Chat Assistant - Protein Analysis");
    ImGui::SetWindowFontScale(1.1f);
    ImGui::PopStyleColor();
    ImGui::Separator();
    ImGui::Spacing();
    ImGui::Spacing();

    // Chat history area
    float inputAreaHeight = 140.0f;
    ImGui::BeginChild("ChatHistory", ImVec2(0, -inputAreaHeight), true);

    {
        std::lock_guard<std::mutex> lock(chatMutex);

        for (const auto& msg : chatHistory) {
            ImVec4 color;
            if (msg.sender == "You") {
                color = ImVec4(0.4f, 0.8f, 1.0f, 1.0f); // Light blue
            }
            else if (msg.sender == "Assistant") {
                color = ImVec4(0.4f, 1.0f, 0.4f, 1.0f); // Light green
            }
            else {
                color = ImVec4(0.8f, 0.8f, 0.4f, 1.0f); // Yellow
            }

            ImGui::PushStyleColor(ImGuiCol_Text, color);
            ImGui::PushTextWrapPos(ImGui::GetContentRegionAvail().x - 15);

            ImGui::Text("%s:", msg.sender.c_str());
            ImGui::SameLine();
            ImGui::TextWrapped("%s", msg.message.c_str());

            ImGui::PopTextWrapPos();
            ImGui::PopStyleColor();
            ImGui::Spacing();
            ImGui::Spacing();
        }

        if (isProcessingMessage) {
            ImGui::PushStyleColor(ImGuiCol_Text, ImVec4(1.0f, 1.0f, 0.4f, 1.0f));
            ImGui::TextWrapped("Assistant is thinking...");
            ImGui::PopStyleColor();
        }

        if (!processingError.empty()) {
            ImGui::PushStyleColor(ImGuiCol_Text, ImVec4(1.0f, 0.4f, 0.4f, 1.0f));
            ImGui::TextWrapped("%s", processingError.c_str());
            ImGui::PopStyleColor();
        }
    }

    if (ImGui::GetScrollY() >= ImGui::GetScrollMaxY())
        ImGui::SetScrollHereY(1.0f);

    ImGui::EndChild();

    ImGui::Separator();
    ImGui::Spacing();
    ImGui::Spacing();

    ImGui::Text("Your message:");
    ImGui::PushItemWidth(-90);

    bool enterPressed = ImGui::InputTextMultiline("##input", inputBuffer,
        IM_ARRAYSIZE(inputBuffer),
        ImVec2(-1, 60),
        ImGuiInputTextFlags_None);
    ImGui::PopItemWidth();

    ImGui::SameLine();
    ImGui::BeginGroup();
    bool sendClicked = ImGui::Button("Send", ImVec2(80, 30));
    if (ImGui::Button("Clear", ImVec2(80, 30))) {
        std::lock_guard<std::mutex> lock(chatMutex);
        chatHistory.clear();
        ChatMessage welcome;
        welcome.sender = "System";
        welcome.message = "Chat history cleared.";
        welcome.timestamp = std::chrono::duration_cast<std::chrono::milliseconds>(
            std::chrono::system_clock::now().time_since_epoch()
        ).count();
        chatHistory.push_back(welcome);
    }
    ImGui::EndGroup();

    if (ImGui::IsKeyPressed(ImGuiKey_Enter) && ImGui::IsKeyDown(ImGuiKey_ModCtrl)) {
        sendClicked = true;
    }

    if (sendClicked && strlen(inputBuffer) > 0) {
        std::string message(inputBuffer);
        sendMessage(message);
        inputBuffer[0] = '\0';
    }

    ImGui::Spacing();
    ImGui::TextDisabled("Tip: Press Ctrl+Enter to send. Ask about protein structures, mutations, or analysis.");

    ImGui::End();
}

void ChatWindow::sendMessage(const std::string& message) {
    std::lock_guard<std::mutex> lock(chatMutex);

    if (isProcessingMessage) {
        return;
    }

    ChatMessage userMsg;
    userMsg.sender = "You";
    userMsg.message = message;
    userMsg.timestamp = std::chrono::duration_cast<std::chrono::milliseconds>(
        std::chrono::system_clock::now().time_since_epoch()
    ).count();
    chatHistory.push_back(userMsg);

    pendingUserMessage = message;
    isProcessingMessage = true;
    processingError.clear();

    std::thread apiThread(sendChatMessageAsync, this, pendingUserMessage);
    apiThread.detach();
}

void ChatWindow::sendChatMessageAsync(ChatWindow* chatWindow, const std::string& userMessage) {
    try {
        std::string contextPrompt = "You are a helpful assistant for protein structure analysis. ";
        contextPrompt += "You help researchers understand protein structures, mutations, and their effects. ";
        contextPrompt += "Previous conversation:\n";

        {
            std::lock_guard<std::mutex> lock(chatWindow->chatMutex);

            size_t historySize = chatWindow->chatHistory.size();
            size_t contextCount = historySize < 5 ? historySize : 5;
            size_t startIndex = historySize - contextCount;

            for (size_t i = startIndex; i < historySize; i++) {
                if (chatWindow->chatHistory[i].sender != "System") {
                    contextPrompt += chatWindow->chatHistory[i].sender + ": " +
                        chatWindow->chatHistory[i].message + "\n";
                }
            }
        }

        contextPrompt += "\nPlease provide a helpful and informative response.\nAssistant:";

        nlohmann::json json_obj;
        json_obj["model"] = "gpt-3.5-turbo-instruct";
        json_obj["prompt"] = contextPrompt;
        json_obj["max_tokens"] = 500;
        json_obj["temperature"] = 0.7;

        auto completion = openai::completion().create(json_obj);
        std::string response = completion["choices"][0]["text"];

        size_t start = response.find_first_not_of(" \n\r\t");
        size_t end = response.find_last_not_of(" \n\r\t");
        if (start != std::string::npos && end != std::string::npos) {
            response = response.substr(start, end - start + 1);
        }

        {
            std::lock_guard<std::mutex> lock(chatWindow->chatMutex);
            ChatMessage assistantMsg;
            assistantMsg.sender = "Assistant";
            assistantMsg.message = response;
            assistantMsg.timestamp = std::chrono::duration_cast<std::chrono::milliseconds>(
                std::chrono::system_clock::now().time_since_epoch()
            ).count();
            chatWindow->chatHistory.push_back(assistantMsg);

            chatWindow->isProcessingMessage = false;
        }
    }
    catch (const std::exception& e) {
        std::lock_guard<std::mutex> lock(chatWindow->chatMutex);
        chatWindow->processingError = std::string("Error: ") + e.what();
        chatWindow->isProcessingMessage = false;
    }
}

void ChatWindow::run() {
    std::cout << "=== Chat window run() ENTRY ===" << std::endl;
    std::cout << "  window ptr: " << window << std::endl;
    std::cout << "  isInitialized: " << isInitialized << std::endl;

    if (!window) {
        std::cerr << "ERROR: No window!" << std::endl;
        return;
    }

    std::cout << "About to call initImGui()..." << std::endl;

    // Add flush to ensure it prints
    std::cout.flush();

    initImGui();

    std::cout << "AFTER initImGui() call" << std::endl;
    std::cout << "  isInitialized: " << isInitialized << std::endl;

    while (!glfwWindowShouldClose(window) && !shouldClose) {
        glfwMakeContextCurrent(window);
        ImGui::SetCurrentContext(imguiContext);
        glfwPollEvents();

        int display_w, display_h;
        glfwGetFramebufferSize(window, &display_w, &display_h);

        if (display_w <= 0 || display_h <= 0) {
            continue;
        }

        ImGui_ImplOpenGL3_NewFrame();
        ImGui_ImplGlfw_NewFrame();

        ImGuiIO& io = ImGui::GetIO();
        io.DisplaySize = ImVec2((float)display_w, (float)display_h);

        ImGui::NewFrame();
        renderUI();
        ImGui::Render();

        glViewport(0, 0, display_w, display_h);
        glClearColor(0.1f, 0.1f, 0.12f, 1.0f);
        glClear(GL_COLOR_BUFFER_BIT);

        ImDrawData* draw_data = ImGui::GetDrawData();
        if (draw_data != nullptr) {
            ImGui_ImplOpenGL3_RenderDrawData(draw_data);
        }

        glfwSwapBuffers(window);
    }

    shouldClose = true;
}

void ChatWindow::addSystemMessage(const std::string& message) {
    std::lock_guard<std::mutex> lock(chatMutex);

    ChatMessage sysMsg;
    sysMsg.sender = "System";
    sysMsg.message = message;
    sysMsg.timestamp = std::chrono::duration_cast<std::chrono::milliseconds>(
        std::chrono::system_clock::now().time_since_epoch()
    ).count();
    chatHistory.push_back(sysMsg);
}