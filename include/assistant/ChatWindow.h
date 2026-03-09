#ifndef CHAT_WINDOW_H
#define CHAT_WINDOW_H

#include <string>
#include <deque>
#include <mutex>
#include <thread>
#include <chrono>
#include "GL/glew.h"
#include "GLFW/glfw3.h"
#include "imgui/imgui.h"
#include "imgui/imgui_impl_glfw.h"
#include "imgui/imgui_impl_opengl3.h"
#include "nlohmann/json.hpp"
#include "openai.hpp"

// Chat message structure
struct ChatMessage {
    std::string sender;
    std::string message;
    long long timestamp;
};

class ChatWindow {
private:
    GLFWwindow* window;
    GLFWwindow* sharedContext; // Reference to main window for context sharing
    bool shouldClose;
    ImFontAtlas* sharedFontAtlas;

    // ImGui context for this window
    ImGuiContext* imguiContext;

    // Thread-safe chat data
    std::mutex chatMutex;
    std::deque<ChatMessage> chatHistory;
    std::string pendingUserMessage;
    bool isProcessingMessage;
    std::string processingError;

    char inputBuffer[1024];

    void initImGui();
    void cleanupImGui();
    void renderUI();
    void sendMessage(const std::string& message);
    static void sendChatMessageAsync(ChatWindow* chatWindow, const std::string& userMessage);
    bool isInitialized;  // Add this

public:
    ChatWindow(int width, int height, GLFWwindow* mainWindow);
    ~ChatWindow();

    void run();
    bool isWindowClosed() const { return shouldClose; }
    void addSystemMessage(const std::string& message);
    GLFWwindow* getWindow() { return window; }
};

#endif