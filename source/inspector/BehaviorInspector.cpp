/*
#include "../include/inspector/BehaviorInspector.h"
#include "../include/graphics/Window.h"
#include "../include/Camera.h"
#include "imgui/imgui.h"
#include "../include/inspector/DrawableModel.h"
#include <vector>

using namespace std;

void ModelBehaviorInspector::render(Window windowObj, Camera& camera,
    std::vector<DrawableModel*>& models) {

    ImGuiIO& io = ImGui::GetIO();
    ImVec2 window_size = io.DisplaySize;
    ImVec2 window_pos = ImVec2(window_size.x - 300, 0.0f);

    ImGui::SetNextWindowPos(ImVec2(window_pos.x, 120));
    auto flags = ImGuiWindowFlags_AlwaysAutoResize | ImGuiWindowFlags_NoMove;

    ImGui::Begin("Model Render", nullptr, flags);

    ImGui::Text("FPS: %f", ImGui::GetIO().Framerate);
    ImGui::Text("Frame Time: %f ms", ImGui::GetIO().DeltaTime);
    ImGui::Separator();

    ImGui::Text("Press ~ to toggle camera mode");
    ImGui::Text("(Orbit / First Person)");
    ImGui::Separator();

    ImGui::Checkbox("Hide Crosshair", &hide_crosshair);

    ImGui::SliderFloat("FOV", &camera.Zoom, 5.0f, 150.0f);

    auto& model = *models[current_model];

    if (ImGui::Button("Reset Camera"))
    {
        camera.Reset(model.avg_pos, -90, -10);
    }

    bool info = false;
    bool instructions = false;
    bool statistics = false;
    bool fileIo = false;

    if (ImGui::BeginMainMenuBar())
    {
        if (ImGui::MenuItem("Quit"))
        {
            glfwSetWindowShouldClose(windowObj.get(), true);
        }

        if (ImGui::MenuItem("File"))
            fileIo = true;

        if (ImGui::MenuItem("Model Info"))
            info = true;

        if (ImGui::MenuItem("Statistics"))
            statistics = true;

        if (ImGui::MenuItem("Instructions"))
            instructions = true;


        // hack to get the menu to the right
        ImGui::SameLine(ImGui::GetWindowWidth() - 126);
        ImGui::Text("Datalens");
        ImGui::EndMainMenuBar();
    }

    if (info)
        ImGui::OpenPopup("info_popup");

    if (instructions)
        ImGui::OpenPopup("instruction_popup");

    if (statistics)
        ImGui::OpenPopup("statistics_popup");

    if (fileIo)
        ImGui::OpenPopup("file_io_popup");

    if (ImGui::BeginPopup("info_popup"))
    {
        ImGui::Text("Meshes: %d", model.mesh_count);
        ImGui::Text("Vertex Count: %d", model.vertex_count);
        ImGui::Text("Material Count: %d", model.material_count);
        ImGui::EndPopup();
    }

    if (ImGui::BeginPopup("instruction_popup"))
    {
        ImGui::Text("Orbit Mode");
        ImGui::BulletText("Click and drag to rotate");
        ImGui::BulletText("Right click and drag to pan");
        ImGui::BulletText("Scroll to move in/out");

        ImGui::Spacing();

        ImGui::Text("First Person Mode");
        ImGui::BulletText("WASD QE to move");
        ImGui::BulletText("Use mouse to aim");
        ImGui::BulletText("Scroll to zoom");
        ImGui::EndPopup();
    }

    if (ImGui::BeginPopup("statistics_popup"))
    {
        ImGui::BulletText("Statistic data 1");
        ImGui::BulletText("Statistic data 2");
        ImGui::BulletText("Statistic data 3");
        ImGui::EndPopup();
    }

    if (ImGui::BeginPopup("file_io_popup"))
    {
        ImGui::Button("Open new file");
        ImGui::Spacing();

        ImGui::Button("Recent files");
        ImGui::Spacing();

        ImGui::Button("Export image");
        ImGui::EndPopup();
    }

    ImGui::ColorEdit3("Background Color", background_color);

    if (ImGui::CollapsingHeader("Model Transform", ImGuiTreeNodeFlags_DefaultOpen))
    {
        ImGui::DragFloat3("Position", &*position, 0.05f);
        ImGui::DragFloat3("Rotation", &*rotation, 0.1f);
        ImGui::DragFloat3("Scale", &*scale, 0.1f);
        ImGui::Checkbox("Animate", &rotatable);
    }

    ImGui::End();

    ImGui::SetNextWindowPos(ImVec2(ImGui::GetIO().DisplaySize.x, 18), 0, ImVec2(1.0f, 0));
    ImGui::Begin("Visual Options", nullptr, flags);
    if (ImGui::Checkbox("Wireframe", &wireframe))
    {
        glPolygonMode(GL_FRONT_AND_BACK, wireframe ? GL_LINE : GL_FILL);
    }
    if (ImGui::Combo("Model", &current_model,
        "1q8i\0"))
    {
        camera.Reset(models[current_model]->avg_pos, -90, -10);
    }
    if (ImGui::Combo("Shader", &current_shader,
        "Flat\0Blinn-Phong\0Toon\0Gradient\0Screen Door\0"))
    {

    }
    ImGui::End();
}
*/