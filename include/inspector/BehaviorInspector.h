#pragma once

#ifndef OPENGL_MODEL_VIEWER_IMGUI_INSPECTOR_H
#define OPENGL_MODEL_VIEWER_IMGUI_INSPECTOR_H

#include "../include/graphics/Window.h"
#include "imgui/imgui.h"
#include "../include/Camera.h"
#include "DrawableModel.h"

/*
 * Hold properties and settings related to visual representation and behavior of models
 */
struct ModelBehaviorInspector {

    float position[3] = { 0, 0, 0 };
    float rotation[3] = { 0, 0, 0 };
    float scale[3] = { 1, 1, 1 };

    float background_color[3] = { 0.5f, 0.5f, 0.5f };

    // First init
    bool wireframe = false; // Indicating whether to render objects in wireframe mode
    bool hide_crosshair = true; // Controlling the visibility of a crosshair in viewport
    bool rotatable = false; // Whether object can be rotated

    int current_model = 0;
    int current_shader = 0;

    void render(Window windowObj, Camera& camera, std::vector<DrawableModel*>& models);
};


#endif //OPENGL_MODEL_VIEWER_IMGUI_INSPECTOR_H
