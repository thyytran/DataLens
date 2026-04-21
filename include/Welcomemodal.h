#pragma once
#pragma once

#include "imgui/imgui.h"

/*
 *  WelcomeModal
 *  ─────────────
 *  Startup splash screen shown once when DataLens first opens.
 *  Introduces the tool and presents a research-use disclaimer
 *  the user must explicitly acknowledge before proceeding.
 *
 *  Usage (in your render loop):
 *
 *      static WelcomeModal welcome;
 *      welcome.draw();   // no-ops automatically after user dismisses
 */
class WelcomeModal {
public:
    WelcomeModal() : dismissed(false), openPending(true) {}

    void draw() {
        // Open the modal on the first frame
        if (openPending) {
            ImGui::OpenPopup("##welcome");
            openPending = false;
        }

        if (dismissed) return;

        // Center the modal on the display
        ImVec2 center = ImGui::GetMainViewport()->GetCenter();
        ImGui::SetNextWindowPos(center, ImGuiCond_Always, ImVec2(0.5f, 0.5f));
        ImGui::SetNextWindowSize(ImVec2(560, 0), ImGuiCond_Always); // auto height

        ImGuiWindowFlags flags =
            ImGuiWindowFlags_NoMove |
            ImGuiWindowFlags_NoResize |
            ImGuiWindowFlags_NoCollapse |
            ImGuiWindowFlags_NoScrollbar;

        if (ImGui::BeginPopupModal("##welcome", nullptr, flags)) {

            // ── Title ─────────────────────────────────────────────────────────
            ImGui::PushStyleColor(ImGuiCol_Text, ImVec4(0.40f, 0.78f, 1.00f, 1.00f));
            ImGui::SetWindowFontScale(1.25f);
            ImGui::Text("DataLens");
            ImGui::SetWindowFontScale(1.00f);
            ImGui::PopStyleColor();

            ImGui::SameLine();
            ImGui::SetWindowFontScale(0.85f);
            ImGui::TextDisabled("AI-Enhanced Protein Mutation Analysis");
            ImGui::SetWindowFontScale(1.00f);

            ImGui::Separator();
            ImGui::Spacing();

            // ── About ─────────────────────────────────────────────────────────
            ImGui::TextWrapped(
                "DataLens is an interactive platform for exploring "
                "the structural and functional effects of protein missense mutations. "
                "It combines real-time 3D visualization with pathogenicity predictions "
                "(AlphaMissense), stability estimates (FoldX), curated evidence from "
                "ClinVar, UniProt, and gnomAD, and AI-synthesized interpretations "
                "to help researchers rapidly triage candidate variants."
            );

            ImGui::Spacing();

            // ── Workflow summary ──────────────────────────────────────────────
            ImGui::TextDisabled("What you can do:");
            ImGui::Bullet(); ImGui::TextWrapped("Load any protein from the RCSB Protein Data Bank by PDB ID.");
            ImGui::Bullet(); ImGui::TextWrapped("Select a residue to instantly retrieve AlphaMissense pathogenicity scores and FoldX DDG stability predictions.");
            ImGui::Bullet(); ImGui::TextWrapped("Ask the AI assistant to interpret the combined evidence in plain language.");
            ImGui::Bullet(); ImGui::TextWrapped("Compare wild-type and mutant conformations side by side.");

            ImGui::Spacing();
            ImGui::Separator();
            ImGui::Spacing();

            // ── Disclaimer box ────────────────────────────────────────────────
            ImGui::PushStyleColor(ImGuiCol_ChildBg, ImVec4(0.35f, 0.20f, 0.05f, 0.40f));
            ImGui::BeginChild("##disclaimer_box", ImVec2(0, 100), true);

            ImGui::PushStyleColor(ImGuiCol_Text, ImVec4(1.00f, 0.75f, 0.30f, 1.00f));
            ImGui::Text("  Research Use Only");
            ImGui::PopStyleColor();

            ImGui::Spacing();
            ImGui::PushStyleColor(ImGuiCol_Text, ImVec4(0.90f, 0.90f, 0.90f, 1.00f));
            ImGui::TextWrapped(
                "DataLens is intended for exploratory research and hypothesis generation "
                "only. All predictions are computational estimates. "
                "Do not use this tool to make clinical, diagnostic, or treatment decisions."
            );
            ImGui::PopStyleColor();

            ImGui::EndChild();
            ImGui::PopStyleColor();

            ImGui::Spacing();

            // ── Dismiss button ────────────────────────────────────────────────
            float buttonWidth = 180.0f;
            ImGui::SetCursorPosX((ImGui::GetContentRegionAvail().x - buttonWidth) * 0.5f);

            ImGui::PushStyleColor(ImGuiCol_Button, ImVec4(0.20f, 0.50f, 0.85f, 1.00f));
            ImGui::PushStyleColor(ImGuiCol_ButtonHovered, ImVec4(0.30f, 0.60f, 1.00f, 1.00f));
            ImGui::PushStyleColor(ImGuiCol_ButtonActive, ImVec4(0.15f, 0.40f, 0.75f, 1.00f));

            if (ImGui::Button("I understand \nGet started", ImVec2(buttonWidth, 0))) {
                dismissed = true;
                ImGui::CloseCurrentPopup();
            }

            ImGui::PopStyleColor(3);
            ImGui::Spacing();

            ImGui::EndPopup();
        }
    }

    bool isDismissed() const { return dismissed; }

private:
    bool dismissed;
    bool openPending;
};