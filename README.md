# DataLens
**AI-Enhanced Interactive Protein Mutation Analysis System**

DataLens unifies fragmented protein mutation analysis tools into a single interactive interface. Select any residue by clicking directly on the 3D structure, and get AlphaMissense pathogenicity scores, FoldX ΔΔG stability predictions, and natural language interpretation in under 3 seconds — compared to 15–30 minutes using traditional workflows.

---

## Features

- **3D protein visualization** — real-time instanced rendering of atoms and backbone connectors via C++/OpenGL
- **Ray-cast residue selection** — click any atom in the viewport to trigger mutation analysis
- **Manual variant input** — select chain, residue, and substitution amino acid from dropdowns
- **AlphaMissense integration** — pathogenicity scores for all 19 possible substitutions at a selected position
- **FoldX ΔΔG prediction** — async stability prediction pipeline with interpretation
- **Structure comparison** — per-residue CA RMSD bar chart and filterable table between wildtype and mutant
- **Natural language interpretation** — Claude API summarizes mutation impact in plain English

---

## Architecture

```
Frontend  (C++)         Backend  (Python)          Data
─────────────────       ──────────────────         ─────────────────────
OpenGL renderer    ←→   FastAPI server        ←→   PostgreSQL
ImGui panels       ←→   FoldX pipeline             AlphaMissense 216M+
RayCaster                SIFTS mapping              PDB / UniProt
```

---

## Requirements

### Frontend
- Visual Studio 2019 or 2022 (MSVC)
- GLFW 3.4
- GLEW 2.1
- curl 7.71.1 (win64-mingw)
- nlohmann/json
- ImGui

All dependencies are bundled in `lib/`. The following DLLs must be copied to your build output directory (`x64/Debug/` or `x64/Release/`):

```
glew32.dll
glfw3.dll
libcurl-x64.dll
libcrypto-1_1-x64.dll
libssl-1_1-x64.dll
```

### Backend
- Python 3.9+
- PostgreSQL with AlphaMissense database loaded
- FoldX executable (place in `backend/foldx/`)

```bash
cd backend
pip install -r requirements.txt
uvicorn main:app --reload
```

---

## Configuration

Create `config.json` in the project root (excluded from git):

```json
{
    "openai_api_key": "your-key-here",
    "backend_url": "http://localhost:8000"
}
```

Set the Visual Studio working directory to `$(SolutionDir)` so the config file is found at runtime:
```
Project Properties → Debugging → Working Directory → $(SolutionDir)
```

---

## Usage

```
fetch pdb <pdb_id>       Load a structure from RCSB (e.g. fetch pdb 7LMK)
showmutations            Load AlphaMissense data for all residues
rotate <x|y|z>           Rotate the model
setdefault               Reset colors and radii to default
leave                    Clear the current structure
exit                     Quit
```

Click any atom in the 3D viewport to select a residue and trigger mutation analysis. Use the **Variant Selector** panel to manually specify a substitution.

---

## Proteins Used for Validation

| Protein | UniProt | PDB |
|---------|---------|-----|
| TP53 | P04637 | 7LMK |
| BRCA1 | P38398 | 1JM7 |
| Histone H4 | P62805 | — |

---

## Project Structure

```
DataLens/
├── include/
│   ├── alphamissense/    AlphaMissenseDB, variant summaries
│   ├── assistant/        ChatWindow, APIHandler
│   ├── bio/              PDB parsing, atom/residue/chain/molecule data
│   ├── graphics/         renderer, shaders, ray caster, tooltip
│   ├── imgui/            ImGui headers and backend bindings
│   ├── inspector/        behavior inspector, drawable mesh/model
│   ├── mapper/           PDB ↔ UniProt coordinate mapping
│   ├── math/             Vec, Mat, MathUtils
│   ├── mutation/         analysis, FoldX pipeline, structure comparison
│   └── openai/           OpenAI API client
├── source/               implementations (.cpp, .inl)
├── lib/                  third-party dependencies (GLFW, GLEW, curl, GLM)
├── config.json           credentials (not committed)
└── README.md
```

---

## Acknowledgements

Developed as a Master's thesis project at Northeastern University (Khoury College of Computer Sciences). Domain validation by protein design scientists at Moderna.
