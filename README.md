# CTPP
This repository contains the implementation of a **Core Tree Labeling** algorithm designed for graph processing and analysis tasks. The key components include a main program (`CoreTreeLabeling.cpp`) and a header file (`CoreTreeLabelling.h`) that implements the core logic.

## Directory Structure
├── CoreTreeLabeling.cpp – Main program entry point
├── CoreTreeLabelling.h – Core algorithm implementation
├── graph/ – Input graph datasets
├── fb/ – Additional graph input files
├── txt/ – Text-based input/output files
├── result*/ – Output result directories
├── run*.sh – Shell scripts for executing the program
├── compare_files.py – Python script for result comparison
├── 1.sh, 2.sh – Example script files
└── readme.txt – Original plain text readme

## 📌 File Descriptions

### 🔹 CoreTreeLabeling.cpp

This file contains the `main()` function and serves as the entry point of the Core Tree Labeling program. Its key responsibilities include:

- Parsing input arguments
- Reading graph data from file
- Calling the labeling functions from `CoreTreeLabelling.h`
- Managing the execution flow
- Saving output results to the specified directory

This file orchestrates the entire labeling process and integrates all core components.

---

### 🔹 CoreTreeLabelling.h

This header file contains the implementation of all key functions used in the core tree labeling algorithm. It is designed as a self-contained module and includes:

- **Graph Preprocessing**:
  - Node degree computation
  - Node sorting and filtering
  - Identification of core vertices
- **Core Tree Construction**:
  - Building tree structures over core nodes
  - Maintaining parent-child relationships
- **Labeling Functions**:
  - Assigning structural or numerical labels to nodes
  - Supporting recursive or BFS-based labeling strategies
- **Utility Functions**:
  - Logging and debugging tools
  - File parsing and result formatting
  - Distance or reachability helper functions

This file is fully templated and intended to be directly included in the main program without separate compilation.

---

## ⚙️ How to Compile and Run

You can compile and run the program using `g++`:

```bash
g++ -std=c++11 -O2 CoreTreeLabeling.cpp -o core_labeling
./core_labeling [input_graph_file] [output_directory]
