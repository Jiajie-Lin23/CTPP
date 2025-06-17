# CTPP
This repository contains the implementation of a **Core Tree Labeling** algorithm designed for graph processing and analysis tasks. The key components include a main program (`CoreTreeLabeling.cpp`) and a header file (`CoreTreeLabelling.h`) that implements the core logic.

## Directory Structure
- ├── `CoreTreeLabelling.cpp` — Main program entry point  
- ├── `CoreTreeLabelling.h` — Core algorithm implementation  
- ├── `graph/` — Input graph datasets  
- ├── `fb/` — Intermediate output file
- ├── `txt/` — Result file
- ├── `result/` — Output result to record the process
- ├── `compare_files.py` — Python script for result comparison  
- ├── `1.sh`, — The script to run the example  

## ⚙️ How to Compile and Run

## 🛠️ 1. Compile

To compile the project:

```bash
g++ -O3 -fopenmp -std=c++11 CoreTreeLabeling.cpp -o run -w
```

This will generate an executable file named `run`.

---

## 🚀 2. Run Instructions

### 2.1 Convert Text Graph to Binary

```bash
./run txt-to-bin graph/ DBLP.txt
```

- `@1` — Path to the input graph directory  
- `@2` — Graph file name (e.g., `DBLP.txt`)  

This step converts a plain-text graph into a binary format for faster processing.

---

### 2.2 Decompose Using Core Tree (BT)

```bash
./run decompose_bt fb/ 7 64
```

- `@1` — Input directory (e.g., preprocessed files or landmarks)  
- `@2` — Allowed tree-width (e.g., 7)  
- `@3` — Number of threads (e.g., 64)

This performs core-tree-based decomposition with specified parameters.

---

### 2.3 Query Distance Between Node Pairs

```bash
./run query-dis fb/ DBLP.txt 7 1
```

- `@1` — Input directory (same as decomposition output)  
- `@2` — Graph file name  
- `@3` — Tree-width allowed  
- `@4` — Number of reference nodes (set to 1 in our experiment)

This runs landmark-based distance queries using the decomposed core tree.

---
