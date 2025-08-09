# Investigating Properties of Scale-Free Networks of Memristors

**Author:** *Mehrdad Hasanpour*  
**Institution:** Institute for Advanced Studies in Basic Sciences (IASBS)  
**Supervisor:** Dr. Ehsan Nedaaee Oskooee  
**Year:** 2016  

---

## 📜 Overview
In 1971, Leon Chua predicted the existence of the **memristor** — the fourth fundamental circuit element. Decades later, researchers developed both **current-controlled memristors (CCM)** and **voltage-controlled memristors (VCM)**.  

This project simulates **scale-free networks of memristors** using a C++ implementation of the **Ariapour–Nedaaee algorithm**. It investigates:
- **Scale-free constant**
- **Clustering coefficient**
- **Average shortest path length**

The networks are generated in various sizes and topologies, and multiple random seeds are used to ensure statistical accuracy. Results provide insight into efficiency, stability, and transmission cost in memristor-based networks.

---

## 🧪 Research Context
- **19 different network sizes** tested
- **7 different random seeds** for each size
- **Voltage variation** studies for CCM and VCM models
- **Key findings**:
  - VCM networks maintain a stable scale-free constant without voltage drop.
  - CCM networks can keep this constant stable with optimized input voltage.
  - All tested networks show **low average path length** → Efficient communication.

---

## 📂 Repository Structure
├── CRT_src/ # Circuit-related simulation components
├── LinSolver/ # Linear equation solvers for network calculations
├── MNA/ # Modified Nodal Analysis implementations
├── NetGen/ # Scale-free network generator (Ariapour–Nedaaee algorithm)
├── NeuronModels/ # Memristor neuron-like models
├── NonLinSolver/ # Non-linear equation solvers
├── Plot/ # Plot generation scripts (e.g., gnuplot, matplotlib)
├── Random/ # Random number generation utilities
├── Simulation/ # Main simulation pipeline
├── Test/ # Unit and integration test codes
├── YounesModel/ # Additional memristor model variant
├── bin/ # Compiled executables
├── include/ # C++ header files
├── lib/ # Libraries or compiled code
├── model/ # XML and configuration model files
├── src/ # General source code
├── doc/ # Documentation and notes
├── Makefile # Build configuration
├── run # Example run script
├── *.xml # Network structure inputs
├── memristor.plt # Gnuplot script for visualizing results



---

## ⚙️ Requirements
- **C++11** or higher
- GNU Make or CMake
- Standard Template Library (STL)
- **Gnuplot** (for visualization)
- (Optional) Python 3 with Matplotlib

---

## 🚀 Build Instructions
From the root directory:
make
or, if using `g++` directly:
g++ -std=c++11 -O2 src/*.cpp -o memristor_network


---

## ▶️ Running the Simulation
After building:
./bin/memristor_network

Optional arguments:

./bin/memristor_network [network_size] [random_seed] [input_voltage]

Example:
./bin/memristor_network 500 42 1.8


---

## 📈 Outputs
- **Text files** in `results/` containing network properties
- **Gnuplot graphs** from `memristor.plt`
- **Log files** with random seed and configuration details for reproducibility

---

## 🔑 Keywords
Scale-free network, memristor, VCM, CCM, complex networks, network analysis, Ariapour–Nedaaee algorithm.

---

## 📚 References
1. Leon O. Chua, “Memristor – The Missing Circuit Element”, IEEE Trans. Circuit Theory, 1971.  
2. Strukov et al., “The Missing Memristor Found”, Nature, 2008.  
3. Wu Hua, “Simulation and Fabrication of Voltage-Controlled Memristors”, 2014.  
4. Ariapour–Nedaaee Network Generation Algorithm.

---

## 📄 License
This project is provided for academic and research purposes.






