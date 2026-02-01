# Automated Active Site Profiler

##  Project Overview
This tool automates the **Structure-Based Drug Design (SBDD)** workflow by identifying the precise amino acid residues that form the binding pocket of a target protein.

Using **Biopython** and **KD-Tree algorithms**, it parses raw PDB/ENT files, locates the ligand (drug), and performs a geometric neighbor search to define the active site micro-environment.

##  Key Features
* **Automated PDB Retrieval:** Fetches structures directly from the RCSB Protein Data Bank.
* **Geometric Neighbor Search:** Uses a KD-Tree (k-dimensional tree) for $O(\log n)$ fast searching of atomic coordinates within a 5.0 Å radius.
* **Residue Profiling:** Outputs a clean list of interacting amino acids (Hydrophobic contacts, Hydrogen bonds) essential for docking grid generation.

##  Tech Stack
* **Language:** Python 3.x
* **Libraries:** `Bio.PDB` (Biopython)
* **Algorithm:** KD-Tree (Nearest Neighbor Search)

##  Example Output (Gleevec bound to Abl Kinase)
```text
BINDING POCKET FOR STI (Gleevec):
Total interacting residues: 29
LEU 248  A
VAL 256  A
ALA 269  A
...