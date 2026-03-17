# fClaHyb_analysis
This repository contains the active analysis workflows for **the F1 Thai hybrid catfish genome** (*Clarias macrocephalus × Clarias gariepinus*), **hatchery resequencing population genomics of North African catfish** (*Clarias gariepinus*), and downstream comparative analyses including **phylogenomics**, **ancestral karyotype reconstruction**, and **polymorphic inversion discovery and interpretation**. It serves as the main technical workflow repository for reproducible downstream analysis built around the public **fClaHyb** genome resource.

---

## 🧭 Repository structure

### 📁 `00_project_docs`
Project notes, planning, workflow overviews, and documentation.

### 📁 `01_assembly`
Genome assembly-related analyses, assembly QC, and assembly preparation scripts.

### 📁 `02_annotation`
Genome annotation workflows and supporting annotation utilities.

### 📁 `03_pop_gen_structure_relatedness`
Population genomics workflows for:
- site discovery
- BEAGLE genotype likelihood generation
- structure / PCA / relatedness-oriented preprocessing
- thinning / SNP panel preparation

### 📁 `04_pop_gen_stats_diversity_scans`
Population-genetic summary statistics and diversity scan workflows, including window-based analyses.

### 📁 `05_phylogenomics`
Phylogenomic workflows, orthology-related analyses, and multi-species comparative sequence processing.

### 📁 `06_ancestral_karyotype`
Macrosynteny, chromosome evolution, and ancestral karyotype reconstruction workflows.

### 📁 `07_polymorphic_inversions`
Detection, validation, interpretation, and summary of polymorphic inversions.

### 📁 `08_final_figures_and_tables`
Final scripts and helpers used to prepare publication figures, tables, and manuscript-linked outputs.

---

## 🧩 General design philosophy

This repository is organized by **analysis module**, not by software tool alone.

Each major folder should contain:
- the scripts for that module
- a small `README_pipeline.txt` explaining run order
- clearly named steps
- only the active/current scripts in the main folder
- older or replaced scripts moved to an `old/` subfolder when needed

The goal is to keep each module understandable on its own while still fitting into the full project workflow.

---
