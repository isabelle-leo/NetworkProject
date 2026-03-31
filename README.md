# NetworkProject

Initiated by Sylvia: summer work with Ioannis, Janne, Isabelle

Data can be accessed from: https://ftp.pride.ebi.ac.uk/pride/data/archive/2022/03/PXD023662/target_psmtable.txt

Script used to make the (first draft/exploratory) object: run__deepmeltome_proteogenomics_IL20211206.R

Script used to probe deeper into the (first draft/exploratory) object: exploration_IL20211206.Rmd

### Overview

This project investigates whether proteogenomics-derived novel peptides exhibit systematically different thermal stability compared to conventional peptides, using mass spectrometry-based Thermal Proteome Profiling (TPP) data from pediatric acute lymphoblastic leukaemia (ALL) cell lines (83 samples across 10 TMT sets, sourced from PRIDE PXD023662).
Novel peptides were identified through a proteogenomics pipeline (Husen et al.) and filtered by BLAST specificity and AI-based validation. We fitted melting curves using the NPARC framework, extracted thermal stability metrics (RSS, AUMC, Tm), and systematically tested whether protein biological features could explain the observed differences between novel and conventional peptides.
