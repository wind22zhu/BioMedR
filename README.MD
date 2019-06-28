# BioMedR

## Introduction

The BioMedR package offers an R package for generating various molecular descriptors for chemicals, proteins, DNAs/RNAs and their interactions. Also, this package realized some similarity searching methods and powerful clustering methods as well as several useful auxiliary tools, which aims at building integrated analysis pipelines from data acquisition, data checking, and descriptor calculation to modeling. See `vignette('BioMedR')` for the comprehensive user guide.

## Installation

To install the BioMedR package in R, simply type

   ```
   install.package('BioMedR')
   ```
or 
   ```
   install.packages('devtools')
   library(devtools)
   install_github('wind22zhu/BioMedR')
   ```
To make the BioMedR package fully functional (especially the Open Babel and rjava related functionalities), we recommend the users some tips: 

[>>Tips 1](http://projects.scbdd.com/si/install_biomedr_tips.pdf)
[>>Tips 2](http://projects.scbdd.com/si/install_rjava.pdf)

Several dependencies of the BioMedR package may require some system-level libraries, check the corresponding manuals of these packages for detailed installation guides.

## Features

Rcpi implemented and integrated the state-of-the-art protein sequence descriptors and molecular descriptors/fingerprints with R. For protein sequences, the Rcpi package could

  * Calculate six protein descriptor groups composed of fourteen types of commonly used structural and physicochemical descriptors that include 9920 descriptors.

  * Calculate six types of generalized scales-based descriptors derived by various dimensionality reduction methods for proteochemometric (PCM) modeling.

  * Parallellized pairwise similarity computation derived by protein sequence alignment and Gene Ontology (GO) semantic similarity measures within a list of proteins.

For small molecules, the Rcpi package could

  * Calculate 307 molecular descriptors (2D/3D), including constitutional, topological, geometrical, and electronic descriptors, etc.

  * Calculate more than ten types of molecular fingerprints, including FP4 keys, E-state fingerprints, MACCS keys, etc., and parallelized chemical similarity search.

  * Calculate three nucleic acid composition features describing the local sequence information by means of kmers (subsequences of DNA sequences).

  * Calculate  six autocorrelation features describing the level of correlation between  two oligonucleotides along a DNA sequence in terms of their specific physicochemical properties.
 
  * Calculate two pseudo nucleotide composition features, which can be used to represent a DNA sequence with a discrete model or vector yet still keep considerable sequence order information, particularly the global or long-range sequence order information, via the physicochemical properties of its constituent oligonucleotides.

  * Parallelized pairwise similarity computation derived by fingerprints and maximum common substructure search within a list of small molecules.

By combining various types of descriptors for drugs, proteins and DNA/RNA in different methods, interaction descriptors representing protein-protein, compound-compound, DNA-DNA,  compound-DNA compound-protein and DNA-protein interactions could be conveniently generated withBioMedR, including:

  * Two types of compound-protein interaction (CPI) descriptors
  * Two types of compound-DNA interaction (CDI) descriptors 
  * Two types of DNA-protein interaction (DPI) descriptors
  * Three types of protein-protein interaction (PPI) descriptors
  * Three types of compound-compound interaction (CCI) descriptors
  * Three types of DNA-DNA interaction (DDI) descriptors

Several useful auxiliary utilities are also shipped with BioMedR:

  * Parallelized molecule and protein sequence retrieval from several online databases, like PubChem, ChEMBL, KEGG, DrugBank, UniProt, RCSB PDB, genBank, etc.

  * Loading molecules stored in SMILES/SDF files and loading protein sequences from FASTA/PDB files

  * Molecular file format conversion

The computed protein sequence descriptors, molecular descriptors/fingerprints, interaction descriptors and pairwise similarities are widely used in various research fields relevant to drug disvery, primarily bioinformatics, chemoinformatics, proteochemometrics and chemogenomics.

## Links

  * Bioconductor Page: http://bioconductor.org/packages/release/bioc/html/BioMedR.html

  * Track Devel: https://github.com/wind22zhu/BioMedR

  * Report Bugs: https://github.com/wind22zhu/BioMedR/issues

## Contact

The BioMedR package is developed by Computational Biology and Drug Design Group, Central South University, China.
  
  * Minfeng Zhu <wind2zhu@163.com> 

  * Dong Jie <biomed@csu.edu.cn>

  * Dongsheng Cao <oriental-cds@163.com>

