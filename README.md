# Paper Title

Intramolecular and water mediated tautomerism of solvated glycine 

[JCIM](https://pubs.acs.org/doi/10.1021/acs.jcim.4c00273)

[arXiv](https://arxiv.org/abs/2311.05917)

# Toturials

[From DFT to DeePKS to DeePMD | DeePKS Basics](https://nb.bohrium.dp.tech/detail/8742877753)

[From DFT to DeePKS to DeePMD | DeePKS Toturials](https://nb.bohrium.dp.tech/detail/7144731675)

[OPES (on-the-fly probability enhanced sampling)](https://bohrium.dp.tech/notebooks/9874998164)

[Voronoi CVs for enhanced sampling autoionization and tautomerism](https://bohrium.dp.tech/notebooks/83327491785)

Note: The above web links have Chinese to English translations.

# Dataset and model

The dataset for training the DP model is uploaded to [AIS Square](https://www.aissquare.com/datasets/detail?pageType=datasets&name=M06-2X_C2H5O2N_H2O&id=238) and [Google Drive](https://drive.google.com/drive/folders/1SLuqSO00_kIsGftYd241bccZFbtmm_S4?usp=drive_link).

The compressed DP model is uploaded to [AIS Square](https://www.aissquare.com/models/detail?pageType=models&name=M06-2X_C2H5O2N_H2O&id=241) and [Google Drive](https://drive.google.com/drive/folders/1SLuqSO00_kIsGftYd241bccZFbtmm_S4?usp=drive_link).

# Packages Used

## 1. plumed_v2.8.1_patch
Requiring the installation of the [OPES](https://www.plumed.org/doc-v2.8/user-doc/html/_o_p_e_s.html) module. 
To use additional Voronoi CVs code, put the three .cpp files above into /your_plumed_path/plumed/src/colvar, and then compile plumed.

## 2. deepmd-kit_v2.1.5
Incorporating lammps and plumed

## 3. deepks-kit_v0.1

## 4. abacus_v3.0.5

## 5. cp2k_v9.1
Incorporating plumed
