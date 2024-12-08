# Paper

Intramolecular and water mediated tautomerism of solvated glycine. [J. Chem. Inf. Model. 2024, 64, 9, 3599–3604.](https://pubs.acs.org/doi/10.1021/acs.jcim.4c00273) [arXiv](https://arxiv.org/abs/2311.05917)

# Toturials

[From DFT to DeePKS to DeePMD | DeePKS Basics](https://nb.bohrium.dp.tech/detail/8742877753)

[From DFT to DeePKS to DeePMD | DeePKS Toturials](https://nb.bohrium.dp.tech/detail/7144731675)

[OPES (on-the-fly probability enhanced sampling)](https://bohrium.dp.tech/notebooks/9874998164)

[Voronoi CVs for enhanced sampling autoionization and tautomerism](https://bohrium.dp.tech/notebooks/83327491785)

Note: The above web links have Chinese to English translations.

# Dataset and model

The dataset for training the DP model is uploaded to [AIS Square](https://www.aissquare.com/datasets/detail?pageType=datasets&name=M06-2X_C2H5O2N_H2O&id=238), [Zenodo](https://zenodo.org/records/14309264) and [Google Drive](https://drive.google.com/drive/folders/1SLuqSO00_kIsGftYd241bccZFbtmm_S4?usp=drive_link).

The compressed DP model is uploaded to [AIS Square](https://www.aissquare.com/models/detail?pageType=models&name=M06-2X_C2H5O2N_H2O&id=241), [Zenodo](https://zenodo.org/records/14309264) and [Google Drive](https://drive.google.com/drive/folders/1SLuqSO00_kIsGftYd241bccZFbtmm_S4?usp=drive_link).

# Packages Used

## 1. plumed_v2.8.1_patch
Requiring the installation of the [OPES](https://www.plumed.org/doc-v2.8/user-doc/html/_o_p_e_s.html) module. 
To use additional Voronoi CVs code, put the three .cpp files above into /your_plumed_path/plumed/src/colvar, and then compile plumed. 
The Voronoi CV VORONOID2.cpp, VORONOIS1.cpp, VORONOIC0.cpp code files are linked to CVs named s_d, s_p, s_a as illustrated in the [paper](https://pubs.acs.org/doi/10.1021/acs.jcim.4c00273) and [SI](https://pubs.acs.org/doi/suppl/10.1021/acs.jcim.4c00273/suppl_file/ci4c00273_si_001.pdf). Other [Voronoi CVs](https://github.com/Zhang-pchao/OilWaterInterface/tree/main/Ion_Diffusion_Coefficient) can be used to calculate the diffusion coefficient for H₃O⁺ or OH⁻ ions.
## 2. deepmd-kit_v2.1.5
Incorporating lammps and plumed, follow [plumed-feedstock](https://github.com/Zhang-pchao/plumed-feedstock/tree/devel) to overlay default plumed version or use [LOAD](https://www.plumed.org/doc-v2.8/user-doc/html/_l_o_a_d.html) command.

## 3. deepks-kit_v0.1

## 4. abacus_v3.0.5

## 5. cp2k_v9.1
Incorporating plumed
