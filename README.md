# Glycine Tautomerism Repository

This repository collects the input files, trained models, and analysis utilities used in the study of intramolecular and water-mediated tautomerism of solvated glycine. It is meant to guide readers from the original **ab initio** data production through DeePKS/DeePMD model training and the enhanced-sampling workflows discussed in the publication linked below.

## Tutorials

- [From DFT to DeePKS to DeePMD \| DeePKS Basics](https://nb.bohrium.dp.tech/detail/8742877753)
- [From DFT to DeePKS to DeePMD \| DeePKS Tutorials](https://nb.bohrium.dp.tech/detail/7144731675)
- [OPES (on-the-fly probability enhanced sampling)](https://bohrium.dp.tech/notebooks/9874998164)
- [Voronoi CVs for enhanced sampling autoionization and tautomerism](https://bohrium.dp.tech/notebooks/83327491785)
- *Note:* The above notebook links offer Chinese-to-English translations.

## Dataset, Model, and MD trajectory Availability

- The dataset used to train the DeePMD model is hosted on [AIS Square](https://www.aissquare.com/datasets/detail?pageType=datasets&name=M06-2X_C2H5O2N_H2O&id=238) and mirrored on [Zenodo](https://zenodo.org/records/14309264).
- The compressed DeePMD model is distributed through [AIS Square](https://www.aissquare.com/models/detail?pageType=models&name=M06-2X_C2H5O2N_H2O&id=241) and [Zenodo](https://zenodo.org/records/14309264).
- glycine_10_12ns.lammpstrj.xz on [Zenodo](https://zenodo.org/records/20265544) contains a 2‑ns segment of the glycine MD trajectory (10–12 ns) with element indices 1=H, 2=O, 3=N, 4=C, dump output every 10 fs.

## Software Packages Used

### PLUMED v2.8.1 Patch

- Requires enabling the [OPES module](https://www.plumed.org/doc-v2.8/user-doc/html/_o_p_e_s.html).
- To activate the additional Voronoi collective variables (CVs), copy the three `.cpp` files provided in this repository to `/your_plumed_path/plumed/src/colvar` and recompile PLUMED.
- The Voronoi CV implementations `VORONOID2.cpp`, `VORONOIS1.cpp`, and `VORONOIC0.cpp` map to the CVs `s_d`, `s_p`, and `s_a` described in the [paper](https://pubs.acs.org/doi/10.1021/acs.jcim.4c00273) and its [supporting information](https://pubs.acs.org/doi/suppl/10.1021/acs.jcim.4c00273/suppl_file/ci4c00273_si_001.pdf). Additional [Voronoi CVs](https://github.com/Zhang-pchao/OilWaterInterface/tree/main/Ion_Diffusion_Coefficient) can be used to compute diffusion coefficients for H₃O⁺ or OH⁻ ions.
- More Voronoi CV examples for water autoionization are available [here](https://github.com/Zhang-pchao/OilWaterInterface/tree/main).

### deepmd-kit v2.1.5

- **Re-compile PLUMED:** Follow the [plumed-feedstock](https://github.com/Zhang-pchao/plumed-feedstock/tree/devel) recipe to overlay the default PLUMED version when building LAMMPS.
- **No Re-compile (quick test):** Use the [`LOAD`](https://www.plumed.org/doc-v2.8/user-doc/html/_l_o_a_d.html) command at runtime if recompilation is not feasible.

### deepks-kit v0.1

### ABACUS v3.0.5

### CP2K v9.1

- Comes with instructions for incorporating PLUMED.

## Repository Structure

- `Analysis_Scripts/`
  - `analysis/`: Post-processing workflows grouped by property (collective variables, free-energy surfaces, hydrogen bonding, geometric descriptors, movies, miscellaneous utilities, radial distribution functions, and dipole analysis via Wannier centers).
  - `dpdata/`: Helper tools for manipulating datasets in `dpdata` format.
  - `kinetics/`: Scripts to extract kinetic information from enhanced-sampling simulations.
  - `reference/`: Reference data and configurations used to validate analysis pipelines.
- `DFT_Calculation/`
  - `cp2k_input/`: CP2K input decks (`geo.xyz`, functional-specific directories, and submission scripts) for generating the **ab initio** training data.
- `DeePKS_Iteration/`
  - `glycine/` and `water_ion/`: Iterative DeePKS training setups, each containing `systems/` inputs, `projector/` definitions, and `iter/` directories for successive refinements.
- `DeePMD_Training/`
  - `partial_dataset/`: Representative DeePMD training subsets for glycine and water-ion environments.
  - `run.json`: Example training configuration.
  - `lcurve.out`: Learning-curve history from a DeePMD run.
  - `frozen_model.pb`: Exported DeePMD model ready for deployment.
- `Enhanced_Sampling/`
  - `config_sample/`: Initial CP2K and PLUMED inputs (`cp2k.inp`, `plumed.dat`, `init.xyz`) for generating configurations.
  - `final_MD/`: LAMMPS and PLUMED inputs for the production simulation (`in.lammps`, `input.plumed`) and the optimized atomic structure.
  - `opes_flooding/`: Assets to perform OPES flooding calculations, including submission scripts, notebooks (e.g., `kstest.ipynb`), transition gathering scripts, and time-rescaling utilities.
- `Voronoi_collective_variables/`
  - Custom PLUMED source files (`VoronoiC0.cpp`, `VoronoiD2.cpp`, `VoronoiS1.cpp`) implementing the Voronoi CVs referenced above.

## Reference Paper

Intramolecular and water mediated tautomerism of solvated glycine. [J. Chem. Inf. Model.](https://pubs.acs.org/doi/10.1021/acs.jcim.4c00273) [arXiv](https://arxiv.org/abs/2311.05917)

```bibtex
@article{Zhang_JChemInfModel_2024_v64_p3599,
  title        = {{Intramolecular and Water Mediated Tautomerism of Solvated Glycine}},
  author       = {Pengchao Zhang and Axel Tosello Gardini and Xuefei Xu and Michele Parrinello},
  year         = 2024,
  journal      = {J. Chem. Inf. Model.},
  volume       = 64,
  number       = 9,
  pages        = {3599--3604},
  doi          = {10.1021/acs.jcim.4c00273},
}
```
