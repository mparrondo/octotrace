<div id="top">

<!-- HEADER STYLE: COMPACT -->
<img src="resources/octotrace_logo.png" width="30%" align="left" style="margin-right: 15px">

# OCTOTRACE
<em></em>

<!-- BADGES -->
<img src="https://img.shields.io/github/license/mparrondo/octotrace?style=for-the-badge&logo=opensourceinitiative&logoColor=white&color=FF8552" alt="license">
<img src="https://img.shields.io/github/last-commit/mparrondo/octotrace?style=for-the-badge&logo=git&logoColor=white&color=FF8552" alt="last-commit">
<img src="https://img.shields.io/github/languages/top/mparrondo/octotrace?style=for-the-badge&color=FF8552" alt="repo-top-language">
<img src="https://img.shields.io/github/languages/count/mparrondo/octotrace?style=for-the-badge&color=FF8552" alt="repo-language-count">

<em>Built with the tools and technologies:</em>


<br clear="left"/>

## Table of Contents

1. [Table of Contents](#table-of-contents)
2. [Overview](#overview)
3. [Project Structure](#project-structure)
4. [Getting Started](#getting-started)
    - [Prerequisites](#prerequisites)
    - [Installation](#installation)
    - [Usage](#usage)
5. [Contributing](#contributing)
6. [License](#license)
7. [Acknowledgments](#acknowledgments)

---

## Overview

<details open>
<summary><b>🇬🇧 English</b></summary>

**octoTrace** contains the analysis code supporting the study *"Genome-wide SNP discovery and a reduced diagnostic panel for geographic assignment of common octopus (Octopus vulgaris) in the Bay of Biscay and adjacent fishing regions."* The repository documents the workflows used to filter genome-wide SNP data, identify informative markers, and evaluate reduced SNP panels for geographic assignment across five operational population units: Northern and Southern Iberian Atlantic, Western and Eastern Mediterranean, and Macaronesia. The selected eight-SNP panel offers a balance between assignment performance and marker number. Its intended resolution is regional; the current data do not support assignment specifically to the Asturian fishery.

</details>

<details>
<summary><b>🇪🇸 Español</b></summary>

**octoTrace** contiene el código de análisis empleado en el estudio *"Genome-wide SNP discovery and a reduced diagnostic panel for geographic assignment of common octopus (Octopus vulgaris) in the Bay of Biscay and adjacent fishing regions."*  El repositorio documenta los flujos de trabajo utilizados para filtrar datos genómicos de SNPs, identificar marcadores informativos y evaluar paneles reducidos para la asignación geográfica en cinco unidades poblacionales operativas: Atlántico ibérico norte y sur, Mediterráneo occidental y oriental, y Macaronesia. El panel seleccionado de ocho SNPs ofrece un equilibrio entre la precisión de asignación y el número de marcadores. Su resolución prevista es regional: los datos actuales no permiten asignar ejemplares específicamente a la pesquería asturiana.

</details>

---

## Project Structure

```sh
└── octotrace/
    ├── LICENSE
    ├── README.md
    └── resources
        └── octotrace_logo.png
    └── scripts
        └── 01_basemaps.R
        └── 02_snpfiltr.R  
        └── 03_mpcrselect.config 
        └── 04_snpaimer.R  
```
---

## Getting Started

### Prerequisites

The analysis was developed in R 4.6.1 on 64-bit Manjaro Linux
(`x86_64-pc-linux-gnu`; kernel `7.1.13-2-MANJARO`). RStudio Desktop is
recommended but not required; the R analyses can also be run from a standard
R session or the command line.

R package dependencies are managed with
[renv](https://rstudio.github.io/renv/). If `renv.lock` is included in this
repository, it records the package versions and installation sources used for
the R-based analyses. An internet connection is generally required for the
initial dependency restore.

The workflow comprises four components:

1. **SNP filtering:** [SNPfiltR](https://cran.r-project.org/package=SNPfiltR)
   and associated R scripts were used to inspect and filter genotype data.
2. **Candidate marker selection:** [mPCRselect](https://github.com/ellieearmstrong/mPCRselect)
   was run as a Nextflow workflow in a bioinformatics high-performance
   computing (HPC) environment. Reproducing this stage as run in the study
   requires access to a suitably configured cluster, the pipeline's external
   software dependencies, and the corresponding configuration files.
   Restoring `renv.lock` alone does not install these HPC dependencies.
3. **Diagnostic panel evaluation:** [snpAIMeR](https://cran.r-project.org/package=snpAIMeR)
   and associated R scripts were used to evaluate marker combinations and
   geographic assignment performance.
4. **Maps and figures:** R scripts were used to generate sampling maps and
   other figures. Recreating the maps may require access to externally
   sourced coastline and bathymetric data if these are not provided with
   the repository.

To restore the R environment from the repository root, run:

```r
install.packages("renv")
renv::restore()
```

The `mPCRselect` stage must be configured separately for the available HPC
environment before the complete workflow can be rerun.

### Installation

Build octotrace from the source and intsall dependencies:

1. **Clone the repository:**

    ```sh
    ❯ git clone https://github.com/mparrondo/octotrace/
    ```

2. **Navigate to the project directory:**

    ```sh
    ❯ cd octotrace
    ```

3. **Install the dependencies:**
This project uses the renv package to provide a reproducible R environment. Package versions and sources are recorded in renv.lock. After downloading or cloning the repository, open octotrace.Rproj in RStudio (or set the working directory to the project root) and restore the required packages by running:
```r
install.packages("renv")  # Run only if renv is not already installed
renv::restore()
```
This step only needs to be performed once when setting up the project.

The `mPCRselect` stage must be configured separately for the available HPC
environment before the complete workflow can be rerun.

---

## Contributing

- **💬 [Join the Discussions](https://github.com/mparrondo/octotrace/discussions)**: Share your insights, provide feedback, or ask questions.
- **🐛 [Report Issues](https://github.com/mparrondo/octotrace/issues)**: Submit bugs found or log feature requests for the `octotrace` project.
- **💡 [Submit Pull Requests](https://github.com/mparrondo/octotrace/blob/main/CONTRIBUTING.md)**: Review open PRs, and submit your own PRs.

<details closed>
<summary>Contributing Guidelines</summary>

1. **Fork the Repository**: Start by forking the project repository to your github account.
2. **Clone Locally**: Clone the forked repository to your local machine using a git client.
   ```sh
   git clone https://github.com/mparrondo/octotrace/
   ```
3. **Create a New Branch**: Always work on a new branch, giving it a descriptive name.
   ```sh
   git checkout -b new-feature-x
   ```
4. **Make Your Changes**: Develop and test your changes locally.
5. **Commit Your Changes**: Commit with a clear message describing your updates.
   ```sh
   git commit -m 'Implemented new feature x.'
   ```
6. **Push to github**: Push the changes to your forked repository.
   ```sh
   git push origin new-feature-x
   ```
7. **Submit a Pull Request**: Create a PR against the original project repository. Clearly describe the changes and their motivations.
8. **Review**: Once your PR is reviewed and approved, it will be merged into the main branch. Congratulations on your contribution!
</details>

<details closed>
<summary>Contributor Graph</summary>
<br>
<p align="left">
   <a href="https://github.com{/mparrondo/octotrace/}graphs/contributors">
      <img src="https://contrib.rocks/image?repo=mparrondo/octotrace">
   </a>
</p>
</details>

---

## License

Octotrace is protected under the [LICENSE](https://choosealicense.com/licenses) License. For more details, refer to the [LICENSE](https://choosealicense.com/licenses/) file.

---

## Acknowledgments

- Credit `contributors`, `inspiration`, `references`, etc.

<div align="right">

[![][back-to-top]](#top)

</div>


[back-to-top]: https://img.shields.io/badge/-BACK_TO_TOP-151515?style=flat-square


---
