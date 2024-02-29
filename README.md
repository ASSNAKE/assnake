# Assnake Metagenomic Analysis System

## Overview
Assnake is a comprehensive metagenomics analysis tool designed for short read sequencing data, focusing on microbiome studies. It combines various bioinformatics tools into an integrated workflow for metagenomic data processing and analysis.

## Key Features
- **Dataset and Sample Management:** Organizes datasets and individual biological samples.
- **Modular Preprocessing and Analysis:** Customizable preprocessing and analysis steps tailored to specific research needs.
- **Preset Management:** Ensures consistent and reproducible analyses through preset configurations.
- **Command Line Interface (CLI):** User-friendly CLI for easy access to functionalities.

## System Architecture
Assnake's architecture is centered around datasets containing unique samples. Each sample can go through several preprocessing stages, captured as `SampleContainers`. Pipelines manage the workflow, applying various `Results` and `Presets` for analysis.

## Installation
1. **Create a Conda Environment:** (Recommended to manage dependencies)
   ```bash
   conda create -n assnake_env python=3.10
   ```

2. **Activate your new Conda Environment:** (Recommended to manage dependencies)
   ```bash
   conda activate assnake_env
   ```

3. **Clone the Repository:**
   ```bash
   git clone https://github.com/ASSNAKE/assnake.git
   ```

4. **Navigate to the Repository and Install:**
   ```bash
   cd assnake
   pip install -e ./
   ```

5. **Initialize Assnake Configuration:**
   ```bash
   assnake config init
   ```

   It will ask you to specify the folder for assnake database. This can grow quite large (10Gb +) due to the conda envs that will be installed into this folder, so choose wisely, contact your system administrator if you are unsure about the directory choice. 


## Quick Start
- To create and manage datasets: `assnake dataset <command>`
- To import sequencing data: `assnake dataset import-reads`
- To execute analysis pipelines: `assnake result <analysis-command>`

## Documentation
Detailed documentation, including usage guidelines and advanced configurations, is available in the `docs` directory. Examples and tutorials are provided for guidance.

## Contributing
Contributions to Assnake are welcome, including documentation improvements, new features, and bug reporting. Please refer to the contributing guidelines for more information.

## License
Assnake is open-source software, available under the MIT License. See `LICENSE.txt` for details.

## Support and Contact
For support or inquiries, please open an issue on the GitHub repository or contact [fedorov.de@gmail.com](mailto:fedorov.de@gmail.com).

---
*Note: This README offers a general overview of the Assnake system. For comprehensive information, refer to the individual documentation files within the repository.*