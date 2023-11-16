## Architecture Overview

### Core Concepts and Classes

1. **Dataset (`Dataset.py`):**
   - Represents a collection of biological samples.
   - Acts as a container for multiple samples, each representing a unique biological entity.

2. **Sample (`Sample.py`):**
   - Represents an individual biological sample within a dataset.
   - Each `Sample` is unique within its `Dataset`.
   - Contains data specific to a single biological entity and is linked to various SampleContainers.

3. **SampleContainer (`SampleContainer.py`):**
   - Encapsulates the sequencing data for a particular preprocessing step of a `Sample`.
   - A `Sample` can have multiple `SampleContainer` instances, each representing the sample processed through different pipelines or preprocessing steps.
   - Holds information such as file paths for read files and metadata about the sample's dataset.

4. **Pipeline (`Pipeline.py`):**
   - Defines a sequence of processing and analytical steps applied to a dataset.
   - Manages the workflow from raw data to final analysis, applying various preprocessing and analysis steps.

5. **Result (`Result.py`):**
   - Specifies a particular computational result or output within a pipeline.
   - Linked to `Pipeline` steps and generates specific outputs from samples or sample sets.

6. **PresetManager (`PresetManager.py`):**
   - Manages preset configurations for different results and pipelines.
   - Facilitates consistent and reproducible processing steps.

### Preprocessing in Assnake
- **Preprocessing Chains:** Integral to the Pipeline, these are sequences of steps for preparing raw data for analysis. For example, a chain might involve raw data going through 'Cutadapt' for primer removal, followed by 'DADA2 Filter and Trim' for quality control.
- **Example Flow:** raw data -> Cutadapt (rmV3V4primers preset) -> DADA2 Filter and Trim (strict preset). Each step transforms the data, leading to new `SampleContainer` instances reflecting these changes.

### System Interconnections and Workflow

```
Dataset
  │
  ├── Sample 1
  │     │
  │     ├── SampleContainer (Preprocessing Step 1)
  │     │
  │     └── SampleContainer (Preprocessing Step 2)
  │
  ├── Sample 2
  │     │
  │     ├── SampleContainer (Preprocessing Step 1)
  │     │
  │     └── SampleContainer (Preprocessing Step 2)
  │
  └── ... (More samples)
```

- **Workflow Processing:**
  - A `Dataset` contains multiple `Sample` instances.
  - Each `Sample` has several `SampleContainer` instances, each representing a different preprocessing stage.
  - The `Pipeline` manages the overall workflow, processing each `Sample` through various stages.

- **Results and Configuration:**
  - `Result` objects define the outputs for each stage in the `Pipeline`.
  - The `PresetManager` ensures consistent parameterization for each `Result`.

- **CLI and Configuration Management:**
  - The system is interacted with through a CLI (`assnake_cli.py`), providing access to all functionalities.
  - Configuration files managed by `config.py` store settings, ensuring consistent behavior and reproducibility.

This architecture ensures that each Sample in a Dataset is processed through a series of steps, with each step captured as a SampleContainer. The modular design allows for flexibility in data processing and analysis, tailored to specific research needs.