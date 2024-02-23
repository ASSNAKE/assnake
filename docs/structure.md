# Assnake Framework Repository Structure

The Assnake Framework repository is organized into several directories and files, each serving a specific purpose in the framework's ecosystem. Below is a detailed breakdown of the structure:

```
├── LICENSE.txt                  # License information for the project
├── MANIFEST.in                  # Manifest file for including non-code files in the distribution
├── README.md                    # Main README file with project overview and instructions
├── assnake                       # Main package directory
│   ├── __init__.py              # Initializes the package
│   ├── cli                      # Command Line Interface directory
│   │   ├── __init__.py          # Initializes the CLI sub-package
│   │   ├── assnake_cli.py       # Main CLI script, entry point for commands
│   │   ├── command_builder.py   # Helper functions to build CLI commands
│   │   └── commands             # Directory containing various CLI command implementations
│   │       ├── __init__.py      # Initializes the commands sub-package
│   │       ├── config_commands.py  # Commands related to configuration management
│   │       ├── dataset_commands.py # Commands for dataset operations
│   │       ├── execute_commands.py # Commands for executing pipelines and results
│   │       └── module_commands.py  # Commands for handling Assnake modules
│   ├── core                     # Core functionalities of Assnake
│   │   ├── Dataset.py           # Dataset class definition and related functions
│   │   ├── Pipeline.py          # Pipeline class for handling analysis workflows
│   │   ├── PresetManager.py     # Manages presets for various analyses
│   │   ├── Result.py            # Defines individual analysis steps (results)
│   │   ├── Sample.py            # Represents biological samples and their data
│   │   ├── SampleContainerSet.py# Manages groups of samples for analysis
│   │   ├── __init__.py          # Initializes the core sub-package
│   │   ├── config.py            # Configuration management functions
│   │   ├── exceptions.py        # Custom exceptions for error handling
│   │   └── snake_module.py      # Handles the integration of external modules
│   ├── snake                    # Snakemake-related files
│   │   ├── config_template.yml  # Template for instance configuration
│   │   ├── snake_base.py        # Base snakemake file for defining workflows
│   │   └── wc_config.yaml       # Configuration for wildcard paths in Snakemake
│   └── utils                    # Utility functions and helpers
│       ├── fs_helpers.py        # File system helper functions
│       └── general.py           # General utility functions
├── docs                         # Documentation directory
│   ├── architecture.md          # Documentation on system architecture
│   ├── params_management.md     # Documentation on parameter management
│   ├── structure.md             # Documentation on project structure
│   └── thoughts.md              # Thoughts and notes on development
├── fastentrypoints.py           # Script for faster CLI entry point loading
├── old_README.md                # Older version of README
├── setup.cfg                    # Configuration for package setup
├── setup.py                     # Script for installing the package
└── tests                        # Test suite for the package
    ├── conftest.py              # Configuration for pytest
    ├── pytest.ini               # INI configuration for pytest
    ├── small_tests.py           # Small unit tests
    ├── test_dataset.py          # Tests related to Dataset functionality
    ├── test_init.py             # Tests for initialization process
    ├── test_snake.py            # Tests for Snakemake integration
    └── util_for_test.py         # Utility functions for tests
    
```
