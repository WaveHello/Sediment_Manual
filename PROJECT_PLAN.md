# SandWave: Sediment Transport Equations Python Library

## Project Overview
SandWave is a comprehensive Python library that extracts and implements mathematical equations from "Sedimentation engineering processes, measurements.pdf" (ASCE Manual 110) for sediment transport calculations.

## Goals
- Parse and extract equations from the sedimentation engineering PDF
- Create a modular Python library with clear function implementations
- Use pixi for dependency management
- Maintain version control with git
- Provide comprehensive documentation and examples

## Project Structure (Chapter-Based Organization)
```
sandwave/
├── src/
│   └── sandwave/
│       ├── __init__.py
│       ├── chapter01_overview/
│       │   └── __init__.py
│       ├── chapter02_transport_morphodynamics/
│       │   └── __init__.py
│       ├── chapter03_gravel_mixtures/
│       │   └── __init__.py
│       ├── chapter04_fine_grained/
│       │   └── __init__.py
│       ├── chapter05_measurements/
│       │   └── __init__.py
│       ├── chapter06_fluvial_geomorphology/
│       │   └── __init__.py
│       ├── chapter07_streambank_erosion/
│       │   └── __init__.py
│       ├── chapter08_river_meandering/
│       │   └── __init__.py
│       ├── chapter09_stream_restoration/
│       │   └── __init__.py
│       ├── chapter10_bridge_scour/
│       │   └── __init__.py
│       ├── chapter11_bridge_countermeasures/
│       │   └── __init__.py
│       ├── chapter12_reservoir/
│       │   └── __init__.py
│       ├── chapter13_ice_effects/
│       │   └── __init__.py
│       ├── chapter14_computational/
│       │   └── __init__.py
│       ├── chapter15_numerical_simulation/
│       │   └── __init__.py
│       ├── chapter16_turbulence/
│       │   └── __init__.py
│       ├── chapter17_watershed/
│       │   └── __init__.py
│       ├── chapter18_engineering_geomorphology/
│       │   └── __init__.py
│       ├── chapter19_sedimentation_hazards/
│       │   └── __init__.py
│       ├── chapter20_sedimentation_law/
│       │   └── __init__.py
│       ├── chapter21_contaminant_processes/
│       │   └── __init__.py
│       ├── chapter22_sediment_oxygen_demand/
│       │   └── __init__.py
│       ├── chapter23_dam_removal/
│       │   └── __init__.py
│       ├── appendices/
│       │   ├── rock_scour.py
│       │   ├── riprap_design.py
│       │   ├── scaling.py
│       │   └── sediment_discharge.py
│       └── utils/
│           ├── __init__.py
│           ├── constants.py
│           ├── validators.py
│           └── unit_conversions.py
├── tests/
├── docs/
├── examples/
├── chapters/ (split PDFs)
├── pixi.toml
├── pyproject.toml
└── README.md
```

## Current Status (Updated 2025-08-09)

### ✅ COMPLETED: Chapter 2 - Sediment Transport and Morphodynamics
- **50 files committed** with comprehensive implementation
- **5 core modules** with 25+ equations implemented:
  - `fluid_mechanics.py` - Reynolds, Froude numbers, shear calculations (Equations 2.1, 2.2, 2.7, 2.8)  
  - `sediment_properties.py` - Settling velocities, particle properties (Equations 2.17-2.19)
  - `bed_load_transport.py` - Meyer-Peter Müller, Einstein-Brown transport (Equations 2.34-2.37)
  - `suspended_load.py` - Rouse profiles, concentration distributions (Equations 2.5, 2.60-2.62)
  - `bed_forms.py` - Dune geometry, roughness, friction factors (Equations 2.40-2.45)
- **Utility framework** with constants, validators, unit conversions
- **Full equation traceability** with ASCE Manual 110 references

### ✅ LATEST UPDATE
- **pyproject.toml added** - Complete Python packaging configuration with hatchling, dev tools, and optional dependencies

### 📋 REMAINING WORK  
- **21 more chapters** to implement (~460+ equations estimated)
- Testing framework and comprehensive test suite
- Documentation generation and examples
- Performance optimization and benchmarking

## Implementation Checklist

### Phase 1: Project Setup
- [x] Initialize git repository
- [x] Create .gitignore (including PDFs)
- [x] Set up pixi project structure
- [x] Create pyproject.toml configuration
- [x] Set up basic Python package structure

### Phase 2: PDF Analysis & Equation Extraction - Chapter 2 COMPLETE
- [x] Read and analyze PDF content (Chapter 2)
- [x] Identify major equation categories:
  - [x] Particle settling velocity equations (Stokes, Dietrich)
  - [x] Bed load transport equations (Meyer-Peter Müller, Einstein-Brown)
  - [x] Suspended sediment transport equations (Rouse profiles, van Rijn)
  - [x] Flow velocity and shear stress equations (logarithmic profiles)
  - [x] Particle size distribution equations (dimensionless diameter)
  - [x] Turbulence and mixing equations (von Kármán constant)
- [x] Extract equations with their context and variable definitions
- [x] Document equation sources with proper equation numbers (2.1-2.62)

### Phase 3: Python Library Development - Chapter 2 COMPLETE
- [x] Implement core equation modules (Chapter 2: Transport & Morphodynamics)
- [x] Add input validation and error handling
- [x] Create utility functions for common calculations
- [x] Define physical constants and default parameters
- [x] Implement unit conversion helpers

### Phase 4: Testing & Documentation
- [ ] Write comprehensive unit tests
- [ ] Create integration tests with realistic scenarios
- [ ] Generate API documentation
- [ ] Create usage examples and tutorials
- [ ] Add docstrings with equation references

### Phase 5: Quality Assurance
- [ ] Code review and refactoring
- [ ] Performance optimization
- [ ] Validate against known solutions/benchmarks
- [ ] Final documentation review

## Technical Dependencies (via pixi)
- Python 3.9+
- NumPy (numerical computations)
- SciPy (scientific functions)
- PyMuPDF/fitz (PDF parsing)
- Matplotlib (plotting/visualization)
- Pytest (testing)
- Sphinx (documentation)

## Next Steps
1. Set up pixi environment
2. Begin PDF content analysis
3. Start with particle settling equations as proof of concept