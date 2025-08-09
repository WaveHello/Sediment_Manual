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

## Chapter 3 Implementation Plan - Transport of Gravel and Sediment Mixtures

### Chapter Structure Analysis (18 chunks, 88 pages)
Based on analysis of first 3 chunks (pages 1-15):

**Foundation Sections:**
- Section 3.1: Fluvial Phenomena (surface armoring, sorting patterns)
- Section 3.2: Engineering Relevance (dam effects, mining impacts)
- Section 3.3: Grain-Size Distributions (18+ equations for statistical analysis)
- Section 3.4: Dimensionless Bank-Full Relations (hydraulic characterization)

### Proposed Module Structure
```
src/sandwave/chapter03_gravel_mixtures/
├── __init__.py
├── grain_size_distributions.py      # Equations 3-1 to 3-12+ 
├── hydraulic_relations.py           # Equations 3-13 to 3-18+
├── sorting_phenomena.py             # Field patterns, armoring processes
├── transport_mechanics.py           # Core transport equations (from later chunks)
├── engineering_applications.py      # Dam effects, mining impacts
└── sampling_methods.py              # Sediment sampling techniques
```

### Implementation Strategy

**Phase 1: Foundation (Chunks 1-6) - ✅ COMPLETE**
1. ✅ `grain_size_distributions.py` - 12 functions, Equations 3-1 to 3-12, 3-29
2. ✅ `hydraulic_relations.py` - 8 functions, Equations 3-13 to 3-22  
3. ✅ `active_layer_mechanics.py` - 10 functions, Equations 3-23 to 3-44
4. ✅ `hiding_functions.py` - 10 functions, Equations 3-60 to 3-74
5. ✅ `flow_hydraulics.py` - 10 functions, Equations 3-75 to 3-76

**Phase 2: Core Transport (Chunks 7-12) - ✅ COMPLETE**
6. ✅ Extract transport equations from chunks 7-12 (30+ equations)
7. ✅ `transport_mechanics.py` - 16 functions with major bed-load transport models

**Phase 3: Applications (Chunks 13-18) - NEXT**  
8. `engineering_applications.py` - Practical engineering tools
9. `sampling_methods.py` - Field measurement techniques

### ✅ COMPLETED: Chapter 3 Core Implementation (Chunks 1-12)
- **66 Python functions implemented** across 6 core modules  
- **~105 equations covered** from grain size statistics to bed-load transport
- **Complete module structure** with proper imports and documentation
- **Input validation and error handling** throughout all functions

**Phase 1 Features (Chunks 1-6):**
- Grain size distributions: ψ scale, lognormal stats, discretization
- Hydraulic relations: bankfull flow, dimensionless parameters, stream classification
- Active layer mechanics: conservation equations, entrainment-deposition, transport rates
- Hiding functions: Egiazaroff, power law, threshold conditions
- Flow hydraulics: shear stress, friction coefficients, 1D flow equations

**Phase 2 Features (Chunks 7-12):**
- Major transport models: Ashida-Michiue, Parker, Wilcock-Crowe, Wu-Wang-Jia
- Hiding functions: Modified Egiazaroff, power law, substrate-based
- Morphologic complexity corrections and transport enhancement
- Equal mobility conditions and static armor predictions
- Generic transport relations with customizable parameters

## Next Steps
1. ✅ Set up pixi environment (COMPLETE)
2. ✅ Begin PDF content analysis (Chapter 3 structure analyzed)
3. ✅ Chapter 3 Phase 1 implementation (COMPLETE - 50 functions, ~75 equations)
4. ✅ Chapter 3 Phase 2 implementation (COMPLETE - 16 transport functions, ~30 equations)
5. **Continue with chunks 13-18 for engineering applications (Phase 3)**