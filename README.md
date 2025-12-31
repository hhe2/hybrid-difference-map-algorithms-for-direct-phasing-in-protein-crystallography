# Crystallographic Phase Reconstruction (Under Construction)

A high-performance computational tool for solving the crystallographic phase problem through iterative density modification algorithms. This software implements hybrid input-output (HIO), difference map, and advanced density modification techniques for ab initio phasing from X-ray diffraction data.

## Overview

Determining macromolecular structures from X-ray crystallography requires both diffraction amplitudes and phases. While amplitudes are directly measurable, phases must be computed. This program employs iterative algorithms combined with real-space constraints to retrieve phases from experimental amplitudes, enabling structure determination without prior molecular replacement models.

## Key Features

- **Multiple Phasing Algorithms**: Classical Hybrid Input-Output (HIO), Continuous HIO, Modified Continuous HIO, Transition HIO, Hybrid Projection Reflection (HPR), Difference Map(DM), Modified DM, Averaged Successive Reflections (ASR), RAAR, Modified ASR/RAAR, and Hybrid Difference Map (HDM) methods
- **MPI Parallelization**: Distributed computing support for rapid convergence
- **Advanced Constraints**: Solvent flattening, histogram matching, and density limiting
- **Genetic Algorithm**: Optional evolutionary optimization for improved phase solutions
- **Flexible Initialization**: Start from random maps, PDB structures, or molecular replacement
- **Real-time Monitoring**: Track R-factors, phase errors, and convergence metrics
- **Standard Formats**: Compatible with MTZ and CIF file formats

## Algorithm Implementation

The software implements three core iterative projection algorithms:

1. **Hybrid Input-Output (HIO)**: Applies negative feedback in solvent regions while maintaining calculated density in protein regions. The feedback parameter β controls convergence behavior.

2. **Difference Map**: Alternates between Fourier-space constraints (observed amplitudes) and real-space constraints (density positivity, solvent flatness) using sophisticated projection operators.

3. **Hybrid Difference Map**: Combines advantages of both HIO and Difference Map approaches for enhanced convergence properties.

All methods incorporate:
- Solvent content requirement:
Solvent content must exceed 60%.

- Environment variable setup for Linux (required):
export LD_LIBRARY_PATH=../ccp4_lib/lib:$LD_LIBRARY_PATH

- Environment variable setup for macOS (optional):
export DYLD_LIBRARY_PATH="/usr/local/Cellar/gcc/14.2.0_1/lib/gcc/14:../ccp4_lib/lib:$DYLD_LIBRARY_PATH"

- MPI initialization (optional, for parallel execution):
source /public/software/profile.d/mpi_intelmpi-2021.3.0.sh

- Phase retrieval procedure:
Step 1: Edit the configuration file parameters.txt.
Step 2: Execute the program with MPI parallelization:
mpirun -np 100 ./program