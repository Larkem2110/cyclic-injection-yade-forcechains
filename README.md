# Cyclic Injection YADE Force-Chain Analysis

This repository contains YADE DEM scripts used for direct shear simulation, cyclic fluid injection, and force-chain tracking in granular fault gouge simulations.

## Repository structure

```text
scripts/
  01_direct_shear_test.py
  02_cyclic_injection.py
  03_force_chain_tracking.py
```

## Scripts

### 01_direct_shear_test.py

Generates a periodic granular gouge sample, applies normal loading, and performs a direct shear test. The final state is saved as `laststate.yade`.

### 02_cyclic_injection.py

Loads `laststate.yade`, relaxes the shear stress to a prescribed fraction of the initial shear stress, and applies cyclic pore-pressure injection using YADE FlowEngine.

### 03_force_chain_tracking.py

Post-processes saved `.yade` state files to detect and track force chains using per-particle stress tensors, a stress threshold, contact topology, and an alpha-angle alignment criterion.

## Requirements

- YADE DEM
- Python
- NumPy

## Usage

The scripts are intended to be used in the following order:

```bash
yade 01_direct_shear_test.py
yade 02_cyclic_injection.py
yade 03_force_chain_tracking.py
```

The force-chain tracking script can also be applied to other saved YADE state files, provided that the required particle, contact, and stress information is available.

## Outputs

The scripts produce YADE state files, text files containing macroscopic quantities, and VTU files for visualization in ParaView.

## Citation

Please cite the archived Zenodo release of this repository if you use these scripts:

https://doi.org/10.5281/zenodo.20534192
