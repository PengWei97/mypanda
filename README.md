mypanda
=====

Fork "mypanda" to create a new MOOSE-based application.

For more information see: [https://mooseframework.org/getting_started/new_users.html#create-an-app](https://mooseframework.org/getting_started/new_users.html#create-an-app)

---

## Basic MOOSE

### Some commonly used commanda lines

```bash
# update moose
mamba update --all --yes
git pull origin next
make -j10
./run_test -j10

# TODO
```
---
## Upcoming & Ongoing works

### Upcoming Work:
1. **Polycrystalline Simulations in Coupled PC-PF Model:**
   - **Objective:** Extend the current bicrystal simulation capabilities to support polycrystalline structures.
   - **Details:** 
     - Implement the framework to handle multiple grains in the coupled crystal plasticity–phase field (CP-FP) model.
     - Aim to simulate more complex microstructure behaviors and interactions within polycrystalline materials.
   - **References:** Based on the methodology outlined in the paper "[Coupled crystal plasticity–phase field simulation of microstructure evolution in dual-phase steels](https://www.sciencedirect.com/science/article/pii/S092702561630605X)".

2. **Incorporation of Grain Rotation Mechanism:**
   - **Objective:** Enhance the simulation model by including the mechanism of grain rotation.
   - **Details:**
     - Develop and integrate models that account for the rotational behavior of grains.
     - Improve the realism and accuracy of microstructure evolution simulations by considering grain rotation.

### Ongoing Work:
1. **Refactoring Code in `prm2_ThermalGNS_AGG2023`:**
   - **Objective:** Improve the existing codebase for better performance, readability, and maintainability.
   - **Details:**
     - Refactor the code related to the thermal grain boundary network structure (GNS) simulations.
     - Ensure that the refactored code is modular, efficient, and adheres to best practices in software development.
   - **Current Progress:**
     - Identified key areas for improvement.
     - Working on breaking down monolithic code into smaller, reusable functions.
     - Conducting tests to ensure that refactoring does not affect the functionality of the simulations.
