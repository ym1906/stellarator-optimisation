# stellarator-project

---

**Table of Contents**

- [Installation](#installation)
- [License](#license)

## Installation

Once you have a bluemira conda environemt setup you can install this project with

```bash
pip install git+https://github.com/ym1906/stellarator-optimisation.git
```
## Running reactor designs

In SIMSOPT, you can import and optimize a boundary shape from another code such as VMEC or generate one from scratch. Similarly, you can import coils or create generic ones and optimize them.

In this example, we build a HELIAS 5B-like stellarator. We start by importing a plasma scenario from VMEC output, then generate and optimize a set of coils to achieve the desired magnetic properties for the plasma. The plasma surface and coils are converted into NURBS surfaces and curves, which are then used to generate CAD data for Bluemira.

This project loads the VMEC file and creates a surface object within the SIMSOPT framework. Next, it generates some generic coils and coil filaments, optimizing them for simplicity (to aid construction), length, and proximity. The coils, which have a finite build comprised of filaments, are then converted to NURBS data, ready for Bluemira.

Bluemira creates the CAD objects from the NURBS data.

```bash
conda activate bluemira
python studies/first/run.py
```

To generate the model there are two run modes which are set in the config json to 'read' to use existing datafiles or 'run' to run the optimisation.

The generated CAD will appear in a display, facilitating further reactor component design and optimisation.

![Alt text](https://github.com/ym1906/stellarator/blob/main/docs/stellarator_coils_cad.png)
