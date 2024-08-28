# stellarator-project

---

**Table of Contents**

- [Installation](#installation)
- [License](#license)

## Running reactor designs

In SIMSOPT, you can import and optimize a boundary shape from another code such as VMEC or generate one from scratch. Similarly, you can import coils or create generic ones and optimize them.

In this example, we build a HELIAS 5B-like stellarator. We start by importing a plasma scenario from VMEC output, then generate and optimize a set of coils to achieve the desired magnetic properties for the plasma. The plasma surface and coils are converted into NURBS surfaces and curves, which are then used to generate CAD data for Bluemira.

```bash
activate simsoptenv
python stellarator_project/helias_plasma_coil_finite_optimisation.py
```

The script above loads the VMEC file and creates a surface object within the SIMSOPT framework. Next, it generates some generic coils and coil filaments, optimizing them for simplicity (to aid construction), length, and proximity. The coils, which have a finite build comprised of filaments, are then converted to NURBS data, ready for Bluemira.

Bluemira can create CAD objects from the NURBS data. First, activate the Bluemira environment and run the example:

```bash
conda activate bluemira
python stellarator_project/stellarator_build.ex.py
```

The generated CAD will appear in a display, facilitating further reactor component design and optimization.

![Alt text](https://github.com/ym1906/stellarator/blob/main/docs/stellarator_coils_cad.png)
