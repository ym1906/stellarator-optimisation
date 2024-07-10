# stellarator-project

---

**Table of Contents**

- [Installation](#installation)
- [License](#license)

## Running reactor designs

In SIMSOPT, you can import and optimize a boundary shape from another code such as VMEC or generate or generate one from scratch.Similarly, you can import coils or generate generic ones and optimise them.

We build a HELIAS 5B like stellarator in the example. This takes a plasma scenario, which is imported from VMEC output, and generates and optimises a set of coil required to achieve the magnetic properites for the plasma. The plasma surface and coils are converted into NURBS surfaces and curves and this data is used to generate CAD for use in Bluemira.

```bash
activate simsoptenv
python stellarator_project/helias_plasma_coil_finite_optimisation.py
```

The above script loads the VMEC file, and create a surface object in the SIMSOPT framework. Then, some generic coils and coil filaments are generated. These are optimised for simplicity (to aid construction), length and closeness. The coils have a finite build, comprised of filaments. This is then converted to NURBS data, which is stored and ready to be read by Bluemira.

Bluemira can created CAD objects, from the NURBS data.
First, active the bluemira enviornment and run the example:

```bash
conda activate bluemira
python stellarator_project/stellarator_build.ex.py
```

The generated CAD will appear in a display. This CAD can facilitate further reactor component design and optimisation.

## Running tests

TODO
