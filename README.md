# stellarator-project

---

**Table of Contents**

- [Installation](#installation)
- [License](#license)

## Setting up and using this template

Start by git cloning this repository (or use the "Use this template" button on GitHub).

```bash
git clone git@github.com:Fusion-Power-Plant-Framework/bluemira-stellarator-project.git
cd bluemira-stellarator-project
```

The following strings should set using a "Find and Replace" tool:

| String to place              | Description                                                                               |
|------------------------------|-------------------------------------------------------------------------------------------|
| stellarator-project             | The project name                                                                          |
| stellarator_project             | The main project directory (usally the project name, using underscores instead of dashes) |
| Graeme Turkington            | The name of the author of the project. If there are multiple, add each to the list.       |
| <graeme.turkington@ukaea.uk>           | The email address of the author                                                           |
| {{ copyright-holder }}       | The copyright holder (may be the same as the author(s))                                   |
| {{ copyright-holder-email }} | The copyright holder contact email (may be the same as the author(s))                     |
| {{ org-name }}               | The GitHub organisation. Only applies if using GitHub.                                    |

***All relevant directories should be set/updated too.**

Once the repository is cloned please run the install scripts as show

```bash
bash scripts/install_bluemira.sh
```

In future we hope to extend this script to improve the setup experience.

## Running reactor designs

The example study can be run as shown:

```
python studies/first/run.py
```

## Running tests

TODO
