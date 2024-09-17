# SPDX-FileCopyrightText: 2024-present {{ copyright-holder }} <{{ copyright-holder-email }}>
#
# SPDX-License-Identifier: MIT
"""First study."""

from pathlib import Path

from stellarator_project.stellarator import main

Path(Path(__file__).parent, Path("data/cad_ouput")).mkdir(parents=True, exist_ok=True)
Path(Path(__file__).parent, Path("data/coils")).mkdir(parents=True, exist_ok=True)
Path(Path(__file__).parent, Path("data/plasma")).mkdir(parents=True, exist_ok=True)
build_config_path = Path(__file__, "../config/config.json").resolve()

reactor = main(build_config_path)
reactor.show_cad()
reactor.save_cad(
    "xyz", filename=Path(__file__, "../data/cad_ouput/plasmastellarator.stp")
)
