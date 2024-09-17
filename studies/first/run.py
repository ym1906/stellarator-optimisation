# SPDX-FileCopyrightText: 2024-present {{ copyright-holder }} <{{ copyright-holder-email }}>
#
# SPDX-License-Identifier: MIT
"""First study."""

from pathlib import Path

from stellarator_project.stellarator import main

Path(Path(__file__).parent, "data/cad_output").mkdir(parents=True, exist_ok=True)
Path(Path(__file__).parent, "data/coils").mkdir(parents=True, exist_ok=True)
Path(Path(__file__).parent, "data/plasma").mkdir(parents=True, exist_ok=True)
build_config_path = Path(__file__, "../config/config.json").resolve()

reactor = main(build_config_path)
reactor.show_cad()
reactor.save_cad(
    "xyz",
    filename=Path(
        Path(__file__).parent, "../data/cad_output/plasmastellarator.stp"
    ).resolve(),
)
