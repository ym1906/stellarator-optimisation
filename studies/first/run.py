# SPDX-FileCopyrightText: 2024-present {{ copyright-holder }} <{{ copyright-holder-email }}>
#
# SPDX-License-Identifier: MIT
"""First study."""

from pathlib import Path

from stellarator_project.stellarator import main

build_config_path = Path(__file__, "config/config.json").resolve()
reactor = main(build_config_path)
