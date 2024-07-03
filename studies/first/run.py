import sys
from pathlib import Path

# Add the bluemira installation path to the system path
bluemira_path = Path("/home/graeme/mambaforge/envs/bluemira/")
sys.path.append(str(bluemira_path))
# Add the project root to the PYTHONPATH
project_root = Path(__file__).parents[2].resolve()
sys.path.append(str(project_root))

from stellarator_project.reactor import main

build_config_path = Path(project_root, "config/config.json").resolve()
reactor = main(build_config_path)
