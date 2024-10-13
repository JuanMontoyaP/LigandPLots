"""
This file generates all the data for plotting
"""
import logging
from pathlib import Path
import docker

from helpers import create_folder

logging.basicConfig(
    format='[%(asctime)s] - %(levelname)s - %(filename)s - %(funcName)s:%(lineno)d - %(message)s',
    level=logging.INFO
)

class GromacsData:
    CONTAINER_WORK_DIR = "/container/data"
    XTC_FILE = "final.xtc"
    XTC_NO_PBC_FILE = "fixed.xtc"

    def __init__(self, path: str) -> None:
        self.path = path
        self.volume: dict = {
            self.path: {
                'bind': self.CONTAINER_WORK_DIR,
                'mode': 'rw'
            }
        }

        create_folder(self.path)
        self.results_folder = Path(f"{self.path}/results")

    def _run_gromacs_container(
        self,
        command,
        working_dir=CONTAINER_WORK_DIR,
        **kwargs
    ) -> docker.models.containers.Container:
        """
        Runs a Gromacs container with the specified command.

        Parameters:
            command (str): The command to be executed inside the container.
            working_dir (str, optional): The working directory inside the container. Defaults
                to CONTAINER_WORK_DIR.
            **kwargs: Additional keyword arguments to be passed to the container run method.

        Returns:
            docker.models.containers.Container: The container object representing the
                running Gromacs container.
        """
        gromacs_image = "jpmontoya19/gromacs:latest"
        client = docker.from_env()
        return client.containers.run(
            gromacs_image,
            command,
            working_dir=working_dir,
            **kwargs
        )

    def generate_minimization_data(self) -> str:
        """
        Generate energy minimization data.
        """
        command = [
            "sh",
            "-c",
            f"echo 10 0 | gmx energy -f em.edr -o {self.results_folder.name}/potential.xvg"
        ]

        logging.info(command[-1])
        energy_minimization = self._run_gromacs_container(
            command,
            volumes=self.volume
        ).decode('utf-8')
        return energy_minimization