"""
This file generates all the data for plotting
"""
import logging
from pathlib import Path
import tarfile
import docker
import numpy as np
import pandas as pd
import MDAnalysis as mda
from MDAnalysis.analysis import rms

from helpers import create_folder, find_files_with_same_pattern

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
    
    def _find_last_tpr_file(self) -> Path:
        """
        Finds the last TPR file in the given path.

        Returns:
            A tuple containing the name and path of the last TPR file.
        """
        last_tpr = find_files_with_same_pattern(self.path, "md_0_*.tpr")[-1]
        return (last_tpr.name, last_tpr)

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
    
    def generate_nvt_data(self):
        """
        Function to generate the temperature data
        """
        command = [
            "sh",
            "-c",
            f"""
                echo 16 0 | \
                gmx energy -f nvt.edr -o {self.results_folder.name}/temperature.xvg
            """
        ]

        logging.info(command[-1])
        temperature = self._run_gromacs_container(
            command,
            volumes=self.volume
        ).decode("utf-8")
        return temperature

    def generate_npt_data(self):
        """
        Function to generate the pressure & density data
        """
        command_pressure = [
            "sh",
            "-c",
            f"""
                echo 18 0 | \
                gmx energy -f npt.edr -o {self.results_folder.name}/pressure.xvg
            """
        ]

        logging.info(command_pressure[-1])
        pressure = self._run_gromacs_container(
            command_pressure,
            volumes=self.volume
        ).decode("utf-8")

        command_density = [
            "sh",
            "-c",
            f"""
                echo 24 0 | \
                gmx energy -f npt.edr -o {self.results_folder.name}/density.xvg
            """
        ]

        logging.info(command_density[-1])
        density = self._run_gromacs_container(
            command_density,
            volumes=self.volume
        ).decode("utf-8")

        return f"{pressure}, \n, {density}"
    
    def generate_final_xtc_file(self):
        """
        Generates the final xvg file joining all the parts of the trajectory
        """
        steps = find_files_with_same_pattern(
            self.path, "my_job.output_*.tar.gz")

        xtc_files = []
        for step in steps:
            with tarfile.open(step, "r:gz") as tar:
                for member in tar.getmembers():
                    if member.name.endswith(".xtc"):
                        tar.extract(member, path=f"{self.path}")
                        xtc_files.append(member.name)

        command = [
            "sh",
            "-c",
            f"gmx trjcat -f {' '.join(xtc_files)} -o {self.results_folder.name}/{self.XTC_FILE}"
        ]

        logging.info(command[-1])
        _ = self._run_gromacs_container(
            command,
            volumes=self.volume
        )
    
    def fix_periodicity(self):
        """
        Function to fix the periodicity of the trajectory.
        """
        filename, _ = self._find_last_tpr_file()
        command = [
            "sh",
            "-c",
            f"""
                echo 1 0 | \
                gmx trjconv \
                    -s {filename} \
                    -f {self.results_folder.name}/{self.XTC_FILE} \
                    -o {self.results_folder.name}/{self.XTC_NO_PBC_FILE} \
                    -pbc mol -center
            """
        ]

        logging.info(command[-1])
        periodicity = self._run_gromacs_container(
            command,
            volumes=self.volume
        ).decode("utf-8")

        return periodicity
    
    def generate_video(self):
        """
        Generate a video from the trajectory
        """
        filename, _ = self._find_last_tpr_file()
        command = [
            "sh",
            "-c",
            f"""
                echo 20 | \
                gmx trjconv \
                    -s {filename} \
                    -f {self.results_folder.name}/{self.XTC_NO_PBC_FILE} \
                    -n index.ndx \
                    -o {self.results_folder.name}/reduced.xtc \
                    -dt 250
            """
        ]

        logging.info(command[-1])
        _ = self._run_gromacs_container(
            command,
            volumes=self.volume
        )

    def get_initial_configuration(self):
        """
        Retrieves the initial configuration of the system.

        Returns:
            None
        """
        filename, _ = self._find_last_tpr_file()

        command = [
            "sh",
            "-c",
            f"""
                echo 20 | \
                gmx trjconv \
                    -s {filename} \
                    -f {self.results_folder.name}/{self.XTC_NO_PBC_FILE} \
                    -n index.ndx \
                    -o {self.results_folder.name}/initial_conf.pdb \
                    -dump 0
            """
        ]

        logging.info(command[-1])

        _ = self._run_gromacs_container(
            command,
            volumes=self.volume
        )

    def generate_com_distance(self):
        """
        Generate the distance between the center of mass of the protein and the ligand
        """
        filename, _ = self._find_last_tpr_file()

        command = [
            "sh",
            "-c",
            f"""
                gmx distance \
                    -s {filename} \
                    -f {self.results_folder.name}/{self.XTC_NO_PBC_FILE} \
                    -oall {self.results_folder.name}/com_dist.xvg \
                    -n index.ndx \
                    -select 'com of group "UNL" plus com of group "Protein"'
            """
        ]

        logging.info(command[-1])
        com = self._run_gromacs_container(
            command,
            volumes=self.volume
        ).decode("utf-8")

        return com

    def generate_sasa_ligand(self):
        """
        Solvent Accessible Surface Area (SASA)
        """
        filename, _ = self._find_last_tpr_file()

        command = [
            "sh",
            "-c",
            f"""
                echo 13  0 | \
                gmx sasa \
                    -s {filename} \
                    -f {self.results_folder.name}/{self.XTC_NO_PBC_FILE} \
                    -surface UNL \
                    -o {self.results_folder.name}/sasa.xvg
            """
        ]

        logging.info(command[-1])
        sasa = self._run_gromacs_container(
            command,
            volumes=self.volume
        ).decode("utf-8")
        return sasa
    
    def generate_interaction_energy(self):
        """
        Generate the Interaction Energy
        """
        filename, _ = self._find_last_tpr_file()

        command = [
            "sh",
            "-c",
            f"""
                echo 20 21 0 | \
                gmx energy \
                    -f ie.edr \
                    -o {self.results_folder.name}/interaction_energy.xvg
            """
        ]

        logging.info(command[-1])
        coulombic_energy = self._run_gromacs_container(
            command,
            volumes=self.volume
        ).decode("utf-8")

        return coulombic_energy
    
    def generate_rmsd(self):
        """
        Generate the RMSD
        """
        u = mda.Universe(
            f"{self.path}/npt.gro",
            f"{self.path}/results/fixed.xtc",
            topology_format="GRO",
            trajectory_format="XTC"
        )

        protein = u.select_atoms('protein')
        ligand = u.select_atoms('resname UNL')

        logging.info(f"Total residues in the protein: {protein.residues.n_residues}")
        logging.info(f"Total atoms in the ligand: {ligand.n_atoms}")

        total_residues = protein.residues.n_residues
        selection = f"backbone and not (resid {total_residues-200} to {total_residues}) and not (resid 1 to 211)"
        rmsd_protein = rms.RMSD(
            protein,
            protein,
            select=selection,
            ref_frame=0
        )
        rmsd_protein.run()
        logging.info(f"RMSD Protein: {np.mean(rmsd_protein.rmsd[:, 2]):.2f} Å")

        rmsd_ligand = rms.RMSD(
            ligand,
            ligand,
            select="all",
            ref_frame=0
        )
        rmsd_ligand.run()
        logging.info(f"RMSD Ligand: {np.mean(rmsd_ligand.rmsd[:, 2]):.2f} Å")

        df = pd.DataFrame({
            "Time (ps)": rmsd_protein.rmsd[:, 1],
            "RMSD Protein (Å)": rmsd_protein.rmsd[:, 2],
            "RMSD Ligand (Å)": rmsd_ligand.rmsd[:, 2]
        })

        df.to_csv(f"{self.results_folder}/rmsd.csv", index=True)

        with open(f"{self.results_folder}/rmsd_summary.txt", "w") as file:
            file.write(df.describe().to_string())

        return df
    
    def generate_radius_gyration(self):
        """
        Generate radius of gyration        
        """
        filename, _ = self._find_last_tpr_file()
        command = [
            "sh",
            "-c",
            f"""
                echo 1 | \
                gmx gyrate \
                -s {filename} \
                -f {self.results_folder.name}/{self.XTC_NO_PBC_FILE} \
                -o {self.results_folder.name}/gyrate.xvg
            """
        ]

        logging.info(command[-1])
        gyr = self._run_gromacs_container(
            command,
            volumes=self.volume
        ).decode('utf-8')

        return gyr
    
    def generate_rmsf(self):
        """
        Calculate the RMSF
        """
        u = mda.Universe(
            f"{self.path}/npt.gro",
            f"{self.path}/results/fixed.xtc",
            topology_format="GRO",
            trajectory_format="XTC"
        )

        protein = u.select_atoms('protein')
        rmsf = rms.RMSF(
            protein,
            ref_frame=0
        ).run()

        resids = protein.residues.resids

        df = pd.DataFrame({
            "Residue ID": resids,
            "RMSF": rmsf.rmsf
        })

        df.to_csv(f"{self.results_folder}/rmsf.csv", index=True)

        return df

         

