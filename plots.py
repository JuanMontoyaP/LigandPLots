"""
This file contains the calls to graph the gromacs data
"""
from pathlib import Path
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from helpers import create_folder, read_xvg_files

class GromacsPlot:
    """
    Class for plotting all the data
    """

    def __init__(self, path: Path, protein: str, ligand: str):
        self.path = path
        self.protein = protein
        self.ligand = ligand
        self.images_path = f"{self.path}/images/"

        create_folder(self.path, "images")

    def plot_data(self, x, y, titles):
        """
        Plot the specified data
        """
        fig, ax = plt.subplots()

        ax.plot(x, y, linewidth=0.75)

        # Add labels and title
        plt.xlabel(titles[0])
        plt.ylabel(titles[1])
        plt.title(titles[2])

        if len(titles) == 4:
            plt.legend(titles[-1])

        plt.xlim(left=x[0], right=x[-1])

        plt.minorticks_on()

        plt.tick_params(axis='both', which='both',
                        direction='in', right=True, top=True)

        return fig, ax
    
    def plot_energy_minimization(self):
        """
        Plot the energy minimization function
        """
        energy_data = read_xvg_files(f"{self.path}/results/potential.xvg")

        _, ax = self.plot_data(
            energy_data[:, 0]/1000,
            energy_data[:, 1],
            [
                'Time $(ns)$',
                'Potential Energy $(kJ/mol)$',
                f"{self.protein} - {self.ligand}"
            ]
        )

        ax.lines[0].set_color('black')

        plt.suptitle('Energy Minimization', fontsize=20, y=1)

        plt.savefig(
            self.images_path + "minimization.png",
            bbox_inches='tight',
            pad_inches=0.1,
            dpi=300
        )

    def plot_temperature(self):
        """
        Plot temperature in the NVT equilibration
        """
        temperature = read_xvg_files(f"{self.path}/results/temperature.xvg")

        _, ax = self.plot_data(
            temperature[:, 0],
            temperature[:, 1],
            [
                'Time $(ps)$',
                'Potential Energy $(K)$',
                f"{self.protein} - {self.ligand}, NVT Equilibration"
            ]
        )

        ax.lines[0].set_color('black')

        plt.ylim(
            temperature[:, 1].min() - 2,    
            temperature[:, 1].max() + 2
        )

        plt.suptitle('Temperature', fontsize=20, y=1)

        plt.savefig(
            self.images_path + "temperature.png",
            bbox_inches='tight',
            pad_inches=0.1,
            dpi=300
        )

    def plot_pressure(self):
        """
        Plot pressure in NPT equilibration
        """
        pressure = read_xvg_files(f"{self.path}/results/pressure.xvg")

        _, ax = self.plot_data(
            pressure[:, 0],
            pressure[:, 1],
            [
                'Time $(ps)$',
                'Pressure $(bar)$',
                f"{self.protein} - {self.ligand}, NPT Equilibration"
            ]
        )

        ax.lines[0].set_color('black')

        plt.ylim(
            pressure[:, 1].min() - 60,    
            pressure[:, 1].max() + 60
        )

        plt.suptitle('Pressure', fontsize=20, y=1)

        plt.savefig(
            self.images_path + "pressure.png",
            bbox_inches='tight',
            pad_inches=0.1,
            dpi=300
        )

    def plot_density(self):
        """
        Plot density in NPT equilibration
        """
        density = read_xvg_files(f"{self.path}/results/density.xvg")

        _, ax = self.plot_data(
            density[:, 0],
            density[:, 1],
            [
                'Time $(ps)$',
                'Density $(kg/m^{3})$',
                f"{self.protein}, NPT Equilibration"
            ]
        )

        ax.lines[0].set_color('black')

        plt.ylim(
            density[:, 1].min() - 10,    
            density[:, 1].max() + 10
        )

        plt.suptitle('Density', fontsize=20, y=1)

        plt.savefig(
            self.images_path + "density.png",
            bbox_inches='tight',
            pad_inches=0.1,
            dpi=300
        )

    def plot_com_distance(self):
        """
        Plot the center of mass distance between the protein and the ligand
        """
        com = read_xvg_files(f"{self.path}/results/com_dist.xvg")

        _, ax = self.plot_data(
            com[:, 0]/1000,
            com[:, 1],
            [
                'Time $(ns)$',
                'Distance $(nm)$',
                f"{self.protein} - {self.ligand}"
            ]
        )

        ax.lines[0].set_color('black')

        plt.ylim(
            com[:, 1].min() - 0.5,    
            com[:, 1].max() + 0.5
        )

        plt.suptitle('Distance between Centers of Mass', fontsize=20, y=1)

        plt.savefig(
            self.images_path + "com.png",
            bbox_inches='tight',
            pad_inches=0.1,
            dpi=300
        )

    def plot_sasa_ligand(self):
        """
        Plot the Solvent Accessible Surface Area of the ligand
        """
        sasa = read_xvg_files(f"{self.path}/results/sasa.xvg")

        _, ax = self.plot_data(
            sasa[:, 0]/1000,
            sasa[:, 1],
            [
                'Time $(ns)$',
                'SASA $(nm^{2})$',
                f"{self.ligand}"
            ]
        )

        ax.lines[0].set_color('black')

        plt.ylim(
            sasa[:, 1].min() - 0.5,    
            sasa[:, 1].max() + 0.5
        )

        plt.suptitle('Solvent Accessible Surface Area', fontsize=20, y=1)

        plt.savefig(
            self.images_path + "sasa.png",
            bbox_inches='tight',
            pad_inches=0.1,
            dpi=300
        )

    def plot_interaction_energy(self):
        """
        Plot the interaction energy
        """
        energy = np.loadtxt(f"{self.path}/results/interaction_energy.xvg", comments=['#', '@'])

        total_energy = energy[:, 1] + energy[:, 2]

        _ = self.plot_data(
            energy[:, 0]/1000,
            list(zip(energy[:, 1], energy[:, 2], total_energy)),
            [
                'Time $(ns)$',
                'Energy $(kJ/mol)$',
                f"{self.ligand}"
            ]
        )

        plt.legend(
            ['Coulombic', 'LJ', 'Total'],
            fancybox=True,
            fontsize='small'
        )

        col_min, col_max = np.min(energy[:, 1]), np.max(energy[:, 1])
        lj_min, lj_max = np.min(energy[:, 2]), np.max(energy[:, 2])
        total_min, total_max = np.min(total_energy), np.max(total_energy)

        plt.ylim(
            min(col_min, lj_min, total_min) - 25,    
            max(col_max, lj_max, total_max) + 25
        )

        plt.suptitle('Interaction energy', fontsize=20, y=1)

        plt.savefig(
            self.images_path + "inter_energy.png",
            bbox_inches='tight',
            pad_inches=0.1,
            dpi=300
        )


    def plot_rmsd(self):
        """
        Plot RMSD
        """
        rmsd = pd.read_csv(f"{self.path}/results/rmsd.csv")

        _ = self.plot_data(
            rmsd['Time (ps)'].values/1000,
            list(zip(
                rmsd['RMSD Protein (Å)'].values/10,
                rmsd['RMSD Ligand (Å)'].values/10
        )),
            [
                'Time $(ns)$',
                'RMSD $(nm)$',
                f"{self.protein}, {self.ligand}",
            ]
        )

        plt.legend(
            [
                self.protein,
                self.ligand
            ],
            fancybox=True,
            fontsize='small'
        )

        plt.ylim(0, (rmsd['RMSD Protein (Å)'].max()+2)/10)

        plt.suptitle('RMSD', fontsize=20, y=1)
        
        plt.tight_layout()
    
        plt.savefig(
            self.images_path + "rmsd.png",
            bbox_inches='tight',
            pad_inches=0.1,
            dpi=300
        )

    def plot_radius_gyration(self):
        """
        Plot radius of gyration
        """
        gyr = read_xvg_files(f"{self.path}/results/gyrate.xvg")

        _, ax = self.plot_data(
            gyr[:, 0]/1000,
            gyr[:, 1],
            [
                'Time $(ns)$',
                '$R_{g}$ $(nm)$',
                f"{self.protein}, Unrestrained MD"
            ]
        )

        ax.lines[0].set_color('black')

        plt.ylim(
            min(gyr[:, 1]) - 0.2,
            max(gyr[:, 1]) + 0.2
        )

        plt.suptitle('Radius of gyration', fontsize=20, y=1)

        plt.savefig(
            self.images_path + "gyration.png",
            bbox_inches='tight',
            pad_inches=0.1,
            dpi=300
        )
