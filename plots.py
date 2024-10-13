"""
This file contains the calls to graph the gromacs data
"""
from pathlib import Path
import matplotlib.pyplot as plt

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

        ax.plot(x, y, linewidth=0.7)

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