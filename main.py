"""
This is file contains orchestrates all the functions to generate Gromacs plots.
"""
#!/usr/bin/env python3 -u
import logging
from pathlib import Path
import click

from gromacs import GromacsData
from plots import GromacsPlot

logging.basicConfig(
    format='[%(asctime)s] - %(levelname)s - %(filename)s - %(funcName)s:%(lineno)d - %(message)s',
    level=logging.DEBUG
)

class LigandPlots:

    def __init__(self, path: Path, protein: str, ligand: str, run_gromacs: bool = False) -> None:
        self.path: Path = path
        self.protein: str = protein
        self.ligand: str = ligand
        self.run_gromacs: bool = run_gromacs

    def generate_gromacs_data(self):
        gromacs: GromacsData = GromacsData(self.path)
        
        logging.info('Generating minimization data')
        energy_minimization: str = gromacs.generate_minimization_data()
        logging.info(energy_minimization)

    def generate_gromacs_plots(self):
        """Generate plots from Gromacs data"""
        plots: GromacsPlot = GromacsPlot(self.path, self.protein, self.ligand)
        plots.plot_energy_minimization()

    def Run(self):
        logging.info('Generating Gromacs data')
        self.generate_gromacs_data()

        logging.info('Generating Gromacs plots')
        self.generate_gromacs_plots()

@click.command()
@click.option(
    '--path',
    type=click.Path(
        exists=True,
        file_okay=True,
        dir_okay=True,
        readable=True
    ),
    help="Path where the data files should be.",
    required=True
)
@click.option(
    '--protein',
    type=str,
    required=True,
    help="Name of the protein to be deployed"
)
@click.option(
    '--ligand',
    type=str,
    required=True,
    help="Name of the ligand to be deployed"
)
@click.option(
    '--run_gromacs',
    is_flag=True,
    help="Run gromacs commands"
)
def cli(
    path: str,
    protein: str,
    ligand: str,
    run_gromacs: bool = False
):
    logging.info(f"Running for protein {protein} and ligand {ligand}")
    ligand = LigandPlots(
        path=Path(path),
        protein=protein,
        ligand=ligand,
        run_gromacs=run_gromacs
    )
    ligand.Run()

if __name__ == '__main__':
    cli()
