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

        logging.info('Generating NVT data')
        temperature_data: str = gromacs.generate_nvt_data()
        logging.info(temperature_data)

        logging.info('Generating NPT data')
        npt_data = gromacs.generate_npt_data()
        logging.info(npt_data)

        if self.run_gromacs:
            logging.info('Generating final xtc file')
            gromacs.generate_final_xtc_file()

            periodicity = gromacs.fix_periodicity()
            logging.info(periodicity)

            gromacs.generate_video()

        logging.info('Generate initial configuration file')
        gromacs.get_initial_configuration()

        logging.info('Calculate COM between ligand and protein')
        com_dist = gromacs.generate_com_distance()
        logging.info(com_dist)

        logging.info('Generate Solvent Accessible Surface Area (SASA)')
        sasa_ligand = gromacs.generate_sasa_ligand()
        logging.info(sasa_ligand)

        logging.info('Generate Coulombic Interaction Energy')
        interaction_energy = gromacs.generate_interaction_energy()
        logging.info(interaction_energy)

        logging.info('Generate RMSD')
        rmsd = gromacs.generate_rmsd()
        logging.info(rmsd)

        logging.info('Generate Radius of Gyration')
        radius_of_gyration = gromacs.generate_radius_gyration()
        logging.info(radius_of_gyration)

        # logging.info('Generate RMSF')
        # rmsf = gromacs.generate_rmsf()
        # logging.info(rmsf)

    def generate_gromacs_plots(self):
        """Generate plots from Gromacs data"""
        plots: GromacsPlot = GromacsPlot(self.path, self.protein, self.ligand)
        
        logging.info('Generating energy minimization plot')
        plots.plot_energy_minimization()

        logging.info('Generating temperature plot')
        plots.plot_temperature()

        logging.info("Plotting pressure")
        plots.plot_pressure()

        logging.info("Plotting density")
        plots.plot_density()

        logging.info("Plotting COM distance")
        plots.plot_com_distance()

        logging.info("Plotting SASA")
        plots.plot_sasa_ligand()

        logging.info("Plotting Interaction Energy")
        plots.plot_interaction_energy()

        logging.info("Plotting RMSD")
        plots.plot_rmsd()

        logging.info("Plotting Radius of Gyration")
        plots.plot_radius_gyration()

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
