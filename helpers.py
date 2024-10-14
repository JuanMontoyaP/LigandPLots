"""
Helper functions
"""
from pathlib import Path
import numpy as np

def create_folder(path: str, folder_name: str = "results") -> None:
    """
    Create a new folder at the specified path.

    Args:
        path (str): The path where the new folder should be created.
        folder_name (str, optional): The name of the new folder. Defaults to "results".

    Returns:
        None
    """
    # Specify the path of the new folder
    folder_path = Path(f'{path}/{folder_name}')

    # Use the mkdir method to create the folder
    folder_path.mkdir(parents=True, exist_ok=True)

def read_xvg_files(file_path: str):
    """
    Read data from xvg files.
    """
    data = []

    with open(file_path, 'r', encoding='utf8') as file:
        for line in file:

            if line.startswith('#') or line.startswith('@'):
                continue

            values = [float(val) for val in line.strip().split()]
            data.append(values)

    data_array = np.array(data, dtype=float)
    return data_array

def find_files_with_same_pattern(path: Path, pattern: str) -> list[str]:
    """
    Find files in the given path that match the specified pattern.

    Args:
        path (Path): The path to search for files.
        pattern (str): The pattern to match against file names.

    Returns:
        list[str]: A list of file paths that match the specified pattern.
    """
    search_directory = Path(path)
    search_pattern = sorted(search_directory.glob(pattern))

    return list(search_pattern)
