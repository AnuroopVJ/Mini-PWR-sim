import pandas as pd
import json
import numpy as np

def extract_and_save_to_json(input_file, output_file):
    """
    Extracts 'DATA (B) 0.1' and 'EN (EV) 1.1' from the CSV file and saves them to a JSON file.

    Parameters:
        input_file (str): Path to the input CSV file.
        output_file (str): Path to the output JSON file.
    """
    try:
        # Load the CSV file
        df = pd.read_csv(input_file)

        # Extract relevant columns
        extracted_df = df[['DATA (B) 0.1', 'EN (EV) 1.1']]

        # Rename columns for better JSON compatibility
        extracted_df.columns = ['data_b', 'energy_ev']

        # Convert to dictionary
        data_dict = extracted_df.to_dict(orient='records')

        # Save to JSON file
        with open(output_file, 'w') as json_file:
            json.dump(data_dict, json_file, indent=4)

        print(f"Data successfully saved to {output_file}")

    except Exception as e:
        print(f"Error during extraction: {e}")


def load_data_from_json(json_file):
    """
    Loads data from the specified JSON file.

    Parameters:
        json_file (str): Path to the JSON file.

    Returns:
        list: List of dictionaries containing 'data_b' and 'energy_ev' pairs.
    """
    try:
        with open(json_file, 'r') as file:
            data = json.load(file)
        return data
    except Exception as e:
        print(f"Error while loading data: {e}")
        return []


def linear_interpolate(data, target_energy):
    """
    Interpolates data to estimate 'data_b' at the specified 'energy_ev'.

    Parameters:
        data (list): List of data points, each containing 'energy_ev' and 'data_b'.
        target_energy (float): The neutron energy in EV units for which we want to interpolate.

    Returns:
        float: Interpolated 'data_b' value.
    """
    try:
        # Extract arrays
        energies = np.array([entry['energy_ev'] for entry in data])
        data_b_values = np.array([entry['data_b'] for entry in data])

        # Interpolation
        interpolated_value = np.interp(target_energy, energies, data_b_values)
        return interpolated_value

    except Exception as e:
        print(f"Error during interpolation: {e}")
        return None


def load_cross_section_data_from_json(json_file, mt_key,cs_type='cross_sections'):
    """
    Loads cross section data for a given MT key from neutron_data.json.

    Parameters:
        json_file (str): Path to the JSON file.
        mt_key (str): The MT key, e.g., "MT_18".

    Returns:
        dict: The cross section data for the specified MT key.
    """
    try:
        with open(json_file, 'r') as file:
            data = json.load(file)
        return data[cs_type][mt_key]
    except KeyError:
        print(f"MT key '{mt_key}' not found in the JSON file.")
        return None
    except FileNotFoundError:
        print(f"File '{json_file}' not found.")
        return None
    except Exception as e:
        print(f"Error loading cross section data: {e}")
        return None


def get_energy_xs_arrays(cross_section_data):
    """
    Extracts energy and cross section arrays from the cross section data.

    Parameters:
        cross_section_data (dict): The cross section data for a specific MT.

    Returns:
        tuple: (energies, xs) as numpy arrays.
    """
    # Handle classic xstable format
    if cross_section_data is None:
        return [], []
    if "xstable" in cross_section_data:
        arr = cross_section_data["xstable"]
        # If xstable is a dict with 'E' and 'xs' keys (your JSON format)
        if isinstance(arr, dict) and "E" in arr and "xs" in arr:
            energies = arr["E"]
            xs = arr["xs"]
            if isinstance(energies, list) and isinstance(xs, list) and len(energies) == len(xs) and len(xs) > 0:
                return energies, xs
            else:
                print("Warning: 'xstable' dict found but 'E' and 'xs' are not valid or xs is empty.")
                return [], []
        # If xstable is a list of [energy, xs] pairs (classic format)
        elif isinstance(arr, list) and all(isinstance(row, (list, tuple)) and len(row) >= 2 for row in arr):
            energies = [row[0] for row in arr]
            xs = [row[1] for row in arr]
            return energies, xs
        else:
            print("Warning: 'xstable' is not in the expected format.")
            return [], []

    # Handle ENDF-like subsection format
    if "subsection" in cross_section_data:
        for subsec in cross_section_data["subsection"].values():
            # Try E/Eint and yi arrays
            # E as dict (keys are str numbers)
            if "E" in subsec and "yi" in subsec:
                E = subsec["E"]
                energies = [E[str(i)] for i in range(1, len(E) + 1)]
                xs = subsec["yi"]
                return energies, xs
            # Eint as list
            if "yields" in subsec and "Eint" in subsec["yields"] and "yi" in subsec["yields"]:
                energies = subsec["yields"]["Eint"]
                xs = subsec["yields"]["yi"]
                return energies, xs
            # Direct Eint/yi at subsection level
            if "Eint" in subsec and "yi" in subsec:
                energies = subsec["Eint"]
                xs = subsec["yi"]
                return energies, xs
    # Not found
    return [], []


def linear_interpolate_xs(energies, xs, target_energy):
    """
    Interpolates the cross section at the specified neutron energy.

    Parameters:
        energies (np.ndarray): Array of energies.
        xs (np.ndarray): Array of cross section values.
        target_energy (float): The neutron energy in eV.

    Returns:
        float: Interpolated cross section value.
    """
    try:
        if len(energies) == 0 or len(xs) == 0:
            return None
        return float(np.interp(target_energy, energies, xs))
    except Exception as e:
        print(f"Error during interpolation: {e}")
        return None


# Example Usage
if __name__ == "__main__":

    json_file = "neutron_data.json"
    mt_key = "MT_102"
    cross_section_data = load_cross_section_data_from_json(json_file, mt_key, cs_type='product_distributions')


    energies, xs = get_energy_xs_arrays(cross_section_data)
    target_energy = 1e-05  # 1 MeV in eV
    interpolated_xs = linear_interpolate_xs(energies, xs, target_energy)
    print(f"Interpolated cross section at {target_energy} eV: {interpolated_xs}")
