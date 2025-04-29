import pandas as pd
import numpy as np
import os

def pre_process_data_set(file_name, nucleus, cytoplasm, resolution_z, resolution_xy):
    data_set = pd.read_csv(file_name)
    if nucleus:
        data_set = data_set[data_set.nuclei_mask > 0.0]
    if cytoplasm:
        data_set = data_set[data_set.nuclei_mask == 0]
    data_set['Z_det'] = data_set['Z_det'] * resolution_z
    data_set['Y_det'] = data_set['Y_det'] * resolution_xy
    data_set['X_det'] = data_set['X_det'] * resolution_xy

    return data_set

def extract_nucleus_data(file_path):
    with open(file_path, 'r') as f:
        lines = f.readlines()

    nucleus_data = []
    current_cell_name = None
    parsing_nucleus = False
    nucleus_section = {"X_POS": [], "Y_POS": [], "Cell_Name": []}

    for line in lines:
        line = line.strip()
        if line.startswith("CELL_START"):
            # Capture the cell name from the line
            current_cell_name = line.split("\t")[1]
        elif line.startswith("Nucleus_START"):
            # Start a new nucleus section
            parsing_nucleus = True
            nucleus_section = {"Cell_Name": [], "X_POS": [], "Y_POS": []}
        elif line.startswith("Nucleus_END"):
            # End the current nucleus section, append to data
            parsing_nucleus = False
            if nucleus_section["X_POS"] and nucleus_section["Y_POS"]:
                # Add the cell name to every row in the nucleus section
                nucleus_section["Cell_Name"] = [current_cell_name] * len(nucleus_section["X_POS"])
                nucleus_data.append(pd.DataFrame(nucleus_section))
        elif parsing_nucleus:
            if line.startswith("X_POS"):
                nucleus_section["X_POS"] = list(map(int, line.split("\t")[1:]))
            elif line.startswith("Y_POS"):
                # Collect Y_POS values until Z_POS is encountered
                nucleus_section["Y_POS"].extend(map(int, line.split("\t")[1:]))
            elif line.startswith("Z_POS"):
                # Stop collecting Y_POS when Z_POS is encountered
                pass

    # Combine all nucleus data into one DataFrame
    if nucleus_data:
        final_df = pd.concat(nucleus_data, ignore_index=True)
    else:
        final_df = pd.DataFrame(columns=["X_POS", "Y_POS", "Cell_Name"])

    return final_df

def extract_cytoplasm_data(file_path):
    with open(file_path, 'r') as f:
        lines = f.readlines()

    nucleus_data = []
    current_cell_name = None
    parsing_nucleus = False
    nucleus_section = {"X_POS": [], "Y_POS": [], "Cell_Name": []}

    for line in lines:
        line = line.strip()
        if line.startswith("CELL_START"):
            # Capture the cell name from the line
            current_cell_name = line.split("\t")[1]
            parsing_nucleus = True #new
            nucleus_section = {"Cell_Name": [], "X_POS": [], "Y_POS": []} #new
        elif line.startswith("CELL_END"):
            # End the current nucleus section, append to data
            parsing_nucleus = False
            if nucleus_section["X_POS"] and nucleus_section["Y_POS"]:
                # Add the cell name to every row in the nucleus section
                nucleus_section["Cell_Name"] = [current_cell_name] * len(nucleus_section["X_POS"])
                nucleus_data.append(pd.DataFrame(nucleus_section))
        elif parsing_nucleus:
            if line.startswith("X_POS"):
                nucleus_section["X_POS"] = list(map(int, line.split("\t")[1:]))
            elif line.startswith("Y_POS"):
                # Collect Y_POS values until Z_POS is encountered
                nucleus_section["Y_POS"].extend(map(int, line.split("\t")[1:]))
            elif line.startswith("Z_POS"):
                # Stop collecting Y_POS when Z_POS is encountered
                pass

    # Combine all nucleus data into one DataFrame
    if nucleus_data:
        final_df = pd.concat(nucleus_data, ignore_index=True)
    else:
        final_df = pd.DataFrame(columns=["X_POS", "Y_POS", "Cell_Name"])

    return final_df


def find_cell_name(sk_3dimg, final_df):
    z_indices, y_indices, x_indices = np.where(sk_3dimg > 0)
    intensities = sk_3dimg[z_indices, y_indices, x_indices]
    positions_df = pd.DataFrame({'Y_POS': y_indices, 'X_POS': x_indices, 'Intensity': intensities})
    merged_df = positions_df.merge(final_df, on=['Y_POS', 'X_POS'], how='inner')
    overlap_counts = merged_df.groupby(['Intensity', 'Cell_Name']).size().reset_index(name='Match_count')
    best_matches = overlap_counts.loc[overlap_counts.groupby('Intensity')['Match_count'].idxmax()]

    return best_matches[['Cell_Name', 'Intensity', 'Match_count']]



def data_prep(file_path_cellpose, file_path_outline_files):
    data = []
    #print(file_path_cellpose)
    for item in os.listdir(file_path_cellpose):
        parts = item.split('_')
        donor = parts[0][:2]  # First 2 characters, e.g., 'D1'
        timepoint = parts[1]  # Second part, e.g., '0H'
        type = 'C' if 'C' in parts[0] else 'HuRKO' #dit kan mis gaan als er meerdere inputs zijn.
        sample = parts[2] if len(parts) > 2 else None
        data.append({'ID': item, 'Donor': donor, 'TimePoint': timepoint, 'Type': type, 'Sample': sample})


    files_mask_df = pd.DataFrame(data)
    #print('files_mask_df', files_mask_df)
    data_2 = []
    for timepoint in os.listdir(file_path_outline_files):
        for item in os.listdir(f"{file_path_outline_files}/{timepoint}/FQ_outline"):
            #TimePoint = timepoint
            parts = item.split('_')
            if len(parts)>1: #dit moet ik even checken
                donor = parts[0][:2]  # First 2 characters, e.g., 'D1'
                type = 'C' if 'C' in parts[0] else 'HuRKO'  # dit kan mis gaan als er meerdere inputs zijn.
                sample = parts[1] if len(parts) > 2 else None #DIT GAAT NOG NIET GOED!
                data_2.append({'ID': item, 'Donor': donor, 'TimePoint': timepoint, 'Type': type, 'Sample': sample})
    files_outline_df = pd.DataFrame(data_2)
    #merge the two dataframes

    name_overlap_df = pd.merge(files_mask_df, files_outline_df, on=["Donor", "TimePoint", "Type", "Sample"], suffixes=('_cp', '_outline'))
    return name_overlap_df

def data_prep_cyto(file_path_cellpose, file_path_outline_files):
    data = []
    #print(file_path_cellpose)
    for folder in os.listdir(file_path_cellpose):
        for item in os.listdir(f'{file_path_cellpose}/{folder}'):
            parts = item.split('_')
            donor = parts[0][:2]  # First 2 characters, e.g., 'D1'
            timepoint = folder  # Second part, e.g., '0H'
            type = 'C' if 'C' in parts[0] else 'HuRKO' #dit kan mis gaan als er meerdere inputs zijn.
            sample = parts[1] if len(parts) > 2 else None
            data.append({'ID': item, 'Donor': donor, 'TimePoint': timepoint, 'Type': type, 'Sample': sample})


    files_mask_df = pd.DataFrame(data)
    #print('files_mask_df', files_mask_df)
    data_2 = []
    for timepoint in os.listdir(file_path_outline_files):
        for item in os.listdir(f"{file_path_outline_files}/{timepoint}/FQ_outline"):
            #TimePoint = timepoint
            parts = item.split('_')
            if len(parts)>1: #dit moet ik even checken
                donor = parts[0][:2]  # First 2 characters, e.g., 'D1'
                type = 'C' if 'C' in parts[0] else 'HuRKO'  # dit kan mis gaan als er meerdere inputs zijn.
                sample = parts[1] if len(parts) > 2 else None #DIT GAAT NOG NIET GOED!
                data_2.append({'ID': item, 'Donor': donor, 'TimePoint': timepoint, 'Type': type, 'Sample': sample})
    files_outline_df = pd.DataFrame(data_2)
    #merge the two dataframes

    name_overlap_df = pd.merge(files_mask_df, files_outline_df, on=["Donor", "TimePoint", "Type", "Sample"], suffixes=('_cp', '_outline'))
    return name_overlap_df

def write_to_excel(resolution_xy, resolution_z, eps, min_samples, repetition, nucleus, cytoplasm, output_filepath, final_data):

    input_params = [
        "Input Parameters:",
        "------------------",
        f"Resolution XY: {resolution_xy}",
        f"Resolution Z: {resolution_z}",
        f"eps: {eps}",
        f"min_samples: {min_samples}",
        f"repetition: {repetition}",
        f"nucleus: {nucleus}",
        f"cytoplasm: {cytoplasm}"
    ]

    with pd.ExcelWriter(output_filepath, engine='xlsxwriter') as writer:
        workbook = writer.book
        worksheet = workbook.add_worksheet("Sheet1")  # Manually add a worksheet
        writer.sheets["Sheet1"] = worksheet

        # Write input parameters at the top
        for i, line in enumerate(input_params):
            worksheet.write(i, 0, line)  # (row, col, value)

        # Write the DataFrame below the input parameters
        final_data.to_excel(writer, sheet_name="Sheet1", startrow=len(input_params) + 2, index=False)
    return
