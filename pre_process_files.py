import pandas as pd
import numpy as np
import os

def pre_process_data_set(file_name, nucleus, cytoplasm, resolution_z, resolution_xy, cytokine_analysis):
    data_set = pd.read_csv(file_name)
    if nucleus:
        data_set = data_set[data_set.nuclei_mask > 0.0]
    if cytoplasm:
        data_set = data_set[data_set.nuclei_mask == 0]

    if cytokine_analysis == 'TNF':
        data_set = data_set[data_set['cytokine'] == 'TNF']
    if cytokine_analysis == 'IFN':
        data_set = data_set[data_set['cytokine'] == 'IFN']

    data_set['Z_det'] = data_set['Z_det'] * resolution_z
    data_set['Y_det'] = data_set['Y_det'] * resolution_xy
    data_set['X_det'] = data_set['X_det'] * resolution_xy

    return data_set

def extract_nucleus_data(file_path):
    #this code rewrite the XYZ positions from the outline_files of the NUCLEUS to an usable dataframe
    with open(file_path, 'r') as f:
        lines = f.readlines()

    nucleus_data = []
    current_cell_name = None
    parsing_nucleus = False
    nucleus_section = {"X_POS": [], "Y_POS": [], "Cell_Name": []}

    for line in lines:
        line = line.strip()
        if line.startswith("CELL_START"):
            current_cell_name = line.split("\t")[1]
        elif line.startswith("Nucleus_START"):
            parsing_nucleus = True
            nucleus_section = {"Cell_Name": [], "X_POS": [], "Y_POS": []}
        elif line.startswith("Nucleus_END"):
            parsing_nucleus = False
            if nucleus_section["X_POS"] and nucleus_section["Y_POS"]:
                nucleus_section["Cell_Name"] = [current_cell_name] * len(nucleus_section["X_POS"])
                nucleus_data.append(pd.DataFrame(nucleus_section))
        elif parsing_nucleus:
            if line.startswith("X_POS"):
                nucleus_section["X_POS"] = list(map(int, line.split("\t")[1:]))
            elif line.startswith("Y_POS"):
                nucleus_section["Y_POS"].extend(map(int, line.split("\t")[1:]))
            elif line.startswith("Z_POS"):
                pass

    if nucleus_data:
        final_df = pd.concat(nucleus_data, ignore_index=True)
    else:
        final_df = pd.DataFrame(columns=["X_POS", "Y_POS", "Cell_Name"])

    return final_df

def extract_cytoplasm_data(file_path):
    #this code rewrite the XYZ positions from the outline_files of the CYTOPLASM to an usable dataframe
    with open(file_path, 'r') as f:
        lines = f.readlines()

    nucleus_data = []
    current_cell_name = None
    parsing_nucleus = False
    nucleus_section = {"X_POS": [], "Y_POS": [], "Cell_Name": []}

    for line in lines:
        line = line.strip()
        if line.startswith("CELL_START"):
            current_cell_name = line.split("\t")[1]
            parsing_nucleus = True #new
            nucleus_section = {"Cell_Name": [], "X_POS": [], "Y_POS": []} #new
        elif line.startswith("CELL_END"):
            parsing_nucleus = False
            if nucleus_section["X_POS"] and nucleus_section["Y_POS"]:
                nucleus_section["Cell_Name"] = [current_cell_name] * len(nucleus_section["X_POS"])
                nucleus_data.append(pd.DataFrame(nucleus_section))
        elif parsing_nucleus:
            if line.startswith("X_POS"):
                nucleus_section["X_POS"] = list(map(int, line.split("\t")[1:]))
            elif line.startswith("Y_POS"):
                nucleus_section["Y_POS"].extend(map(int, line.split("\t")[1:]))
            elif line.startswith("Z_POS"):
                pass

    if nucleus_data:
        final_df = pd.concat(nucleus_data, ignore_index=True)
    else:
        final_df = pd.DataFrame(columns=["X_POS", "Y_POS", "Cell_Name"])

    return final_df


def find_cell_name(sk_3dimg, final_df):
    #This code matches the intensities from the cellpose mask to the outline coordinates from the outline_files

    #select where there are cells and find intensities
    z_indices, y_indices, x_indices = np.where(sk_3dimg > 0)
    intensities = sk_3dimg[z_indices, y_indices, x_indices]

    #make a datafram
    positions_df = pd.DataFrame({'Y_POS': y_indices, 'X_POS': x_indices, 'Intensity': intensities})

    #merge the dataframe from the outline with the intensity coordinates
    merged_df = positions_df.merge(final_df, on=['Y_POS', 'X_POS'], how='inner')

    #group the cells and intensities, match_count is the number of overlaps between the two in pixels
    overlap_counts = merged_df.groupby(['Intensity', 'Cell_Name']).size().reset_index(name='Match_count')

    #find the maximal overlap count
    best_matches = overlap_counts.loc[overlap_counts.groupby('Intensity')['Match_count'].idxmax()]

    return best_matches[['Cell_Name', 'Intensity', 'Match_count']]



def data_prep(file_path_cellpose, file_path_outline_files):
    #this code is made to find the overlap between the nucleus cellpose file and the outline file path (only the names thus!!!)
    #if files are differently saved than explained in the GITHUB, change this code (probably change the for loop structure)

    #cellpose
    data = []
    for item in os.listdir(file_path_cellpose):
        parts = item.split('_')
        if len(parts) > 1:
            donor = parts[0][:2]  #IF DONOR>9 THAN THIS NEEDS TO BE CHANGED
            timepoint = parts[1]
            type = parts[0][2:]
            # sample = parts[2] if len(parts) > 2 else None #This was for the dataset of Valeria
            sample = parts[3]  # This is for dataset Zan
            data.append({'ID': item, 'Donor': donor, 'TimePoint': timepoint, 'Type': type, 'Sample': sample})
    files_mask_df = pd.DataFrame(data)

    #outline files
    data_2 = []
    for timepoint in os.listdir(file_path_outline_files):
        for item in os.listdir(f"{file_path_outline_files}/{timepoint}/FQ_outline"):
            parts = item.split('_')
            if len(parts)>1:
                donor = parts[0][:2] #IF DONOR>9 THAN THIS NEEDS TO BE CHANGED
                type = parts[0][2:]
                # sample = parts[1] if len(parts) > 2 else None #This was for the dataset of Valeria
                sample = parts[3]  # This is for data set Zan
                data_2.append({'ID': item, 'Donor': donor, 'TimePoint': timepoint, 'Type': type, 'Sample': sample})
    files_outline_df = pd.DataFrame(data_2)

    #two are merged:
    name_overlap_df = pd.merge(files_mask_df, files_outline_df, on=["Donor", "TimePoint", "Type", "Sample"], suffixes=('_cp', '_outline'))

    return name_overlap_df

def data_prep_cyto(file_path_cellpose, file_path_outline_files):
    #this code is made to find the overlap between the cytoplasm cellpose file and the outline file path (only the names thus!!!)
    #if files are differently saved than explained in the GITHUB, change this code (probably change the for loop structure)

    #cellpose
    data = []
    for folder in os.listdir(file_path_cellpose):
        for item in os.listdir(f'{file_path_cellpose}/{folder}'):
            parts = item.split('_')
            if len(parts) > 1:
                donor = parts[0][:2] #IF DONOR>9 THAN THIS NEEDS TO BE CHANGED
                timepoint = folder
                type = parts[0][2:]
                # sample = parts[1] if len(parts) > 2 else None #This was for the dataset of Valeria
                sample = parts[3]  # This is for dataset Zan
                data.append({'ID': item, 'Donor': donor, 'TimePoint': timepoint, 'Type': type, 'Sample': sample})
    files_mask_df = pd.DataFrame(data)

    #outline_files
    data_2 = []
    for timepoint in os.listdir(file_path_outline_files):
        for item in os.listdir(f"{file_path_outline_files}/{timepoint}/FQ_outline"):
            parts = item.split('_')
            if len(parts)>1:
                donor = parts[0][:2]
                type = parts[0][2:]
                # sample = parts[1] if len(parts) > 2 else None #This was for the dataset of Valeria
                sample = parts[3]  # This is for data set Zan
                data_2.append({'ID': item, 'Donor': donor, 'TimePoint': timepoint, 'Type': type, 'Sample': sample})
    files_outline_df = pd.DataFrame(data_2)

    # merge the two dataframes
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
        worksheet = workbook.add_worksheet("Sheet1")
        writer.sheets["Sheet1"] = worksheet
        worksheet_2 = workbook.add_worksheet("Sheet2")
        writer.sheets["Sheet2"] = worksheet_2

        for i, line in enumerate(input_params):
            worksheet_2.write(i, 0, line)

        final_data.to_excel(writer, sheet_name="Sheet1", startrow=len(input_params) + 2, index=False)
    return
