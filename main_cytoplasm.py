import numpy as np
import pandas as pd
from sklearn.cluster import DBSCAN
from skimage import io
import re
from tqdm import tqdm
from pre_process_files import pre_process_data_set, data_prep, data_prep_cyto, extract_cytoplasm_data, extract_nucleus_data, find_cell_name, write_to_excel
from post_process_dbscan import process_dbscan
from plotting import plot_dbscan

def dbscan (eps, min_samples, df) :
    #code to perform dbscan and retract number of clusters
    dbscan = DBSCAN(eps=eps, min_samples = min_samples)
    labels = dbscan.fit_predict(df)
    arr = np.unique(labels)
    N_cluster = len(arr) - 1 if -1 in arr else len(arr)
    return N_cluster

def main(file_path_cellpose_n, file_path_cellpose_c, file_path_outline_files, eps, min_samples, resolution_z, resolution_xy, repetition, cytokine_analysis, nucleus, cytoplasm):
    data_spots = pre_process_data_set(file_name, nucleus, cytoplasm, resolution_z, resolution_xy, cytokine_analysis)

    #find the nucleus and the cluster dataset
    name_overlap_n = data_prep(file_path_cellpose_n, file_path_outline_files)
    name_overlap_c = data_prep_cyto(file_path_cellpose_c, file_path_outline_files)
    name_overlap_df = pd.merge(name_overlap_n, name_overlap_c, on=["Donor", "TimePoint", "Type", "Sample"],
                               suffixes=('_n', '_c'))


    ### Run Code
    final_data = pd.DataFrame()
    error_data = pd.DataFrame()

    for i in tqdm(range(len(name_overlap_df))):
        try:
            timepoint = name_overlap_df.loc[i, 'TimePoint']
            name_cytoplasm = name_overlap_df.loc[i, 'ID_cp_c']
            name_nucleus = name_overlap_df.loc[i, 'ID_cp_n']

            outline_file_path = f"{file_path_outline_files}/{timepoint}/FQ_outline/{name_overlap_df.loc[i, 'ID_outline_n']}"
            file_path_tiff_n = f"{file_path_cellpose_n}/{name_nucleus}"
            file_path_tiff_c = f"{file_path_cellpose_c}/{timepoint}/{name_cytoplasm}"

            #open images
            segmented_array_cytoplasm = io.imread(file_path_tiff_c)
            segmented_array_nucleus = io.imread(file_path_tiff_n)

            #find outline file path
            cytoplasm_df = extract_cytoplasm_data(outline_file_path)
            nucleus_df = extract_nucleus_data(outline_file_path)

            #match intensity and cell together
            cell_name_cytoplasm_df = find_cell_name(segmented_array_cytoplasm, cytoplasm_df)
            cell_name_nucleus_df = find_cell_name(segmented_array_nucleus, nucleus_df)
            cell_name_cytoplasm_df = cell_name_cytoplasm_df.groupby('Cell_Name').agg({
                'Intensity': lambda x: ', '.join(map(str, sorted(x)))}).reset_index()
            cell_name_nucleus_df = cell_name_nucleus_df.groupby('Cell_Name').agg({
                'Intensity': lambda x: ', '.join(map(str, sorted(x)))}).reset_index()

            merged = pd.merge(cell_name_cytoplasm_df,cell_name_nucleus_df, on=['Cell_Name'], how = 'inner', suffixes = ['_cytoplasm', '_nucleus'])
            part_name = re.split(r'[_,\s]+', name_nucleus)
            merged['ID'] = merged['Cell_Name'].astype(str).apply(lambda x: f"{part_name[0]}_{part_name[1]}__{part_name[2]}_{x}")

            #find the stacks that are cutted and use this to find the cytoplasm volume
            first_stack_nucleus = int(re.split(r'[_,\s]+', name_nucleus)[4])
            last_stack_nucleus = int(re.split(r'[_, \s]+', name_nucleus)[5])
            first_stack_cytoplasm = int(re.split(r'[_,\s]+', name_cytoplasm)[4])
            last_stack_cytoplasm = int(re.split(r'[_, \s]+', name_cytoplasm)[5])

            stack_difference = first_stack_nucleus-first_stack_cytoplasm

            first_stack_nucleus -= first_stack_cytoplasm
            last_stack_nucleus -= first_stack_cytoplasm
            last_stack_cytoplasm -= first_stack_cytoplasm
            first_stack_cytoplasm -= first_stack_cytoplasm

            nucleus_volume = np.argwhere(segmented_array_nucleus > 0)
            nucleus_volume[:,0] = nucleus_volume[:, 0]+stack_difference

            if first_stack_cytoplasm > first_stack_nucleus and last_stack_cytoplasm < last_stack_nucleus:
                nucleus_volume = nucleus_volume[(nucleus_volume[:, 0] <= last_stack_cytoplasm) & (nucleus_volume[:, 0] >= 0)]

            elif first_stack_cytoplasm > first_stack_nucleus:
                nucleus_volume = nucleus_volume[nucleus_volume[:,0] >= 0]

            elif last_stack_cytoplasm < last_stack_nucleus:
                nucleus_volume = nucleus_volume[nucleus_volume[:, 0] <= last_stack_cytoplasm]

            #make the nucleus volume 0, so background. After this it is the real cytoplasm volume
            segmented_array_cytoplasm[nucleus_volume[:,0], nucleus_volume[:,1], nucleus_volume[:,2]] = 0

            #loop to find simulated clusters using DBSCAN
            final_df = pd.merge(data_spots, merged, on=['ID'], how = 'inner')
            for id in final_df['ID'].unique():
                if len(final_df[final_df['ID'] == id]) >= min_samples:
                    #Select data
                    data = final_df[final_df.ID == id]
                    df = data.iloc[:, [3, 4, 5]].values

                    #Find real number of clusters
                    N_cluster_real = dbscan(eps, min_samples, df)

                    #Find intensities, number of particles and the cytokine type(s)
                    data_df = merged[merged['ID'] == id]
                    intensities = [int(num.strip()) for num in data_df["Intensity_cytoplasm"].iloc[0].split(',')]
                    number_of_particles = len(data)
                    cytokine = data['cytokine'].unique()

                    #Find the real volume
                    positions_cytoplasm = np.argwhere(np.isin(segmented_array_cytoplasm, intensities))
                    volume = len(positions_cytoplasm)

                    try:
                        N_cluster_simulated = []
                        rng = np.random.default_rng()

                        #Randomly make repetition x random indices from the number of particle
                        random_indices = np.stack([rng.choice(positions_cytoplasm.shape[0], size=number_of_particles,  replace=False) for _ in range(repetition)])

                        #Select and alter distance based on resolution
                        random_positions = positions_cytoplasm[random_indices]
                        random_positions[:, :, 0] = random_positions[:, :, 0]*resolution_z
                        random_positions[:, :, 1] = random_positions[:, :, 1]*resolution_xy
                        random_positions[:, :, 2] = random_positions[:, :, 2]*resolution_xy

                        #Repetition over all the random sets and perform DBSCAN (THIS IS A SLOW STEP)
                        for j in range(repetition):
                            N_cluster = dbscan(eps, min_samples, random_positions[j])
                            N_cluster_simulated.append(N_cluster)

                        #Calculate P-value and save
                        pvalue = 1-((sum((n > N_cluster_real) + (n < N_cluster_real) for n in N_cluster_simulated) + 1) / (repetition + 1))
                        new_data_row = {'ID': id, 'Cytokine': cytokine, 'Volume': volume, 'Production': number_of_particles,
                                                'N_cluster_real': N_cluster_real, 'N_cluster_simulated': np.mean(N_cluster_simulated),
                                                'P_value': pvalue, 'Intensities': intensities}
                        final_data = pd.concat([final_data, pd.DataFrame([new_data_row])], ignore_index=True)

                    except:
                        #IF THIS ERROR IS PRINTED, THE CELL VOLUME IS TOO SMALL TO MAKE NUMBER_OF_PARTICLES TIMES REPETITION SETS OF COORDINATES
                        print('For this cell there is not enough data')
                        new_error_row = {'ID': id, 'Cytokine': cytokine, 'Volume': volume, 'Production': number_of_particles}
                        error_data = pd.concat([error_data, pd.DataFrame([new_error_row])],
                                                                   ignore_index=True)

        except:
            #IF THIS ERROR IS PRINTED SOMETHING WENT WRONG WITH THE IMAGE
            print('Error in code, for image:', name_overlap_df.loc[i, 'ID_cp_n'])
            print('Code is continued without this image')

    return final_data, error_data

if __name__ == '__main__':
    #input parameters
    eps = 600
    min_samples = 3
    resolution_z = 200
    resolution_xy = 65
    repetition = 3000
    nucleus = False
    cytoplasm = True
    p_value = 0.05

    #Fill in BOTH if interested in both TNF and IFN, or fill in 'TNF' or 'IFN' when interested in a single analysis
    cytokine_analysis ='TNF'

    #input files
    # file_name = "F:/Eva/TryOuts/Spots_with_localization.csv"
    # file_path_cellpose_n = "F:/Eva/HuRKO_mask/Cellpose_DONE"
    # file_path_cellpose_c = "f:/Eva/CY5_HuRKO/CellPose2"
    # file_path_outline_files = "F:/Eva/Outline_files"
    # output_filepath = f"C:/Users/Eva Heemskerk/Documents/Delft/TU Delft/MEP/4. Results/Control_clustering_code/Cytoplasm_cluster_control_data_minsamp3.xlsx"

    #filepaths ZAN
    file_name = "F:/Eva/Data_Zan/Spots_with_localization.csv"
    file_path_cellpose_n = "F:/Eva/Data_Zan/Zan_nucleus_data"
    file_path_cellpose_c = "F:/Eva/Data_Zan/Cytoplasm_images/CellPose"
    file_path_outline_files = "F:/Eva/Data_Zan/Outline_files"
    output_filepath = f'F:/Eva/Data_Zan/Results/RDI/CC_cytoplasm_{cytokine_analysis}.xlsx'


    final_data, error_data = main(file_path_cellpose_n, file_path_cellpose_c, file_path_outline_files, eps, min_samples, resolution_z, resolution_xy, repetition, cytokine_analysis, nucleus, cytoplasm)
    write_to_excel(resolution_xy, resolution_z, eps, min_samples, repetition, nucleus, cytoplasm, output_filepath, final_data)
    #error_data.to_excel("error_data_cytoplasm.xlsx")

    particles_data = process_dbscan(output_filepath, file_name, resolution_xy, resolution_z, eps, min_samples,
                                    nucleus, cytoplasm, p_value)

    #plot data
    plot_dbscan(file_name, particles_data, nucleus, cytoplasm)
