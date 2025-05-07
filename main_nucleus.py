import numpy as np
import pandas as pd
from skimage import io
from sklearn.cluster import DBSCAN
import os
from tqdm import tqdm
from pre_process_files import pre_process_data_set, data_prep, extract_nucleus_data, find_cell_name, write_to_excel
from post_process_dbscan import process_dbscan
from plotting import plot_dbscan

def dbscan (eps, min_samples, df):
    #code to perform dbscan and retract number of clusters
    dbscan = DBSCAN(eps=eps, min_samples = min_samples)
    labels = dbscan.fit_predict(df)
    arr = np.unique(labels)
    N_cluster = len(arr) - 1 if -1 in arr else len(arr)
    return N_cluster

def real_vs_simulated_clusters(name, cell_name_df, data_set, eps, min_samples, sk_3dimg, repetition, final_data):
    for cell_name in cell_name_df['Cell_Name'].unique():
        #Find data, intensities, number of particles and the cytokine type(s)
        intensity = cell_name_df[cell_name_df.Cell_Name == cell_name]['Intensity'].values
        pre_name = '_'.join(name.split('_')[:4])
        NAME = pre_name + '_' + cell_name
        data = data_set[data_set.ID == NAME]
        number_of_particles = len(data)
        cytokine = data['cytokine'].unique()

        #Find real number of clusters
        if number_of_particles >= min_samples:
            df = data.iloc[:, [3, 4, 5]].values
            N_cluster_real = dbscan(eps, min_samples, df)

            # Find the real volume
            positions = np.argwhere(np.isin(sk_3dimg, intensity))
            volume = len(positions)

            try:
                N_cluster_simulated = []
                rng = np.random.default_rng()

                #Randomly make repetition x random indices from the number of particle
                random_indices = np.stack([rng.choice(positions.shape[0], size=number_of_particles,  replace=False) for _ in range(repetition)])

                #Select and alter distance based on resolution
                random_positions = positions[random_indices]
                random_positions[:, :, 0] = random_positions[:, :, 0] * resolution_z
                random_positions[:, :, 1] = random_positions[:, :, 1] * resolution_xy
                random_positions[:, :, 2] = random_positions[:, :, 2] * resolution_xy

                #Repetition over all the random sets and perform DBSCAN (THIS IS A SLOW STEP)
                for j in range(repetition):
                    N_cluster = dbscan(eps, min_samples, random_positions[j])
                    N_cluster_simulated.append(N_cluster)

                #Calculate P-value and save
                pvalue = 1 - ((sum((n > N_cluster_real) + (n < N_cluster_real) for n in N_cluster_simulated) + 1) / (
                            repetition + 1))
                new_data_row = {'ID': id, 'Cytokine': cytokine, 'Volume': volume, 'Production': number_of_particles,
                                'N_cluster_real': N_cluster_real, 'N_cluster_simulated': np.mean(N_cluster_simulated),
                                'P_value': pvalue, 'Intensities': intensity}
                final_data = pd.concat([final_data, pd.DataFrame([new_data_row])], ignore_index=True)

            except:
                #IF THIS ERROR IS PRINTED, THE CELL VOLUME IS TOO SMALL TO MAKE NUMBER_OF_PARTICLES TIMES REPETITION SETS OF COORDINATES
                print('For this cell there is not enough data')

    return final_data

def main(overlap_df, file_path_outline_files, file_path_cellpose, data_set, eps, min_samples, repetition):
    final_data = pd.DataFrame()
    for i in tqdm(range(len(overlap_df))):
        try:
            outline_file = overlap_df['ID_outline'].iloc[i]
            timepoint = overlap_df['TimePoint'].iloc[i]

            #find nucleus coordinared
            file_path = f"{file_path_outline_files}/{timepoint}/FQ_outline/{outline_file}"
            nucleus_df = extract_nucleus_data(file_path)

            #CELLPOSE MASK, RETURNS A DATAFRAMA THAT LINKS THE CELL NAMES TO THE RIGHT INTENSITIES (MORE THAN ONE INTENSITY!)
            tif_file = overlap_df['ID_cp'].iloc[i]
            sk_3dimg = io.imread(f"{file_path_cellpose}/{tif_file}")
            cell_name_df = find_cell_name(sk_3dimg, nucleus_df)

            #Find real versus the simulated number of clusters
            final_data = real_vs_simulated_clusters(tif_file, cell_name_df, data_set, eps, min_samples, sk_3dimg, repetition, final_data)

        except:
            #IF THIS ERROR IS PRINTED SOMETHING WENT WRONG WITH THE IMAGE
            print('For this cell there is not enough data', overlap_df['ID_cp'].iloc[i])
            print('Code is continued without this image')

    return final_data



if __name__ == '__main__':
    # parameters:
    resolution_xy = 65
    resolution_z = 200
    eps = 600
    min_samples = 3
    repetition = 3000
    nucleus = True
    cytoplasm = False
    p_value = 0.05

    #Fill in BOTH if interested in both TNF and IFN, or fill in 'TNF' or 'IFN' when interested in a single analysis
    cytokine_analysis ='BOTH'

    #filepaths val:
    # file_name = 'C:/Users/Eva Heemskerk/Documents/Delft/TU Delft/MEP/4. Results/DATA/Spots_with_localization.csv'
    # file_path_cellpose = "F:/Eva/HuRKO_mask/Cellpose_DONE"
    # file_path_outline_files = "F:/Eva/Outline_files"
    # output_filepath = f"C:/Users/Eva Heemskerk/Documents/Delft/TU Delft/MEP/4. Results/Control_clustering_code/Nucleus_cluster_control_data_IFN_poging2.xlsx"

    #filepaths zan:
    file_name = "F:/Eva/Data_Zan/Spots_with_localization.csv"
    file_path_cellpose = "F:/Eva/Data_Zan/Zan_nucleus_data"
    file_path_outline_files = "F:/Eva/Data_Zan/Outline_files"
    output_filepath = f'F:/Eva/Data_Zan/Results/Clusters/CC_nucleus_{cytokine_analysis}.xlsx'


    data_set = pre_process_data_set(file_name, nucleus, cytoplasm, resolution_z, resolution_xy, cytokine_analysis)
    overlap_df = data_prep(file_path_cellpose, file_path_outline_files)
    final_data = main(overlap_df, file_path_outline_files, file_path_cellpose, data_set, eps, min_samples, repetition)

    write_to_excel(resolution_xy, resolution_z, eps, min_samples, repetition, nucleus, cytoplasm, output_filepath, final_data)

    particles_data = process_dbscan(output_filepath, file_name, resolution_xy, resolution_z, eps, min_samples,
                                    nucleus, cytoplasm, p_value)

    #plot data
    plot_dbscan(file_name, particles_data, nucleus, cytoplasm, cytokine_analysis)
