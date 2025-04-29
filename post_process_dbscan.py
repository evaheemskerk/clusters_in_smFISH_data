import pandas as pd
from sklearn.cluster import DBSCAN
from tqdm import tqdm
import numpy as np
from pre_process_files import pre_process_data_set

def dbscan_details (eps, min_samples, particles) :
    particles_df = pd.DataFrame()
    dbscan = DBSCAN(eps=eps, min_samples = min_samples)
    labels = dbscan.fit_predict(particles)
    for i in range(len(particles)):
        new_row = {'label': labels[i], 'Z':particles[i,0], 'Y': particles[i,1], 'X':particles[i, 2]}
        particles_df = pd.concat([particles_df, pd.DataFrame([new_row])], ignore_index=True)
    return particles_df

def calculate_distance(particles):
    # distances
    distance = []
    for j in range(len(particles)):
        a = j + 1
        while a < len(particles):
            dist = np.sqrt((particles[a, 0] - particles[j, 0]) ** 2 + (
                    particles[a, 1] - particles[j, 1]) ** 2 + (
                                   particles[a, 2] - particles[j, 2]) ** 2)
            distance.append(dist)
            a += 1
    return distance


def process_dbscan (file_name_dbscan, file_name_spots, resolution_xy, resolution_z, eps, min_samples, nucleus, cytoplasm, p_value):

    data = pd.read_excel(file_name_dbscan, sheet_name="Sheet1")
    data['Type'] = data['ID'].str.split('_').str[0].str[2:]
    data['TimePoint'] = data['ID'].str.split('_').str[1]


    data_spots = pre_process_data_set(file_name_spots, nucleus, cytoplasm, resolution_z, resolution_xy)
    data = data[data.P_value <= p_value]

    particles_data = pd.DataFrame()
    for timepoint in tqdm(data['TimePoint'].unique()):
        data_timepoint = data[data.TimePoint == timepoint]
        for type in data_timepoint['Type'].unique():
            data_type = data_timepoint[data_timepoint['Type'] == type]
            for id in data_type['ID'].unique():
                spots = data_spots[data_spots['ID'] == id]
                #print(spots)
                cytokines = spots['cytokine'].unique()
                particles = spots.iloc[:, [3, 4, 5]].values
                particles_df = dbscan_details(eps, min_samples, particles)
                arr = particles_df['label'].unique()
                N_clusters = len(arr) - 1 if -1 in arr else len(arr)
                for cluster in particles_df['label'].unique():
                    if cluster == -1:
                        continue
                    particles = particles_df[particles_df['label'] == cluster].iloc[:, [1, 2, 3]].values
                    n_particles_per_cluster = len(particles)
                    distances = calculate_distance(particles)
                    new_row = {'ID': id, 'Cytokine': str(cytokines), 'Type': type, 'Timepoint': timepoint, 'N_clusters': N_clusters, 'Cluster_label': cluster , 'N_particles_per_cluster': n_particles_per_cluster, 'Distances': distances}
                    particles_data = pd.concat([particles_data, pd.DataFrame([new_row])], ignore_index=True)

    return particles_data

if __name__ == '__main__':
    resolution_xy = 65
    resolution_z = 200
    eps = 600
    min_samples = 3
    nucleus = True
    cytoplasm = False
    p_value = 0.05

    ##files
    file_name_dbscan = "C:/Users/Eva Heemskerk/Documents/Delft/TU Delft/MEP/4. Results/Control_clustering_code/Nucleus_cluster_control_data_IFN.xlsx"
    file_name_spots = 'C:/Users/Eva Heemskerk/Documents/Delft/TU Delft/MEP/4. Results/DATA/Spots_with_localization.csv'

    particles_data = process_dbscan (file_name_dbscan, file_name_spots, resolution_xy, resolution_z, eps, min_samples, nucleus, cytoplasm, p_value)
    particles_data.to_excel('particles_data_nucleus_IFN.xlsx', index=False)
