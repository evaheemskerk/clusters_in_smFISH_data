import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
import pandas as pd


def plot_dbscan(file_name_spots, cluster_data, nucleus, cytoplasm, cytokine_analysis):
    df = pd.read_csv(file_name_spots)
    if cytoplasm == True:
        df = df[df['nuclei_mask'] == 0]
    if nucleus == True:
        df = df[df['nuclei_mask'] > 0]

    if cytokine_analysis == 'TNF':
        df = df[df['cytokine'] == 'TNF']
    if cytokine_analysis == 'IFN':
        df = df[df['cytokine'] == 'IFN']

    # Part 1
    #find general_data
    cell_data = pd.DataFrame()
    cell_data['ID'] = df['ID'].unique()
    cell_data['Type'] = cell_data['ID'].str.split('_').str[0].str[2:]
    cell_data['Timepoint'] = cell_data['ID'].str.split('_').str[1]

    general_data = cell_data.groupby(['Type', 'Timepoint']).size().reset_index(name='N_cells')
    print(general_data)

    ### Which cytokines do the cells produce?
    df_unique = cluster_data.drop_duplicates(subset=['ID'])
    cluster_data_CT = df_unique[['ID', 'Cytokine', 'Type', 'Timepoint']]
    cytokine_type = cluster_data_CT.groupby(['Type', 'Timepoint', 'Cytokine']).size().reset_index(name='N_clustered_cells')
    cytokine_type = pd.merge(cytokine_type, general_data, on=['Timepoint', 'Type'], suffixes=('_clustered', '_all'))
    cytokine_type['percentage'] = cytokine_type['N_clustered_cells']/cytokine_type['N_cells']
    cytokine_type = cytokine_type.sort_values(by=['Type', 'Timepoint'])
    cytokine_type['cytokine'] = cytokine_type['Cytokine'].apply(lambda x: ', '.join(x.replace("'", "").split()) if isinstance(x, str) else str(x))
    cytokine_type['cytokine'] = cytokine_type['cytokine'].apply(lambda x: x.strip("[]").replace("', '", ", ") if isinstance(x, str) else str(x))


    for timepoint in cytokine_type['Timepoint'].unique():
        data = cytokine_type[cytokine_type['Timepoint'] == timepoint]
        sns.barplot(data = data, x = 'Cytokine', y = 'percentage', hue = 'Type')
        plt.ylim(0,0.50)
        plt.title(f'Which cytokine do the cells produce? Timepoint: {timepoint}')
        plt.show()


    ### How many clusters are there per cell?
    df_unique = cluster_data.drop_duplicates(subset=['ID'])
    cluster_data_NC = df_unique[['ID', 'Type', 'Timepoint', 'N_clusters', 'Cytokine']]
    print(cluster_data_NC)
    cytokine_Nclusters = cluster_data_NC.groupby(['Type', 'Timepoint', 'N_clusters', 'Cytokine']).size().reset_index(name='cluster_N_cells')
    cytokine_Nclusters = cytokine_Nclusters.merge(general_data[['Type', 'Timepoint', 'N_cells']],
                                         on=['Type', 'Timepoint'], how='left')

    #Calculate the percentage
    cytokine_Nclusters['percentage'] = cytokine_Nclusters['cluster_N_cells'] / cytokine_Nclusters['N_cells']
    cytokine_Nclusters = cytokine_Nclusters.sort_values(by=['Type', 'Timepoint'])
    cytokine_Nclusters['N_clusters'] = pd.to_numeric(cytokine_Nclusters['N_clusters'])

    for timepoint in cytokine_Nclusters['Timepoint'].unique():
        data = cytokine_Nclusters[cytokine_Nclusters['Timepoint'] == timepoint]
        sns.pointplot(data=data, x='N_clusters', y='percentage', hue='Type', dodge=True)
        plt.title(f'How many clusters per cell? Timepoint: {timepoint}')
        plt.show()

    ### Number of particles per cluster
    print(cytokine_Nclusters)
    g = sns.relplot(data = cytokine_Nclusters, kind ='line', x = 'Timepoint', y = 'N_clusters', hue = 'Type')
    g.fig.subplots_adjust(top=0.9)  # Extra ruimte voor de titel
    g.fig.suptitle(f'Numbers of clusters per cell per timepoint', fontsize=14)
    g.set_ylabels(label = 'Average number of clusters')
    g.set(ylim=(0, 25))
    g._legend.set_bbox_to_anchor((1.0, 0.8)) 
    g._legend.set_title("Type")  
    plt.show()

    ## N_particles per cluster
    for timepoint in cluster_data['Timepoint'].unique():
        data = cluster_data[cluster_data['Timepoint'] == timepoint]
        sns.histplot(data=data, x='N_particles_per_cluster', binwidth = 1, hue = 'Type', multiple= 'dodge')
        plt.title(f'Number of particles per cluster. TP: {timepoint}')
        plt.xlim(0,40)
        plt.show()

    ###Number of particles per cluster
    g = sns.relplot(data=cluster_data, kind='line', x='Timepoint', y='N_particles_per_cluster', hue='Type')
    g.fig.subplots_adjust(top=0.9)
    g.fig.suptitle('Average number of particles per cluster per timepoint', fontsize=14)
    g.set_ylabels(label = 'Number of particles per cluster')
    g.set(ylim=(2.5, 7))
    g._legend.set_bbox_to_anchor((1.0, 0.8))
    g._legend.set_title("Type")
    plt.show()


    # ## Part 2
    df = pd.read_csv(file_name_spots)
    if cytoplasm == True:
        df = df[df['nuclei_mask'] == 0]
    if nucleus == True:
        df = df[df['nuclei_mask'] > 0]

    if cytokine_analysis == 'TNF':
        df = df[df['cytokine'] == 'TNF']
    if cytokine_analysis == 'IFN':
        df = df[df['cytokine'] == 'IFN']

    df = df.groupby('ID')['cytokine'].unique().reset_index()
    df['cytokine'] = df['cytokine'].apply(lambda x: ', '.join(x) if isinstance(x, (list, np.ndarray)) else str(x))
    df['Type'] = df['ID'].str.split('_').str[0].str[2:]
    df['Timepoint'] = df['ID'].str.split('_').str[1]

    producer_data = df.groupby(['Type', 'Timepoint', 'cytokine']).size().reset_index(name='N_cells')

    producer_data = producer_data.merge(cytokine_type[['Type', 'Timepoint', 'cytokine', 'N_clustered_cells']], on = ['Type', 'Timepoint', 'cytokine'])
    producer_data['percentage'] = producer_data['N_clustered_cells']/producer_data['N_cells']
    producer_data = producer_data.sort_values(by=['Type', 'Timepoint'], ascending=[True, True])
    producer_data['Timepoint'] = pd.Categorical(
        producer_data['Timepoint'],
        categories=sorted(producer_data['Timepoint'].unique()),
        ordered=True
    )
    print(producer_data)

    for cytokine in producer_data['cytokine'].unique():
        data = producer_data[producer_data['cytokine'] == cytokine]
        sns.barplot(data=data, x = 'Timepoint', y = 'percentage', hue = 'Type')
        plt.ylim(0,0.55)
        plt.ylabel('Proportion')
        plt.xlabel('Time point')
        plt.title(f'Proportion of cells that have clusters for producers: {cytokine}')
        plt.show()

    return

if __name__ == '__main__':
    file_name_spots = "C:/Users/Eva Heemskerk/Documents/Delft/TU Delft/MEP/4. Results/DATA/Spots_with_localization.csv"
    cluster_data = pd.read_excel("C:/Users/Eva Heemskerk/PycharmProjects/MEP/particles_data_nucleus_both.xlsx")
    nucleus = False
    cytoplasm = True

    # Fill in BOTH if interested in both TNF and IFN, or fill in 'TNF' or 'IFN' when interested in a single analysis
    cytokine_analysis = 'BOTH'

    plot_dbscan(file_name_spots, cluster_data, nucleus, cytoplasm, cytokine_analysis)
