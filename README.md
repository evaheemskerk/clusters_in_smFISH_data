# clusters_in_smFISH_data
This algorithm allows you to find the clusters in smFISH data. 

# Content
* main_nucleus: run when interested in the nucleus
* main_cytoplasm: run when interested in the cytoplasm 

## Installation
Before using the python algorithm, ensure you have the following libraries installed in your local environment: numpy pandas seaborn tqdm scikit.learn scikit.image matplotlib.pyplot
If not; use the package manager [pip](https://pip.pypa.io/en/stable/) to install these.

```bash
pip install numpy pandas seaborn tqdm scikit.learn scikit.image matplotlib.pyplot
```
Make sure to use matplotlib.pyplot version 3.9.2. 
Download all the python files in this repository. 

## Inputs 
Change the inputs at the bottom of the script. 

Input parameters:
* **nucleus** = Boolean set on True when interested in nucleus 
* **cytoplasm** = Boolean set on True when interested in cytoplasm 
* **resolution_z** = resolutioin in Z-axis
* **resolution_xy** = resolution in XY-direction
* **eps** = input value for DBSCAN
* **min_samples** = input value for DBSCAN
* **p_value** = p_value for statistic test 
  
Input files:
* **file_name** = File path to the Spots_with_localization.csv file
* **file_path_cellpose_n** = File path to the segmentation images of the nucleus
* **file_path_cellpose_c** = File path to the segmentation images of the cytoplasm (does not exist when main_cytoplasm is run)
* **file_path_outline_files** =  File path to the outline files
* **file_path_output** = Output file path

# File handling
Be sure that all outline files are saved as follows: 
```python
file_path_outline_files\{timepoint}\FQ_outline
```
And the files are saved as:
```python
{Donor}{Type}_{Timepoint}__{Sample}
```

Be sure that the nucleus segmentation images files are all saved in one folder and the images are saved as follows: 
```python
{Donor}{Type}_{Timepoint}__{Sample}_{stain}_{firstcut}_{second_cut}
```

Be sure that the cytoplasm segmentation images files are all saved as follows: 
```python
file_path_cellpose_c\{timepoint}
```
And the files are saved as:
```python
{Donor}{Type}_{Timepoint}__{Sample}_{stain}_{firstcut}_{second_cut}
```

# Run
Run the script. 
