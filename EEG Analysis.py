# -*- coding: utf-8 -*-
"""
Created on Wed Apr  3 17:53:14 2024

@author: lukeb
"""
#%% Imports
import pandas as pd
import math 
import numpy as np
import matplotlib.pyplot as plt
from scipy.signal import find_peaks
from scipy.stats import rankdata
import itertools
from matplotlib import cm
from matplotlib.colors import Normalize

#%% Data import
file_paths = [
    r"C:\Users\lukeb\OneDrive\Documents\research\s01_ex01_s01.csv",
    r"C:\Users\lukeb\OneDrive\Documents\research\s01_ex01_s02.csv",
    r"C:\Users\lukeb\OneDrive\Documents\research\s01_ex01_s03.csv",
    # Add more paths as needed
]

# Each time series is 24000 long 
chunk_size = 12000

#%% Load each dataset into a DataFrame and collect them into a list
datasets = []
for file_path in file_paths:
    data = pd.read_csv(file_path)
    datasets.append(data)
    
#%% define time to be length of timeseries    
for data in datasets:
    time = np.arange(1, len(data) + 1)
    
#%% def analyze_multiple_datasets
def analyze_multiple_datasets(datasets):
    results = []
    for data in datasets:
        result = analyze_single_dataset(data)
        results.append(result)
    return results

#%% def analyze_single_dataset
def analyze_single_dataset(data):
    # Extract time series columns
    node_columns = ['T7', 'F8', 'Cz', 'P4']
    time_series_dict = {col: data[col].values for col in node_columns}
    
    # Calculate ordinal patterns
    ordinal_patterns_dict = calculate_ordinal_patterns(time_series_dict)
    
    # Calculate PE for ordinal patterns in each chunk
    pe_chunk_OPdict = calculate_pe_for_chunks(ordinal_patterns_dict)
    
    # Calculate PE for periodics in each chunk
    PE_p111 = calculate_periodics_pe(ordinal_patterns_dict)  # <-- Here is the issue
    
    # Combine PE for ordinal patterns and periodics
    PE = combine_pe_values(pe_chunk_OPdict, PE_p111)
    
    # Calculate Phi values
    phi_values = calculate_phi_values(ordinal_patterns_dict, PE_p111)
    
    return {
        'Time': np.arange(1, len(data) + 1),
        'Permutation Entropy': PE,
        'Phi values': phi_values
    }

#%% def peaks_time_diff
def peaks_time_diff(time_series):
    peaks, _ = find_peaks(time_series)
    t_at_peak = time[peaks]
    time_diff = np.diff(t_at_peak)
    return time_diff

#%% def ordinal_distribution
def ordinal_distribution(data, dx=3, dy=1, taux=1, tauy=1, return_missing=True, tie_precision=8):
    try:
        ny, nx = np.shape(data) 
        data = np.array(data)
    except:
        nx = np.shape(data)[0]
        ny = 1
        data = np.array([data])

    if tie_precision is not None:
        data = np.round(data, tie_precision)

    partitions = np.concatenate(
        [
            [np.concatenate(data[j:j+dy*tauy:tauy,i:i+dx*taux:taux]) for i in range(nx-(dx-1)*taux)] 
            for j in range(ny-(dy-1)*tauy)
        ]
    )

    symbols = np.apply_along_axis(rankdata, 1, partitions, method='min') - 1
    symbols, symbols_count = np.unique(symbols, return_counts=True, axis=0)

    probabilities = symbols_count / len(partitions)

    if return_missing == False:
        return symbols, probabilities
    else:
        all_symbols = list(map(list, list(itertools.permutations(np.arange(dx*dy)))))
        all_probs = [0] * 6
        for i in range(len(all_symbols)):
            for j in range(len(symbols)):
                if np.array_equal(all_symbols[i], symbols[j]):
                    all_probs[i] += probabilities[j]
    
    return all_symbols, all_probs

#%% def calculate_ordinal_patterns
def calculate_ordinal_patterns(time_series_dict):
    ordinal_patterns_dict = {}
    for key, time_series in time_series_dict.items():
        ordinal_patterns_list = []
        for i in range(0, len(time_series), chunk_size):
            chunk = time_series[i:i+chunk_size]
            peak_time_diff = peaks_time_diff(chunk)
            ordinal_patterns = ordinal_distribution(peak_time_diff)
            ordinal_patterns_list.append(ordinal_patterns)
        ordinal_patterns_dict[key] = ordinal_patterns_list
    return ordinal_patterns_dict

#%% def calculate_periodics_pe
def calculate_periodics_pe(ordinal_patterns_dict):
    p111_chunk_dict = {}
    for key, patterns_list in ordinal_patterns_dict.items():
        p111_chunk_dict[key] = []
        for patterns in patterns_list:
            _, probabilities = patterns
            summed_probability = (1 - sum(probabilities))
            pe = permutation_entropy({key: (None, [summed_probability])})[key]
            p111_chunk_dict[key].append(pe)  # Calculate and append PE directly
    return p111_chunk_dict

#%% def permutation_entropy
def permutation_entropy(ordinal_distributions, probs=True):
   pe_dict = {}
   for series_name, (patterns, probabilities) in ordinal_distributions.items():
       pe = 0.0
       for p in probabilities:
           if p > 0:
               pe -= (p * math.log(p))/math.log(7)
       pe_dict[series_name] = pe
       
   return pe_dict

#%% def calculate_pe_for_chunks
def calculate_pe_for_chunks(ordinal_patterns_dict, probs=True):
    pe_dict = {}
    for key, patterns_list in ordinal_patterns_dict.items():
        pe_list = []
        for patterns in patterns_list:
            # Calculate PE for each chunk of ordinal patterns
            pe_chunk = permutation_entropy({key: patterns})
            pe_list.append(pe_chunk[key])  # Append PE for the current chunk
        pe_dict[key] = pe_list
    return pe_dict

#%% def combine_pe_values
def combine_pe_values(pe_chunk_OPdict, PE_p111):
    PE = {}
    for key, chunk_values in pe_chunk_OPdict.items():
        combined_values = []
        for i, (chunk_value, pe_value) in enumerate(zip(chunk_values, PE_p111[key])):
            combined_value = chunk_value + pe_value
            combined_values.append(combined_value)
        PE[key] = combined_values
    return PE

#%% def calculate_phi_values
def calculate_phi_values(ordinal_patterns_dict, PE_p111):
    phi1 = {}
    phi2 = {}
    phi3 = {}
    for key, patterns_list in ordinal_patterns_dict.items():
        phi1[key] = []
        phi2[key] = []
        phi3[key] = []
        for i, patterns in enumerate(patterns_list):
            prob_012 = patterns[1][patterns[0].index([0, 1, 2])]
            prob_210 = patterns[1][patterns[0].index([2, 1, 0])]
            prob_102 = patterns[1][patterns[0].index([1, 0, 2])]
            prob_201 = patterns[1][patterns[0].index([2, 0, 1])]
            prob_021 = patterns[1][patterns[0].index([0, 2, 1])]
            prob_120 = patterns[1][patterns[0].index([1, 2, 0])]
            prob_p111 = PE_p111[key][i]
            phi1_val = math.sqrt(prob_012**2 + prob_210**2 + (prob_p111**2 / math.sqrt(3)))
            phi2_val = math.sqrt(prob_102**2 + prob_201**2 + (prob_p111**2 / math.sqrt(3)))
            phi3_val = math.sqrt(prob_021**2 + prob_120**2 + (prob_p111**2 / math.sqrt(3)))
            phi1[key].append(phi1_val)
            phi2[key].append(phi2_val)
            phi3[key].append(phi3_val)
            
    return phi1, phi2, phi3

#%% results
results = analyze_multiple_datasets(datasets)
for i, result in enumerate(results):
    print(f"Results for Dataset {i+1}:")
    for node, pe_values in result['Permutation Entropy'].items():
        print(f"Node: {node}")
        for j, pe_value in enumerate(pe_values):
            print(f"Chunk {j+1} - Permutation Entropy: {pe_value}")
            #print(f"Chunk {j+1} - Phi1 value: {result['Phi values'][0][node][j]}")
            #print(f"Chunk {j+1} - Phi2 value: {result['Phi values'][1][node][j]}")
            #print(f"Chunk {j+1} - Phi3 value: {result['Phi values'][2][node][j]}")

#%% Plotting PE

# Initialize a figure and 3D axes
fig = plt.figure(figsize=(10, 6))
ax = fig.add_subplot(111, projection='3d')

# Get unique nodes across all datasets
all_nodes = []
for result in results:
    all_nodes.extend(result['Permutation Entropy'].keys())
all_nodes = list(set(all_nodes))

# Set up color map and normalization
cmap = plt.get_cmap('tab10')  # You can choose any colormap you prefer
norm = Normalize(vmin=0, vmax=len(all_nodes) - 1)

# Plot PE vs. time for each dataset
for i, result in enumerate(results):
    time = result['Time']
    for j, node in enumerate(all_nodes):
        pe_values = result['Permutation Entropy'].get(node, [])
        if pe_values:
            color = cmap(norm(j))
            ax.plot(time[:len(pe_values)], [i+1]*len(pe_values), pe_values, color=color, label=f'Dataset {i+1}, Node {node}')

# Set labels and title
ax.set_xlabel('Time')
ax.set_ylabel('Dataset')
ax.set_zlabel('Permutation Entropy')
ax.set_title('Permutation Entropy vs. Time for Each Dataset')
ax.legend(loc='upper left', bbox_to_anchor=(1, 1))
plt.tight_layout()
plt.show()

#%% Citations
# Ordpy for ordianl patterns and PE:A. A. B. Pessa, H. V. Ribeiro, ordpy: A Python package for data analysis with permutation entropy and ordinal network methods, Chaos 31, 063110 (2021). arXiv:2102.06786
# people who ran tests: Abo Alzahab, N., Di Iorio, A., Apollonio, L., Alshalak, M., Gravina, A., Antognoli, L., Baldi, M., Scalise, L., & Alchalabi, B. (2021). Auditory evoked potential EEG-Biometric dataset (version 1.0.0). PhysioNet. https://doi.org/10.13026/ps31-fc50.
# Physio.net (Filtered data): Goldberger, A., Amaral, L., Glass, L., Hausdorff, J., Ivanov, P. C., Mark, R., ... & Stanley, H. E. (2000). PhysioBank, PhysioToolkit, and PhysioNet: Components of a new research resource for complex physiologic signals. Circulation [Online]. 101 (23), pp. e215–e220.
# PE on graph signal: J.S. Fabila-Carrasco, C. Tan, and J. Escudero,“Permutation Entropy for Graph Signal”, arXiv:2110.00628, 2021.     
 
