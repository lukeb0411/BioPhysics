# -*- coding: utf-8 -*-
"""
Created on Fri Mar 29 17:22:48 2024

@author: lukeb
"""
#%% Notes
# Filtered data from Physio.net
# Calculates Peak_time_diff then breaks into chunks
#%% imports
import pandas as pd
import math 
import numpy as np
import matplotlib.pyplot as plt
from scipy.signal import find_peaks
from scipy.stats import rankdata
import itertools
from matplotlib import cm
#%% Data:
data = pd.read_csv(r"C:\Users\lukeb\OneDrive\Documents\research\s01_ex01_s03.csv")

# Extract time series columns
node_columns = ['T7', 'F8', 'Cz', 'P4']
time = data['Time'].tolist()
time_series_dict = {col: data[col].tolist() for col in node_columns}

#%% def peaks_time_diff
def peaks_time_diff(time_series_dict):
    peaks, _ = find_peaks(time_series_dict)
    tAtpeak = []

    for i in range(len(peaks)):
        tAtpeak.append(time[peaks[i]])
    timeDiff = []
    
    for j in range(len(tAtpeak)-1):
        timeDiff.append(np.abs(tAtpeak[j] - tAtpeak[j+1]))
    
    return timeDiff
#%% def ordinal_distribution
def ordinal_distribution(data, dx=3, dy=1, taux=1, tauy=1, return_missing=True, tie_precision=8):

# Updated version of ordinal distribution that considers periodic words in the total population caluclation
# but doesn't treat them as individual words in the return    

    try:
        ny, nx = np.shape(data) 
        data   = np.array(data)
    except:
        nx     = np.shape(data)[0]
        ny     = 1
        data   = np.array([data])

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

    probabilities = symbols_count/len(partitions)

    if return_missing==False:
        return symbols, probabilities
    
    else:
        all_symbols   = list(map(list,list(itertools.permutations(np.arange(dx*dy)))))
        all_probs     = [0]*6
        for i in range(len(all_symbols)):
            for j in range(len(symbols)):
                if np.array_equal(all_symbols[i], symbols[j]):
                    all_probs[i] += probabilities[j]
                
                
    
    return all_symbols, all_probs

#%% Calc ordinal patterns in Chunks
#peak_time_diff array  is 6402 long
chunk_size = 3201
ordinal_patterns_dict = {}

# Step 1: Calculate peak time differences for each time series
peak_time_diff_dict = {}
for key, array in time_series_dict.items():
    peak_time_diff = peaks_time_diff(array)
    peak_time_diff_dict[key] = peak_time_diff
    
# Iterate over each time series in peak_time_diff_dict
for key, peak_time_diff_array in peak_time_diff_dict.items():
    ordinal_patterns_list = []
    
    # Split the peak time differences array into chunks
    for i in range(0, len(peak_time_diff_array), chunk_size):
        chunk = peak_time_diff_array[i:i+chunk_size]
        
        # Calculate ordinal patterns for the chunk
        ordinal_patterns = ordinal_distribution(chunk)
        ordinal_patterns_list.append(ordinal_patterns)
    
    # Store ordinal patterns for the current time series in the dictionary
    ordinal_patterns_dict[key] = ordinal_patterns_list

# Print ordinal pattern probabilities for each time series
#for key, patterns_list in ordinal_patterns_dict.items():
    #print(f"{key}:")
    #for i, patterns in enumerate(patterns_list):
        #print(f"Chunk {i+1}: Ordinal Patterns: {patterns}")

#%% Periodics
p111_chunk_dict = {}

# Iterate over each time series in peak_time_diff_dict
for key, peak_time_diff_array in peak_time_diff_dict.items():
    p111_chunk_dict[key] = []  # Initialize a list to store summed probabilities for each chunk
    
    # Split the peak time differences array into chunks
    for i in range(0, len(peak_time_diff_array), chunk_size):
        chunk = peak_time_diff_array[i:i+chunk_size]
        
        # Calculate ordinal patterns for the chunk
        _, probabilities = ordinal_distribution(chunk)
        
        # Sum the probabilities in the chunk and store it
        summed_probability = (1- sum(probabilities))
        p111_chunk_dict[key].append(summed_probability)

# Print p111 dictionary
#for key, summed_probabilities in p111_chunk_dict.items():
    #print(f"{key}: Periodics: {summed_probabilities}")       
#%% Define calculate_pe_for_chunks
def permutation_entropy(ordinal_distributions, probs=True):
   pe_dict = {}
   for series_name, (patterns, probabilities) in ordinal_distributions.items():
       pe = 0.0
       for p in probabilities:
           if p > 0:
               pe -= (p * math.log(p))/math.log(6)
       pe_dict[series_name] = pe
       
   return pe_dict

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

#%% Calculate PE for each chunk
pe_chunk_OPdict = calculate_pe_for_chunks(ordinal_patterns_dict)

# Print the calculated PE for each chunk in each time series
for key, pe_list in pe_chunk_OPdict.items():
    print(f"Time Series: {key}")
    for i, pe_chunk in enumerate(pe_list):
       print(f"Chunk {i+1}: PE = {pe_chunk}")
       
#%% PE Graph
# Create a PE vs. Time graph with a line chart for each time series
plt.figure(figsize=(10, 6))

for key, pe_list in pe_chunk_OPdict.items():
    plt.plot(time[:len(pe_list)], pe_list, label=key)

plt.title('Permutation Entropy vs. Time')
plt.xlabel('Time')
plt.ylabel('Permutation Entropy')
plt.legend()
plt.grid(True)
plt.show()

#%% Phi's 
phi1 = {}
phi2 = {}
phi3 = {}

# Iterate over each time series in ordinal patterns dictionary
for key, patterns_list in ordinal_patterns_dict.items():
    phi1[key] = []  # Initialize list for storing phi1 values for each chunk
    phi2[key] = []  # Initialize list for storing phi2 values for each chunk
    phi3[key] = []  # Initialize list for storing phi3 values for each chunk
    
    # Iterate over each chunk of ordinal patterns
    for i, patterns in enumerate(patterns_list):
        # Extract probabilities of ordinal patterns 012 and 210
        prob_012 = patterns[1][patterns[0].index([0, 1, 2])]
        prob_210 = patterns[1][patterns[0].index([2, 1, 0])]
        
        # Extract probabilities of ordinal patterns 102 and 201
        prob_102 = patterns[1][patterns[0].index([1, 0, 2])]
        prob_201 = patterns[1][patterns[0].index([2, 0, 1])]
        
        # Extract probabilities of ordinal patterns 021 and 120
        prob_021 = patterns[1][patterns[0].index([0, 2, 1])]
        prob_120 = patterns[1][patterns[0].index([1, 2, 0])]
        
        # Extract probability of p111
        prob_p111 = p111_chunk_dict[key][i]
        
        # Calculate Euclidean magnitudes for phi1, phi2, and phi3
        phi1_val = math.sqrt(prob_012**2 + prob_210**2 + (prob_p111**2 / math.sqrt(3)))
        phi2_val = math.sqrt(prob_102**2 + prob_201**2 + (prob_p111**2 / math.sqrt(3)))
        phi3_val = math.sqrt(prob_021**2 + prob_120**2 + (prob_p111**2 / math.sqrt(3)))
        
        # Append phi values to respective dictionaries
        phi1[key].append(phi1_val)
        phi2[key].append(phi2_val)
        phi3[key].append(phi3_val)

# Print phi dictionaries
#for key, values in phi1.items():
    #print(f"{key}: Phi1 values: {values}")

#for key, values in phi2.items():
    #print(f"{key}: Phi2 values: {values}")

#for key, values in phi3.items():
    #print(f"{key}: Phi3 values: {values}")

#%% Phi plotting 
keys = list(phi1.keys())  # Assuming keys are the same for all dictionaries
fig = plt.figure(figsize=(10, 10))
ax = fig.add_subplot(111, projection='3d')

# Define a color map
colors = cm.rainbow([i / len(keys) for i in range(len(keys))])

for i, key in enumerate(keys):
    ax.scatter(phi2[key], phi3[key], phi1[key], color=colors[i], label=key)

ax.set_title('Phi space')
ax.set_xlabel('Phi2', labelpad=20)
ax.set_ylabel('Phi3', labelpad=20)
ax.set_zlabel('Phi1', labelpad=20)
plt.legend()
plt.show()
#%% Citations
#Ordpy for ordianl patterns and PE:A. A. B. Pessa, H. V. Ribeiro, ordpy: A Python package for data analysis with permutation entropy and ordinal network methods, Chaos 31, 063110 (2021). arXiv:2102.06786
#people who ran tests: Abo Alzahab, N., Di Iorio, A., Apollonio, L., Alshalak, M., Gravina, A., Antognoli, L., Baldi, M., Scalise, L., & Alchalabi, B. (2021). Auditory evoked potential EEG-Biometric dataset (version 1.0.0). PhysioNet. https://doi.org/10.13026/ps31-fc50.
#Physio.net (Filtered data): Goldberger, A., Amaral, L., Glass, L., Hausdorff, J., Ivanov, P. C., Mark, R., ... & Stanley, H. E. (2000). PhysioBank, PhysioToolkit, and PhysioNet: Components of a new research resource for complex physiologic signals. Circulation [Online]. 101 (23), pp. e215–e220.
#PE on graph signal: J.S. Fabila-Carrasco, C. Tan, and J. Escudero,“Permutation Entropy for Graph Signal”, arXiv:2110.00628, 2021.     