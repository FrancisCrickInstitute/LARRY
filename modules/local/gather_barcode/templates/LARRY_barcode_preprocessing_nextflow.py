#!/opt/conda/envs/pymc-dev/bin/python



##Libraries
import os

os.mkdir('./mplconfigdir')
os.mkdir('./compiledir')
os.environ['MPLCONFIGDIR'] = './mplconfigdir'

import pytensor as pt
pt.config.compiledir = './compiledir'

from typing import NamedTuple, List, Dict
import pyfastx
import re
from collections import defaultdict, Counter
from umi_tools import UMIClusterer
import numpy as np
import pandas as pd
from scipy.stats import hypergeom
import matplotlib.pyplot as plt

# Set a random seed for reproducibility
np.random.seed(42)


##Declare object classes
class Read(NamedTuple):
    cell_barcode: str
    umi: str
    LARRY_barcode: str

class Molecule(NamedTuple):
    cell_barcode: str
    LARRY_barcode: str
    count: int 
    umi_set: set

class ClusteredSeq(NamedTuple):
    main_seq: str
    counts: int
    associated_seqs: tuple

class ClusteredUMI(NamedTuple):
    umi_seq: set
    counts: int

class Cell(NamedTuple):
    cell_id: str
    LARRY_barcode: str
    count: int 
        
class clonal_group(NamedTuple):
    clonal_id: str
    cell_id: str
    LARRY_barcode: str
    count: int

##Declare functions
def UMIclusterer(seq_dict: dict, ham_dist: int = 1) -> list[ClusteredSeq]:

    """
    cluster sequences by similarity using umi-tools
    The clustering is performed using the network-based method 'cluster'. 
    For more information please refer to https://umi-tools.readthedocs.io/en/latest/index.html
    
    Parameters
    ----------
    seq_dict: dict
        Dictionary with keys being sequences and values being counts
    threshold: int
        hamming distance threshold. Default: 1
        
    Returns
    ----------
    List of the NamedTuple class ClusteredSeq, which consists of the main sequence of the group ('main_seq'), 
    the number of sequences in the group ('counts') and all sequences in the group ('associated_seqs')
    """

    clusterer=UMIClusterer( cluster_method="cluster")
    seqs = {bytes(seq, encoding="ascii"):count for seq, count in seq_dict.items()}
    clustered_seqs =clusterer(seqs, threshold = ham_dist)
    seq_groups = []
    for seq_group in clustered_seqs:
        seq_group_bc = [bc.decode() for bc in seq_group]
        seq_groups.append(ClusteredSeq(main_seq = seq_group_bc[0], 
                                       counts = len(seq_group_bc), 
                                       associated_seqs = tuple(seq_group_bc)))
    return(seq_groups)
    
     
def getUMIset(umi_seqs, ham_dist: int = 1):

    """
    cluster UMI sequences by similarity using umi-tools
    The clustering is performed using the network-based method 'cluster'. 
    For more information please refer to https://umi-tools.readthedocs.io/en/latest/index.html
    
    Parameters
    ----------
    umi_seqs: dict
        Dictionary with keys being UMI sequences and values being counts
        
    Returns
    ----------
    List of the NamedTuple class ClusteredUMI, which consists of the set of unique UMI sequences and the  
    the number of UMIs ('counts') in the set
    """

    clustered_umis = UMIclusterer(seq_dict = umi_seqs , ham_dist = ham_dist)
    umi_counts = len(clustered_umis)
    umis = []
    for seq_clust in clustered_umis:
        umis.append([seq for seq in seq_clust.associated_seqs])
    umis = set(np.concatenate(umis).flat)
    return(ClusteredUMI(umi_seq = umis,
                       counts = umi_counts))

##Script
#Extract barcodes and UMI from fastq file
filtered_seqs = []

for name, seq, qual in pyfastx.Fastq("${fastq}", build_index=False):
    split_name = name.split("_")
    cell_barcode = split_name[-2]
    umi = split_name[-1]
    if (len(cell_barcode) == 16) and  (len(umi) == 12):
        filtered_seqs.append(Read(cell_barcode = cell_barcode, umi = umi, LARRY_barcode = seq))

#Remove duplicates
dedup_seqs = set(filtered_seqs)

#Group all the LARRY barcodes that have the same cell_barcode and same umi
CB_UMI_group = defaultdict(list)
for read in dedup_seqs:
    CB_UMI_group[(read.cell_barcode, read.umi)].append(read.LARRY_barcode)

#Only larry barcodes are kept that have a unique cell_barcode + umi + larry. In case of more thatn one we account for larry sequencing error.
CB_UMI_group_filtered = []
more_larry_than_expect = []
for (CB, UMI), LARRY in CB_UMI_group.items():
    larry_ids = dict(Counter(LARRY))
    if len(larry_ids) == 1:
        larry_consensus = list(larry_ids.keys())[0]
        CB_UMI_group_filtered.append(Read(cell_barcode = CB, umi = UMI, LARRY_barcode = larry_consensus))
    elif len(larry_ids) > 1:
        larcs = UMIclusterer(seq_dict = larry_ids , ham_dist = int(${params.ham_larry}))
        if len(larcs) == 1:
            CB_UMI_group_filtered.append(Read(cell_barcode = CB, umi = UMI, LARRY_barcode = larcs[0].main_seq))
        elif len(larcs) != 1:
            more_larry_than_expect.append(Read(cell_barcode = CB, umi = UMI, LARRY_barcode = larcs))
            
#Group the cell_barcodes and the LARRY_barcodes together to measure the umi expression
CB_LARRY_group = defaultdict(list)
for read in CB_UMI_group_filtered:
    CB_LARRY_group[(read.cell_barcode, read.LARRY_barcode)].append(read.umi)

#Count the number of umis for every CB LARRY combination.
cb_larry_counts = []
all_larry_counts = []
for (CB, LARRY), UMI in CB_LARRY_group.items():
    umi_counts = dict(Counter(UMI))
    #Here you perform clustering with hamming distance 1
    umi_seqs = getUMIset(umi_counts , ham_dist = ${params.ham_umi})
    all_larry_counts.append(umi_seqs.counts)
    cb_larry_counts.append(Molecule(cell_barcode = CB, 
                                    count = umi_seqs.counts, 
                                    LARRY_barcode = LARRY, 
                                    umi_set = tuple(umi_seqs.umi_seq)))

# Convert list to a DataFrame with a single row
df = pd.DataFrame([all_larry_counts])

# Write DataFrame to CSV
df.to_csv("${meta.id}_" + "numbers_single_row.csv", index=False, header=False)


#Determine the count cutoff
def determine_cutoff(counts , epsilon = 0.01):

    #Initialise comparison variable, we start with -1 because first comparison is with initialised mn_new (0)
    comparison = - 1

    #Initialise p_new variable
    mn_new = 0

    #Initialise rate difference variable
    rate_difference = 1

    while rate_difference > epsilon:

        #MLE of parameter
        mn = np.mean(counts)

        #Update old and new
        mn_old = mn_new
        mn_new = mn

        mn_mean = np.mean([mn_old , mn_new])
        difference = abs(mn_new - mn_old)
        rate_difference = difference / mn_mean

        print(rate_difference)

        #Updata count data
        counts = [count - 1 for count in counts if count > 1]

        #Add comparison value
        comparison += 1
        
    return comparison

cutoff = determine_cutoff(all_larry_counts , epsilon = 0.01)

#Make figure of cutoffs and mean of shifted count data
#Calculate the mean of shifted count data

count_data = np.array(all_larry_counts)
mean_shift_data = []

for ct in range(0,50):
    count_data_shift = [count - ct for count in count_data if count > ct]
    
    mean_shift_data.append([ct +0.5 , np.mean(count_data_shift), len(count_data_shift)])

#Mean shift data to dataframe
mean_plot = pd.DataFrame(mean_shift_data , columns = ["cutoff" ,"mean" ,"amount"])

# Get unique values and their counts
unique_values, counts = np.unique(count_data, return_counts=True)

# Display the frequency table
frequency_table = dict(zip(unique_values, counts))
plot_count_data = pd.DataFrame(list(frequency_table.items()), columns=['count', 'freq'])
plot_count_data["rel_freq"] = plot_count_data["freq"] / len(count_data)
plot_count_data["log_freq"] = np.log(plot_count_data["freq"])

# Create a figure and a set of subplots
fig, ax1 = plt.subplots(figsize=(10, 6))

# Bar plot
ax1.bar(plot_count_data["count"], plot_count_data["log_freq"], label='Count data', color='black')
ax1.set_xlabel('Count')
ax1.set_ylabel('Log Frequency', color='black')
ax1.tick_params(axis='y', labelcolor='black')

# Create a second y-axis
ax2 = ax1.twinx()

# Line plot
ax2.plot(mean_plot["cutoff"], mean_plot["mean"], marker='o', linestyle='-', color='red', label='Count mean')
ax2.set_ylabel('Mean', color='red')
ax2.set_ylim(0, np.max(mean_plot["mean"]) + 0.1 * (np.max(mean_plot["mean"])))
ax2.tick_params(axis='y', labelcolor='red')

# Add vertical line at x = 3
plt.axvline(x=cutoff - 0.5, color='blue', linestyle='--', linewidth=3, label=f'Cutoff x={cutoff - 0.5}')

# Set x-axis limits
ax1.set_xlim(0, 50)

# Add legends
ax1.legend(loc='upper left')
ax2.legend(loc='upper right')

# Title and grid
plt.title("Mean of shifted dataset")
ax1.grid(True)  # Show grid for the bar plot

# Show the plot
plt.tight_layout()  # Adjust layout to prevent clipping

plt.savefig("${meta.id}_" + str(cutoff) + "_mean_shifted_count_data.png")

#Remove LARRY barcodes that are below the treshhold.
cb_larry_counts_filtered = [el for el in cb_larry_counts if el.count >= cutoff]

#Filter based on minimum larry umi parameter
cell_filtered_df = pd.DataFrame(cb_larry_counts_filtered)

#Add sample ids
cell_filtered_df = cell_filtered_df.assign(sample_id = "${meta.id}")

#Select columns of interest
barcode_df = cell_filtered_df[["sample_id", "LARRY_barcode", "cell_barcode", "count"]]

#Rename the cells just to make sure the naming is definitely unique
barcode_df["cell_id"] = barcode_df[["sample_id", "cell_barcode"]].apply(".".join, axis = 1)
barcode_df = barcode_df.drop(["sample_id", "cell_barcode"], axis = 1)

#Convert pandas dataframe to tuples
barcodes = list(barcode_df.itertuples(index = False, name = "Cell"))

#Get list of all the LARRY and cell combinations
all_cells = []
for cell in barcodes:
    all_cells.append(Cell(LARRY_barcode = cell.LARRY_barcode, cell_id = cell.cell_id, count = cell.count))

# LARRY_groups will group cells with exactly the same LARRY barcode sequences together
LARRY_groups = defaultdict(list)
for cell in all_cells:
    LARRY_groups[cell.LARRY_barcode].append({cell.cell_id: cell.count})

#Sum the total LARRY counts for a particular larry barcode. Sum counts of multiple cells together.
larrys = {}
for bc, cell_counts in LARRY_groups.items():
    larc = bytes(bc, encoding="ascii")
    umi_count_sums = sum([sum(list(x.values())) for x in cell_counts])
    larrys[larc] = umi_count_sums

#Here we perform larry clustering across cells.
clusterer = UMIClusterer(cluster_method="cluster")
clustered_bcs = clusterer(larrys, threshold = int(${params.ham_larry}))

#Make the clonal groups
clonal_groups = []
i = 1
for bc_group in clustered_bcs:
    clonal_group_id = "clonal_group_" + str(i)
    for lbc in bc_group:
        for cb_count in LARRY_groups[lbc.decode()]:
            for cb, count in cb_count.items():
                clonal_groups.append(clonal_group(
                    clonal_id= clonal_group_id,
                    cell_id = cb,
                    LARRY_barcode = lbc.decode(),
                    count = count))
    i += 1


#Calculate subgroup averages
def sub_group_filtering(dataset , group_list , sub_group_list):
    dataset.loc[:,"UMI_avg"] = dataset.groupby(by = group_list)["count"].transform(lambda x : x.mean()).to_frame()
    dataset.loc[:,"UMI_avg_subgroup"] = dataset.groupby(by = sub_group_list)["count"].transform(lambda x : x.mean()).to_frame() 
    
    return dataset

#Calculate for the clones with different LARRY labels
clonal_group_info_1 = pd.DataFrame(clonal_groups)
clonal_group_info_1 = sub_group_filtering(clonal_group_info_1 , ["clonal_id"] , ["clonal_id","LARRY_barcode"])
clonal_group_final_1 = clonal_group_info_1.drop(columns = ['UMI_avg', 'UMI_avg_subgroup'])

def cluster_merge(dataset , cutoff_value = 0.5) :

    jaccard_output = {}

    #Create clone dictionaries with keys being the clone name and the values the cells.
    list_of_clones = dataset.groupby('clonal_id')['cell_id'].apply(list).apply(list).tolist()

    dict_clones = {str(number): clone for number , clone in enumerate(list_of_clones)}

    jaccard_array = np.ones((len(list_of_clones), len(list_of_clones)))

    jaccard_final_dict = {}

    while np.sum(jaccard_array) > 0:
        i = 0
        jaccard_array = np.ones((len(list_of_clones), len(list_of_clones)))

        while i < len(list_of_clones):

            jaccard_array[i,i] = 0

            j = i + 1
            while j < len(list_of_clones):

                group1_set = set(list_of_clones[i])
                group2_set = set(list_of_clones[j])
                overlap_size = len(group1_set.intersection(group2_set))
                group1_size = len(group1_set)
                group2_size = len(group2_set)
                jaccard_similarity =  (overlap_size / (group1_size + group2_size - overlap_size))

                if jaccard_similarity >= cutoff_value:
                    jaccard_array[i,j] = jaccard_similarity
                    jaccard_array[j,i] = jaccard_similarity

                else:
                    jaccard_array[i,j] = 0
                    jaccard_array[j,i] = 0


                j += 1
            i += 1
        #Make into dataframe
        df = pd.DataFrame(jaccard_array, index=dict_clones.keys(), columns=dict_clones.keys())

        clone_1 = df.stack().idxmax()[0]
        clone_2 = df.stack().idxmax()[1]

        zero_clones = df.columns[np.sum(df, axis = 1) == 0]

        jaccard_final_dict.update({key : value for key, value in dict_clones.items() if key in zero_clones})
        dict_clones = {key : value for key, value in dict_clones.items() if key not in zero_clones}

        if np.sum(jaccard_array) == 0:
            break

        dict_clones[clone_1 + "_" + clone_2] = list(set(dict_clones[clone_1] + dict_clones[clone_2]))
        dict_clones.pop(clone_1)
        dict_clones.pop(clone_2)

        list_of_clones = list(dict_clones.values())
        
    new_clones = {"clone_" + str(number) : el for number , el in enumerate(jaccard_final_dict.values())}
    new_clones_flat = [(key, value) for key, values in new_clones.items() for value in values]
    new_clones_flat_df = pd.DataFrame(new_clones_flat, columns=['Clone', 'Cell'])
    cells_in_one_clone_idx = new_clones_flat_df.groupby("Cell").filter(lambda x: len(x) == 1).index
    new_clones_flat_df = new_clones_flat_df.loc[cells_in_one_clone_idx]
    new_clones = pd.merge(dataset , new_clones_flat_df , how = "right" , left_on = "cell_id", right_on = "Cell")
    new_clones = new_clones.drop(['clonal_id',"cell_id","LARRY_barcode","count"], axis=1)
    new_clones = new_clones.drop_duplicates()

    return new_clones

clone_group_merged = cluster_merge(clonal_group_final_1 , cutoff_value = 0.5)

output_name_jaccard = "${meta.id}_" + str(cutoff) + "_clone_output.csv"

#Write out the dataframe
clone_group_merged.to_csv(output_name_jaccard, index = False)