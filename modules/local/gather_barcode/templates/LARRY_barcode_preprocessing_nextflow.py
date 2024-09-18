#!/opt/conda/bin/python
##/usr/local/bin/python

##Libraries
import os

#os.mkdir('./mplconfigdir')
os.environ['MPLCONFIGDIR'] = './mplconfigdir'

from typing import NamedTuple, List, Dict
import pyfastx
import re
from collections import defaultdict, Counter
from umi_tools import UMIClusterer
import numpy as np
import pandas as pd

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
        larcs = UMIclusterer(seq_dict = larry_ids , ham_dist = int(${ham_larry}))
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
for (CB, LARRY), UMI in CB_LARRY_group.items():
    umi_counts = dict(Counter(UMI))
    #Here you perform clustering with hamming distance 1
    umi_seqs = getUMIset(umi_counts , ham_dist = int(${params.ham_umi}))
    cb_larry_counts.append(Molecule(cell_barcode = CB, 
                                    count = umi_seqs.counts, 
                                    LARRY_barcode = LARRY, 
                                    umi_set = tuple(umi_seqs.umi_seq)))

#For every cell get the LARRY and UMI information
cell_barcode_grouped = defaultdict(list)
for molecule in cb_larry_counts:
    cell_barcode_grouped[molecule.cell_barcode].append((molecule.LARRY_barcode, molecule.count, molecule.umi_set))

#This is within cell filtering.
LARRY_filtered = []
for cb, counting in cell_barcode_grouped.items():
    lbc_umi_combi = {lbc: umi_seqs for lbc, umi_cts, umi_seqs in counting}
    
    def getUMIinClust(clustered_seqs):
        umi_seqs = {}
        for clust_seq in clustered_seqs:
            umi_bcs = list(np.concatenate(list(lbc_umi_combi[seq] for seq in clust_seq.associated_seqs)).flat)
            umi_bcs = dict(Counter(umi_bcs))
            umis = getUMIset(umi_bcs , ham_dist = int(${params.ham_umi}))
            umi_seqs[clust_seq.main_seq] = umis
        return(umi_seqs)
        
    if len(counting) == 1:
        LARRY_filtered.append(Molecule(cell_barcode = cb, 
                                       count = counting[0][1],
                                       LARRY_barcode = counting[0][0],
                                       umi_set = counting[0][2]))
    elif len(counting) != 1:
        larry_bcs = {lbc: umi_cts for lbc, umi_cts, umi_seqs in counting}
        clustered_lbcs = UMIclusterer(seq_dict = larry_bcs , ham_dist = int(${ham_larry}))
        # if all LARRY barcodes within a cell likely stem from sequencing errors,
        # we select the most common one and subsequently re-run umi-tools on the umis
        # to ensure we do not count any umi more than once
        
        # if a cell is associated with only on LARRY barcode group, i.e. all detected LARRY barcodes 
        # originate from the same lineage tracing barcode we only select the group-defining one
        # (as specified by umi-tools clusterer)
        if len(clustered_lbcs) == 1:
            umi_clu = getUMIinClust(clustered_lbcs)
            for lbc, umis_clus in umi_clu.items():
                LARRY_filtered.append(Molecule(cell_barcode = cb, 
                                               count = umis_clus.counts,
                                               LARRY_barcode = lbc,
                                               umi_set = tuple(umis_clus.umi_seq)))
#         # If a cell is associated with more than one LARRY barcode group (i.e. a cell is associated 
#         # with more than one clonal group), we check that none of the detected barcodes is detected
#         # due to contamination
        elif len(clustered_lbcs) > 1:
            umi_clu = getUMIinClust(clustered_lbcs)
            larrys = {lbc:umis_clus.counts for lbc, umis_clus in umi_clu.items()}
            for larc, reads in larrys.items():
                # we filter out LARRY barcodes within a cell that occur less than half often in comprarison to
                # the most commonly detected LARRY barcode. This assumption is based on the code UMI-tools uses 
                # (https://umi-tools.readthedocs.io/en/latest/the_methods.html). 
                # We thereby assume that LARRY barcodes that are detected at that lower level are due to 
                # contamination, e.g., from ambient RNA. This assumption needs further testing, but may be a good 
                # first approximation.
                if reads >= (max(larrys.values()) * float(${within_cell_cutoff})):
                    LARRY_filtered.append(Molecule(cell_barcode = cb, 
                                           count = reads,
                                           LARRY_barcode = larc,
                                           umi_set = tuple(umi_clu[larc].umi_seq)))

#Group LARRY barcodes together to look at expression level within the same clone.
larry_barcode_grouped = defaultdict(list)
for molecule in LARRY_filtered:
    larry_barcode_grouped[molecule.LARRY_barcode].append((molecule.cell_barcode, 
                                                          molecule.count, 
                                                          molecule.umi_set))

#Perform a within clone filtering.
cell_filtered = []
for larry, counting in larry_barcode_grouped.items():
    if len(counting) == 1:
        cell_filtered.append(Molecule(cell_barcode = counting[0][0], 
                                      count = counting[0][1], 
                                      LARRY_barcode = larry,
                                      umi_set = counting[0][2]))
    elif len(counting) > 1:
        cells =  {cb:umi_nb for cb, umi_nb, umi_set in counting}
        cell_umis =  {cb:umi_set for cb, umi_nb, umi_set in counting}
        # We filter out all cells in which the LARRY barcode was detected at less than a quarter of the expression 
        # of the most commonly expressed barcode across all cells within the sample
        # NOTE: this is a first approximation and the filtering may need to be adjusted based on more objective
        # filtering parameter estimation
        exp_filt = dict(filter(lambda elem: elem[1] >= np.mean(list(cells.values())) * float(${within_clone_cutoff}), cells.items()))
        for cb, count in exp_filt.items():
            cell_filtered.append(Molecule(cell_barcode = cb, 
                                          count = count, 
                                          LARRY_barcode = larry, 
                                          umi_set = cell_umis[cb]))

#Filter based on minimum larry umi parameter
cell_filtered_df = pd.DataFrame(cell_filtered)
cell_filtered_umi_df = cell_filtered_df[cell_filtered_df["count"] > int(${min_larry_umi})]

#Add sample ids
cell_filtered_df = cell_filtered_df.assign(sample_id = "${meta.id}")
cell_filtered_umi_df = cell_filtered_umi_df.assign(sample_id = "${meta.id}")

#Select columns of interest
cell_filtered_df = cell_filtered_df[["sample_id", "LARRY_barcode", "cell_barcode", "count"]]
barcode_df = cell_filtered_umi_df[["sample_id", "LARRY_barcode", "cell_barcode", "count"]]

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
clustered_bcs = clusterer(larrys, threshold = int(${ham_larry}))

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

#Combine overlapping lists
def merge_overlapping_lists(lists):
    merged_lists = []
    
    while lists:
        # Start with the first list in the remaining lists
        first, *rest = lists
        first = set(first)
        
        # Initialize a list of lists that will be merged in this iteration
        merged = False
        for i, current_list in enumerate(rest):
            # If there's an overlap, merge the lists
            if first & set(current_list):
                lists = rest[:i] + rest[i+1:]  # Remove the current list from the remaining lists
                lists.append(list(first | set(current_list)))  # Add the merged list back to the list
                merged = True
                break

        # If no lists were merged in this iteration, move first list to merged_lists
        if not merged:
            merged_lists.append(list(first))
            lists = rest

    return merged_lists

#Calculate for the clones with different LARRY labels
clonal_group_info_1 = pd.DataFrame(clonal_groups)
clonal_group_info_1 = sub_group_filtering(clonal_group_info_1 , ["clonal_id"] , ["clonal_id","LARRY_barcode"])
filt_out_1 = clonal_group_info_1.query('UMI_avg_subgroup < UMI_avg*${within_clone_cutoff}').index
clonal_group_final_1 = clonal_group_info_1.drop(index = filt_out_1, columns = ['UMI_avg', 'UMI_avg_subgroup'])

'''
#Take all the clones with overlapping LARRY together
multiple_clone_cell_list = {cell : [clone_id for clone_id  in clonal_group_final_1[clonal_group_final_1.cell_id == cell]["clonal_id"].tolist()] for cell in set(clonal_group_final_1.cell_id)}
clone_merge_dict = {clonal_group : "clonal_supergroup_" + str(number)  for number , clonal_groups in enumerate(merge_overlapping_lists(list(multiple_clone_cell_list.values()))) for clonal_group in clonal_groups}
clone_merge = pd.DataFrame({"clonal_id" : list(clone_merge_dict.keys()), "clone_merge_id" : list(clone_merge_dict.values())})
clonal_group_info_2 = sub_group_filtering(pd.merge(clonal_group_final_1 , clone_merge),["clone_merge_id"] , ["clonal_id"])
filt_out_2 = clonal_group_info_2.query('UMI_avg_subgroup < UMI_avg*${within_clone_cutoff}').index
clonal_group_final_2 = clonal_group_info_2.drop(index = filt_out_2, columns = ['UMI_avg', 'UMI_avg_subgroup'])

#Keep clones with multiple LARRY label separate
multiple_clone_cell_list = {cell : [clone_id for clone_id  in clonal_group_final_2[clonal_group_final_2.cell_id == cell]["clonal_id"].tolist()] for cell in set(clonal_group_final_2.cell_id)}
multiple_clone_cell_list_sort_join = {key : "|".join(sorted(value)) for key , value in multiple_clone_cell_list.items()}
separate_clusters = {combined_clonal : "clonal_sepgroup_" + str(number) for number , combined_clonal in enumerate(set(multiple_clone_cell_list_sort_join.values()))}
separate_cluster_cells = {cell : separate_clusters[combine_clone] for cell , combine_clone in multiple_clone_cell_list_sort_join.items()}
clone_sep_merge = pd.DataFrame({"cell_id" : list(separate_cluster_cells.keys()) , "clone_sep_id" : list(separate_cluster_cells.values())})
clonal_group_info_3 = sub_group_filtering(pd.merge(clonal_group_final_2 , clone_sep_merge) ,["clone_sep_id"] , ["clone_sep_id","LARRY_barcode"])
filt_out_3 = clonal_group_info_3.query('UMI_avg_subgroup < UMI_avg*${within_clone_cutoff}').index
clonal_group_final_3 = clonal_group_info_3.drop(index = filt_out_3, columns = ['UMI_avg', 'UMI_avg_subgroup'])
'''


#Write out the dataframe
clonal_group_final_1.to_csv("${meta.id}_${ham_larry}_${within_cell_cutoff}_${within_clone_cutoff}_${min_larry_umi}_clone_output.csv")