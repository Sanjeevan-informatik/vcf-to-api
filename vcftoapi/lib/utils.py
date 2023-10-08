from simplejson import JSONEncoder

import numpy as np


from scipy.cluster.hierarchy import ClusterWarning, ClusterNode
from typing import  List
from warnings import simplefilter


from Bio import SeqIO,  Phylo, AlignIO
from numpy import where
import numpy as np
import pandas as pd
from numpy import fill_diagonal
from pandas import DataFrame
from itertools import combinations


simplefilter("ignore", ClusterWarning)

RED = '\033[31m'
RESET = '\033[0m'

def print_error(msg):
    print(f"{RED}ERROR: {msg}{RESET}")
    pass

class ApiError(Exception):
    status_code = 200
    def __init__(self, message, status_code=None, payload=None):
        Exception.__init__(self)
        self.message = message
        if status_code is not None:
            self.status_code = status_code
        self.payload = payload

    def to_dict(self):
        rv = dict(self.payload or ())
        rv['status'] = 'error'
        rv['message'] = self.message
        return rv


class StrictEncoder(JSONEncoder):
    def __init__(self, *args, **kwargs):
        kwargs["allow_nan"] = False
        kwargs["ignore_nan"] = True
        super().__init__(*args, **kwargs)


def _scipy_tree_to_newick_list(node: ClusterNode, newick: List[str], parentdist: float, leaf_names: List[str]) -> List[str]:
    """Construct Newick tree from SciPy hierarchical clustering ClusterNode

    This is a recursive function to help build a Newick output string from a scipy.cluster.hierarchy.to_tree input with
    user specified leaf node names.

    Notes:
        This function is meant to be used with `to_newick`

    Args:
        node (scipy.cluster.hierarchy.ClusterNode): Root node is output of scipy.cluster.hierarchy.to_tree from hierarchical clustering linkage matrix
        parentdist (float): Distance of parent node of `node`
        newick (list of string): Newick string output accumulator list which needs to be reversed and concatenated (i.e. `''.join(newick)`) for final output
        leaf_names (list of string): Leaf node names

    Returns:
        (list of string): Returns `newick` list of Newick output strings
    """
    if node.is_leaf():
        return newick + [f'{leaf_names[node.id]}:{parentdist - node.dist}']

    if len(newick) > 0:
        newick.append(f'):{parentdist - node.dist}')
    else:
        newick.append(');')
    newick = _scipy_tree_to_newick_list(node.get_left(), newick, node.dist, leaf_names)
    newick.append(',')
    newick = _scipy_tree_to_newick_list(node.get_right(), newick, node.dist, leaf_names)
    newick.append('(')
    return newick


def to_newick(tree: ClusterNode, leaf_names: List[str]) -> str:
    """Newick tree output string from SciPy hierarchical clustering tree

    Convert a SciPy ClusterNode tree to a Newick format string.
    Use scipy.cluster.hierarchy.to_tree on a hierarchical clustering linkage matrix to create the root ClusterNode for the `tree` input of this function.

    Args:
        tree (scipy.cluster.hierarchy.ClusterNode): Output of scipy.cluster.hierarchy.to_tree from hierarchical clustering linkage matrix
        leaf_names (list of string): Leaf node names

    Returns:
        (string): Newick output string
    """
    newick_list = _scipy_tree_to_newick_list(tree, [], tree.dist, leaf_names)
    return ''.join(newick_list[::-1])


def sequence_vcf(_result , samples, vcf_to_seq, consensus):

    number_of_sample= len(samples)


    with open(vcf_to_seq, "w") as outfile:

        for index in range(number_of_sample):


             # Get the reference nucleotides (as letters ATCG)
            sliced_reference_1 = _result.reference_allele[_result.slice_variant_calls]

            # Get the alternate nucleotides (as letters ATCG)
            sliced_alternates_1 = _result.alternate_alleles[_result.slice_variant_calls]
            sliced_alternates_1 =  [row[0] for row in sliced_alternates_1]


            ref_seq = _result.sliced_variant_calls[:, index, 0]
            alter_seq = _result.sliced_variant_calls[:, index , 1]
            ref_seq = np.where(ref_seq == -1, 0, ref_seq)
            alter_seq = np.where(alter_seq == -1, 0, alter_seq)




            itemindex = np.where(ref_seq == 1)

            if len(itemindex[0]) != 0:

                # Get the reference nucleotides (as letters ATCG)
                sliced_reference_1 = _result.reference_allele[_result.slice_variant_calls]


                # Get the alternate nucleotides (as letters ATCG)
                sliced_alternates_1 = _result.alternate_alleles[_result.slice_variant_calls]
                sliced_alternates_1 =  [row[0] for row in sliced_alternates_1]

                ref_seq_data_1 = sliced_reference_1


                for data_index in itemindex[0]:

                    ref_seq_data_1[data_index] = sliced_alternates_1[data_index]
            else :

                # Get the reference nucleotides (as letters ATCG)
                sliced_reference_1 = _result.reference_allele[_result.slice_variant_calls]
                ref_seq_data_1 = sliced_reference_1


            itemindex = np.where(alter_seq == 1)

            if len(itemindex[0]) != 0:


                  # Get the reference nucleotides (as letters ATCG)
                sliced_reference_1 = _result.reference_allele[_result.slice_variant_calls]


                # Get the alternate nucleotides (as letters ATCG)
                sliced_alternates_1 = _result.alternate_alleles[_result.slice_variant_calls]
                sliced_alternates_1 =  [row[0] for row in sliced_alternates_1]

                ref_seq_data_2 = sliced_reference_1

                for data_index in itemindex[0]:

                    ref_seq_data_2[data_index] = sliced_alternates_1[data_index]
            
            else :
                # Get the reference nucleotides (as letters ATCG)
                sliced_reference_1 = _result.reference_allele[_result.slice_variant_calls]
                ref_seq_data_2 = sliced_reference_1


            if consensus == "R":
                sequence_data = ''.join(ref_seq_data_1)

            elif consensus == "A":
                sequence_data = ''.join(ref_seq_data_2)

            else:
                """
                sequence_data_1 = ''.join(ref_seq_data_1)
                sequence_data_2 = ''.join(ref_seq_data_2)
                sequence_data = sequence_data_1 + sequence_data_2
                """

                #https://www.bioinformatics.org/sms/iupac.html
                sequence_data = ""
                for x in range(len(ref_seq_data_1)):
                    
                    if ref_seq_data_1[x] == ref_seq_data_2[x]:
                        sequence_data += ref_seq_data_1[x]

                    elif  ref_seq_data_1[x] == "A" and ref_seq_data_2[x] == "G":
                        sequence_data += "R"

                    elif  ref_seq_data_1[x] == "C" and ref_seq_data_2[x] == "T":
                        sequence_data += "Y"

                    elif  ref_seq_data_1[x] == "G" and ref_seq_data_2[x] == "C":
                        sequence_data += "S"

                    elif  ref_seq_data_1[x] == "A" and ref_seq_data_2[x] == "T":
                        sequence_data += "W"

                    elif ref_seq_data_1[x] == "G" and ref_seq_data_2[x] == "T":
                        sequence_data += "K"
                    
                    elif ref_seq_data_1[x] == "G" and ref_seq_data_2[x] == "T":
                        sequence_data += "M"

                    elif  ref_seq_data_1[x] == "A" and ref_seq_data_2[x] == "T":
                        sequence_data += "W"

                    elif ref_seq_data_1[x] == "G" and ref_seq_data_2[x] == "T":
                        sequence_data += "K"
                    
                    elif ref_seq_data_1[x] == "G" and ref_seq_data_2[x] == "T":
                        sequence_data += "M"

                    elif  ref_seq_data_2[x] == "A" and ref_seq_data_1[x] == "G":
                        sequence_data += "R"

                    elif  ref_seq_data_2[x] == "C" and ref_seq_data_1[x] == "T":
                        sequence_data += "Y"

                    elif  ref_seq_data_2[x] == "G" and ref_seq_data_1[x] == "C":
                        sequence_data += "S"

                    elif  ref_seq_data_2[x] == "A" and ref_seq_data_1[x] == "T":
                        sequence_data += "W"

                    elif ref_seq_data_2[x] == "G" and ref_seq_data_1[x] == "T":
                        sequence_data += "K"
                    
                    elif ref_seq_data_2[x] == "G" and ref_seq_data_1[x] == "T":
                        sequence_data += "M"


                    elif  ref_seq_data_1[x] in [ "C", "G", "T"] and   ref_seq_data_2[x] in [ "C", "G", "T"]:
                        sequence_data += "B"

                    elif  ref_seq_data_1[x] in [ "A", "G", "T"] and   ref_seq_data_2[x] in  [ "A", "G", "T"]:
                        sequence_data += "D"

                    elif  ref_seq_data_1[x] in [ "A", "C", "T"] and   ref_seq_data_2[x] in  [ "A", "C", "T"]:
                        sequence_data += "H"

                    elif  ref_seq_data_1[x] in [ "A", "C", "G"] and   ref_seq_data_2[x] in  [ "A", "C", "G"]:
                        sequence_data += "V"

                    else:
                        sequence_data += "N"

            """

            my_seq = ""
            for x in range(len(ref_seq_data_1)):
                    
                print(ref_seq_data_1[x] +"  "+alter_seq_data[x])
                if ref_seq_data_1[x] == alter_seq_data[x]:
                    my_seq += ref_seq_data_1[x]

                elif ref_seq_data_1[x] in [  "A", "C"]  and  alter_seq_data[x] in [  "A", "C"]:
                    my_seq += "M"
                elif  ref_seq_data_1[x] in [  "A", "G"] and  alter_seq_data[x] in [  "A", "C"]:
                    my_seq += "R"
                elif  ref_seq_data_1[x] in [  "A", "T"] and  alter_seq_data[x] in [  "A", "C"]:
                    my_seq += "W"
                elif ref_seq_data_1[x] in [  "C", "G"] and  alter_seq_data[x] in [  "A", "C"]:
                    my_seq += "S"
                elif ref_seq_data_1[x] in [  "C", "T"] and  alter_seq_data[x] in [  "A", "C"]:
                    my_seq += "Y"
                elif ref_seq_data_1[x] in [  "G", "T"] and  alter_seq_data[x] in [  "A", "C"]:
                    my_seq += "K"
                elif ref_seq_data_1[x] in [  "A", "C", "G"] and  alter_seq_data[x] in [  "A", "C"]:
                    my_seq += "V"
                elif ref_seq_data_1[x] in [  "A", "C", "T"] and  alter_seq_data[x] in [  "A", "C"]:
                    my_seq += "H"
                elif ref_seq_data_1[x] in  [  "A", "G", "T"] and  alter_seq_data[x] in [  "A", "C"]:
                    my_seq += "D"
                elif ref_seq_data_1[x] in [ "C", "G", "T"] and  alter_seq_data[x] in [  "A", "C"]:
                    my_seq += "B"
                elif ref_seq_data_1[x] in [ "G" "A", "T" "C"] and  alter_seq_data[x] in [  "A", "C"]:
                    my_seq += "N"

            """
                
            outfile.write(">" + samples[index] +"\n")

            outfile.write(sequence_data+"\n")


'''
this code will explain how to compute the distance matrix for multiple sequence aligen data.
'''

def pairwise_dist(first_seq, second_seq):
    '''
    compute the pairwise distances between two sequences
    '''
    total_distance = 0
    
    if len(first_seq) > len(second_seq):
        total_distance =  len(first_seq) - len(second_seq)
    elif len(second_seq) > len(first_seq):
        total_distance = len(second_seq) - len(first_seq)
    
    
    for char_1, char_2 in zip(first_seq, second_seq):
        if char_1 != char_2:
            total_distance += 1
    return total_distance


def distance_matrix(dist_dict):
    '''
   Create a blank dataframe and label its column and index headings
   with all of the creatures.
    '''
    dist_matrix = DataFrame(index=dist_dict.keys(), columns=dist_dict.keys())
    fill_diagonal(dist_matrix.values, 0)
    return dist_matrix


def generate_distance_matrix(dist_matrix, pairwise_dict, start, counter, end, sequence_dict, sequence_df):
    '''
    Utilize vectorized assignments to fill in the complete distance matrix.
    Imagine allocating lists for each organism that are at right angles and have a vertex of zero.
    Setting the start position to the end position and adding a decrementing counter value to the
    end value update the start and end positions for the right angles.
    '''
    for i, j in zip(range(len(sequence_dict.keys()) - 1), range(1, len(sequence_dict.keys()))):
        sequence_df.iloc[i, j:] = list(pairwise_dict.values())[start:end]
        sequence_df.iloc[j:, i] = sequence_df.iloc[i, j:]
        start = end
        end += counter
        counter -= 1
    return dist_matrix


def generate_distance_matrix_for_genotype(_result, number_of_sample):
    '''
     input : this block of code will take the JSON dataset which contains the sample name and their genotype data 
     which is present in array format 

     output : this block would generate the distance matrix of multiple sequence 

    '''



    df = pd.DataFrame(_result.numbers_of_alternate_alleles, index= _result.samples_selected_mapped)
    df = df[:number_of_sample]
    df = df.replace(-1, 0)
    df = df.astype(str)
    columns=df.index
    sequence_dict = {}

    for x in range(len(columns)):

        sequence_dict[columns[x]] = str(''.join(df.iloc[x]))

    sequence_df = distance_matrix(sequence_dict)
    Seq_pairwise_diff = {
        key: pairwise_dist(sequence_dict[key[0]], sequence_dict[key[1]])
        for key in combinations(sequence_dict, 2)
    }

    Seq_df = generate_distance_matrix(sequence_df, Seq_pairwise_diff, 0,
                            len(sequence_dict.keys()) - 2,
                            len(sequence_dict.keys()) - 1,sequence_dict,sequence_df)
    
    return Seq_df

def generate_distance_matrix_for_sequence(all_fasta_file):
    '''
        input : this block of code will take the fattest dataset which contains the sample name and their sequence data 
            which is present in array format

        output : this block would generate the distance matrix of multiple sequence 

    '''

    sequence_list = [] # To keep order of sequence
    sequence_dict = {}
    for record in SeqIO.parse(open(all_fasta_file, "r"), "fasta"):
        sequence_list.append(record.id)
        sequence_dict[record.id] = str(record.seq)

    sequence_df = distance_matrix(sequence_dict)
    Seq_pairwise_diff = {
        key: pairwise_dist(sequence_dict[key[0]], sequence_dict[key[1]])
        for key in combinations(sequence_dict, 2)
    }

    Seq_df = generate_distance_matrix(sequence_df, Seq_pairwise_diff, 0,
                                len(sequence_dict.keys()) - 2,
                                len(sequence_dict.keys()) - 1, sequence_dict,sequence_df)
    
    return Seq_df


