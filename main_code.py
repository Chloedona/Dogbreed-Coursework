#IMPORTS#
import re

import matplotlib.pyplot as plt
import os
from Bio import Phylo
from Bio.Phylo.TreeConstruction import DistanceTreeConstructor, DistanceMatrix
import numpy as np

#I/O FUNCTIONS#
def read_fasta(file):
    """
    Reads a FASTA file and returns a dictionary
    with header as key and sequence as a value
    """
    sequences = {}

    with open(file, 'r') as f:
        name = None
        seq = []

        for line in f:
            line = line.strip()

            if line.startswith('>'):
                if name:
                    sequences[name] = ''.join(seq)
                name = line[1:]
                seq = []
            else:
                seq.append(line)
        if name:
            sequences[name] = ''.join(seq)

    return sequences

def breed_name(header):
    """
    Extracts breed name from FASTA header e.g. [breed = boxer]
    Uses re from the standard library to search for the breed name and returns it.
    """
    match = re.search(r'\[breed\s*=\s*([^\]]+)\]', header)
    if match:
        return match.group(1).strip()


#SEQUENCE ANALYSIS#
def percent_identity(seq1, seq2):
    """
    Calculates the percent identity between two sequences. Takes into consideration different lenths of sequences
    by using None for mising bases. Returns the percentage of matching bases divided by the length of the longer sequence.
    """
    s1 = str(seq1).upper()
    s2 = str(seq2).upper()
    max_len = max(len(s1), len(s2))

    if max_len == 0:
        return 0.0
    
    matches = 0
    for i in range (max_len):
        base1 = s1[i] if i < len(s1) else None
        base2 = s2[i] if i < len(s2) else None
        if base1 == base2:
            matches += 1
        
    return (matches / max_len) * 100

def percent_difference(seq1, seq2):
    """
    Returns the percentage difference between sequences by subtracting percent identity from 100.
    """
    return 100 - percent_identity(seq1, seq2)

def find_closest(mystery_seq, database):
    """
    Finds the closest breed match from the database using the mystery sequence
    Returns the breed name and similarity.
    """
    best_breed = None
    best_similarity = -1

    for breed, seq in database.items():
        similarity = percent_identity(mystery_seq, seq)
        if similarity > best_similarity:
            best_similarity = similarity
            best_breed = breed
    return best_breed, best_similarity

#STATISTICS#
def p_values(scores):
    """
    Calculates the p-value for the similarity between the mystery sequence and the closest breed match.
    """
    identities = list(scores.values())
    num_breeds = len(identities)

    if num_breeds == 0:
        return {}
    
    pvals = {}
    for breed, score in scores.items():
        num_same_or_higher = sum(1 for s in identities if s >= score)
        pvals[breed] = num_same_or_higher / num_breeds

    return pvals

#PHLOGENY#
def build_tree(sequences):
    """
    Construct phylogenetic tree from mystery DNA and database
    """
    names = list(sequences.keys())
    matrix = []

    for i, name1 in enumerate(names):
        row = []
        for j in range(i + 1):
            name2 = names[j]
            if i == j:
                row.append(0.0)
            else:
                dist = percent_difference(sequences[name1], sequences[name2]) / 100
                row.append(dist)
        matrix.append(row)

    distance_matrix = DistanceMatrix(names, matrix)
    constructor = DistanceTreeConstructor()
    return constructor.nj(distance_matrix)

def plot_tree_image(tree, output_dir):
    """
    Displays the tree to user
    """
    image_path = os.path.join(output_dir, "phylogenetic_tree.png")

    # Scale figure height with the number of tips to reduce label overlap.
    num_tips = tree.count_terminals()
    fig_height = max(12, num_tips * 0.35)
    fig, ax = plt.subplots(figsize=(24, fig_height))

    # Sorting clades to make the tree easier to visualise.
    tree.ladderize()
    Phylo.draw(tree, axes=ax, do_show=False, show_confidence=False)

    label_fontsize = 10 if num_tips > 80 else 12
    for text in ax.texts:
        text.set_fontsize(label_fontsize)

    # Adding title and adjusting layout to prevent overlap.
    ax.set_title("Phylogenetic Tree", fontsize=24, pad=12)
    plt.tight_layout()
    plt.savefig(image_path, dpi=300, bbox_inches='tight')
    plt.close(fig)

    return image_path

def write_report(output_dir, best_breed, best_similarity, pval, tree_image_path):
    """
    Writes a report summarising the results of the analysis and saves it to the output directory.
    """
    report_path = os.path.join(output_dir, "report.txt")
    with open(report_path, 'w') as out:
        out.write("\nChloe's Dog Breed Analysis Report\n")
        out.write(f"Closest breed match: {best_breed}\n")
        out.write(f"Percent Identity: {best_similarity:.2f}%\n")
        out.write(f"P-value (closest breed): {pval:.4f}\n")
        out.write("\nPhylogenetic tree outputs:\n")
        out.write(f"Phylogenetic tree image saved at: {tree_image_path}\n")

#MAIN PROGRAM#
"""
This is the main program that runs the dog breed analysis. The program will read the database and mystery sequence 
and make a new directory called "Results" saving a report and the tree image to display to the user.


"""
if __name__ == "__main__":

    output_dir = "Results"
    os.makedirs(output_dir, exist_ok=True)
    
    database = read_fasta("Data/dog_breeds.fa")
    mystery = read_fasta("Data/mystery.fa")

    mystery_name = list(mystery.keys())[0]
    mystery_seq = mystery[mystery_name]
    breed_display = {header: breed_name(header) for header in database}

    print("Finding closest breed match...")
    best_breed, best_similarity = find_closest(mystery_seq, database)
    best_breed_name = breed_name(best_breed) or best_breed

    all_results = []
    observed_scores = {}

    for db_breed, db_seq in database.items():
        identity = percent_identity(mystery_seq, db_seq)
        all_results.append((db_breed, identity))
        observed_scores[db_breed] = identity

    print("Calculating p-values...")
    pvals = p_values(observed_scores)

    print("Building phylogenetic tree...")
    sequences_for_tree = {breed_display[header]: seq for header, seq in database.items()}
    sequences_for_tree[mystery_name] = mystery_seq
    tree = build_tree(sequences_for_tree)
    tree_image_path = plot_tree_image(tree, output_dir)

    print("Writing report...")
    report_path = os.path.join(output_dir, "report.txt")
    write_report(output_dir, best_breed_name, best_similarity, pvals[best_breed], tree_image_path)

    print(f"Analysis complete. Report saved at: {report_path}")
    


