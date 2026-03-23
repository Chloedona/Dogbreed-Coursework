#IMPORTS#
# This section imports the libraries needed for the program.
# re is used for regular expression matching to extract breed names from FASTA headers.
# matplotlib is used for plotting the phylogenetic tree.
# os is used for file and directory operations.
# Bio is used for phylogenetic tree construction and manipulation.
import re
import matplotlib.pyplot as plt
import os
from Bio import Phylo
from Bio.Phylo.TreeConstruction import DistanceTreeConstructor, DistanceMatrix

#This I/O section contains functions for reading FASTA files and extracting breed names from headers.
#I/O FUNCTIONS#
def read_fasta(file):
    """
    Reads a FASTA file and returns a dictionary
    with header as key and sequence as a value
    """
    # Initialize an empty dictionary to store sequences
    sequences = {}
    # Open the FASTA file for reading, iterate through each line and append sequences to the dictionary
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
    # Use regular expression to find the breed name in the header and return it if found stripping any blanks.
    match = re.search(r'\[breed\s*=\s*([^\]]+)\]', header)
    if match:
        return match.group(1).strip()

#This section contains functions for sequence analysis, including percent identity, difference 
#and finding the closest breed match.
#SEQUENCE ANALYSIS#
def percent_identity(seq1, seq2):
    """
    Calculates the percent identity between two sequences. Takes into consideration different lenths of sequences
    by using None for mising bases. Returns the percentage of matching bases divided by the length of the longer sequence.
    """
    # Convert sequences to uppercase and calculate the max length
    s1 = str(seq1).upper()
    s2 = str(seq2).upper()
    max_len = max(len(s1), len(s2))
    # If both sequences are empty, return 0% identity to avoid division by zero.
    if max_len == 0:
        return 0.0
    # Iterate through both sequences up to the longer sequence, counting matches and adding to the counter.
    matches = 0
    for i in range (max_len):
        base1 = s1[i] if i < len(s1) else None
        base2 = s2[i] if i < len(s2) else None
        if base1 == base2:
            matches += 1
    # Return the value by dividing the number of matches by the length of the longer sequence and multiplying by 100 to get a percentage.
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
    # Set initial values for the best breed and similarity score. 
    # Iterate through the database, calculating the percent identity using previous function.
    # If the similarity is higher than the current best similarity, update the best breed and similarity score.
    best_breed = None
    best_similarity = -1

    for breed, seq in database.items():
        similarity = percent_identity(mystery_seq, seq)
        if similarity > best_similarity:
            best_similarity = similarity
            best_breed = breed
    # Return the best breed and similarity score found in the database.
    return best_breed, best_similarity

#This section contains functions for statistical analysis
#STATISTICS#
def p_values(scores):
    """
    Calculates the p-value for the similarity between the mystery sequence and the closest breed match.
    """
    # initialise the list of identities and number of breeds in the database. If there are no breeds, return an empty dictionary.
    identities = list(scores.values())
    num_breeds = len(identities)

    if num_breeds == 0:
        return {}
    # Calculate the p-value for each breed by counting how many breeds have the same or higher similarity score
    # and dividing by the total number of breeds.
    pvals = {}
    for breed, score in scores.items():
        num_same_or_higher = sum(1 for s in identities if s >= score)
        pvals[breed] = num_same_or_higher / num_breeds

    return pvals

#This section contains functions for phylogenetic tree construction and plotting, as well as writing the final report.
#PHLOGENY#
def build_tree(sequences):
    """
    Construct phylogenetic tree from mystery DNA and database
    """
    # Initialise the names and distance matrix for the tree construction.
    names = list(sequences.keys())
    matrix = []
    # Iterate through the names of the sequences, calculating the percent difference between each pair of sequences
    # and filling in the distance matrix. If the sequences are the same, the distance is 0.0. 
    # Otherwise, the distance is calculated as the percent difference divided by 100 to convert it to a value between 0 and 1.
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
    # Use neighbour joining method from Bio.Phylo to construct the tree from the distance matrix and return it.
    distance_matrix = DistanceMatrix(names, matrix)
    constructor = DistanceTreeConstructor()
    return constructor.nj(distance_matrix)

def plot_tree_image(tree, output_dir):
    """
    Displays the tree to user
    """
    # Define the path to save the tree image in the output directory.
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
    # Define the path to save the report in the output directory and write the report with the closest breed match, percent identity,
    #  p-value and tree image path.
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
# Use if __name__ == "__main__" to ensure that the code only runs when the script is executed directly.
if __name__ == "__main__":
# Make an output directory to display the report and tree image
# If the directory already exists, it will not raise an error due to exist_ok=True.
    output_dir = "Results"
    os.makedirs(output_dir, exist_ok=True)
# Read the database and mystery sequence from the specified FASTA files.
    database = read_fasta("Data/dog_breeds.fa")
    mystery = read_fasta("Data/mystery.fa")
# Extract the mystery sequence and its name from the mystery dictionary.
# Create a display dictionary to map headers to breed names for easier interpretation of results.
    mystery_name = list(mystery.keys())[0]
    mystery_seq = mystery[mystery_name]
    breed_display = {header: breed_name(header) for header in database}
# Find the closest breed match to the mystery sequence and calculate the percent identity.
    print("Finding closest breed match...") # Prints as calculating to give user feedback on progress.
    best_breed, best_similarity = find_closest(mystery_seq, database)
    best_breed_name = breed_name(best_breed) or best_breed
# initialise lists and dictionaries to store results for p-value calculation and tree construction.
    all_results = []
    observed_scores = {}
# Iterate through the database, calculating the percent identity between the mystery sequence and each breed sequence,
# and storing the results in the all_results list and observed_scores dictionary for later use in p-value calculation and tree construction.
    for db_breed, db_seq in database.items():
        identity = percent_identity(mystery_seq, db_seq)
        all_results.append((db_breed, identity))
        observed_scores[db_breed] = identity

    print("Calculating p-values...")
    pvals = p_values(observed_scores)
# Build the phylogenetic tree using the sequences from the database and the mystery sequence
# and save the tree image to the output directory.
    print("Building phylogenetic tree...")
    sequences_for_tree = {breed_display[header]: seq for header, seq in database.items()}
    sequences_for_tree[mystery_name] = mystery_seq
    tree = build_tree(sequences_for_tree)
    tree_image_path = plot_tree_image(tree, output_dir)
# Output the report to the results folder with the closest breed match, percent identity, p-value and tree image path.
    print("Writing report...")
    report_path = os.path.join(output_dir, "report.txt")
    write_report(output_dir, best_breed_name, best_similarity, pvals[best_breed], tree_image_path)
# Display message to the user that analysis is complete.
    print(f"Analysis complete. Report saved at: {report_path}")
    


