#IMPORTS#
import matplotlib.pylot as plt




#I/O FUNCTIONS#
def read_fasta(file)
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

#SEQUENCE ANALYSIS#
def percent_identity(seq1, seq2):
    """
    Calculate the percent identity between two sequences
    """
    return
def percent_difference(seq1, seq2):
    """
    Return the percentage difference between sequences
    """
    return
def find_closest(mystery_seq, database):
    """
    Find the closest breed match from the database using the myserty sequence.
    Return the breed name and similarity.
    """
    return 

#STATISTICS#
def p_value():
    return

#PHLOGENY#
def build_tree():
    """
    Construct phylogenetic tree from mystery DNA and database
    """

def plot_tree():
    """
    Displays the tree to user
    """




#MAIN PROGRAM#
if __name__ == "__main__":

    database = read_fasta("dog_breeds.fa")
    mystery = read_fasta("mystery.fa")


