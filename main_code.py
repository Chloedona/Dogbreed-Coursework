#IMPORTS#
import matplotlib.pylot as plt




#I/O FUNCTIONS#
def read_fasta(file)
    """
    Reads a FASTA file
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





#PHLOGENY#




#MAIN PROGRAM#
if __name__ == "__main__":

    database = read_fasta("dog_breeds.fa")
    mystery = read_fasta("mystery.fa")


