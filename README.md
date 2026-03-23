# Dogbreed-Coursework
This repository is for my Biocomputing coursework (MSc Bioinformatics)
This code is designed to:

Identify the most similar sequence
• Context: you want to develop a DNA identification service for dog breeds
• Identify the closest sequence in the database to the provided sequence 
• Input (Data folder): sequence database (dog_breeds.fa), test sequence ('mystery.fa')
• Outputs (Results folder): the closest sequence, and the difference with stretch goals:
Stretch goal 1: Probabilities across database, p-value 
Atretch goal 2: reconstructed phylogeny

There are multiple folders within this repository as follows:
- Data provides the input data (mystery sequence and database for sequence analysis)
- Results. This is where the output of the code will save as 'report.txt' and 'pyhylogenetic_tree.png'
- Tests. This is internal testing I used for testing my code
- Main_code.py is outside of any folder, but is the main code used for running the analysis pipeline.
- Each subfolder had a readme file to explain the contents of each folder.

This code uses packages:
BioPython for sequence analysis