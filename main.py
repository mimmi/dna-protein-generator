from dna_toolkit import *
from utilities import *

virus_name = "SARS-CoV-2"

f = open(f'{virus_name}.dna',"r")
dna_seq = f.read().upper()

if not validate_seq(dna_seq):
    raise Exception('Invalid DNA Seq')

print(f'Analyzing Genome of Severe acute respiratory syndrome coronavirus 2 isolate Wuhan-Hu-1')
print(f'[1] + Sequence Length: {len(dna_seq)}')
print(f'[2] + Nucleotide Frequency: {nucleotide_frequency(dna_seq)}')
rna = transcription(dna_seq)
print(f'[3] + RNA Transcribed to Memory')
dna = dna_seq
print(f'[4] + Stored DNA (5\' 3\') to Memory. Length({len(dna)})')
reverse_complement_dna = reverse_complement(dna)
print(f'[5] + Stored Reverse Complement DNA (3\' 5\') to Memory. Length({len(reverse_complement_dna)})')
print(f'[6] + GC Content: {gc_content(dna)}%')
print(f'[7] + GC Content: {gc_content_subsec(dna, 2000)}')
amino_acid_sequences = gen_reading_frames(dna)
print(f'[8] + Translated and Stored Amino Acid Sequence in Memory.')

proteins = []
for amino_acid_sequence in amino_acid_sequences:
    proteins_in_sequence = proteins_from_rf(amino_acid_sequence)
    for protein in proteins_in_sequence:
        if len(protein) > 44:
            proteins.append(protein)

with open(f'{virus_name}.proteins', mode='wt', encoding='utf-8') as protein_store:
    protein_store.write('\n'.join(proteins))