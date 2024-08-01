from Bio import SeqIO

fasta_sequences = SeqIO.parse(open('/home/vit/Projects/cryptobench/data/G-download-uniprot-sequences-and-create-annotations/ahoj-v2/sequences.fasta'),'fasta')
for fasta in fasta_sequences:
    name, sequence = fasta.id.split('|')[1], str(fasta.seq)
    with open(f'/home/vit/Projects/cryptobench/data/G-download-uniprot-sequences-and-create-annotations/ahoj-v2/sequences/{name}.txt', 'w') as out_file:
        out_file.write(sequence)