from Bio import SeqIO

fasta_file = "JoinFullClean.fa"
result_file = "organelle.fa"

wanted = set(["cpDNA","mtDNA"])

fasta_sequences = SeqIO.parse(open(fasta_file),'fasta')
with open(result_file, "w") as f:
    for seq in fasta_sequences:
        if seq.id in wanted:
            SeqIO.write([seq], f, "fasta")
