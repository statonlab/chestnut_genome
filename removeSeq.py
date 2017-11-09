from Bio import SeqIO

fasta_file = "MaskedJoinFullClean.fa"
result_file = "ContamMaskedContigs.fa"

wanted = set(["contig0011799","cpDNA","mtDNA"])

fasta_sequences = SeqIO.parse(open(fasta_file),'fasta')
with open(result_file, "w") as f:
    for seq in fasta_sequences:
        if seq.id not in wanted:
            SeqIO.write([seq], f, "fasta")
