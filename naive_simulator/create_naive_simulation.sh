curl -O ftp://ftp.ensembl.org/pub/release-90/fasta/homo_sapiens/cdna/Homo_sapiens.GRCh38.cdna.all.fa.gz
gunzip Homo_sapiens.GRCh38.cdna.all.fa.gz
mkdir rawInput
mv Homo_sapiens.GRCh38.cdna.all.fa rawInput
python separate_ensembl_transcripts_into_genes.py
python simulate.py