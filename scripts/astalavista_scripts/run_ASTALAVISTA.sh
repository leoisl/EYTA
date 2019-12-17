#see http://sammeth.net/confluence/pages/viewpage.action?pageId=5177667
#-t asta - use astalavista tool: asta - AStalavista event retriever
#-c - used to get the sequence of the events
#-e [ASI] - used to get only Alternative Splicing Internal - AS events that are comprised between 2 constitutive exons
#-d 2 - dimension of the events - 2 = pairwise events
#-i is the annotation
#-a [SEQ] SEQ: output splice site sequences of event flanks
set -eu
FILES=Homo_sapiens.GRCh38.90.gtf_divided_by_gene_id/*
for f in $FILES
do
  echo "Running ASTALAVIST on $f ..."
  /data2/leandro/EYTA/ASTALAVISTA/astalavista-4.0/bin/astalavista -t asta -c /data2/leandro/EYTA/ASTALAVISTA/ASTALAVISTA_run_on_GRCh38.p10_Ensembl90/GRCh38/ -e [ASI] -d 2 -i $f -a [SEQ] --threads 1
done
