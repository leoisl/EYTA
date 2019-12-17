#map several stuff with STAR in order to help the debug 
set -eux
files=(Homo_sapiens.GRCh38.90.gtf_astalavista.gtf.all_pairwise_events.fa Homo_sapiens.GRCh38.90.gtf.shallow_transcriptome_0.1.fa ASTALAVISTABubblesNotFoundByEYTA EYTABubblesChimeric EYTABubblesThatIReallyCantExplain ASTALAVISTABubblesFoundByEYTA ASTALAVISTABubblesToBeFound)
nbThreads=8

for i in "${files[@]}"
do
   :
   python fix_headers.py $i
   /data2/leandro/repeats_on_transcriptome/STAR_2.5.3a/STAR-2.5.3a/bin/Linux_x86_64_static/STARlong --runThreadN $nbThreads --genomeDir Star_index_chr1/ --readFilesIn ${i}.fixed_header.fa --outFileNamePrefix ${i}_mapped --outSAMtype BAM SortedByCoordinate
   samtools index ${i}_mappedAligned.sortedByCoord.out.bam
done

rm *mappedLog.out
rm *mappedLog.progress.out
rm *mappedSJ.out.tab
