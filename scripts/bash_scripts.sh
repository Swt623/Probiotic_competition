sample_list=$(ls *.gz | awk -F\_ '{print $1}' | uniq)

for i in $sample_list
do 
 	cat /projects/p31414/$i\_S*_L00{5,6,7,8}_R1_001.fastq.gz > /projects/p31414/RNA/$i\.fastq.gz
done  

# bowtie version 1.2.2    
cat *gz > neg_c.fastq.gz
gunzip neg_c.fastq.gz
fastq_to_fasta -i neg_c.fastq -o neg_c.fasta
bowtie2-build neg_c.fasta neg_ref

# kneaddata

# trim-galore
trim_galore -q 20  --phred33 --fastqc  --gzip *.fastq -j 12    

#sortmerna
sortmerna -ref /home/jhk7793/Application/sortmerna_db/silva-bac-16s-id90.fasta \
			-ref /home/jhk7793/Application/sortmerna_db/silva-bac-23s-id98.fasta \
			-ref /home/jhk7793/Application/sortmerna_db/silva-euk-18s-id95.fasta \
			-ref /home/jhk7793/Application/sortmerna_db/silva-euk-28s-id98.fasta \
			-ref /home/jhk7793/Application/sortmerna_db/silva-arc-16s-id95.fasta \
			-ref /home/jhk7793/Application/sortmerna_db/silva-arc-23s-id98.fasta \
	-reads ./AA24-1_kneaddata_trimmed.fq.gz \
	--threads 12 \
	--fastx \
	-workdir ./sortmerna/AA24-1 \
	--other AA24-1_kneaddata_trimmed_nonrna \
	--aligned AA24-1_kneaddata_trimmed_rna

for i in $(ls *fq.gz| awk -F\. '{print $1}')
do
	sortmerna -ref /projects/p31414/BG/silva/smr_v4.3_default_db.fasta \
				-ref /projects/p31414/BG/silva/smr_v4.3_sensitive_db.fasta \
				-ref /projects/p31414/BG/silva/smr_v4.3_fast_db.fasta \
		-reads ./$i --threads 12 --fastx -workdir ./sortmerna/$i/ \
		--other $i\_nonrna --aligned $i\_rna         
done                                                

#fastqc

#number of reads and bases
for i in $(ls *.fq.gz)   
do
	echo $i    
	awk 'BEGIN{sum=0;}{if(NR%4==2){sum+=length($0);}}END{print sum;}' $i                                                  
	done                       

for F in *.fastq ; do echo -n "$F :" && awk 'NR%4 == 2 {N+=length($0);} END { printf("%d\n",N);}' $F ; done              


echo $(cat *.fastq | wc -l)/4|bc  


# DNA sequence assemble
module load spades
spades.py -1 /projects/b1042/HartmannLab/probiotic_competition_2021/ABBL_CRE/combined/kneaddata-CRE231/CRE231_new_R1_kneaddata_paired_1.fastq -2 /projects/b1042/HartmannLab/probiotic_competition_2021/ABBL_CRE/combined/kneaddata-CRE231/CRE231_new_R1_kneaddata_paired_2.fastq -o CRE231-assembl

# Quast: check assembly quality
module load python-anaconda3/2019.10
module load blast/2.7.1
export PATH=$PATH:/software/blast/ncbi-blast-2.7.1+/bin
python /projects/p31421/software/quast-5.0.2/quast.py -o ./scaffold_quast -t 12 ./CRE231-scaffolds.fasta

# prokka annotation
module load python-anaconda3/2019.10
source activate prokka-env
prokka /projects/b1042/HartmannLab/probiotic_competition_2021/assembly/CRE231-scaffolds.fasta --outdir /projects/b1042/HartmannLab/probiotic_competition_2021/Annotation/prokka_1_14_6 --prefix CRE231_scaffolds --cpus 12

# mapping
bowtie2 --threads 24 -x ./$BASENAME\_mapping.contigs \ 
		-1 ../../kneaddata/$BASENAME/$BASENAME\_*_kneaddata_paired_1.fastq \ 
		-2 ../../kneaddata/$BASENAME/$BASENAME\_*_kneaddata_paired_2.fastq \ 
		-U ../../kneaddata/$BASENAME/$BASENAME\_*_kneaddata_unmatched_1.fastq \ 
		-U ../../kneaddata/$BASENAME/$BASENAME\_*_kneaddata_unmatched_2.fastq \ 
		-S $BASENAME.bowtie2.sam 

# sort and index sam files 
for i in $(ls *.sam | awk -F\_ '{print $1}' | uniq);  do samtools view -@ 24 -bS $i\_hisat2_ABBL18.sam | samtools sort -@ 24 -o ./sorted_bam/$i\.sorted.bam; done  

# convert gff from prokka to gtf
grep -v "#" $infile | grep "ID=" | cut -f1 -d ';' | sed 's/ID=//g' | cut -f1,4,5,7,9 |  awk -v OFS='\t' '{print $1,"PROKKA","CDS",$2,$3,".",$4,".","gene_id " $5}'

#hisat2
module load hisat2
hisat2-build #use fna file from prokka 
hisat2 -x $INDEX -p 24 -U $READS -S $SAMPLE.sam

#featureCounts
source activate prokka
featureCounts --verbose -a ../PROKKA_06012021.gtf -t CDS -o final_featureCounts.txt ./*.bam -F GTF -g gene_id   

# CARD (RGI main)
module load python-anaconda3/2019.10
source activate rgi

rgi main --input_sequence /projects/b1042/HartmannLab/probiotic_competition_2021/Annotation/prokka_1_14_6/CRE231_scaffolds.ffn --output_file /projects/b1042/HartmannLab/probiotic_competition_2021/CARD/CRE231_CARD --input_type contig --num_threads 24 --local --clean
