#!/bin/bash

help()
{
	cat <<- EOF
	Desc: humann3 analysis metagenome
	Usage: bash run_humann.sh -h
	Author: wnse
	Parameters:
		
		-h --help       	Print the help info and exit.

		-v --version    	Print the version info.

		-i1 --infastq    	The input fastq R1.

		-i2 --infastq    	The input fastq R2.

		-o --out			The output files prefix.

		-t --threads    	The threads for processing.

		-n --nt_db			The humann nucleotide-database directory.

		-p --pro_db			The humann protein-database directory.

		-m --mpa_db			The metaphlan bowtie2db directory.

		-r --rgi_db_json	The rgi card-db json file.

		-vf --vf_db_pro		The virus factor db protein diamond(dmnd) file.

	License: none
	
EOF
	exit 0
}

version_info()
{
	cat <<- EOF

		humann3 metagenome 20230303

	EOF
		exit 0
}


while true;do
	case "$1" in
		-i1|--infq1)
			infq1=$2;shift 2;;
		-i2|--infq2)
			infq2=$2;shift 2;;
		-o|--out)
			out=$2;shift 2;;
		-t|--threads)
			threads=$2;shift 2;;
		-n|--nt_db)
			nt_db=$2;shift 2;;
		-p|--pro_db)
			pro_db=$2;shift 2;;
		-m|--mpa_db)
			mpa_db=$2;shift 2;;
		-r|--rgi_db_json)
			rgi_db_json=$2;shift 2;;
		-vf|--vf_db_pro)
			vf_db_pro=$2;shift 2;;
		-h|--help)
			help;;
		-v|--version)
			version_info;;
		*)
			break;
	esac
done


source /opt/conda/bin/activate 
conda activate humann

fastp -w ${threads} -i ${infq1} -I ${infq2} -o ${out}/clean_1.fq.gz -O ${out}/clean_2.fq.gz -j ${out}/fastp.json -h ${out}/fastp.html
cat ${out}/clean_1.fq.gz ${out}/clean_2.fq.gz > ${out}/clean.fq.gz
humann --threads ${threads} --bypass-translated-search --metaphlan-options "--bowtie2db ${mpa_db}" --nucleotide-database ${nt_db} --protein-database {pro_db} --input ${out}/clean.fq.gz --output ${out}/humann
humann_split_stratified_table --input ${out}/humann/clean_pathabundance.tsv --output ${out}/humann/

conda activate megahit
megahit -t ${threads} -1 ${out}/clean_1.fq.gz -2 ${out}/clean_2.fq.gz -o ${out}/megahit_out --min-contig-len 1500


conda activate rgi
rgi load -i ${rgi_db_json} --local
rgi main -n ${threads} --input_sequence ${out}/megahit_out/final.contigs.fa --output_file ${out}/card --local -a DIAMOND

bowtie2-build --threads ${threads} ${out}/final.contigs.fa.temp.contigToORF.fsa ${out}/ORF-index
bowtie2 -p ${threads} -x ${out}/ORF-index -1 ${out}/clean_1.fq.gz -2 ${out}/clean_2.fq.gz -S ${out}/raw2orf.sam
samtools view -bS ${out}/raw2orf.sam 
samtools sort -@ ${threads} -o ${out}/raw2orf.bam ${out}/raw2orf.sam 
samtools index ${out}/raw2orf.bam
samtools idxstats ${out}/raw2orf.bam > ${out}/raw2orf.bam.idxstats

conda activate humann
diamond blastx --threads ${threads} --db ${vf_db_pro} --query ${out}/final.contigs.fa.temp.contigToORF.fsa --out ${out}/vfdb.diamond.txt --outfmt 6 --evalue 1e-5



