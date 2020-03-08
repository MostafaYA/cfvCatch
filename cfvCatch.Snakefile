"""----------------------------
Friedrich-Loeffler-Institut (https://www.fli.de/), IBIZ
date: Febraury, 20, 2020
Author: Mostafa Abdel-Glil (mostafa.abdel-glil@fli.de)
-------------------------------

# ToDo:
make abricate database for the iscfe1
Assembly
* assembler_options =config["assembler_options"] ---snake variable and rule
* -opts {params.assembler} options does not wotŕk in the bash --- bash
 * shovill --outdir $OUTDIR --cpus $THREADS --R1 $READ1 --R2 $READ2 --force --ram 10 $assembler_options $assembler_options_additional

"""
import os
import tempfile
#configfile: "config.yaml"
#working_dir=config['working_dir']
snakemake_folder=config['snakemake_folder']
raw_data_dir=config['raw_data_dir']
fasta_dir=config['fasta_dir']
results_dir=config['results_dir']
temporary_todelete= results_dir + "temporary_todelete/"

#samples
SAMPLES, = glob_wildcards( raw_data_dir + "{sample}_R1.fastq.gz")
GENOMES, = glob_wildcards( fasta_dir + "{genome}.fasta")
SAMPLES = set(SAMPLES)
GENOMES = set(GENOMES)
DATASET =  SAMPLES | GENOMES
#scripts_paths
bin_dir= snakemake_folder + "bin/"
dbs_dir= snakemake_folder + "dbs/"
envs_folder= snakemake_folder + "envs/"
#mlst_parameters
if "MLST_options" in config:
    MLST_options= config["MLST_options"]
else:
    MLST_options= ""
#abricate_parameters
if "abricate_options" in config:
    abricate_options= config["abricate_options"]
else:
    abricate_options= ""
if "snippy_options" in config:
    snippy_options= config["snippy_options"]
else:
    snippy_options= ""
if "assembler_options" in config:
    assembler_options= config["assembler_options"]
else:
    assembler_options= ""

assembler= "shovill",
#assembler_options= ""
reference= snakemake_folder + "dbs/CP000487.gbk"
iscfe1= snakemake_folder + "dbs"
snippyFolder= snakemake_folder + "dbs/wgsSnps/"
pcr_primers= snakemake_folder + "dbs/primers_cfetus.txt"
in_silico_pcr= snakemake_folder + "bin/in_silico_pcr/"
snp_markers= snakemake_folder + "dbs/snp_markers.bed"

rule all:
    input:
        Assembly =expand( results_dir + "{sample}/contigs.fa", sample=DATASET),
        iscfe1_results =expand( results_dir + "{genome}/iscfe1.tab", genome=DATASET),
        contig = expand(  results_dir + "{genome}/contigs.fa", genome=DATASET),
        snippy_snps= expand( results_dir + "{sample}/snippy/snps.tab", sample=DATASET),
        pcr = expand(results_dir + "{sample}/pcr.tab", sample=DATASET),
        eleven_snps= expand(results_dir + "{sample}/cfv_snp_markers.tab", sample=DATASET),
        iscfe1_summary = results_dir + "iscfe1_summary.tab",
        #snippy_snps_assemblies= expand( results_dir + "{genome}/snippy/snps.tab", genome=DATASET),
        snippycore= results_dir + "core.aln",
        phylogeny_tree= results_dir + "core.newick",
        mlstresults= results_dir + "mlst.tab",
        snpdistsresults= results_dir + "distances.tab",
        summary=results_dir + "Summary.xls"

"""
collect assemblies
"""
rule collect_assemblies:
    input:
        assemblies = fasta_dir + "{genome}.fasta",
    output:
        contig = results_dir + "{genome}/contigs.fa",
    conda:
        envs_folder + "bioawk.yaml" #spades, sickle
    shell:
        "bash {bin_dir}fa_rename.sh {input.assemblies} {output.contig} " #ln -s -f ##* Prokka does not like contigs ID > 37, '--centre X --compliant' is not the way always
"""
assembly
"""
rule assembly:
    input:
        r1 = raw_data_dir + "{sample}_R1.fastq.gz",
        r2 = raw_data_dir + "{sample}_R2.fastq.gz",
    output:
        contig = results_dir + "{sample}/contigs.fa",
    threads: 16 #increasing threads produces errors with spades
    conda:
        envs_folder + "shovill.yaml" #spades, sickle
    params:
        assembler = assembler_options,
        assembly_dir = directory (results_dir + "{sample}")
    shell:
        "bash {bin_dir}assembly.sh -a {assembler} --R1 {input.r1} --R2 {input.r2} -o {output.contig} -d {params.assembly_dir}/{assembler}  -p {threads} {assembler_options} " #-t trimmomatic -f true

"""
reference
"""
rule reference:
    output:
        ref_fasta= results_dir + "ref.fa",
    conda:
        envs_folder + "reference.yaml"
    shell:
        "seqret -auto -filter -osformat2 fasta {reference} > {output.ref_fasta} && \
        samtools faidx {output.ref_fasta}"

"""
iscfe1
"""
rule iscfe1:
    input:
        contig = results_dir + "{sample}/contigs.fa",
    output:
        iscfe1 = results_dir + "{sample}/iscfe1.tab",
    threads: 16
    conda:
        envs_folder + "abricate.yaml"
    params:
        iscfe1 = iscfe1,
    shell:
        "abricate --nopath --db iscfe1 --datadir {params.iscfe1} --threads {threads} {abricate_options} {input.contig} > {output.iscfe1}"

"""
iscfe1_summary
"""
rule iscfe1_summary:
    input:
        iscfe1_all = expand( results_dir + "{sample}/iscfe1.tab", sample=DATASET),
    output:
        iscfe1_summary = results_dir + "iscfe1_summary.tab",
    threads: 16
    conda:
        envs_folder + "abricate.yaml"
    shell:
        """var={results_dir} && abricate --nopath --summary {input.iscfe1_all}  | sed "s#$var##g" | sed "s#/iscfe1.tab##g" | tee -a {output.iscfe1_summary}"""

"""
pcr
"""
rule pcr:
    input:
        contig = results_dir + "{sample}/contigs.fa",
    output:
        pcr = results_dir + "{sample}/pcr.tab",
        pcr_amplicon = results_dir + "{sample}/pcr_amplicon.fasta",
    threads: 16
    #conda:
    #    envs_folder + "abricate.yaml"
    params:
        primers = pcr_primers,
    shell:
        "perl {in_silico_pcr}in_silico_PCR.pl -s {input.contig} -p {params.primers} > {output.pcr} 2> {output.pcr_amplicon}"


"""
snippy
"""
#export LD_LIBRARY_PATH=/usr/local/lib:/usr/lib:/usr/local/lib64:/usr/lib64
rule snippy:
    input:
        r1 = raw_data_dir + "{sample}_R1.fastq.gz",
        r2 = raw_data_dir + "{sample}_R2.fastq.gz",
    output:
        snippy_snps= results_dir + "{sample}/snippy/snps.tab",
        eleven_snps= results_dir + "{sample}/cfv_snp_markers.tab",
        #snippy_folder=  results_dir + "{sample}",
    threads: 32
    conda: envs_folder + "snippy.yaml"
    params:
        snippy_outdir= results_dir + "{sample}/snippy",
        #snippy_outdir_cp= results_dir + "{sample}/",
    shell:
        """
        snippy --force  --cpus {threads} --ram {threads} --outdir {params.snippy_outdir} --ref {reference} --R1 {input.r1} --R2 {input.r2} {snippy_options} 2>&1 | sed 's/^/[snippy] /' && \
        samtools mpileup -ugB -I -l {snp_markers}  -f {params.snippy_outdir}/ref.fa {params.snippy_outdir}/snps.bam | bcftools call -c  | bcftools query -f '%CHROM\\t%POS\\t%REF\\t%ALT\\n' --print-header | tee -a {output.eleven_snps}
        """
rule snippy_ln:
    input:
        snippy_snps= results_dir + "{sample}/snippy/snps.tab",
    output:
        snippy_snps_tab=  results_dir + "{sample}/snps.tab",
        snippy_snps_aligned= results_dir + "{sample}/snps.aligned.fa",
        #snippy_snps_rawvcf= results_dir + "{sample}/snps.raw.vcf",
        snippy_snps_vcf= results_dir + "{sample}/snps.vcf",
        #snippy_snps_bam= results_dir + "{sample}/snps.bam",
        #snippy_snps_bai= results_dir + "{sample}/snps.bam.bai",
        snippy_snps_log= results_dir + "{sample}/snps.log",
    conda: envs_folder + "snippy.yaml"
    params:
        snippy_outdir= results_dir + "{sample}/snippy",
    shell:
        "ln -s -f $(realpath {params.snippy_outdir}/snps.tab) {output.snippy_snps_tab}; ln -s -f $(realpath {params.snippy_outdir}/snps.aligned.fa) {output.snippy_snps_aligned} ;  ln -s -f $(realpath {params.snippy_outdir}/snps.vcf) {output.snippy_snps_vcf};  ln -s -f $(realpath {params.snippy_outdir}/snps.log) {output.snippy_snps_log}" #ln -s -f $(realpath {params.snippy_outdir}/snps.bam) {output.snippy_snps_bam}; ln -s -f $(realpath {params.snippy_outdir}/snps.bam.bai) {output.snippy_snps_bai}; ln -s -f $(realpath {params.snippy_outdir}/snps.raw.vcf) {output.snippy_snps_rawvcf};
rule snippy_assemblies:
    input:
        contig = fasta_dir + "{genome}.fasta"
    output:
        snippy_snps= results_dir + "{genome}/snippy/snps.tab",
        eleven_snps= results_dir + "{genome}/cfv_snp_markers.tab",
    threads: 32
    conda: envs_folder + "snippy.yaml"
    params:
        snippy_outdir= results_dir + "{genome}/snippy",
    shell:
        """
        snippy --force  --cpus {threads} --ram {threads} --outdir {params.snippy_outdir} --ref {reference} --ctgs {input.contig} {snippy_options} 2>&1 | sed 's/^/[snippy] /' && \
        samtools mpileup -ugB -I -l {snp_markers}  -f {params.snippy_outdir}/ref.fa {params.snippy_outdir}/snps.bam | bcftools call -c  | bcftools query -f '%CHROM\\t%POS\\t%REF\\t%ALT\\n' --print-header | tee -a {output.eleven_snps}
        """

def snippy_folders(wildcards):
    return expand(results_dir + "{sample}", sample=SAMPLES)
def snippy_folders_assemblies(wildcards):
    return expand(results_dir + "{genome}", genome=GENOMES)

rule snippy_core:
    input:
        #snippy_folders=  expand( results_dir + "{sample}", sample=SAMPLES),
        #snippy_folders_assemblies=  expand( results_dir + "{genome}", genome=GENOMES),
        a=expand( results_dir + "{sample}/snps.tab", sample=SAMPLES),
        b=expand( results_dir + "{sample}/snps.aligned.fa", sample=SAMPLES),
        #c=expand( results_dir + "{sample}/snps.raw.vcf", sample=SAMPLES),
        d=expand( results_dir + "{sample}/snps.vcf", sample=SAMPLES),
        #e=expand( results_dir + "{sample}/snps.bam", sample=SAMPLES),
        #f=expand( results_dir + "{sample}/snps.bam.bai", sample=SAMPLES),
        g=expand( results_dir + "{sample}/snps.log", sample=SAMPLES),
        h=expand( results_dir + "{genome}/snps.tab", genome=GENOMES),
        i=expand( results_dir + "{genome}/snps.aligned.fa", genome=GENOMES),
        #j=expand( results_dir + "{genome}/snps.raw.vcf", genome=GENOMES),
        k=expand( results_dir + "{genome}/snps.vcf", genome=GENOMES),
        #l=expand( results_dir + "{genome}/snps.bam", genome=GENOMES),
        #m=expand( results_dir + "{genome}/snps.bam.bai", genome=GENOMES),
        n=expand( results_dir + "{genome}/snps.log", genome=GENOMES),
    output:
        snippycore= results_dir + "core.aln",
    conda: envs_folder + "snippy.yaml"
    params:
        snippycoreoutdir= results_dir + "core",
        snippy_folders = snippy_folders,
        snippy_folders_assemblies = snippy_folders_assemblies,
        snippyFolder = snippyFolder,
        #snippy_outdir_cp= results_dir + "{sample}/",
    shell:
        """ var={results_dir}Reference && snippy-core  --ref {reference} {params.snippy_folders} {params.snippyFolder}/* $(echo {params.snippy_folders_assemblies} | sed "s#$var##g") --prefix {params.snippycoreoutdir}  2>&1 | sed 's/^/[snippy-core] /'"""
"""
Core genome phylogeny using FastTree
"""
rule FastTree: #run fasttree
    input:
        snippycore= results_dir + "core.aln",
    output:
        phylogeny_tree= results_dir + "core.newick",
        phylogeny_svg= results_dir + "core.svg"
    #threads: 64
    #benchmark: benchmarks_folder + "FastTree.txt"
    conda: envs_folder + "FastTree.yaml" #needs revision
    shell:
        "FastTree -nt -gtr {input.snippycore} > {output.phylogeny_tree} && \
        nw_display -S -s -w 1024 -l 'font-size:12;font-family:sans-serif;' -i 'opacity:0' -b 'opacity:0' -v 16 {output.phylogeny_tree}  > {output.phylogeny_svg}"

"""
mlst
"""
rule MLST:
    input:
        contigs= expand(results_dir + "{sample}/contigs.fa", sample=SAMPLES),
        contigs_assemblies= expand(results_dir + "{genome}/contigs.fa", genome=GENOMES),
        #ref_fasta= results_dir + "ref.fa",
    output:
        mlstresults= results_dir + "mlst.tab",
    conda:
        envs_folder + "mlst.yaml"
    shell:
        """var={results_dir} && mlst --quiet {MLST_options} {input.contigs} {input.contigs_assemblies} | sed "s#$var##g" | sed "s#/contigs.fa##g" | tee -a {output.mlstresults}"""   #{input.ref_fasta}
"""
snpdists
"""
rule snpdists:
    input:
        snippycore= results_dir + "core.aln",
    output:
        snpdistsresults= results_dir + "distances.tab",
    conda:
        envs_folder + "snpdists.yaml"
    shell:
        "snp-dists -b {input.snippycore} > {output.snpdistsresults}"

"""
combines results
outpur CSV, XLS
"""
rule combine_Results:
    input:
        mlstresults= results_dir + "mlst.tab",
        pcr = expand(results_dir + "{sample}/pcr.tab", sample=SAMPLES),
        iscfe1_summary = results_dir + "iscfe1_summary.tab",
        eleven_snps= expand(results_dir + "{sample}/cfv_snp_markers.tab", sample=SAMPLES),
    output:
        allxls= results_dir + "Summary.xls",
        allcsv= results_dir + "Summary.csv",
    conda:
        envs_folder + "rxls.yaml"
    script:
        bin_dir +"/summary.R"
