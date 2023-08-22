import os
import sys
import shutil
import re
import csv
import Bio.SeqIO
import subprocess
from snakemake.utils import min_version
from Bio.Seq import Seq
import pandas as pd

min_version("5.4.1")

configfile: "config.yaml"

rule all:
    input:
        "transcriptome.fasta",
        "transcriptome.gene_trans_map",
        "transcriptome_stats.txt",
        "transcriptome_exN50.tsv",
        "transcriptome_expression_isoform.tsv",
        "ExN50_plot.pdf",
        "fastqc_before_trim",
        "multiqc_before_trim.txt",
        "samples_trimmed.txt",
        "fastqc_after_trim",
        "multiqc_after_trim.txt",
        "trinity_out_dir/recursive_trinity.cmds.completed",
        "busco_report.txt"
        # "raw_fastqc_outputs/summary.txt",
        # "trimmed_fastqc_outputs/summary.txt",
        # "multiqc_reports/multiqc_raw_data_report.html",
        # "multiqc_reports/multiqc_trimmed_data_report.html"

rule fastqc_before_trim:
  """
  Runs fastQC on individual input reads files.
  """
  input:
    samples=config["samples_file"]
  output:
    directory("fastqc_before_trim")
  log:
    "logs/fastqc_before_trim.log"
  conda:
    "envs/qc.yaml"
  params:
    memory="4"
  threads:
    1
  shell:
    """
    # run fastqc on input files
    # making directory for the fastqc summary files
    mkdir {output} &> {log}
    FILES=$(awk '{{ printf "%s\\n%s\\n", $3,$4}}' {input}) &>> {log}
    # running fasqc on the files
    for file in $FILES; do fastqc -f fastq -o {output} $file; done &>> {log}
    """


rule multiqc_before_trim:
  """
  Creates multiQC report from individual fastQC reports.
  """
  input:
    "fastqc_before_trim"
  output:
    out_dir=directory("multiqc_before_trim"),
    report="multiqc_before_trim.txt"
  log:
    "logs/multiqc_before_trim.log"
  conda:
    "envs/qc.yaml"
  params:
    memory="4"
  threads:
    1
  shell:
    """
    # common error resolved by those two export commands
    export LC_ALL=C.UTF-8
    export LANG=C.UTF-8
    mkdir {output.out_dir} &> {log}
    multiqc -o {output.out_dir} {input} &>> {log}
    cp multiqc_before_trim/multiqc_data/multiqc_fastqc.txt {output.report} &>> {log}
    """

checkpoint trimmomatic_split:
  """
  Splits file with information about all reads files into separate files
  so they can be processed with trimmomatic in parallel.
  """   
  input:
    samples=config["samples_file"]
  output:
    directory("trimmomatic")
  log:
    "logs/trimmomatic_split.log"
  shell:
    """
    mkdir -p {output} &> {log}
    split -d -l 1 {input} {output}/sample_ &>> {log}
    """

rule trimmomatic_parallel:
  """
  Processes individual read files with trimmomatic.
  """
  input:
    "trimmomatic/sample_{job_index}"
  output:
    "trimmomatic/completed_{job_index}"
  log:
    "logs/trimmomatic_parallel{job_index}.log"
  conda:
    "envs/trimmomatic.yaml"
  params:
    memory="10"
  threads:
    4
  shell:
    """
    # note: trimmomatic can use gzipped files directly
    read SAMPLE REPLICATE F_READS R_READS < {input}
    # If the sample line is empty, ignore it
    if [ -z "$REPLICATE" ]; then
      touch {output}
      exit 0
    fi
    if [ ! -z "$R_READS" ]; then
      trimmomatic PE -threads {threads} $F_READS $R_READS trimmomatic/{wildcards[job_index]}.R1-P.qtrim.fastq.gz trimmomatic/{wildcards[job_index]}.R1-U.qtrim.fastq.gz trimmomatic/{wildcards[job_index]}.R2-P.qtrim.fastq.gz trimmomatic/{wildcards[job_index]}.R2-U.qtrim.fastq.gz {config[trimmomatic_parameters]} &> {log}
      echo $SAMPLE	$REPLICATE	trimmomatic/{wildcards[job_index]}.R1-P.qtrim.fastq.gz	trimmomatic/{wildcards[job_index]}.R2-P.qtrim.fastq.gz > {output} 2>> {log}
    else
      trimmomatic SE -threads {threads} $F_READS trimmomatic/{wildcards[job_index]}.U.qtrim.fastq.gz {config[trimmomatic_parameters]} &> {log}
      echo $SAMPLE      $REPLICATE      trimmomatic/{wildcards[job_index]}.U.qtrim.fastq.gz > {output} 2>> {log}
    fi
    """

def trimmomatic_completed_parallel_jobs(wildcards):
  """
  Returns names of files with information about files processed with trimmomatic.
  """
  parallel_dir = checkpoints.trimmomatic_split.get(**wildcards).output[0]
  job_ids = glob_wildcards(os.path.join(parallel_dir, "sample_{job_index}")).job_index
  completed_ids = expand(os.path.join(parallel_dir,"completed_{job_index}"), job_index=job_ids)
  return completed_ids

rule trimmomatic_merge:
  """
  Creates a file with information about all reads files processed with trimmomatic.
  """
  input:
    trimmomatic_completed_parallel_jobs
  output:
    samples_trimmed="samples_trimmed.txt"
  log:
    "logs/trimmomatic_merge.log"
  shell:
    """
    cat {input} > {output.samples_trimmed} 2> {log}
    """

rule fastqc_after_trim:
  """
  Runs fastQC on individual trimmed reads files.
  """
  input:
    samples="samples_trimmed.txt"
  output:
    directory("fastqc_after_trim")
  log:
    "logs/fastqc_after_trim.log"
  conda:
    "envs/qc.yaml"
  params:
    memory="4"
  threads:
    1
  shell:
    """
    # run fastqc on input files
    # making directory for the fastqc summary files
    mkdir {output} &> {log}
    FILES=$(awk '{{ printf "%s\\n%s\\n", $3,$4}}' {input}) &>> {log}
    # running fasqc on the files
    for file in $FILES; do fastqc -f fastq -o {output} $file; done &>> {log}
    """

rule multiqc_after_trim:
  """
  Creates multiQC report from individual fastQC reports after the trimming.
  """
  input:
    "fastqc_after_trim"
  output:
    out_directory=directory("multiqc_after_trim"),
    report="multiqc_after_trim.txt"
  log:
    "logs/multiqc_after_trim.log"
  conda:
    "envs/qc.yaml"
  params:
    memory="4"
  threads:
    1
  shell:
    """
    # common error resolved by those two export commands
    export LC_ALL=C.UTF-8
    export LANG=C.UTF-8
    mkdir {output.out_directory} &> {log}
    multiqc -o {output.out_directory} {input} &>> {log}
    cp multiqc_after_trim/multiqc_data/multiqc_fastqc.txt {output.report} &>> {log}
    """

rule trinity_inchworm_chrysalis:
  """
  Runs first two stages of Trinity assembly (Inchworm and Chrysalis).
  """
  input:
    samples="samples_trimmed.txt",
  output:
    "trinity_out_dir/recursive_trinity.cmds"
  log:
    "logs/trinity_inchworm_chrysalis.log"
  conda:
    "envs/trinity_utils.yaml"
  params:
    memory="50"
  threads:
    16
  shell:
    """
    Trinity --no_distributed_trinity_exec --max_memory {params.memory}G --CPU {threads} --samples_file {input} {config[trinity_parameters]} {config[strand_specific]} &> {log}
    """


checkpoint trinity_butterfly_split:
  """
  Preparation for last stage of Trinity (Butterfly) parallelization by splitting
  the commands into independent parts.
  """
  input:
    "trinity_out_dir/recursive_trinity.cmds"
  output:
    directory("trinity_out_dir/parallel_jobs")
  log:
    "logs/trinity_split.log"
  shell:
    """
    mkdir -p {output} &> {log}
    # Note: it is important to use -d and not --numeric-suffixes, see https://github.com/transXpress/transXpress-snakemake/issues/12
    # Note #2: we split into 1000 chunks to avoid running over the command line limit when too many parallel jobs are created, see https://bitbucket.org/snakemake/snakemake/issues/878/errno-7-argument-list-too-long-path-to-bin 
    split -n l/1000 -e -d {input} {output}/job_ &>> {log}
    """

rule trinity_butterfly_parallel:
  """
  Runs Trinity Butterfly commands (which were split into parts) in parallel.
  """
  input:
    "trinity_out_dir/parallel_jobs/job_{job_index}"
  output:
    "trinity_out_dir/parallel_jobs/completed_{job_index}"
  log:
    "logs/trinity_parallel{job_index}.log"
  conda:
    "envs/trinity_utils.yaml"
  params:
    memory="10"
  threads:
    1
  shell:
    """
    bash {input} &> {log}
    cp -p {input} {output} &>> {log}
    """

def trinity_completed_parallel_jobs(wildcards):
  """
  Returns filenames of the files processed in parallel.
  """
  parallel_dir = checkpoints.trinity_butterfly_split.get(**wildcards).output[0]
  job_ids = glob_wildcards(os.path.join(parallel_dir, "job_{job_index}")).job_index
  completed_ids = expand(os.path.join(parallel_dir,"completed_{job_index}"), job_index=job_ids)
  return completed_ids

rule trinity_butterfly_parallel_merge:
  input:
    jobs=trinity_completed_parallel_jobs,
    cmds="trinity_out_dir/recursive_trinity.cmds"
  output:
    cmds_completed="trinity_out_dir/recursive_trinity.cmds.completed"
  log:
    "logs/trinity_butterfly_parallel_merge.log"
  params:
    memory="10"
  threads:
    1
  shell:
    """
    # Can crash if there are too many parallel jobs
    # See https://bitbucket.org/snakemake/snakemake/issues/878/errno-7-argument-list-too-long-path-to-bin 
    cat {input.jobs} > {output.cmds_completed} 2> {log}
    """

rule trinity_final:
  """
  Runs the final Trinity assembly.
  """
  input:
    cmds_completed="trinity_out_dir/recursive_trinity.cmds.completed",
    samples="samples_trimmed.txt"
  output:
    transcriptome="trinity_out_dir.Trinity.fasta", #new output files live in base directory now, not trinity_out_dir.
    gene_trans_map="trinity_out_dir.Trinity.fasta.gene_trans_map" #according to trinity dev, this is a feature, not a bug
  log:
    "logs/trinity_final.log"
  conda:
    "envs/trinity_utils.yaml"
  params:
    memory="200"
  threads:
    16
  shell:
    """
    Trinity --max_memory {params.memory}G --CPU {threads} --samples_file {input.samples} {config[trinity_parameters]} {config[strand_specific]} &>> {log}
    """

rule transcriptome_copy:
  """
  Copies the assembled transcriptome to the transxpress directory.
  """
  input:
    transcriptome=rules.rnaspades.output.transcriptome if config["assembler"]=="rnaspades" else rules.trinity_final.output.transcriptome,
    gene_trans_map=rules.rnaspades.output.gene_trans_map if config["assembler"]=="rnaspades" else rules.trinity_final.output.gene_trans_map,
  output:
    transcriptome="transcriptome.fasta",
    gene_trans_map="transcriptome.gene_trans_map",
    #create copies the same way older Trinity versions did, for parity, in case other dependencies on that path exist?
    redundant_transcriptome="trinity_out_dir/Trinity.fasta",
    redundant_gene_trans_map="trinity_out_dir/Trinity.fasta.gene_trans_map"
  log:
    "logs/transcriptome_copy.log"
  shell:
    """
    cp -p {input.transcriptome} {output.transcriptome} &> {log}
    cp -p {input.gene_trans_map} {output.gene_trans_map} &>> {log}
    #create copies the same way older Trinity versions did, for parity, in case other dependencies on that path exist?
    cp -p {input.transcriptome} {output.redundant_transcriptome} &> {log}
    cp -p {input.gene_trans_map} {output.redundant_gene_trans_map} &>> {log}
    """

rule trinity_stats:
  """
  Runs Trinity script to get statistics about the assembled transcriptome
  (number of transcripts, genes, GC content, EXN50).
  """
  input:
    transcriptome="transcriptome.fasta",
    expression="transcriptome_expression_isoform.tsv"
  output:
    stats="transcriptome_stats.txt",
    exN50="transcriptome_exN50.tsv",
    exN50plot="ExN50_plot.pdf"
  log:
    "logs/trinity_exN50.log"
  conda:
    "envs/trinity_utils.yaml"
  params:
    memory="2"
  threads:
    1
  shell:
    """
    TRINITY_HOME=$(python -c 'import os;import shutil;TRINITY_EXECUTABLE_PATH=shutil.which("Trinity");print(os.path.dirname(os.path.join(os.path.dirname(TRINITY_EXECUTABLE_PATH), os.readlink(TRINITY_EXECUTABLE_PATH))))')

    assembler={config[assembler]}
    if [ "$assembler" = 'trinity' ]; then
        $TRINITY_HOME/util/TrinityStats.pl {input.transcriptome} > {output.stats} 2> {log}
    else
        $TRINITY_HOME/util/TrinityStats.pl {input.transcriptome} | sed -e 's/trinity/rnaspades/g' > {output.stats} 2> {log}
    fi
    $TRINITY_HOME/util/misc/contig_ExN50_statistic.pl {input.expression} {input.transcriptome} > {output.exN50} 2>> {log}
    $TRINITY_HOME/util/misc/plot_ExN50_statistic.Rscript {output.exN50} &>> {log}
    """

rule busco:
  """
  Runs BUSCO to assess the completeness of the transcriptome.
  """
  input:
    transcriptome="transcriptome.fasta"
  output:
    out_directory=directory("busco"),
    report="busco_report.txt"
  log:
    "logs/busco.log"
  params:
    memory="10",
    docker_image="ezlabgva/busco:v5.5.0_cv1"
  threads:
    4
  shell:
    """
    lineage={config[lineage]}
    docker run -u $(id -u) -v $(pwd):/busco_wd {params.docker_image} bash -c "
        echo $lineage &> /busco_wd/{log};
        if [ -z \"$lineage\" ]; then
            busco -m transcriptome -i /busco_wd/{input.transcriptome} -o /busco_wd/{output.out_directory} --auto-lineage -c {threads} &>> /busco_wd/{log};
        else
            busco -m transcriptome -i /busco_wd/{input.transcriptome} -o /busco_wd/{output.out_directory} -l $lineage -c {threads} &>> /busco_wd/{log};
        fi;

        status=$?;

        if [ $status -eq 0 ]; then
            echo 'BUSCO run completed successfully' &>> /busco_wd/{log};
            cp /busco_wd/busco/short_summary*.txt /busco_wd/{output.report} &>> /busco_wd/{log};
        else
            echo 'BUSCO run failed' &>> /busco_wd/{log};
            exit 1 &>> /busco_wd/{log};
        fi
    "
    """

rule transdecoder_longorfs:
  """
  Runs first stage of Transdecoder extracting the long open reading frames.
  """
  input:
    transcriptome="transcriptome.fasta",
  output:
    orfs="transcriptome.orfs"
  log:
    "logs/transdecoder_longorfs.log"
  conda:
    "envs/transdecoder.yaml"
  params:
    memory="2"
  threads:
    1
  shell:
    """
    rm -rf {input.transcriptome}.transdecoder_dir &> {log}
    TransDecoder.LongOrfs -t {input.transcriptome} --output_dir transdecoder &>> {log} 
    cp -p transdecoder/longest_orfs.pep {output.orfs} &>> {log}
    """

rule kallisto:
  """
  Runs Trinity script to perform transcript expression quantification 
  using Kallisto.
  """
  input:
    samples="samples_trimmed.txt",
    transcriptome="transcriptome.fasta",
    gene_trans_map="transcriptome.gene_trans_map"
  output:
    "transcriptome_expression_isoform.tsv",
    "transcriptome_expression_gene.tsv",
    "kallisto.gene.counts.matrix"
  log:
    "logs/kallisto.log"
  conda:
    "envs/trinity_utils.yaml"
  params:
    memory="16" # increased memory from 2 to 8 since it was not sufficient
  threads:
    8
  shell:
    """
    TRINITY_HOME=$(python -c 'import os;import shutil;TRINITY_EXECUTABLE_PATH=shutil.which("Trinity");print(os.path.dirname(os.path.join(os.path.dirname(TRINITY_EXECUTABLE_PATH), os.readlink(TRINITY_EXECUTABLE_PATH))))')

    assembler="{config[assembler]}"
    strand_specific="{config[strand_specific]}"

    if [ $assembler = "rnaspades" ]
    then
      if [[ $strand_specific = "--ss rf" ]]
      then
        $TRINITY_HOME/util/align_and_estimate_abundance.pl --transcripts {input.transcriptome} --SS_lib_type RF --seqType fq --samples_file {input.samples} --prep_reference --thread_count {threads} --est_method kallisto --gene_trans_map {input.gene_trans_map} &> {log}
      elif [[ $strand_specific = "--ss fr" ]]
      then
        $TRINITY_HOME/util/align_and_estimate_abundance.pl --transcripts {input.transcriptome} --SS_lib_type FR --seqType fq --samples_file {input.samples} --prep_reference --thread_count {threads} --est_method kallisto --gene_trans_map {input.gene_trans_map} &> {log}
      else
        $TRINITY_HOME/util/align_and_estimate_abundance.pl --transcripts {input.transcriptome} --seqType fq --samples_file {input.samples} --prep_reference --thread_count {threads} --est_method kallisto --gene_trans_map {input.gene_trans_map} &> {log}
      fi
    else
      $TRINITY_HOME/util/align_and_estimate_abundance.pl --transcripts {input.transcriptome} {config[strand_specific]} --seqType fq --samples_file {input.samples} --prep_reference --thread_count {threads} --est_method kallisto --gene_trans_map {input.gene_trans_map} &> {log}
    fi
    
    $TRINITY_HOME/util/abundance_estimates_to_matrix.pl --est_method kallisto --name_sample_by_basedir --gene_trans_map {input.gene_trans_map} */abundance.tsv &>> {log}
    if [ -f kallisto.isoform.TMM.EXPR.matrix ]; then
      cp -p kallisto.isoform.TMM.EXPR.matrix {output[0]} &>> {log}
    elif [ -f kallisto.isoform.TPM.not_cross_norm ]; then
      cp -p kallisto.isoform.TPM.not_cross_norm {output[0]} &>> {log}
    else
      echo Neither kallisto.isoform.TMM.EXPR.matrix or kallisto.isoform.TPM.not_cross_norm were produced
      exit 1
    fi
    if [ -f kallisto.gene.TMM.EXPR.matrix ]; then
      cp -p kallisto.gene.TMM.EXPR.matrix {output[1]} &>> {log}
    elif [ -f kallisto.gene.TPM.not_cross_norm ]; then
      cp -p kallisto.gene.TPM.not_cross_norm {output[1]} &>> {log}
    else
      echo Neither kallisto.gene.TMM.EXPR.matrix or kallisto.gene.TPM.not_cross_norm were produced
      exit 1
    fi
    """