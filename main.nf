#!/usr/bin/env nextflow
params.name             = "RNA-seq"
params.email            = "sidhantnp@yahoo.com"
params.reads            = "/home/sidhantpuntambekar/workshop/data/fastq/*{*_R1,*_R2}.fastq.gz"



log.info "RNA-seq Pipeline"
log.info "====================================="
log.info "name         : ${params.name}"
log.info "email        : ${params.email}"
log.info "reads        : ${params.reads}"
log.info "\n"



annotation = Channel.fromPath("/home/sidhantpuntambekar/workshop/data/annotation/*")
genome = Channel.fromPath("/home/sidhantpuntambekar/workshop/data/genome/*")
sample_info_for_de = Channel.fromPath("/home/sidhantpuntambekar/workshop/data/sample_info/*")
index = Channel.fromPath("/home/sidhantpuntambekar/workshop/data/index")
reads = Channel.fromFilePairs(params.reads, size: -1)
  .ifEmpty { error "Can't find any reads matching: ${params.reads}" }
  .into {
    reads_for_fastqc;
    reads_for_mapping
  }



annotation.into {
  annotation_for_count;
  annotation_for_transcriptome;
  annotation_for_de;
}



process make_transcriptome {

  maxForks 1
  publishDir 'results/genome'

  input:
  file annotation from annotation_for_transcriptome
  file genome from genome

  output:
  file "transcriptome.fa" into transcriptome
  file "gencode.vM21.annotation.gtf"

  script:
  """
  gffread -w transcriptome.fa -g ${genome} ${annotation}
  """
}



process fastqc {

  maxForks 1
  publishDir 'results/fastqc'

  input:
  set sample_id, file(fastqz) from reads_for_fastqc

  output:
  file "*.zip" into fastqc

  script:
  """
  fastqc --threads 4 -f fastq -q ${fastqz}
  """
}



process map {

  maxForks 1
  publishDir 'results/bam'

  input:
  set sample_id, file(reads), file(index) from reads_for_mapping.combine(index)

  output:
  set sample_id, file("*Aligned.out.bam") into mapped_genome
  set sample_id, file("*toTranscriptome.out.bam") into mapped_transcriptome
  file '*' into star

  script:
  """
  STAR  --runThreadN 5 \
  --genomeDir ${index} \
  --readFilesIn ${reads.findAll{ it =~ /\_R1\./ }.join(',')} \
                ${reads.findAll{ it =~ /\_R2\./ }.join(',')} \
  --readFilesCommand zcat \
  --outSAMtype BAM Unsorted \
  --outSAMunmapped Within \
  --outSAMattributes NH HI NM MD AS \
  --outReadsUnmapped Fastx \
  --quantMode TranscriptomeSAM \
  --outFileNamePrefix ${sample_id}_
  """
}



mapped_genome.into {
  mapped_for_count;
  mapped_for_igv;
  mapped_for_markduplicates
}



process count {

  publishDir 'results/feature_counts'

  input:
  set sample_id, file(bam), file(annotation) from mapped_for_count.combine(annotation_for_count)

  output:
  file '*.fCounts'
  file '*.fCounts*' into counts

  script:
  """
  featureCounts  -C \
    -p \
    -T 1 \
    -g gene_id \
    -a ${annotation} \
    -o ${sample_id}.fCounts \
    ${bam}
  """
}



process salmon {

  publishDir 'results/salmon'

  input:
  set sample_id, file(bam), file(transcript_fasta) from mapped_transcriptome.combine(transcriptome)

  output:
  file("*") into salmon
  file("*") into salmon_for_de

  script:
  """
  salmon quant -l A \
    -p 1 \
    -t ${transcript_fasta} \
    -o ${sample_id} \
    -a ${bam} \
    --numBootstraps 30
  """
}



process sort_bam {

  publishDir 'results/igv'

  input:
  set sample_id, file(bam_file) from mapped_for_igv

  output:
  set sample_id, file("*.bam"), file('*.bai') into sorted_bam

  script:
  """
  samtools sort --threads 1 \
    -m 4G \
    -o ${sample_id}.bam \
    ${bam_file}
  samtools index ${sample_id}.bam
  """
}



process collect_fastqc {

  input:
  file fastqc from fastqc.collect()

  output:
  file "fastqc" into fastqc_collected

  script:
  """
  mkdir fastqc
  mv ${fastqc} fastqc/.
  """
}



process collect_star {

  input:
  file star from star.collect()

  output:
  file "star" into star_collected

  script:
  """
  mkdir star
  mv ${star} star/.
  """
}



process collect_salmon {

  input:
  file salmon from salmon.collect()

  output:
  file "salmon" into salmon_collected

  script:
  """
  mkdir salmon
  mv ${salmon} salmon/.
  """
}



process collect_counts {

  input:
  file counts from counts.collect()

  output:
  file "counts" into counts_collected

  script:
  """
  mkdir counts
  mv ${counts} counts/.
  """
}



process multiqc {

  publishDir "results/multiqc"

  input:
  file fastqc from fastqc_collected
  file star from star_collected
  file salmon from salmon_collected
  file counts from counts_collected

  output:
  set file('*_multiqc_report.html'), file('*_data/*')

  script:
  """
  multiqc ${fastqc} \
    ${star} \
    ${salmon} \
    ${counts} \
    --title '${params.name}' \
    --cl_config "extra_fn_clean_exts: [ '_1', '_2' ]"
  """
}

process differential_expression {

  publishDir "reports"

  input:
  file annotation from annotation_for_de
  file salmon from salmon_for_de.collect()
  file sample_info from sample_info_for_de

  output:
  file "*.html"

  script:
  """
  cp ${baseDir}/bin/*.R* .
  Rscript -e 'rmarkdown::render("differential_expression.Rmd", params = list(annotation_file = "${annotation}"))'
  """
}
