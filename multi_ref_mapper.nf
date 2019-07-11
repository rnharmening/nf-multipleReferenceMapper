

params.fastas  = "demo/references/*.fa"
params.reads       = "demo/reads/*_R{1,2}.fastq"
params.outdir     = "./results"
params.singleEnd   = false
params.keepOnlyMapped = false


def helpMessage() {
    log.info """
    nextflow run multi_ref_mapper.nf --reads '*_{1,2}.fastq' --fastas  'references/*.fa'
    Mandatory arguments:
      --reads                       Path to input data (surrounded with quotes)
      --fastas

    Options:
      --outdir                     Specifies the output directory [./results/]
      --singleEnd                   Set iff input reads are single end [false]
    """.stripIndent()
}

// Show help emssage
if (params.help){
    helpMessage()
    exit 0
}

ref_ch = Channel.fromPath( params.fastas )

Channel.fromFilePairs( params.reads, size: params.singleEnd ? 1 : 2 )
        .ifEmpty { exit 1, "Cannot find any reads matching: ${params.reads}\n" +\
            "NB: Path needs to be enclosed in quotes!\nNB: Path requires at least one * wildcard!\n" +\
            "If this is single-end data, please specify --singleEnd on the command line." }
        .into { reads_ch; reads_fastqc_ch }

process fastqc {
  tag "$sample_id"
  publishDir "${params.outdir}/FastQC", mode: 'copy'

  input:
    set sample_id, file(reads) from reads_fastqc_ch

  output:
    file fastqc_out into fastqc_results_ch
    file(fastqc_out)
  
  script:
  fastqc_out = "*_fastqc.{zip,html}"
  """
    fastqc -q $reads
  """
}

process adapterremoval {
  tag "$sample_id"

  input:
    set sample_id, file(reads) from reads_ch

  output:
    set sample_id, file(adapter_removed_reads) into adapter_removed_ch

  script:
  base = reads[0].baseName
  adapter_removed_reads = "output/*.truncated"
  """
  mkdir -p output
  mkdir -p output
  AdapterRemoval --file1 ${reads[0]} --file2 ${reads[1]} --basename ${base} --threads ${task.cpus}
  mv ${base}.pair*.truncated output/
  """
}

process prinseq_derep {
  tag "$read_id"

  input:
    set read_id, file(reads) from adapter_removed_ch
  
  output:
    set read_id, file(derep_reads) into derep_reads_ch

  script:
  out="derep_fastqs"
  derep_reads = "${out}/*.fastq"
  """
    prinseq-lite.pl -out_format 3 -derep 23 -fastq ${reads[0]} -fastq2 ${reads[1]} -out_good $read_id -out_bad null
    mkdir -p ${out}
    mv *.fastq ${out}
  """
}

process bwa_index {
  input:
    file fasta from ref_ch
  output:
    set fasta_id, file(index) into bwa_index_ch 
  
  script:
  index = "BWA_index"
  fasta_id = fasta.baseName
  """
    bwa index $fasta
    mkdir ${index} && mv ${fasta}* ${index}
  """
}

process bwamem {
  tag "$sample_id"

  publishDir "${params.outdir}/BWA_mapping", mode: 'copy'

  input:
    set read_id, file(reads), val(ref_id), file(bwa_index) from derep_reads_ch.combine(bwa_index_ch)  

  output:
    set sample_id, file(bam_file) into mapping_pair_ch 
    file bam_file
    file "*.bai"
  
  script:
    sample_id = "${read_id}_v_${ref_id}"
    bam_file  = "${sample_id}.sorted.bam"
    filter = params.keepOnlyMapped ? "-F 4" : "" 
  """
  bwa mem -t ${task.cpus} -o tmp.sam ${bwa_index}/${ref_id}.fa ${reads[0]} ${reads[1]} 
  samtools view $filter -b -@ ${task.cpus} -o tmp.bam tmp.sam
  samtools sort -@ ${task.cpus} -o ${bam_file} tmp.bam
  samtools index -@ ${task.cpus}  ${bam_file}
  rm tmp.[sb]am
  """
}

/*
process dedup {
  input:
     set sample_id, file(bam_file) from mapping_pair_ch 

  output:

  
  script:
  """
  
  """
}
*/

process qualimap {
  tag "$sample_id"
  publishDir "${params.outdir}/qualimap", mode: 'copy'

  echo=false

  input:
     set sample_id, file(bam_file) from mapping_pair_ch 

  output:
    file "$out" into qualimap_results_ch
    file "$out"

  script:
  out="${sample_id}/"
  """
    qualimap bamqc -bam $bam_file -nt ${task.cpus} -outdir $out
  """
}

process multiqc {
  tag "MultiQC"

  echo =false
  publishDir "${params.outdir}/MultiQC", mode: 'copy'

  input:
    file fastqc from fastqc_results_ch.collect().ifEmpty([])
    file qualimap from qualimap_results_ch.collect().ifEmpty([])
  
  output:
    file "*multiqc_report.html" into multiqc_report
    file "*.html"
  
  script:
  """
    multiqc -f ./
  """
}