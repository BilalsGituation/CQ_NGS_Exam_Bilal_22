nextflow.enable.dsl = 2

// After creating repo, it absolutely did not work to push the input data,
// nor the output data, to GitHub. Please use the rawdata.tar.gz file provided by our Teacher in order to run this Nextflow pipeline!
/*
Rules:
-One fastq folder as input
-All fastQ files should be analysed, in parallel when possible
-Single output file
-Resistance prediction Database should be called as a cmd line parameter

Steps of this workflow:

-fastp trim
-Fast QC report of fastp output
-Resistance prediction of trimmed data (using Srst2)

*/

// Starting material: Zipped rawdata file, with rawdata folder in it.
// In my workflow, this will be extracted in and into the directory you run this
// script from.

// Extraction of raw data:
// Used "tar -xvf archive.tar.gz" method
// It countains 3 fastq files in a directory

// So this workflow needs the processes fastp, fastqc and Srst2

// More instructions from the paper (or me to user?):
// RUN THIS IN A NEWLY-CREATED DIRECTORY WITH:
// 1) A NEXTFLOW CONFIG FILE
// 2) THIS SCRIPT
// 3) THE UNZIPPED rawdata.tar.gz ARCHIVE
// 4) NO OTHER FILES (maybe the exam answer (.md), and/or code licence)


params.outdir = "ExamResults"
//params.db = null



process fastp {
  container "https://depot.galaxyproject.org/singularity/fastp:0.20.1--h8b12597_0"
  publishDir "${params.outdir}", mode: "copy", overwrite: true
  input:
    path fastq
  output:
    path "${fastq.getSimpleName()}_fastp.fastq", emit: trim
    path "${fastq.getSimpleName()}_fastp.html"
    path "${fastq.getSimpleName()}_fastp.json", emit: fastpreport
  script:
    """
    fastp -i ${fastq} -o ${fastq.getSimpleName()}_fastp.fastq -h ${fastq.getSimpleName()}_fastp.html -j ${fastq.getSimpleName()}_fastp.json
    """
    }

process fastqc {
  publishDir "${params.outdir}", mode: "copy", overwrite: true
  container "https://depot.galaxyproject.org/singularity/fastqc%3A0.11.9--hdfd78af_1"
  input:
    path fastP_plusRaw
  output:
    path "${fastP_plusRaw}"
  script:
    """
    fastqc ${fastP_plusRaw}
    """
}

process srst2 {
  publishDir "${params.outdir}", mode: "copy", overwrite: true
  container "https://depot.galaxyproject.org/singularity/srst2%3A0.2.0--py27_2"
  input:
    path s2_in1
    val dbChannel
  output:
    path "patient*txt", emit: report
    path "*"
  script:
    """
    srst2 --input_se ${s2_in1} --output ${s2_in1.getSimpleName()}.txt --log --gene_db ~/Documents/${dbChannel}
    """
}

// I will leave the commented-out workflow steps as they are.
// The fastqc step was working for me but you can check it if you want
workflow {
  infolder_fastq = "/rawdata"
  //database = Channel.fromPath("${params.db}")
  fastqChannel = Channel.fromPath('rawdata/*.fastq').flatten()
  dbChannel = Channel.value("${params.db}")
  fastpout = fastp(fastqChannel.flatten())
  s2_in1 = (fastpout.trim.flatten())
  //s2_in2 = s2_in1.combine(dbChannel)
  fastP_plusRaw = fastqChannel.concat(fastpout.trim)
  qcd = fastqc(fastP_plusRaw.flatten())
  texts = srst2(s2_in1, dbChannel)
  texts.report.view()
  //resistDB.view()
}

// cmd line:
// nextflow NGS_Exam_Bilal.nf --db CARD_v3.0.8_SRST2.fasta --outdir AlonaResults -profile singularity
