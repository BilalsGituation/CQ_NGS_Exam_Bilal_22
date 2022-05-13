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

// More instructions:
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
    path "*"
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

workflow {
  infolder_fastq = "/rawdata"
  fastqChannel = Channel.fromPath('rawdata/*.fastq').flatten()
  dbChannel = Channel.value("${params.db}")
  fastpout = fastp(fastqChannel.flatten())
  fastpout.trim.view()
  s2_in1 = (fastpout.trim.flatten())
  // Edit: these files are too big to run fastqc on the raw data
  // so I am going to omit the raw data/fastp data channel
  qcd = fastqc(fastpout.trim.flatten())
  qcd.view()
  texts = srst2(s2_in1, dbChannel)
  texts.report.view()
}

// cmd line:
// nextflow NGS_Exam_Bilal.nf --db CARD_v3.0.8_SRST2.fasta --outdir AlonaResults -profile singularity

/* cmd line output before submission:
N E X T F L O W  ~  version 21.10.6
Launching `NGS_Exam_Bilal.nf` [focused_bernard] - revision: ef73ffcb00
executor >  local (8)
executor >  local (9)
[f0/671709] process > fastp (3)  [100%] 3 of 3 ✔
executor >  local (9)
[f0/671709] process > fastp (3)  [100%] 3 of 3 ✔
[f2/f32361] process > fastqc (3) [100%] 3 of 3 ✔
executor >  local (9)
[f0/671709] process > fastp (3)  [100%] 3 of 3 ✔
[f2/f32361] process > fastqc (3) [100%] 3 of 3 ✔
[93/794a60] process > srst2 (1)  [ 33%] 1 of 3
[/home/cq/Documents/work/d1/8c3bcc618ed1ecfe0a25f93968bfe2/patient1_fastp_fastqc.html,executor >  local (9)
[f0/671709] process > fastp (3)  [100%] 3 of 3 ✔
[f2/f32361] process > fastqc (3) [100%] 3 of 3 ✔
[50/f7f28c] process > srst2 (2)  [ 66%] 2 of 3
[/home/cq/Documents/work/a0/717f033a4dd40e050fffd4805c100f/patient2_fastp_fastqc.html, /home/cq/Documents/work/a0/717f033a4dd40e050fffd4805c100f/patient2_fastp_fastqc.zip]
executor >  local (9)
[f0/671709] process > fastp (3)  [100%] 3 of 3 ✔
[f2/f32361] process > fastqc (3) [100%] 3 of 3 ✔
[49/20b7e0] process > srst2 (3)  [100%] 3 of 3 ✔
[/home/cq/Documents/work/f2/f32361d02eebc46dde360a3166b150/patient3_fastp_fastqc.html, /home/cq/Documents/work/f2/f32361d02eebc46dde360a3166b150/patient3_fastp_fastqc.zip]
[/home/cq/Documents/work/93/794a60ab0d92d375ac47f176e50259/patient1_fastp.txt__fullgenexecutor >  local (9)
[f0/671709] process > fastp (3)  [100%] 3 of 3 ✔
[f2/f32361] process > fastqc (3) [100%] 3 of 3 ✔
[49/20b7e0] process > srst2 (3)  [100%] 3 of 3 ✔
[/home/cq/Documents/work/f2/f32361d02eebc46dde360a3166b150/patient3_fastp_fastqc.html, /home/cq/Documents/work/f2/f32361d02eebc46dde360a3166b150/patient3_fastp_fastqc.zip]
[/home/cq/Documents/work/93/794a60ab0d92d375ac47f176e50259/patient1_fastp.txt__fullgenes__CARD_v3.0.8_SRST2__results.txt, /home/cq/Documents/work/93/794a60ab0d92d375ac47f176e50259/patient1_fastp.txt__genes__CARD_v3.0.8_SRST2__results.txt]
[/home/cq/Documents/work/50/f7f28c0e5476b398ef74ed1c88b6df/patient2_fastp.txt__fullgenes__CARD_v3.0.8_SRST2__results.txt, /home/cq/Documents/work/50/f7f28c0e5476b398ef74ed1c88b6df/patient2_fastp.txt__genes__CARD_v3.0.8_SRST2__results.txt]
[/home/cq/Documents/work/49/20b7e0a689c852fe171364efebccd5/patient3_fastp.txt__fullgenes__CARD_v3.0.8_SRST2__results.txt, /home/cq/Documents/work/49/20b7e0a689c852fe171364efebccd5/patient3_fastp.txt__genes__CARD_v3.0.8_SRST2__results.txt]

Completed at: 12-May-2022 12:11:46
Duration    : 1m 24s
CPU hours   : 0.1
Succeeded   : 9

(base) cq@bioinfobox:PATH_CENSORED ls AlonaResults/
patient1_fastp.fastq
patient1_fastp_fastqc.html
patient1_fastp_fastqc.zip
patient1_fastp.html
patient1_fastp.json
patient1_fastp.txt__fullgenes__CARD_v3.0.8_SRST2__results.txt
patient1_fastp.txt__genes__CARD_v3.0.8_SRST2__results.txt
patient1_fastp.txt.log
patient1_fastp.txt__patient1_fastp.CARD_v3.0.8_SRST2.pileup
patient1_fastp.txt__patient1_fastp.CARD_v3.0.8_SRST2.sorted.bam
patient2_fastp.fastq
patient2_fastp_fastqc.html
patient2_fastp_fastqc.zip
patient2_fastp.html
patient2_fastp.json
patient2_fastp.txt__fullgenes__CARD_v3.0.8_SRST2__results.txt
patient2_fastp.txt__genes__CARD_v3.0.8_SRST2__results.txt
patient2_fastp.txt.log
patient2_fastp.txt__patient2_fastp.CARD_v3.0.8_SRST2.pileup
patient2_fastp.txt__patient2_fastp.CARD_v3.0.8_SRST2.sorted.bam
patient3_fastp.fastq
patient3_fastp_fastqc.html
patient3_fastp_fastqc.zip
patient3_fastp.html
patient3_fastp.json
patient3_fastp.txt__fullgenes__CARD_v3.0.8_SRST2__results.txt
patient3_fastp.txt__genes__CARD_v3.0.8_SRST2__results.txt
patient3_fastp.txt.log
patient3_fastp.txt__patient3_fastp.CARD_v3.0.8_SRST2.pileup
patient3_fastp.txt__patient3_fastp.CARD_v3.0.8_SRST2.sorted.bam
*/
