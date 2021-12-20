#!/usr/bin/env nextflow



// defind parameter

params.fasta = ( "/home/labs/aharoni/yedidyab/rna_anotation_pipeline/genome/Chr02.fa" )
params.transcriptomeFile = ("/home/labs/aharoni/yedidyab/rna_anotation_pipeline//home/labs/aharoni/yedidyab/rna_anotation_pipeline/transc/*" )
params.spesies = ("WT")
params.outDir = ("/home/labs/aharoni/yedidyab/rna_anotation_pipeline/output_pipeline")

// defind channels

fasta_ch = Channel.fromPath(params.fasta)
trans = Channel.fromPath(params.transcriptomeFile)


    
println """\


         R N A   ANNOTATION   P I P E L I N E
         ====================================
         fasta file    : ${params.fasta}
         """
         .stripIndent()



/*
 * STEP 1 - mask repetative elements (EDTA)
 */
 

process EDTA {
  publishDir "${params.outDir}/EDTA", mode: 'copy'
  conda '/home/labs/aharoni/yedidyab/anaconda3/envs/EDTA'


  input:
  file fasta from fasta_ch

  output:
  path('*fa.mod.EDTA.TEanno.gff3') into (masked_genome_ch,masked_genome_ch2)

  script:
  """
  EDTA.pl --genome ${fasta} --anno 1 -force 1
  """
}

/*
 * STEP 2 - soft mask (bed tools)
 */
 

process softMask {
  publishDir "${params.outDir}/softMask", mode: 'copy'

  input:
  file gff from masked_genome_ch

  output:
  path('*.soft.gff') into (soft_mask_ch,soft_mask_ch2)

  script:
  """
  name=`echo $gff | awk '{split(\$0,a,"."); print a[1]}'`
  bedtools maskfasta -soft -fi ${params.fasta} -bed ${gff} -fo \${name}.soft.gff
  """
}

/*
 * STEP 3 - preparing to trinity (combine R1&R2)
 */
 

process preTrinity {
  publishDir "${params.outDir}/preTrinity", mode: 'copy'


  input:
  path transc from trans.collect()

  output:
  tuple file('left_1.gz'), file('right_2.gz') into (readsLeftRight1,readsLeftRight2)

  script:
  """
  cat *1.*gz > left_1.gz
  cat *2.*gz > right_2.gz
  """
}

/*
 * STEP 4 - Genome-independent transcriptome assembly
 */
 
process independenTtrinity {
  publishDir "${params.outDir}/independenTtrinity", mode: 'copy'


  input:
  tuple file(R1), file(R2) from readsLeftRight1

  output:
  file '*' into independent

  script:
  """
  Trinity --seqType fq \
  --max_memory 20G \
  --CPU 16 \
  --normalize_by_read_set \
  --output TrinityOut \
  --left $R1 \
  --right $R2 \
  --trimmomatic --jaccard_clip
  """
}


/*
 * STEP 5 - STAR DB 
 */
 
process starDb {
  publishDir "${params.outDir}/starDb", mode: 'copy'
  
  output:
  file 'db_file' into starDB_ch
  
  script:
  """
  mkdir db_file
  STAR --runThreadN 20 --genomeDir db_file --runMode genomeGenerate --genomeSAindexNbases 10 --genomeFastaFiles ${params.fasta} 
  """
}

/*
 * STEP 6 - STAR mapping
 */
 
process mappingStar {
  publishDir "${params.outDir}/mappingStar", mode: 'copy'

  input:
  tuple file(R1), file(R2) from readsLeftRight2
  file starDB from starDB_ch

  output:
  file "*.sortedByCoord.out.bam" into (bam_ch_to_count,bam_ch_to_umi2)
  
  script:
  """
  STAR --genomeDir $starDB \
  --runThreadN 20 \
  --readFilesIn $R1 $R2 \
  --outFileNamePrefix "STAR_" \
  --outSAMtype BAM SortedByCoordinate\
  --outSAMunmapped Within\
  --outSAMattributes Standard\
  --twopassMode Basic\
  --readFilesCommand zcat   
  """
  
}

/*
 * STEP 7 - Genome-guided transcriptome assembly
 */
 
process GenomeTtrinity {
  publishDir "${params.outDir}/genomeTtrinity", mode: 'copy'


  input:
  file bam from bam_ch_to_count

  output:
  file '*' into genomeTrinity

  script:
  """
   Trinity --genome_guided_bam $bam --genome_guided_max_intron 10000 --max_memory 20G --CPU 20 --jaccard_clip
  """
}

/*
 * STEP 8 - passa - 
 */
 
process tarscriptPASA {
  publishDir "${params.outDir}/pasa", mode: 'copy'


  input:
  file gg_fasta from genomeTrinity
  file independent_fasta from independent

  output:

  script:
  """
  cat $independent_fasta $gg_fasta > transcripts.fasta
  
  /usr/local/src/PASApipeline/misc_utilities/accession_extractor.pl $independent_fasta tdn.accs
  
  /usr/local/src/PASApipeline/Launch_PASA_pipeline.pl \
  -c alignAssembly.config \
  -C \
  -R -g $gg_fasta \
  --ALIGNERS blat,gmap\
  -t $independent_fasta \
  --transcribed_is_aligned_orient \
  --TDN tdn.accs
  
  
  /usr/local/src/PASApipeline/scripts/build_comprehensive_transcriptome.dbi \
           -c alignAssembly.config \
           -t transcripts.fasta \
           --min_per_ID 95 \
           --min_per_aligned 30
  
  """
}