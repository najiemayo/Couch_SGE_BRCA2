import "../tasks/SamToFastq.wdl" as SamToFastq
import "../tasks/Cutadapt.wdl" as Cutadapt
import "../tasks/catfastq.wdl" as catfastq
import "../tasks/bwa.wdl" as bwa
import "../tasks/SortSam.wdl" as SortSam
import "../tasks/SamtoBam.wdl" as SamtoBam


workflow SingleReadToSam {
    File FASTQ1
    File FASTQ2
    String OUTDIR
    File PROFILE
    String JAVA
    String PICARD

    String CUTADAPT
    String CUTADAPT_OPTIONS = "-g AGCATACCCAGGGCTTCATAATCAC -B GTAAAATAAAGTGATTAT "

    String GZIP
    String CAT

    String REF_GENOME
    String BWA
    String BWA_OPTIONS = "-T 20 "

#    call SamToFastq.SamToFastq{ input:
#            BAMFILE = BAMFILE,
#            PICARD = PICARD,
#            OUTDIR = OUTDIR,
#            PROFILE = PROFILE
#    }

    call Cutadapt.cutadapt{ input:
         LeftFQ = FASTQ1,
         RightFQ = FASTQ2,
         CUTADAPT = CUTADAPT,
         CUTADAPT_OPTIONS = CUTADAPT_OPTIONS,
         OUTDIR = OUTDIR,
         PROFILE = PROFILE
    }

    call catfastq.catfastq{ input:
         LeftFQ_cut = cutadapt.LeftFQ_cut,
	 RightFQ_cut = cutadapt.RightFQ_cut,
	 PROFILE = PROFILE,
	 CAT = CAT,
	 GZIP = GZIP,
	 OUTDIR = OUTDIR
    }

    call bwa.bwa{ input:
        MergedFQ_sp = catfastq.MergedFQ_sp,
    	REF_GENOME = REF_GENOME,
        BWA = BWA,
        BWA_OPTIONS = BWA_OPTIONS,
        OUTDIR = OUTDIR,
        PROFILE = PROFILE
    }

    call SortSam.sortsam{ input:
    	 SAMPLE_SAM = bwa.SamFile,
	 REF_GENOME = REF_GENOME,
	 JAVA = JAVA,
	 PICARD = PICARD
    }

    output {
        File SamOUT = sortsam.SortedSamFile
    }
}