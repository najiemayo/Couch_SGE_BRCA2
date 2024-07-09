task samtobam {

    File SORTED_SAM
    String REF_GENOME

    String SAMTOOLS
    String SAMTOOLS_OPTION= " -b -S "
    String OUTDIR

    File PROFILE
    String SAMPLE_NAME=basename(SORTED_SAM, ".sorted.sam")

    command {
        source ${PROFILE}
	${SAMTOOLS} view ${SAMTOOLS_OPTION} ${SORTED_SAM} -o ${OUTDIR}/${SAMPLE_NAME}.bam
    }

    output {
        File BamFile = "${OUTDIR}/${SAMPLE_NAME}.bam"
    }
}