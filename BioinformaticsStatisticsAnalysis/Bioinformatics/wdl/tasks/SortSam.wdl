task sortsam {

    File SAMPLE_SAM
    String REF_GENOME

    String JAVA
    String PICARD
    String PICARD_OPTION= "SORT_ORDER=coordinate --CREATE_INDEX"
    String OUTDIR

    File PROFILE
    String SAMPLE_NAME=basename(SAMPLE_SAM, ".sam")

    command {
        source ${PROFILE}
	${JAVA} -jar ${PICARD} SortSam I=${SAMPLE_SAM} O=${OUTDIR}/${SAMPLE_NAME}.sorted.bam ${PICARD_OPTION}
    }

    output {
        File SortedSamFile = "${OUTDIR}/${SAMPLE_NAME}.sorted.bam"
    }
}