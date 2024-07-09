task SamToFastq {
    File BAMFILE
    String PICARD
    String OUTDIR
    File PROFILE
    String SAMPLE_NAME=basename(BAMFILE, ".bam")

    command {
        source ${PROFILE}
        ${PICARD} SamToFastq I=${BAMFILE} F=${OUTDIR}/${SAMPLE_NAME}_1.fq F2=${OUTDIR}/${SAMPLE_NAME}_2.fq
    }

    output {
        File LeftFQ = "${OUTDIR}/${SAMPLE_NAME}_1.fq"
        File RightFQ ="${OUTDIR}/${SAMPLE_NAME}_2.fq"
    }
}