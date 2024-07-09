task bwa {

    File MergedFQ_sp
    String REF_GENOME

    String BWA
    String BWA_OPTIONS = "-C -T 20 "
    String OUTDIR

    File PROFILE
    String SAMPLE_NAME=basename(MergedFQ_sp, ".merged.fq.gz")

    command {
        source ${PROFILE}
        ${BWA} mem ${REF_GENOME} ${MergedFQ_sp} ${BWA_OPTIONS} -o ${OUTDIR}/${SAMPLE_NAME}.sam
        rm ${MergedFQ_sp}
    }

    output {
        File SamFile = "${OUTDIR}/${SAMPLE_NAME}.sam"
    }
}