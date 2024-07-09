task cutadapt {

    File LeftFQ
    File RightFQ

    String CUTADAPT
    String CUTADAPT_OPTIONS = "-g AGCATACCCAGGGCTTCATAATCAC -B GTAAAATAAAGTGATTAT "
    String OUTDIR

    File PROFILE
    String SAMPLE_NAME=basename(LeftFQ, "_1.fq")

    command {
        ${CUTADAPT} ${CUTADAPT_OPTIONS} -o ${OUTDIR}/${SAMPLE_NAME}_1.cut.fq -p ${OUTDIR}/${SAMPLE_NAME}_2.cut.fq ${LeftFQ} ${RightFQ}
    }

    output {
        File LeftFQ_cut = "${OUTDIR}/${SAMPLE_NAME}_1.cut.fq"
        File RightFQ_cut = "${OUTDIR}/${SAMPLE_NAME}_2.cut.fq"
    }
}