task seqprep {

    File LeftFQ_cut
    File RightFQ_cut

    String SEQPREP
    String SEQPREP_OPTIONS
    String OUTDIR

    File PROFILE
    String SAMPLE_NAME=basename(LeftFQ_cut, "_1.cut.fq")

    command {
        source ${PROFILE}
        ${SEQPREP} ${SEQPREP_OPTIONS} -f ${LeftFQ_cut} -r ${RightFQ_cut} -1 ${OUTDIR}/${SAMPLE_NAME}_1.sp.fq.gz -2 ${OUTDIR}/${SAMPLE_NAME}_2.sp.fq.gz -s ${OUTDIR}/${SAMPLE_NAME}.merged.fq.gz
        rm ${LeftFQ_cut} ${RightFQ_cut}
    }

    output {
        File MergedFQ_sp = "${OUTDIR}/${SAMPLE_NAME}.merged.fq.gz"
    }
}