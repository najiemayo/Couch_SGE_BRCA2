task catfastq {

    File LeftFQ_cut
    File RightFQ_cut

    String OUTDIR
    String CAT
    String GZIP
    File PROFILE
    String SAMPLE_NAME=basename(LeftFQ_cut, "_1.cut.fq")

    command {
        source ${PROFILE}
        ${CAT} ${LeftFQ_cut} ${RightFQ_cut} | ${GZIP} -c >  ${OUTDIR}/${SAMPLE_NAME}.merged.fq.gz
#        rm ${LeftFQ_cut} ${RightFQ_cut}
    }

    output {
        File MergedFQ_sp = "${OUTDIR}/${SAMPLE_NAME}.merged.fq.gz"
    }
}