task makeVCF {

    File PROFILE
    String AllPossibleScript

    String SEQ
    String CHR
    Int START
    String OUTNAME
    String LOGLEVEL = "INFO"

    command {
        source ${PROFILE}
        python ${AllPossibleScript} -s ${SEQ} -c ${CHR} -f ${START} -V ${LOGLEVEL} -o ${OUTNAME}.raw.vcf
    }

    output {
        File VCFout = "${OUTNAME}.raw.vcf"
    }
}
