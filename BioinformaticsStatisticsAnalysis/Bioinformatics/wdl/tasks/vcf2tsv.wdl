task vcf2tsv{
        String VCF2TSV
        File  InputVcf
	String SampleName
	String PYTHON

        command{
                ${PYTHON} ${VCF2TSV} -v ${InputVcf}  -o ${SampleName}.SpliceAI.parsed.tsv
        }

        output{
                File OutputVcf = "${SampleName}.SpliceAI.parsed.tsv"
        }
}
