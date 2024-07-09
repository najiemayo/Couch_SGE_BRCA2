import "../tasks/spliceai.wdl" as SpliceAI
import "../tasks/spliceai_parse.wdl" as spliceai_parse
import "../tasks/vcf2tsv.wdl" as vcf2tsv

workflow SeqToAnno {
    File PROFILE
    File INPUT_VCF
    String SAMPLENAME

    String SpliceaiEnvProfile
    String ReferenceGenome
    String SpliceaiScript
    String SpliceaiParseScript
    
    String VCF2TSV
    String PYTHON

    call SpliceAI.spliceai{ input:
    	SpliceaiEnvProfile = SpliceaiEnvProfile,
	    SpliceaiScript = SpliceaiScript,
	    ReferenceGenome = ReferenceGenome,
	    InputVcf = INPUT_VCF,
	    OutputVcfName = SAMPLENAME + ".SpliceAI.vcf"
    }

    call spliceai_parse.spliceai_parse{ input:
            PYTHON = PYTHON,
            SpliceaiEnvProfile = SpliceaiEnvProfile,
	    SpliceaiParseScript = SpliceaiParseScript,
	    InputVcf = spliceai.OutputVcf,
	    OutputVcfName = SAMPLENAME + ".SpliceAI.parsed.vcf"
    }

    call vcf2tsv.vcf2tsv{ input:
    	 SampleName = SAMPLENAME,
	 VCF2TSV = VCF2TSV,
	 InputVcf = spliceai_parse.OutputVcf,
	 PYTHON = PYTHON
    }

    output {
        File OutputVcf = vcf2tsv.OutputVcf
    }
}
