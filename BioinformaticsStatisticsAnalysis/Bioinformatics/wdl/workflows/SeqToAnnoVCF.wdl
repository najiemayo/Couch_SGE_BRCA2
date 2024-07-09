import "../tasks/MakeVCFFromSeqAndPos.wdl" as MakeVCF
import "../tasks/Cava.wdl" as Cava
import "../tasks/spliceai.wdl" as SpliceAI
import "../tasks/spliceai_parse.wdl" as spliceai_parse

workflow SeqToAnno {

    File PROFILE
    String AllPossibleScript

    String SEQ
    String CHR
    Int START
    String SAMPLENAME
    String PYTHON
    String CAVA_PYTHON
    String CavaScript
    File CavaConfig

    String SpliceaiEnvProfile
    String ReferenceGenome
    String SpliceaiScript
    String SpliceaiParseScript
    

    call MakeVCF.makeVCF{ input:
        PROFILE = PROFILE,
        AllPossibleScript = AllPossibleScript,
        SEQ = SEQ,
        CHR = CHR,
        START = START,
        OUTNAME = SAMPLENAME
    }

    call Cava.cava{ input:
        PROFILE = PROFILE,
	CAVA_PYTHON = CAVA_PYTHON,
        CavaScript = CavaScript,
        CavaConfig = CavaConfig,
        InputVcf = makeVCF.VCFout,
        OutputVcfName = SAMPLENAME
    }

    call SpliceAI.spliceai{ input:
    	SpliceaiEnvProfile = SpliceaiEnvProfile,
	SpliceaiScript = SpliceaiScript,
	ReferenceGenome = ReferenceGenome,
	InputVcf = cava.OutputVcf,
	OutputVcfName = SAMPLENAME + ".SpliceAI.vcf"
    }
    
    call spliceai_parse.spliceai_parse{ input:
        PYTHON = PYTHON, 
       	SpliceaiEnvProfile = SpliceaiEnvProfile,
	SpliceaiParseScript = SpliceaiParseScript,
	InputVcf = spliceai.OutputVcf,
	OutputVcfName = SAMPLENAME
    }

    output {
        File VCFout = spliceai_parse.OutputVcf
    }
}
