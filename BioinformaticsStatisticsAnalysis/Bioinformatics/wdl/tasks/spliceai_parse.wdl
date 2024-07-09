
task spliceai_parse {
	String SpliceaiEnvProfile
	String SpliceaiParseScript
	File InputVcf
	String OutputVcfName
	String PYTHON

	Int SpliceaiParseThreads = 1
	String SpliceaiParseSoftMemLimit = "7936M"
	String SpliceaiParseHardMemLimit = "8G"


	command {
		source ${SpliceaiEnvProfile}

		echo "running spliceai_parse"
		${PYTHON} ${SpliceaiParseScript} -v ${InputVcf} -o ${OutputVcfName}

	}

	output {
		File OutputVcf = "${OutputVcfName}"
	}

	runtime{
		cpu: "${SpliceaiParseThreads}"
		s_vmem: "${SpliceaiParseSoftMemLimit}"
		h_vmem: "${SpliceaiParseHardMemLimit}"
	}

}

