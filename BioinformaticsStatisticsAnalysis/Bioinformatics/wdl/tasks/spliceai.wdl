task spliceai {
	String SpliceaiEnvProfile
	String SpliceaiScript
	File InputVcf
	String OutputVcfName
	String ReferenceGenome
	String CliOptions = " -A grch38 "

	Int SpliceaiThreads = 1
	String SpliceaiSoftMemLimit = "7936M"
	String SpliceaiHardMemLimit = "8G"


	command {
		source ${SpliceaiEnvProfile}

		echo "running spliceai"
		${SpliceaiScript} -I ${InputVcf} -O ${OutputVcfName} -R ${ReferenceGenome} ${CliOptions}
	}

	output {
		File OutputVcf = "${OutputVcfName}"
	}

	runtime{
		cpu: "${SpliceaiThreads}"
		s_vmem: "${SpliceaiSoftMemLimit}"
		h_vmem: "${SpliceaiHardMemLimit}"
	}

}

