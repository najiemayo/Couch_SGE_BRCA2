task cava {

        File PROFILE
    	String CavaScript
	File CavaConfig
    	File InputVcf
	String OutputVcfName
	
	Int CavaThreads = 1
	String CavaSoftMemLimit = "7936M"
	String CavaHardMemLimit = "8G"
	String CAVA_PYTHON


	command {
		source ${PROFILE}

		echo "running cava"
		${CAVA_PYTHON} ${CavaScript} -i ${InputVcf} -c ${CavaConfig} -o "${OutputVcfName}"

	}

	output {
		File OutputVcf = "${OutputVcfName}.vcf"
	}

	runtime{
		cpu: "${CavaThreads}"
		s_vmem: "${CavaSoftMemLimit}"
		h_vmem: "${CavaHardMemLimit}"
	}

}

