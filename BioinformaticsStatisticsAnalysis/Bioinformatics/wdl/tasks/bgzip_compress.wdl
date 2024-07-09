###############################################################################
# Name: bgzipCompress
# Purpose: Compress the input file with bgzip.
# @param: File BashPreamble - Bash setup parameters.
# @param: String Bgzip - File path to bgzip.
# @param: File CompressVcf - Input VCF to be compressed.
# @param: String CompressedVcf - Name of input VCF with ".gz" suffix.
# @return: File OutputVcf - Compressed VCF.
###############################################################################
task bgzipCompress{
        File BashPreamble
        String Bgzip
        File CompressVcf
        String CompressedVcf=basename(CompressVcf)+".gz"

        Int BgzipCompressThreads = 1
        String BgzipCompressSoftMemLimit = "7936M"
        String BgzipCompressHardMemLimit = "8G"

        command{
                source ${BashPreamble}
                ${Bgzip} -c ${CompressVcf} > "${CompressedVcf}"
        }
        output{
                File OutputVcf = "${CompressedVcf}"
        }

        runtime{
                cpu: "${BgzipCompressThreads}"
                s_vmem: "${BgzipCompressSoftMemLimit}"
                h_vmem: "${BgzipCompressHardMemLimit}"
        }

}
