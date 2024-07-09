#################################################################################
# Name: vcfSort
# Purpose: Sort the VCF by CHROM -> POS -> REF -> ALT
# @param: File BashPreamble - Bash setup parameters.
# @param: File InputVcf - VCF to be sorted.
# @param: File VcfSort - Path to vcf-sort script in the vcf-tools toolkit.
# @return File OutputVcf - Sorted VCF.
################################################################################
task vcfSort{
        File BashPreamble
        File? InputVcf
        String VcfSort
        String SortedVcf = "sorted.vcf"

        Int SortThreads = 1
        String SortSoftMemLimit = "7936M"
        String SortHardMemLimit = "8G"

        command{
                source ${BashPreamble}
                ${VcfSort} -c ${InputVcf} > "${SortedVcf}"
        }

        output{
                File OutputVcf = "${SortedVcf}"
        }

        runtime{
                cpu: "${SortThreads}"
                s_vmem: "${SortSoftMemLimit}"
                h_vmem: "${SortHardMemLimit}"
        }

}
