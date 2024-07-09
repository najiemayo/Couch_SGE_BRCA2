###############################################################################
# Name: tabixIndexing
# Purpose: Create a tabix index file of a bgzip compressed VCF.
# @param: File BashPreamble - Bash setup parameters.
# @param: String Tabix - Path to tabix.
# @param: File CompressedVcf - Bgzip compressed VCF to make an index of.
# @return: File IndexFile - Tabix index file.
###############################################################################
task tabixIndex{
        File BashPreamble
        String Tabix
        String CompressedVcf

        command{
                source ${BashPreamble}
                ${Tabix} -f -p vcf "${CompressedVcf}"
        }

        output{
                File IndexFile = "${CompressedVcf}" + ".tbi"
        }
}
