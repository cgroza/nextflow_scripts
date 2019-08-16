params.xg = "$workflow.launchDir/EU_AF_index.xg"
params.gbwt = "$workflow.launchDir/EU_AF_index.gbwt"

bams = Channel.fromPath("bams/*.bam")
file("$workflow.launchDir/gams").mkdir()

process alignFastq {
time '2d'
cpus 40
memory '100 GB'
publishDir "$workflow.launchDir/gams"

input:
file bam from bams

output:
file "${bam_name}.gam" into aln

script:
bam_name = bam.getBaseName()
"""
vg map -k 18 --threads 40 --xg-name $params.xg --gbwt-name $params.gbwt -b $bam > ${bam_name}.gam
"""
}
