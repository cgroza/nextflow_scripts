params.xg = "EU_AF_index.xg"
params.gbwt = "EU_AF_index.gbwt"
params.gcsa = "EU_AF_index.gcsa"
params.bams = "bams/*.bam"

Channel.fromPath(params.xg, params.gbwt, params.gcsa).set{index_ch}
bams = Channel.fromPath(params.bams)
file("gams").mkdir()

process alignFastq {
    time '2d'
    cpus 40
    memory '100 GB'
    publishDir "$workflow.launchDir/gams"

    input:
    file bam from bams
    set file(xg), file(gbwt), file(gcsa) from index_ch

    output:
    file "${bam_name}.gam" into aln

    script:
    bam_name = bam.getSimpleName()
    """
vg map -k 18 --threads 40 --gcsa-name $gcsa --xg-name $xg --gbwt-name $gbwt -b $bam > ${bam_name}.gam
"""
}
