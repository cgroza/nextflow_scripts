params.xg = "EU_AF_index.xg"
params.gbwt = "EU_AF_index.gbwt"
params.gcsa = "EU_AF_index.gcsa"
params.gcsa_lcp = "${params.gcsa}.lcp"
params.bams = "bams/*.bam"

xg_ch = Channel.fromPath(params.xg)
gbwt_ch = Channel.fromPath(params.gbwt)
gcsa_ch = Channel.fromPath(params.gcsa)
gcsa_lcp_ch = Channel.fromPath(params.gcsa_lcp)

bams = Channel.fromPath(params.bams)
file("gams").mkdir()

process alignFastq {
    time '2d'
    cpus 40
    memory '100 GB'
    publishDir "$workflow.launchDir/gams"

    input:
    file bam from bams
    file("index.xg") from xg_ch
    file("index.gbwt") from gbwt_ch
    file('index.gcsa') from gcsa_ch
    file('index.gcsa.lcp') from gcsa_lcp_ch

    output:
    file "${bam_name}.gam" into aln

    script:
    bam_name = bam.getSimpleName()
    """
vg map -k 18 --threads 40 --gcsa-name $gcsa --xg-name $xg --gbwt-name $gbwt -b $bam > ${bam_name}.gam
"""
}
