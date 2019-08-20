params.xg = "EU_AF_index.xg"
params.gbwt = "EU_AF_index.gbwt"
params.gcsa = "EU_AF_index.gcsa"
params.gcsa_lcp = "${params.gcsa}.lcp"
params.bams = "bams"

xg_ch = Channel.fromPath(params.xg)
gbwt_ch = Channel.fromPath(params.gbwt)
gcsa_ch = Channel.fromPath(params.gcsa)
gcsa_lcp_ch = Channel.fromPath(params.gcsa_lcp)

bams = Channel.fromPath("${params.bams}/*.bam")
file("gams").mkdir()

process alignFastq {
    time '2d'
    cpus 40
    memory '100 GB'
    publishDir "$workflow.launchDir/gams"

    input:
    set file("index.xg"), file("index.gbwt"), file('index.gcsa'), file('index.gcsa.lcp'), file(bam) from xg_ch.combine(gbwt_ch).combine(gcsa_ch).combine(gcsa_lcp_ch).combine(bams)

    output:
    file "${bam_name}.gam" into aln

    script:
    bam_name = bam.getSimpleName()
    """
vg map -k 18 --threads 40 --gcsa-name index.gcsa --xg-name index.xg --gbwt-name index.gbwt -b $bam > ${bam_name}.gam
"""
}
