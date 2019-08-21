params.xg = "EU_AF_index.xg"
params.gbwt = "EU_AF_index.gbwt"
params.gcsa = "EU_AF_index.gcsa"
params.gcsa_lcp = "${params.gcsa}.lcp"
params.reads = "reads"
params.paired = true

xg_ch = Channel.fromPath(params.xg)
gbwt_ch = Channel.fromPath(params.gbwt)
gcsa_ch = Channel.fromPath(params.gcsa)
gcsa_lcp_ch = Channel.fromPath(params.gcsa_lcp)

reads_ch = Channel.fromPath("${params.reads}/*")
file("gams").mkdir()

process alignFastq {
    time '2d'
    cpus 40
    memory '100 GB'
    publishDir "$workflow.launchDir/gams"

    input:
    set file("index.xg"), file("index.gbwt"), file('index.gcsa'), file('index.gcsa.lcp'), file(reads) from xg_ch.combine(gbwt_ch).combine(gcsa_ch).combine(gcsa_lcp_ch).combine(reads_ch)

    output:
    file "${reads_name}.gam" into aln

    script:
    extension = reads.getExtension()
    reads_name = reads.getSimpleName()
    vg_flag = ""
    switch(extension) {
        case "fastq.gz": case "fq.gz":
            vg_flag = "-f"
            break
        case "bam":
            vg_flag = "-b"
            break
        case "gam":
            vg_flag = "-G"
            break
    }
    if(params.paired) {
        vg_flag = "-i ${vg_flag}"
    }
    """
vg map -k 18 --threads 40 --gcsa-name index.gcsa --xg-name index.xg --gbwt-name index.gbwt ${vg_flag} ${reads} > ${reads_name}.gam
"""
}
