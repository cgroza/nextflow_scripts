params.ref = "/home/cgroza/Homo_sapiens.hg19.noalts.fa"
params.genome = "EU_AF"
params.outdir = workflow.launchDir
params.vcf = "$params.outdir/renamed_Epi_EU_AF_phased.vcf.gz"

vcf_ch = Channel.fromPath(params.vcf)

process makeVg {
time '6h'
memory '100 GB'
cpus 6

input:
file vcf from vcf_ch

output:
file "graphs/*.vg" into vgs_ch_gbwt
file "graphs/*.vg" into vgs_ch_xg
file "graphs/*.vg" into vgs_ch_gcsa

script:
"""
module load tabix
tabix $vcf
mkdir graphs
(seq 1 22; echo X; echo Y) | parallel -j 6  "vg construct -a -p -C -R chr{} -v $vcf -r $params.ref -t 1 -m 32 > graphs/chr{}.vg"
vg ids -j \$(for i in \$(seq 1 22; echo X; echo Y); do echo graphs/chr\$i.vg; done)
"""
}

process indexGBWT {
cpus 40
time '1d'
memory '150 GB'
publishDir "$params.outdir", mode: 'copy'

input:
file "*.vg" from vgs_ch_gbwt.collect()

output:
file "${params.genome}_threads.db" into db_thread
file "${params.genome}_thread_seqs.fa" into hap_seqs
file "${params.genome}_index.gbwt" into hap_index

script:
"""
vg index -F ${params.genome}_threads.db -H ${params.genome}_thread_seqs.fa -G ${params.genome}_index.gbwt -t 40 -b ~/scratch/temp  -v $params.vcf *.vg
"""
}

process indexXg {
cpus 40
time '1d'
memory '100 GB'
publishDir "$params.outdir", mode: 'copy'

input:
file threads from db_thread
file "*.vg" from vgs_ch_xg.collect()

output:
file "${params.genome}_index.xg" into xg_ch

script:
"""
vg index -x "${params.genome}_index.xg" -t 40 -b ~/scratch/temp  -F $threads *.vg
"""
}

process indexGCSA {
cpus 40
time '1d'
memory '180 GB'
publishDir "$params.outdir", mode: 'copy'

input:
file "*.vg" from vgs_ch_gcsa.collect()
file gbwt from hap_index

output:
file "${params.genome}_index.gcsa"
file "${params.genome}_index.gcsa.lcp"
file "graphs/*.vg"

script:
"""
mkdir graphs
ls *.vg | parallel -j 12 "vg prune -u -g $gbwt {} > graphs/pruned.{}"
vg index -g ${params.genome}_index.gcsa pruned.*.vg
"""
}
