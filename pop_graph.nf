params.ref = "/home/cgroza/Homo_sapiens.hg19.noalts.fa"
params.vcf = "$workflow.launchDir/renamed_Epi_EU_AF_phased.vcf.gz"
params.genome = "EU_AF"

process makeVg {
time '6h'
memory '100 GB'
cpus 6

output:
file "graphs/*.vg" into vgs

script:
"""
mkdir graphs
(seq 1 22; echo X; echo Y) | parallel -j 6  "vg construct -a -p -C -R chr{} -v $params.vcf -r $params.ref -t 1 -m 32 > graphs/chr{}.vg"
vg ids -j \$(for i in \$(seq 1 22; echo X; echo Y); do echo graphs/chr\$i.vg; done)
"""
}

process indexGBWT {
cpus 40
time '1d'
memory '150 GB'
publishDir "$workflow.launchDir"

input:
file "*.vg" from vgs.collect()

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
publishDir "$workflow.launchDir"

input:
file threads from db_thread
file "*.vg" from vgs.collect()

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
memory '150 GB'
publishDir "$workflow.launchDir"

input:
file "*.vg" from vgs.collect()
file gbwt from hap_index

output:
file "${params.genome}_index.gcsa"
file "${params.genome}_index.gcsa.lcp"

script:
"""
ls *.vg | parallel -j "vg prune -u -g $gbwt {} > pruned.{}"
vg index -g ${params.genome}_index.gcsa pruned.*.vg
"""
}
