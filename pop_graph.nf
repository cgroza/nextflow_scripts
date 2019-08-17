params.ref = "/home/cgroza/Homo_sapiens.hg19.noalts.fa"
params.genome = "EU_AF"
params.outdir = workflow.launchDir
params.vcf = "$params.outdir/renamed_Epi_EU_AF_phased.vcf.gz"

Channel.fromPath(params.vcf).into{vcf_con_ch; vcf_gbwt_ch}

process makeVg {
time '6h'
memory '100 GB'
cpus 6

input:
file vcf from vcf_con_ch

output:
file "graphs/*.vg" into vgs_ch_gbwt
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
file vcf from vcf_gbwt_ch

output:
file "${params.genome}_index.gbwt" into gbwt_ch
file "${params.genome}_index.xg" into xg_ch

script:
"""
TMPDIR=/home/cgroza/scratch/temp
vg index -x ${params.genome}_index.xg -G ${params.genome}_index.gbwt -v $vcf *.vg
"""
}


process indexGCSA {
cpus 40
time '1d'
memory '180 GB'
publishDir "$params.outdir", mode: 'copy'

input:
file "*.vg" from vgs_ch_gcsa.collect()
file gbwt from gbwt_ch

output:
file "${params.genome}_index.gcsa"
file "${params.genome}_index.gcsa.lcp"
file "graphs/*.vg"

script:
"""
mkdir graphs
TMPDIR=/home/cgroza/scratch/temp
ls *.vg | parallel -j 3 "vg prune -u -g $gbwt {} > graphs/pruned.{}"
vg index -g ${params.genome}_index.gcsa pruned.*.vg
"""
}
