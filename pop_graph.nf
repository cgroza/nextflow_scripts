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
file "*.tbi" into vcf_index_ch
file "*.vcf.gz" into vcf_ch
file "mapping" into mapping_ch

script:
"""
module load tabix
tabix -p vcf $vcf
(seq 1 22; echo X; echo Y) | parallel -j 24 "tabix -h $vcf chr{} > chr{}.vcf ; bgzip chr{}.vcf ; tabix chr{}.vcf.gz"

mkdir graphs
(seq 1 22; echo X; echo Y) | parallel -j 6  "vg construct -a -p -C -R chr{} -v chr{}.vcf.gz -r $params.ref -t 1 -m 32 > graphs/chr{}.vg"
vg ids -m mapping -j \$(for i in \$(seq 1 22; echo X; echo Y); do echo graphs/chr\$i.vg; done)
"""
}

process indexGBWT_XG {
cpus 40
time '1d'
memory '150 GB'
publishDir "$params.outdir", mode: 'copy'

input:
file "*.vg" from vgs_ch_gbwt.collect()
file "*.vcf.gz" from vcf_ch
file "*.vcf.gz.tbi" from vcf_index_ch

output:
file "${params.genome}_index.gbwt" into gbwt_ch
file "${params.genome}_index.xg" into xg_ch

script:
"""
TMPDIR=/home/cgroza/scratch/temp
(seq 1 22; echo X; echo Y) | parallel -j 8 "vg index -G chr{}.gbwt -v {}.vcf.gz {}.vg"
vg gbwt -m -f -o ${params.genome}_index.gbwt chr*.gbwt
vg index -x ${params.genome}_index.xg *.vg
"""
}


process indexGCSA {
cpus 40
time '2d'
memory '180 GB'
publishDir "$params.outdir", mode: 'copy'

input:
file "*.vg" from vgs_ch_gcsa.collect()
file "*.gbwt" from gbwt_ch
file mapping from mapping_ch

output:
file "${params.genome}_index.gcsa"
file "${params.genome}_index.gcsa.lcp"

script:
"""
mkdir graphs
TMPDIR=/home/cgroza/scratch/temp
for i in \$(seq 1 22; echo X; echo Y); do
    vg prune -a -m ${mapping} -u -g \${i}.gbwt \${i}.vg > graphs/chr\${i}.pruned.vg
done

vg index -g ${params.genome}_index.gcsa graphs/*.vg
"""
}
