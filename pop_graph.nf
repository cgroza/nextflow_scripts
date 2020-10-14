params.ref = "/home/cgroza/Homo_sapiens.hg19.noalts.fa"
params.genome = "EU_AF"
params.index_gcsa = true
params.index_xg = true
params.index_gbwt = true
params.outdir = workflow.launchDir
params.vcf = "$params.outdir/renamed_Epi_EU_AF_phased.vcf.gz"
params.tmp = "/home/cgroza/scratch/temp"

Channel.fromPath(params.vcf).into{vcf_con_ch; vcf_gbwt_ch}

process makeVg {
    time '12h'
    memory '100 GB'
    cpus 40

    publishDir "$params.outdir/graphs", mode: 'copy', pattern: "*.vg"

    input:
    file vcf from vcf_con_ch

    output:
    file "graphs/*.vg" into vgs_ch_gbwt, vgs_ch_gcsa, vgs_ch_xg
    file "*.tbi" into vcf_index_ch
    file "*.vcf.gz" into vcf_ch
    file "mapping" into mapping_ch

    script:
    """
module load tabix
tabix -p vcf $vcf
(seq 1 22; echo X; echo Y) | parallel -j 24 "tabix -h $vcf {} > chr{}.vcf ; bgzip chr{}.vcf ; tabix chr{}.vcf.gz"

mkdir graphs
(seq 1 22; echo X; echo Y) | parallel -j 6  "vg construct -f -S -a -C -R {} -v chr{}.vcf.gz -r $params.ref -t 1 -m 32 > graphs/chr{}.vg"
vg ids -m mapping -j \$(for i in \$(seq 1 22; echo X; echo Y); do echo graphs/chr\$i.vg; done)
"""
}

if(params.index_gbwt) {
    process indexGBWT {
        cpus 40
        time '2d'
        memory '180 GB'
        publishDir "$params.outdir", mode: 'copy', pattern: "${params.genome}_index.gbwt"

        input:
        file "*" from vgs_ch_gbwt.collect()
        file "*" from vcf_ch.collect()
        file "*" from vcf_index_ch.collect()

        output:
        file "${params.genome}_index.gbwt"
        file "${params.genome}_index.xg" into xg_ch
        file "*.gbwt" into gbwt_ch

        script:
        """
        TMPDIR=${params.tmp}
        (seq 1 22; echo X; echo Y) | parallel -j 8 "touch -h chr{}.vcf.gz.tbi ; vg index -G chr{}.gbwt -v chr{}.vcf.gz chr{}.vg"
        vg gbwt -m -f -o ${params.genome}_index.gbwt chr*.gbwt
        """
    }
}

if(params.index_xg) {
    process index_XG {
        publishDir "$params.outdir", mode: 'copy', pattern: "${params.genome}_index.xg"
        cpus 40
        time '6h'
        memory '180GB'

        input:
        file "*" from vgs_ch_xg.collect()

        output:
        file "*.xg"

        script:
        """
        TMPDIR=${params.tmp}
        vg index -b \${TMPDIR} -L -x ${params.genome}_index.xg *.vg
        """
    }
}

if(params.index_gcsa) {
    process indexGCSA {
        cpus 40
        time '2d'
        memory '180 GB'
        publishDir "$params.outdir", mode: 'copy'

        input:
        file "*" from vgs_ch_gcsa.collect()
        file "*" from gbwt_ch.collect()
        file mapping from mapping_ch

        output:
        file "${params.genome}_index.gcsa"
        file "${params.genome}_index.gcsa.lcp"

        script:
        """
        mkdir graphs
        TMPDIR=${params.tmp}
        cp ${mapping} mapping.backup
        for i in \$(seq 1 22; echo X; echo Y); do
        vg prune -a -m mapping.backup -u -g chr\${i}.gbwt chr\${i}.vg > graphs/chr\${i}.pruned.vg
        done

        vg index -b \${TMPDIR} -g ${params.genome}_index.gcsa -f mapping.backup graphs/*.pruned.vg
        """
    }
}
