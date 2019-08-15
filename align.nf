params.xg = "$baseDir/EU_AF_index.xg"
params.gbwt = "$baseDir/EU_AF_index.gbwt"

bams = Channel.fromPath("bams/*.bam")

process bamToFastqGz {
time '12h'
cpus 12
memory "20 GB"

input:
file bam from bams

output:
set val(bam_name), file("fastq.fq.gz") into fastq_ch
stdout conversion

script:
bam_name = bam.getBaseName()
"""
module load java picard
echo Converting ${bam_name} to FASTQ
java -jar \$PICARD_HOME/picard.jar SamToFastq I=$bam INTERLEAVE=true FASTQ=fastq.fq.gz
"""

}

process alignFastq {
time '2d'
cpus 40
publishDir "$baseDir/gams"

input:
set val(bam_name), fqgz from fastq_ch

output:
file "${bam_name}.gam" into aln
stdout aligning

script:
"""
mkdir gams
echo Aligning ${bam_name}
vg map -k 18 --threads 40 --xg-name $params.xg --gbwt-name $params.gbwt -i $fqgz > ${bam_name}.gam
"""
}

conversion.subscribe { println it }
aligning.subscribe { println it }
