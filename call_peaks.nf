params.pop_gams = "pop_gams.tsv"

params.pop_graph = "pop/"
params.pop_name = "pop"

params.genome_size = 3100000000
params.fragment_length = 200
params.qvalue = 0.05
params.paired = true
params.sort = false
params.storeDir = "store"

params.outDir = workflow.launchDir

chromosomes = "chr1,chr2,chr3,chr4,chr5,chr6,chr7,chr8,chr9,chr10,chr11,chr12,chr13,chr14,chr15,chr16,chr17,chr18,chr19,chr20,chr21,chr22,chrX,chrY"


Channel.fromPath("${params.pop_graph}/graphs/*.vg").set{linear_vg_ch}

def pop_gams = []
new File(params.pop_gams).eachLine { line ->
    pop_gams << line
}

println("Pop gams: ")
print(pop_gams)

Channel.fromPath(pop_gams).set{gam_ch}

Channel.fromPath(
    ["${params.pop_graph}/${params.pop_name}.xg",
     "${params.pop_graph}/${params.pop_name}.gbwt",
     "${params.pop_graph}/${params.pop_name}.gcsa",
     "${params.pop_graph}/${params.pop_name}.gcsa.lcp"]).into{pop_index_treatment_ch; pop_index_control_ch}


process vgToJsonPop {
    cpus = 40
    memory '120 GB'
    time '24h'
    storeDir "${params.storeDir}/json"

    input:
    file "graphs/*" from linear_vg_ch.collect()

    output:
    file "graphs/*" into vg2json_ch

    script:
    """
    (seq 1 22; echo X; echo Y) | parallel -j 3 'vg view -Vj graphs/chr{}.vg > graphs/chr{}.json ; graph_peak_caller create_ob_graph graphs/chr{}.json ; vg stats -r graphs/chr{}.vg  | cut -f 2 > graphs/node_range_chr{}.txt'
"""
}

process linearPathsPop {
    cpus = 40
    memory '120 GB'
    time '24 h'
    storeDir "${params.storeDir}/linear"

    input:
    file "graphs/*" from vg2json_ch.collect()

    output:
    file "graphs" into control_linear_ch, treatment_linear_ch, peak_linear_ch

    script:
    """
   (seq 1 22; echo X; echo Y) | parallel -j 3 graph_peak_caller find_linear_path -g graphs/chr{}.nobg graphs/chr{}.json chr{} graphs/chr{}_linear_pathv2.interval
"""
}


process processGamPop {
    cpus = 1
    memory "30 GB"
    time = "12 h"

    input:
    set file(gam), file("graphs") from gam_ch.combine(treatment_linear_ch).view()

    output:
    set file(gam), file("json") into treatment_json_ch

    script:
    name = gam.getSimpleName()

    """
    mkdir json
    vg view -aj $gam > json/${name}_pop.json
    graph_peak_caller split_vg_json_reads_into_chromosomes ${chromosomes} json/${name}_pop.json graphs/
    rm json/${name}_pop.json
"""
}

process callPeaksPop{
    cpus = 40
    memory '170 GB'
    time '24h'

    publishDir "$params.outDir/peaks", pattern: "${name}_peaks.narrowPeak", mode: "copy"

    input:
    set file(gam), file("json"), file("graphs") from treatment_json_ch.combine(peak_linear_ch).map{ it.flatten()}.view()

    output:
    set val(name), file("${name}_peaks.narrowPeak") into pop_peaks_ch

    script:
    name = gam.getSimpleName()

    """
    (seq 1 22; echo X; echo Y) | parallel -j 3 'graph_peak_caller count_unique_reads chr{} graphs/ json/${name}_pop_ | tail -n 1 > counted_unique_reads_chr{}.txt'
    read_length=\$(vg view -X $gam | head -2 | tail -1 | wc -c)
    unique_reads=\$(awk 'BEGIN{i=0}{i = i + \$1}END{print i}' counted_unique_reads_chr*.txt)
    (seq 1 22; echo X; echo Y) | parallel -j 3 "graph_peak_caller callpeaks -q ${params.qvalue} -g graphs/chr{}.nobg -s json/${name}_pop_chr{}.json -G ${params.genome_size} -p True -f ${params.fragment_length} -r \$read_length -u \$unique_reads -n chr{}"

    rename 'touched' '_touched' *touched*
    rename 'background' '_background' *background*
    rename 'direct' '_direct' *direct*
    rename 'fragment' '_fragment' *fragment*
    rename 'pvalues' '_pvalues' *pvalues*

    (seq 1 22; echo X; echo Y) | parallel -j 3 "graph_peak_caller callpeaks_whole_genome_from_p_values -q ${params.qvalue} -d graphs/ -n '' -f ${params.fragment_length} -r \${read_length} chr{}"
    (seq 1 22; echo X; echo Y) | parallel -j 3 'graph_peak_caller peaks_to_linear chr{}_max_paths.intervalcollection graphs/chr{}_linear_pathv2.interval chr{} chr{}_linear_peaks.bed'
    cat *_linear_peaks.bed | awk '\$2<\$3' | sort-bed - > ${name}_peaks.narrowPeak
"""
}
