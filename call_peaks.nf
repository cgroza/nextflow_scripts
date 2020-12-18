params.pop_gams = "pop_gams.tsv"

params.pop_graph = "pop/"
params.pop_name = "pop"

params.genome_size = 3100000000
params.qvalue = 0.05
params.paired = true
params.sort = false
params.storeDir = "store"

params.outDir = workflow.launchDir

chromosomes = "1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,X,Y"


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
    (seq 1 22; echo X; echo Y) | parallel -j 3 'vg view -j graphs/chr{}.vg > graphs/{}.json ; graph_peak_caller create_ob_graph graphs/{}.json ; vg stats -r graphs/chr{}.vg  | cut -f 2 > graphs/node_range_{}.txt'
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
   (seq 1 22; echo X; echo Y) | parallel -j 3 graph_peak_caller find_linear_path -g graphs/{}.nobg graphs/{}.json {} graphs/{}_linear_pathv2.interval
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
    vg view -aj $gam > json/${name}.json
    graph_peak_caller split_vg_json_reads_into_chromosomes ${chromosomes} json/${name}.json graphs/
    rm json/${name}.json
"""
}

process callPeaksPop{
    cpus = 40
    memory '170 GB'
    time '12h'

    publishDir "$params.outDir/peaks", mode: "copy"

    input:
    set file(gam), file("json"), file("graphs") from treatment_json_ch.combine(peak_linear_ch).map{ it.flatten()}.view()

    output:
    file("${name}_peaks.narrowPeak")
    file("${name}.tar.gz")

    script:
    name = gam.getSimpleName()

    """
    (seq 1 22; echo X; echo Y) | parallel -j 10 'graph_peak_caller count_unique_reads {} graphs/ json/${name}_ | tail -n 1 > counted_unique_reads_{}.txt'
    (seq 1 22; echo X; echo Y) | parallel -j 10 'graph_peak_caller estimate_shift {} graphs/ json/${name}_ 2 100 | tail -n 1 | sed -n -e "s/.*Found shift: \\([:digit:]*\\)/\\1/p" > fragment_size_{}.txt'
    read_length=\$(vg view -X $gam | head -2 | tail -1 | wc -c)
    fragment_length=\$(awk 'BEGIN{i=0}{i = i + \$1}END{print int(i/NR)}' fragment_size_*.txt)
    if [ \$fragment_length -le \$read_length ]
    then
        fragment_length=200
    fi
    unique_reads=\$(awk 'BEGIN{i=0}{i = i + \$1}END{print i}' counted_unique_reads_*.txt)
    (seq 1 22; echo X; echo Y) | parallel -j 10 "graph_peak_caller callpeaks -q ${params.qvalue} -g graphs/{}.nobg -s json/${name}_{}.json -G ${params.genome_size} -p True -f \$fragment_length -r \$read_length -u \$unique_reads -n {}"

    rename 'touched' '_touched' *touched*
    rename 'background' '_background' *background*
    rename 'direct' '_direct' *direct*
    rename 'fragment' '_fragment' *fragment*
    rename 'pvalues' '_pvalues' *pvalues*

    (seq 1 22; echo X; echo Y) | parallel -j 3 "graph_peak_caller callpeaks_whole_genome_from_p_values -q ${params.qvalue} -d graphs/ -n '' -f \$fragment_length -r \$read_length {}"
    (seq 1 22; echo X; echo Y) | parallel -j 10 'graph_peak_caller peaks_to_linear {}_max_paths.intervalcollection graphs/{}_linear_pathv2.interval {} {}_linear_peaks.bed'
    cat *_linear_peaks.bed | awk '\$2<\$3' | sort-bed - > ${name}_peaks.narrowPeak
    tar cvzf ${name}.tar.gz *.intervalcollection *.fasta
"""
}
