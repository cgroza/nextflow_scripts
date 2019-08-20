params.fastq_dir = "fastq"
params.control_fastq = "control.fastq.gz"

params.ref_graph = "ref/"
params.pop_graph = "pop/"
params.ref_name = "ref"
params.pop_name = "pop"

params.genome_size = 2480000000
params.fragment_length = 200
params.read_length = 36

chromosomes = "chr1,chr2,chr3,chr4,chr5,chr6,chr7,chr8,chr9,chr10,chr11,chr12,chr13,chr14,chr15,chr16,chr17,chr18,chr19,chr20,chr21,chr22,chrX,chrY"

Channel.fromPath("${params.pop_graph}/graphs/*.vg").set{linear_vg_ch}
Channel.fromPath("${params.ref_graph}/graphs/*.vg").set{ref_linear_vg_ch}
Channel.fromPath("${params.fastq_dir}/*").into{fastq_ch; ref_fastq_ch}
Channel.fromPath(params.control_fastq).into{control_fastq_ch; ref_control_fastq_ch}

Channel.fromPath(
    ["${params.ref_graph}/${params.ref_name}.xg",
     "${params.ref_graph}/${params.ref_name}.gcsa",
     "${params.ref_graph}/${params.ref_name}.gcsa.lcp"]).into{ref_index_treatment_ch; ref_index_control_ch}

Channel.fromPath(
    ["${params.pop_graph}/${params.pop_name}.xg",
     "${params.pop_graph}/${params.pop_name}.gbwt",
     "${params.pop_graph}/${params.pop_name}.gcsa",
     "${params.pop_graph}/${params.pop_name}.gcsa.lcp"]).into{pop_index_treatment_ch; pop_index_control_ch}



process vgToJsonPop {
    cpus = 40
    memory '100 GB'
    time '12h'

    input:
    file "graphs/*" from linear_vg_ch.collect()

    output:
    file "graphs/*" into vg2json_ch

    script:
    """
    (seq 1 22; echo X; echo Y) | parallel -j 3 'vg view -Vj graphs/chr{}.vg > graphs/chr{}.json ; graph_peak_caller create_ob_graph graphs/chr{}.json ; vg stats -r graphs/chr{}.vg  | cut -f 2 > graphs/node_range_chr{}.txt'
"""
}

process vgToJsonRef {
    cpus = 40
    memory '100 GB'
    time '12h'

    input:
    file "graphs/*" from ref_linear_vg_ch.collect()

    output:
    file "graphs/*" into ref_vg2json_ch

    script:
    """
    (seq 1 22; echo X; echo Y) | parallel -j 3 'vg view -Vj graphs/chr{}.vg > graphs/chr{}.json ; graph_peak_caller create_ob_graph graphs/chr{}.json ; vg stats -r graphs/chr{}.vg  | cut -f 2 > graphs/node_range_chr{}.txt'
"""
}

process linearPathsPop {
    cpus = 40
    memory '100 GB'
    time '12 h'

    input:
    file "graphs/*" from vg2json_ch.collect()

    output:
    file "graphs/*" into control_linear_ch, treatment_linear_ch, peak_linear_ch

    script:
    """
    (seq 1 22; echo X; echo Y) | parallel -j 6 graph_peak_caller find_linear_path -g graphs/chr{}.nobg graphs/chr{}.json chr{} graphs/chr{}_linear_pathv2.interval
"""
}


process linearPathsRef {
    cpus = 40
    memory '100 GB'
    time '12 h'

    input:
    file "graphs/*" from ref_vg2json_ch.collect()

    output:
    file "graphs/*" into ref_control_linear_ch, ref_treatment_linear_ch, ref_peak_linear_ch

    script:
    """
    (seq 1 22; echo X; echo Y) | parallel -j 6 graph_peak_caller find_linear_path -g graphs/chr{}.nobg graphs/chr{}.json chr{} graphs/chr{}_linear_pathv2.interval
"""
}

process alignControlRef {
    cpus = 40
    memory '100 GB'
    time '12h'

    input:
    set file(xg), file(gcsa), file(gcsa_lcp), file(fastq), file("graphs/*") from ref_index_control_ch.collect().combine(ref_control_fastq_ch).combine(ref_control_linear_ch.collect()).view()

    output:
    file("control_json/*") into ref_control_json_ch

    script:
    name = fastq.getSimpleName()
    """
    mkdir control_gam
    vg map -f $fastq -x $xg -g $gcsa -t 40 -u 1 -m 1 > "control_gam/${name}_ref.gam"

    mkdir control_json
    vg view -aj control_gam/${name}_ref.gam > control_json/${name}_ref.json
    graph_peak_caller split_vg_json_reads_into_chromosomes ${chromosomes} control_json/${name}_ref.json graphs/
    rm control_json/${name}_ref.json
"""
}


process alignControlPop {
    cpus = 40
    memory '100 GB'
    time '12h'

    input:
    set file(xg), file(gbwt), file(gcsa), file(gcsa_lcp), file(fastq), file("graphs/*") from pop_index_control_ch.collect().combine(control_fastq_ch).combine(control_linear_ch.collect()).view()
    output:
    file("control_json/*") into control_json_ch

    script:
    name = fastq.getSimpleName()
    """
    mkdir control_gam
    vg map -f $fastq -1 $gbwt -x $xg -g $gcsa -t 40 -u 1 -m 1 > "control_gam/${name}_pop.gam"

    mkdir control_json
    vg view -aj control_gam/${name}_pop.gam > control_json/${name}_pop.json
    graph_peak_caller split_vg_json_reads_into_chromosomes ${chromosomes} control_json/${name}_pop.json graphs/
    rm control_json/${name}_pop.json
"""
}

process alignSampleRef {
    cpus = 40
    memory '100 GB'
    time '12h'

    input:
    set file(xg), file(gcsa), file(gcsa_lcp), file(fastq), file("graphs/*") from ref_index_treatment_ch.collect().combine(ref_fastq_ch).combine(ref_treatment_linear_ch.collect()).view()

    output:
    set val(name), file("json/*") into ref_treatment_json_ch
    set val(name), file("gam/${name}_ref.gam") into ref_treatment_gam_ch

    script:
    name = fastq.getSimpleName()
    """
    mkdir gam
    vg map -f $fastq -x $xg -g $gcsa -t 40 -u 1 -m 1 > gam/${name}_ref.gam

    mkdir json
    vg view -aj gam/${name}_ref.gam > json/${name}_ref.json
    graph_peak_caller split_vg_json_reads_into_chromosomes ${chromosomes} json/${name}_ref.json graphs/
    rm json/${name}_ref.json
"""
}

process sortSampleRef {
    cpus = 40
    memory '100 GB'
    time = '12h'

    publishDir "$params.outDir", pattern: "ref_${name}.sorted.gam"

    input:
    set val(name), file(gam) from ref_treatment_gam_ch

    output:
    file "ref_${name}.sorted.gam"

    script:
    """
    vg gamsort ${gam} -i ${name}.sorted.gam.gai -t 40 > ref_${name}.sorted.gam
"""
}

process alignSamplePop {
    cpus = 40
    memory '100 GB'
    time '12h'

    input:
    set file(xg), file(gbwt), file(gcsa), file(gcsa_lcp), file(fastq), file("graphs/*") from pop_index_treatment_ch.collect().combine(fastq_ch).combine(treatment_linear_ch.collect()).view()

    output:
    set val(name), file("json/*") into treatment_json_ch
    set val(name), file("gam/${name}_pop.gam") into treatment_gam_ch

    script:
    name = fastq.getSimpleName()
    """
    mkdir gam
    vg map -f $fastq -1 $gbwt -x $xg -g $gcsa -t 40 -u 1 -m 1 > gam/${name}_pop.gam

    mkdir json
    vg view -aj gam/${name}_pop.gam > json/${name}_pop.json
    graph_peak_caller split_vg_json_reads_into_chromosomes ${chromosomes} json/${name}_pop.json graphs/
    rm json/${name}_pop.json
"""
}

process sortSamplePop {
    cpus = 40
    memory '100 GB'
    time = '12h'

    publishDir "$params.outDir", pattern: "${name}.sorted.gam"

    input:
    set val(name), file(gam) from treatment_gam_ch

    output:
    file "${name}.sorted.gam"

    script:
    """
    vg gamsort ${gam} -i ${name}.sorted.gam.gai -t 40 > ${name}.sorted.gam
"""
}

process callPeaksPop{
    cpus = 40
    memory '100 GB'
    time '12h'

    publishDir "$params.outDir", pattern: "ref_${sample}_peaks.narrowPeak"

    input:
    set val(name), file("json/*"), file("control_json/*"), file("graphs/*") from treatment_json_ch.combine(control_json_ch.collect()).combine(peak_linear_ch.collect()).view()

    output:
    file "ref_${sample}_peaks.narrowPeak"

    script:
    """
    graph_peak_caller count_unique_reads ${chromosomes} graphs/ json/${name}_pop | tail -n 1 > counted_unique_reads.txt
    unique_reads=\$(cat counted_unique_reads.txt)
    (seq 1 22; echo X; echo Y) | parallel -j 6 'graph_peak_caller callpeaks -g graphs/chr{}.nobg -s json/${name}_pop{}.json -c control_json/${name}_pop.json -G ${params.genome_size} -p True -f ${params.fragment_length} -r ${params.read_length} -u \$unique_reads -n chr{}'
    (seq 1 22; echo X; echo Y) | parallel -j 6 'graph_peak_caller callpeaks_whole_genome_from_p_values -d graphs/ -n "" -f ${params.fragment_length} -r ${params.read_length} chr{}'
    (seq 1 22; echo X; echo Y) | graph_peak_caller peaks_to_linear chr{}_max_paths.intervalcollection graphs/chr{}_linear_pathv2.interval chr{} chr{}_linear_peaks.bed
    cat *_linear_peaks.bed | sort-bed - > ${sample}_peaks.narrowPeak
"""
}

process callPeaksRef{
    cpus = 40
    memory '100 GB'
    time '12h'

    publishDir "$params.outDir", pattern: "ref_${sample}_peaks.narrowPeak"

    input:
    set val(name), file("json/*"), file("control_json/*"), file("graphs/*") from ref_treatment_json_ch.combine(ref_control_json_ch.collect()).combine(ref_peak_linear_ch.collect()).view()

    output:
    file "ref_${sample}_peaks.narrowPeak"

    script:
    """
    graph_peak_caller count_unique_reads ${chromosomes} graphs/ json/${name}_pop | tail -n 1 > counted_unique_reads.txt
    unique_reads=\$(cat counted_unique_reads.txt)
    (seq 1 22; echo X; echo Y) | parallel -j 6 'graph_peak_caller callpeaks -g graphs/chr{}.nobg -s json/${name}_pop{}.json -c control_json/${name}_pop.json -G ${params.genome_size} -p True -f ${params.fragment_length} -r ${params.read_length} -u \$unique_reads -n chr{}'
    (seq 1 22; echo X; echo Y) | parallel -j 6 'graph_peak_caller callpeaks_whole_genome_from_p_values -d graphs/ -n "" -f ${params.fragment_length} -r ${params.read_length} chr{}'
    (seq 1 22; echo X; echo Y) | graph_peak_caller peaks_to_linear chr{}_max_paths.intervalcollection graphs/chr{}_linear_pathv2.interval chr{} chr{}_linear_peaks.bed
    cat *_linear_peaks.bed | sort-bed - > ref_${sample}_peaks.narrowPeak
"""
}
