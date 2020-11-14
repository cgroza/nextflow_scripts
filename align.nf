params.xg = "EU_AF_index.xg"
params.gbwt = false
params.gcsa = "EU_AF_index.gcsa"
params.gcsa_lcp = "${params.gcsa}.lcp"
params.reads = "reads"
params.paired = true
params.prefix = "aligned_"

xg_ch = Channel.fromPath(params.xg)
gcsa_ch = Channel.fromPath(params.gcsa)
gcsa_lcp_ch = Channel.fromPath(params.gcsa_lcp)

Channel.fromPath(params.reads).splitCsv(header:true).map{row -> file(row.path)}.view().set{reads_ch}

vg_flag = ""
if(params.paired) {
	vg_flag = "-i"
}

file("gams").mkdir()

if(params.gbwt) {
	gbwt_ch = Channel.fromPath(params.gbwt)
	process alignFastq {
		time '1d'
		cpus 40
		memory '100 GB'
		publishDir "$workflow.launchDir/gams", mode: 'move'

		input:
		set file("index.xg"), file("index.gbwt"), file('index.gcsa'), file('index.gcsa.lcp'), file(reads) from xg_ch.combine(gbwt_ch).combine(gcsa_ch).combine(gcsa_lcp_ch).combine(reads_ch)

		output:
		file "${reads_name}.gam" into aln

		script:
		extension = reads.getExtension()
		reads_name = params.prefix + reads.getSimpleName()
		switch(extension) {
			case "gz":
				vg_flag = "${vg_flag} -f"
				break
			case "bam":
				vg_flag = "${vg_flag} -b"
				break
			case "gam":
				vg_flag = "${vg_flag} -G"
				break
		}
	"""
		vg map --threads 40 --gcsa-name index.gcsa --xg-name index.xg --gbwt-name index.gbwt ${vg_flag} ${reads} > ${reads_name}.gam
		"""
	}
}
else {
	process alignFastqNoGBWT {
		time '1d'
		cpus 40
		memory '100 GB'
		publishDir "$workflow.launchDir/gams", mode: 'move'

		input:
		set file("index.xg"), file('index.gcsa'), file('index.gcsa.lcp'), file(reads) from xg_ch.combine(gcsa_ch).combine(gcsa_lcp_ch).combine(reads_ch)

		output:
		file "${reads_name}.gam" into aln

		script:
		extension = reads.getExtension()
		reads_name = params.prefix + reads.getSimpleName()
		switch(extension) {
			case "gz":
				vg_flag = "${vg_flag} -f"
				break
			case "bam":
				vg_flag = "${vg_flag} -b"
				break
			case "gam":
				vg_flag = "${vg_flag} -G"
				break
		}
	"""
		vg map --threads 40 --gcsa-name index.gcsa --xg-name index.xg ${vg_flag} ${reads} > ${reads_name}.gam
		"""
	}

}
