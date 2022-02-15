xpore_img = "docker://quay.io/biocontainers/xpore:2.1--pyh5e36f6f_0"
module_name = "xpore"

# we have to do an extra eventalign because xpore needs read_index instead of read_name
rule xpore_eventalign:
    input:
        fastq=rules.merge_fastq.output.fastq,
        index=rules.f5c_index.output.index,
        bam=rules.alignmemt_postfilter.output.bam,
        fasta=rules.get_transcriptome.output.fasta,
        kmer_model="resources/f5c/r9.4_70bps.u_to_t_rna.5mer.template.model",
    output:
        tsv=join("results", module_name, rule_name, "{cond}_{rep}_data.tsv"),
        summary=join("results", module_name, rule_name, "{cond}_{rep}_summary.tsv"),
    log:
        join("logs", module_name, rule_name, "{cond}_{rep}.log"),
    threads: get_threads(config, rule_name)
    params:
        opt=get_opt(config, rule_name),
    resources:
        mem_mb=lambda wildcards, attempt, mem=get_mem(config, rule_name): attempt * mem,
    container:
        f5c_container
    shell:
        "f5c eventalign {params.opt} -t {threads} --kmer-model {input.kmer_model} -r {input.fastq} -b {input.bam} -g {input.fasta} --summary {output.summary}  > {output.tsv} 2> {log}"


rule_name = "xpore_dataprep"


rule xpore_dataprep:
    input:
        eventalign=rules.xpore_eventalign.output.tsv,
    output:
        outdir=directory(join("results", module_name, rule_name, "{cond}_{rep}")),
    log:
        join("logs", module_name, rule_name, "{cond}_{rep}.log"),
    threads: get_threads(config, rule_name)
    params:
        opt=get_opt(config, rule_name),
    resources:
        mem_mb=lambda wildcards, attempt, mem=get_mem(config, rule_name): attempt * mem,
    container:
        xpore_img
    shell:
        """
        xpore dataprep --eventalign {input.eventalign} --n_processes {threads} \
            --out_dir {output.outdir} 2> {log}
        """