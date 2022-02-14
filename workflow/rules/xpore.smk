xpore_img = "docker://quay.io/biocontainers/xpore:2.1--pyh5e36f6f_0"
module_name = "xpore"
rule_name = "xpore_dataprep"


rule xpore_dataprep:
    input:
        eventalign=rules.f5c_eventalign.output.tsv,
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
        xpore dataprep --eventalign {input.eventalign} –-n_processes {threads} \
            --out_dir {output.outdir} 2> {log}
        """