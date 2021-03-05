rule busco:
    input:
        lambda wc: pep.get_sample(wc.sample).assembly_file_path
    output:
        directory("results/busco/{sample}/")
    params:
        lin=lambda wc: pep.get_sample(wc.sample).busco_lineage,
    conda:
        "../envs/busco.yaml"
    log:
        "results/busco/{sample}.log"
    resources:
        time=240,
        mem=32000,
        cpus=24
    threads:
        24
    shell:
        """
        busco -i {input} -l {params.lin} -o {wildcards.sample} --out_path {output} --download_path {output} -m genome --cpu {threads} 2> {log}
        """
