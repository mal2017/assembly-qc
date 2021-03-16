def get_fqs(wc):
    s= pep.get_sample(wc.sample)
    sub_ix =  pep.get_sample(wc.sample).subsample_name.index(wc.sub)
    fq_r1s = s.coverage_path[sub_ix] if len(s.subsample_name) > 1 else s.coverage_path
    #fq_r2s = s.coverage_r2_path[sub_ix] if len(s.subsample_name) > 1 else s.coverage_r2_path
    return([fq_r1s])

rule minimap2_index:
    input:
        target=lambda wc: pep.get_sample(wc.sample).assembly_file_path
    output:
        "results/idx/{sample}.mmi"
    log:
        "results/logs/minimap2_index/{sample}.log"
    resources:
        time=60,
        mem=20000,
        cpus=4
    conda:
        "../envs/minimap2.yaml"
    shell:
        """
        minimap2 -x map-ont -d {output} {input.target} 2> {log}
        """

rule minimap2_aln:
    input:
        target=rules.minimap2_index.output,
        query=get_fqs
    output:
        "results/aln/{sample}/{sample}.{sub}.srt.bam"
    log:
        "results/logs/minimap2/{sample}.{sub}.log"
    params:
        extra="-ax map-ont"  # optional
    threads:
        24
    resources:
        time=240,
        mem=20000,
        cpus=24
    conda:
        "../envs/minimap2.yaml"
    shell:
        """
        minimap2 {params.extra} -t {threads} {input.target} {input.query} 2> {log} | \
            samtools sort -O BAM -o {output}
        """

rule bam_index:
    input:
        "results/aln/{sample}/{sample}.{sub}.srt.bam"
    output:
        "results/aln/{sample}/{sample}.{sub}.srt.bam.bai"
    threads:
        8
    resources:
        time=20,
        mem=10000,
        cpus=8
    conda:
        "../envs/minimap2.yaml"
    shell:
        "samtools index -@ {threads} {input}"

rule taxdump:
    output:
        directory("results/blobtools/taxdump")
    resources:
        time=240,
        mem=10000,
        cpus=2
    shell:
        """
        mkdir {output}
        cd {output}
        curl -L ftp://ftp.ncbi.nih.gov/pub/taxonomy/new_taxdump/new_taxdump.tar.gz | tar xzf -;
        """

rule nt:
    output:
        directory("results/blobtools/nt")
    resources:
        time=240,
        mem=10000,
        cpus=2
    shell:
        """
        mkdir {output}
        wget "ftp://ftp.ncbi.nlm.nih.gov/blast/db/nt.??.tar.gz" -P {output} && \
        for file in {output}/*.tar.gz; \
            do tar xf $file -C {output} && rm $file; \
        done
        """

rule uniprot_dl:
    input:
        rules.nt.output,
        rules.taxdump.output
    output:
        "results/blobtools/uniprot-fa/reference_proteomes.tar.gz"
    conda:
        "../envs/aria2.yaml"
    shell:
        """
        aria2c -x4 -o {output} \
            ftp://ftp.ebi.ac.uk/pub/databases/uniprot/current_release/knowledgebase/reference_proteomes/$(curl \
                -vs ftp.ebi.ac.uk/pub/databases/uniprot/current_release/knowledgebase/reference_proteomes/ 2>&1 | \
                awk '/tar.gz/ {{print $9}}')
        """

rule uniprot_db:
    input:
        rules.uniprot_dl.output
    output:
        directory("results/blobtools/uniprot/")
    resources:
        time= 480,
        cpus=24,
        mem=20000,
    threads:
        24
    conda:
        "../envs/diamond.yaml"
    shell:
        """
        mkdir -p {output}
        cp {input} {output}/
        cd {output}
        tar xf reference_proteomes.tar.gz
        touch reference_proteomes.fasta.gz
        find . -mindepth 2 | grep "fasta.gz" | grep -v 'DNA' | grep -v 'additional' | xargs cat >> reference_proteomes.fasta.gz

        echo "accession\taccession.version\ttaxid\tgi" > reference_proteomes.taxid_map
        zcat */*/*.idmapping.gz | grep "NCBI_TaxID" | awk '{{print $1 "\t" $1 "\t" $3 "\t" 0}}' >> reference_proteomes.taxid_map

        diamond makedb -p {threads} --in reference_proteomes.fasta.gz --taxonmap reference_proteomes.taxid_map --taxonnodes ../taxdump/nodes.dmp -d reference_proteomes.dmnd
        """


rule blastn_tax:
    input:
        fa = lambda wc: pep.get_sample(wc.sample).assembly_file_path,
        db = rules.nt.output,
    output:
        "results/blobtools/blast/{sample}.blastn_tax.out"
    conda:
        "../envs/blast.yaml"
    resources:
        time=480,
        mem=10000,
        cpus=32
    threads:
        32
    shell:
        """
        blastn -db {input.db}/nt \
           -query {input.fa} \
           -outfmt '6 qseqid staxids bitscore std' \
           -max_target_seqs 10 \
           -max_hsps 1 \
           -evalue 1e-25 \
           -num_threads {threads} \
           -out {output}
        """

rule diamond_blastx_tax:
    input:
        uniprot = rules.uniprot_db.output,
        fa = lambda wc: pep.get_sample(wc.sample).assembly_file_path,
    output:
        o = "results/blobtools/blast/{sample}.diamond_blastx_tax.out",
	    t = temp(directory("results/blobtools/blast/tmp-diamond-{sample}"))
    conda:
        "../envs/diamond.yaml"
    resources:
        time=8000,
        mem=128000,
        cpus=16
    threads:
        16
    shell:
        """
	mkdir -p {output.t} &&
        diamond blastx \
            --query {input.fa} \
            --db {input.uniprot}/reference_proteomes.dmnd \
	    --tmpdir {output.t} \
            --outfmt 6 qseqid staxids bitscore \
            -b 12 -c1 \
            --verbose \
            --max-target-seqs 1 \
            --evalue 1e-25 \
            --threads {threads} \
            > {output.o}
        """

rule blobtools_create:
    input:
        diamond = rules.diamond_blastx_tax.output.o,
        blastn = rules.blastn_tax.output,
        bams = lambda wc: expand("results/aln/{s}/{s}.{sub}.srt.bam",s=wc.sample, sub=pep.get_sample(wc.sample).subsample_name),
        taxdump = rules.taxdump.output,
        fa = lambda wc: pep.get_sample(wc.sample).assembly_file_path,
        bais = lambda wc: expand("results/aln/{s}/{s}.{sub}.srt.bam.bai",s=wc.sample, sub=pep.get_sample(wc.sample).subsample_name),
    output:
        "results/blobtools/work/{sample}.blobDB.json",
    conda:
         "../envs/blobtools.yaml"
    resources:
        time=240,
        mem=32000,
        cpus=24
    threads:
        24
    params:
        cov = lambda wc: ["--bam results/aln/{s}/{s}.{sub}.srt.bam".format(s=wc.sample, sub=x) for x in pep.get_sample(wc.sample).subsample_name],
        o = "results/blobtools/work/{sample}"
    shell:
        """
        blobtools create -i {input.fa} \
            -t {input.diamond} \
            -t {input.blastn} \
            --nodes {input.taxdump}/nodes.dmp \
            --names {input.taxdump}/names.dmp \
            {params.cov} \
            -o {params.o}
        """


# rule get_taxid_mapping:
#     output:
#         "results/blobtools/taxids/uniprot.taxid.tsv"
#     shell:
#         """
#         wget -O - https://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/idmapping/idmapping.dat.gz | \
#          grep NCBI_TaxID | \
#          gunzip -c > {output}
#         """
#
# rule blobtools_taxify:
#     input:
#         hits = "results/blobtools/blast/{sample}.{blast_type}.out",
#         mappings = rules.get_taxid_mapping.output,
#     output:
#         "results/blobtools/blast/{sample}.{blast_type}.hits.tsv"
#     conda:
#         "../envs/blobtools.yaml"
#     shell:
#         """
#         blobtools taxify -f {input.hits} \
#             -m {input.mappings} \
#             -s 0 -t2
#         """
#
# rule blobtools2_create:
#     input:
#         lambda wc: pep.get_sample(wc.sample).assembly_file_path
#     output:
#         directory("results/blobtools/work/{sample}")
#     shell:
#         """
#         ~/blobtoolkit/blobtools2/blobtools create --fasta {input} {output}
#         """
#
# rule blobtools2_add_blastn:
#     input:
#         blast = rules.blastn_tax.output,
#         taxdump = rules.taxdump.output,
#         dir = rules.blobtools2_create.output,
#     output:
#         touch('results/blobtools/{sample}.blastn_added')
#     shell:
#         """
#         ~/blobtoolkit/blobtools2/blobtools add \
#             --hits {input.blast} \
#             --taxrule bestsumorder \
#             --taxdump {input.taxdump} \
#             {output}
#         """
#
# rule blobtools2_add_blastx:
#     input:
#         blast = rules.diamond_blastx_tax.output,
#         taxdump = rules.taxdump.output,
#         dir = rules.blobtools2_create.output,
#     output:
#         touch('results/blobtools/{sample}.blastx_added')
#     shell:
#         """
#         ~/blobtoolkit/blobtools2/blobtools add \
#             --hits {input.blast} \
#             --taxrule bestsumorder \
#             --taxdump {input.taxdump} \
#             {output}
#         """
#
# rule blobtools2_add_coverage:
#     input:
#         bams = lambda wc: expand("results/aln/{s}/{s}.{sub}.srt.bam",s=wc.sample, sub=pep.get_sample(wc.sample).subsample_name),
#         dir = rules.blobtools2_create.output,
#     params:
#         cov = lambda wc: ["--cov results/aln/{s}/{s}.{sub}.srt.bam".format(s=wc.sample, sub=x) for x in pep.get_sample(wc.sample).subsample_name]
#     output:
#         touch('results/blobtools/{sample}.cov_added')
#     shell:
#         """
#         ~/blobtoolkit/blobtools2/blobtools add \
#             {params.cov} \
#             {input.dir}
#         """
#
# rule blobtools2_add_busco:
#     input:
#         busco = rules.busco.output,
#         dir = rules.blobtools2_create.output,
#     params:
#         lin = lambda wc: pep.get_sample(wc.sample).busco_lineage
#     output:
#         touch('results/blobtools/{sample}.busco_added')
#     shell:
#         """
#         ~/blobtoolkit/blobtools2/blobtools add \
#             --busco {input.busco}/{wildcards.sample}/run_{params.lin}/full_table.tsv \
#             {input.dir}
#         """
#
# rule blobtools2_add_fasta:
#     input:
#         fa = lambda wc: pep.get_sample(wc.sample).assembly_file_path,
#         dir = rules.blobtools2_create.output,
#     output:
#         touch('results/blobtools/{sample}.fasta_added')
#     shell:
#         """
#         ~/blobtoolkit/blobtools2/blobtools add \
#             --fasta {input.fa} \
#             {input.dir}
#         """
