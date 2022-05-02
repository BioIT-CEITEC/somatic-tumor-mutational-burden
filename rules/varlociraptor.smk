

rule varlociraptor_alignment_properties:
    input:
        ref=expand("{ref_dir}/seq/{ref_name}.fa",ref_dir=reference_directory,ref_name=config["reference"])[0],
        ref_idx=expand("{ref_dir}/seq/{ref_name}.fa.fai",ref_dir=reference_directory,ref_name=config["reference"])[0],
        bam="mapped/{sample_name}.bam",
        # bam=lambda wildcards: expand("mapped/{sample_name}_{type}.bam",sample_name = sample_tab.sample_name,type=["tumor","normal"])[0],
    output:
        "varlociraptor/{sample_name}/varlociraptor/alignment-properties.json",
    log:
        "logs/{sample_name}/callers/varlociraptor_estimate_alignment_properties.log",
    conda: "../wrappers/varlociraptor/env.yaml"
    shell:
        "varlociraptor estimate alignment-properties {input.ref} --bam {input.bam} > {output} 2> {log}"

rule varlociraptor_preprocess:
    input:
        ref=expand("{ref_dir}/seq/{ref_name}.fa",ref_dir=reference_directory,ref_name=config["reference"])[0],
        ref_idx=expand("{ref_dir}/seq/{ref_name}.fa.fai",ref_dir=reference_directory,ref_name=config["reference"])[0],
        bam="mapped/{sample_name}.bam",
        bai="mapped/{sample_name}.bam.bai",
        alignment_props="varlociraptor/{sample_name}/varlociraptor/alignment-properties.json",
        candidates="variant_calls/{sample_name}/mutect2/MuTect2.vcf",
    output:
        "varlociraptor/{sample_name}/varlociraptor_observations.bcf",
    params:
        extra=config["varlociraptor_max_depth"],#"--max-depth 200", #overit jaky mame pokryti a jestli to dava smysl
    log:
        "logs/{sample_name}/callers/varlociraptor_preprocess.log",
    conda: "../wrappers/varlociraptor/env.yaml"
    shell:
        "varlociraptor preprocess variants {params.extra} --candidates {input.candidates} "
        "--alignment-properties {input.alignment_props} {input.ref} --bam {input.bam} --output {output} "
        "2> {log}"

rule varlociraptor_call_paired:
    input:
        normal=lambda wildcards: expand("varlociraptor/{val}/observations.bcf",val=sample_tab.loc[sample_tab.sample_name == wildcards.sample_name, "sample_name_normal"])[0],
        tumor=lambda wildcards: expand("varlociraptor/{val}/observations.bcf",val=sample_tab.loc[sample_tab.sample_name == wildcards.sample_name, "sample_name_tumor"])[0],
    output:
       "varlociraptor/{sample_name}/varlociraptor_calls.bcf",
    log:
        "logs/{sample_name}/callers/varlociraptor_calling.log",
    params:
        extra="",#config["params"]["varlociraptor"]["call"],
    conda: "../wrappers/varlociraptor/env.yaml"
    shell:
        "(varlociraptor call variants tumor-normal --tumor {input.tumor} --normal {input.normal} "
        " > {output}) 2> {log}"

# varlociraptor call variants tumor-normal --purity 0.75 --tumor tumor.bcf --normal normal.bcf > calls.bcf

rule sort_calls:
    input:
        "varlociraptor/{sample_name}/varlociraptor_calls.bcf",
    output:
        "varlociraptor/{sample_name}/varlociraptor_calls.sorted.bcf",
    log:
        "logs/{sample_name}/callers/bcftools_sort.log",
    conda: "../wrappers/varlociraptor/env.yaml"
    shell:
        "bcftools sort --temp-dir `mktemp -d` "
        "-Ob {input} > {output} 2> {log}"

rule bcf2vcf:
    input:
        "varlociraptor/{sample_name}/varlociraptor_calls.sorted.bcf",
    output:
        "varlociraptor/{sample_name}/varlociraptor_calls.vcf",
    log:
         "logs/{sample_name}/callers/bcftools_bcf2vcf.log",
    conda: "../wrappers/varlociraptor/env.yaml"
    shell:
        "bcftools view -O v {input} > {output} 2>{log} "


# ANNOTATE VARIANTS
# snpEFF
rule snpeff_varlociraptor:
    input:  "varlociraptor/{sample_name}/varlociraptor_calls.vcf",
    output: "varlociraptor/{sample_name}/varlociraptor_calls.annotated.vcf",
    log:    "logs/{sample_name}/callers/snpeff_varlociraptor.log",
    threads: 10
    conda:  "../wrappers/snpeff/env.yaml"
    shell: 
        """   
        bcftools view bcf | gzip -c > vcf
     
        snpEff -Xms16G GRCh37.75 -canon {input} > {output} 2>{log}
        """
# -fi intervals.bed

rule varlociraptor_TMB:
    input:  "varlociraptor/{sample_name}/varlociraptor_calls.annotated.vcf",
    output: "TMB_by_varlociraptor/{sample_name}_TMB.pdf",
    log:    "logs/{sample_name}/callers/varlociraptor_TMB.log",
    params:
        coding_genome_size= config["exome_genome_size"], #3.5e7,   
        events="SOMATIC_TUMOR",  #get_mutational_burden_events,
        sample="tumor",#nazev sloupce v VCF
    conda: "../wrappers/varlociraptor/env.yaml"
    shell:
        "(varlociraptor estimate mutational-burden "
        "--plot-mode hist "
        "--coding-genome-size {params.coding_genome_size} "
        "--events {params.events} "
        "--sample {params.sample} "
        "< {input} | vl2pdf > {output}) 2> {log}"

# "< {input} | vl2svg > {output}) 2> {log}"
        # "< {input} > {output}) 2> {log}"
