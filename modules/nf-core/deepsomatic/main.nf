process DEEPSOMATIC {
    tag "$meta.id"
    label 'process_high'

    container "docker.io/google/deepsomatic:1.10.0" // Keep in sync with versions_deepsomatic.
    input:
    tuple val(meta), path(input_normal), path(index_normal), path(input_tumor), path(index_tumor), path(intervals)
    tuple val(meta2), path(fasta)
    tuple val(meta3), path(fai)
    tuple val(meta4), path(gzi)

    output:
    tuple val(meta), path("${prefix}.vcf.gz")      ,  emit: vcf
    tuple val(meta), path("${prefix}.vcf.gz.tbi")  ,  emit: vcf_tbi
    tuple val(meta), path("${prefix}.g.vcf.gz")    ,  emit: gvcf, optional: true
    tuple val(meta), path("${prefix}.g.vcf.gz.tbi"),  emit: gvcf_tbi, optional: true
    tuple val("${task.process}"), val("deepsomatic"), val("1.10.0"), emit: versions_deepsomatic, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    if (workflow.profile.tokenize(',').intersect(['conda', 'mamba']).size() >= 1) {
        error "DEEPSOMATIC module does not support Conda. Please use Docker / Singularity / Podman instead."
    }
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.id}"
    def regions = intervals ? "--regions=${intervals}" : ""
    def gvcf = task.ext.output_gvcf ? "--output_gvcf=${prefix}.g.vcf.gz" : ""

    """
    python3 - \\
        --ref=${fasta} \\
        --reads_normal=${input_normal} \\
        --reads_tumor=${input_tumor} \\
        --output_vcf=${prefix}.vcf.gz \\
        ${gvcf} \\
        ${args} \\
        ${regions} \\
        --intermediate_results_dir=tmp \\
        --num_shards=${task.cpus} <<'PY'
    import argparse
    import os
    import re
    import sys

    import pysam

    parser = argparse.ArgumentParser(allow_abbrev=False)
    parser.add_argument('--reads_tumor')
    parser.add_argument('--sample_name_tumor')
    parser.add_argument('--make_examples_extra_args')
    parser.add_argument('--postprocess_variants_extra_args', action='append', default=[])
    options, remaining = parser.parse_known_args()

    def extra_args(value):
        return re.findall(r"[^,]+=[\\"'][^\\"']*[\\"']|[^,]+", value or '')

    sample_name = options.sample_name_tumor
    for argument in extra_args(options.make_examples_extra_args):
        key, value = argument.split('=', 1)
        if key.lstrip('-') == 'sample_name_tumor':
            sample_name = value.strip("\\"'")
    if not sample_name:
        with pysam.AlignmentFile(options.reads_tumor, check_sq=False) as alignment:
            read_groups = alignment.header.to_dict().get('RG', [])
        sample_names = {group['SM'] for group in read_groups if group.get('SM')}
        if len(sample_names) != 1 or any(not group.get('SM') for group in read_groups):
            sys.exit('DEEPSOMATIC requires one tumour read-group sample name; '
                     'supply --sample_name_tumor for missing or conflicting SM values.')
        sample_name = sample_names.pop()

    # Empty VCF-only calls have no records from which to recover the sample name.
    postprocess_args = [
        argument
        for value in options.postprocess_variants_extra_args
        for argument in extra_args(value)
        if argument.split('=', 1)[0].lstrip('-') != 'sample_name'
    ]
    postprocess_args.append('sample_name="' + sample_name + '"')
    remaining.extend(['--reads_tumor', options.reads_tumor])
    if options.sample_name_tumor is not None:
        remaining.extend(['--sample_name_tumor', options.sample_name_tumor])
    if options.make_examples_extra_args is not None:
        remaining.extend(['--make_examples_extra_args', options.make_examples_extra_args])
    remaining.extend(['--postprocess_variants_extra_args', ','.join(postprocess_args)])
    os.execvp('run_deepsomatic', ['run_deepsomatic', *remaining])
    PY
    """

    stub:
    if (workflow.profile.tokenize(',').intersect(['conda', 'mamba']).size() >= 1) {
        error "DEEPSOMATIC module does not support Conda. Please use Docker / Singularity / Podman instead."
    }
    prefix = task.ext.prefix ?: "${meta.id}"
    def gvcf = task.ext.output_gvcf ? "echo \"\" | gzip > ${prefix}.g.vcf.gz; touch ${prefix}.g.vcf.gz.tbi" : ""
    """
    echo "" | gzip > ${prefix}.vcf.gz
    touch ${prefix}.vcf.gz.tbi
    ${gvcf}
    """
}
