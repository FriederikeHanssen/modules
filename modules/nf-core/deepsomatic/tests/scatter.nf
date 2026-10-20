include { DEEPSOMATIC } from '../main'
include { GATK4_MERGEVCFS } from '../../gatk4/mergevcfs/main'

workflow DEEPSOMATIC_SCATTER {
    take:
    samples
    fasta
    fai
    gzi

    main:
    DEEPSOMATIC(samples, fasta, fai, gzi)
    vcfs = DEEPSOMATIC.out.vcf.map { meta, vcf -> vcf }.collect()
        .map { files -> [[id: 'pair'], files] }
    GATK4_MERGEVCFS(vcfs, [[:], []])

    emit:
    shards = DEEPSOMATIC.out.vcf
    shard_indexes = DEEPSOMATIC.out.vcf_tbi
    vcf = GATK4_MERGEVCFS.out.vcf
    tbi = GATK4_MERGEVCFS.out.tbi
}
