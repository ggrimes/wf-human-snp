#!/usr/bin/env nextflow

// Developer notes
// 
// This template workflow provides a basic structure to copy in order
// to create a new workflow. Current recommended pratices are:
//     i) create a simple command-line interface.
//    ii) include an abstract workflow scope named "pipeline" to be used
//        in a module fashion.
//   iii) a second concreate, but anonymous, workflow scope to be used
//        as an entry point when using this workflow in isolation.

nextflow.enable.dsl = 2

params.help = ""
params.bam = ""
params.ref = ""
params.bed = ""
params.output = "output"
params.threads = 8
params.download = ""
params.vcf = ""
params.vcf_bed = ""
params.python="python3"
params.pypy="pypy3"
params.ref_pct_full = 0.3
params.var_pct_full = 0.3
params.GVCF = "False"
params.snp_min_af = 0.0
params.indel_min_af = 0.0
params.sample_name = "SAMPLE"
params.vcf_fn = "EMPTY"
params.ctg_name = "EMPTY"
params.include_all_ctgs = "False"


def helpMessage(){
    log.info """
Workflow template'

Usage:
    nextflow run epi2melabs/wf-clair [options]

Script Options:
    --bam       FILE    Path to bam file
    --bed       FILE    Path to bed file
    --ref       FiLE    Path to ref file
    --output    DIR     Output path
"""
}

process run_clair3_stage_1 {
    label "clair3"
    input:
        path bam
        path bai
        path ref
        path fai
        path bed
    output:
        path "clair_output"
        path "clair_output/pileup.vcf.gz"
        path "clair_output/pileup.vcf.gz.tbi"
        path "clair_output/log/parallel_1_call_var_bam_pileup.log"
        path "clair_output/log/1_call_var_bam_pileup.log"
        path "clair_output/tmp/CONTIGS"
        path "clair_output/tmp/CHUNK_LIST"
        path "clair_output/tmp/split_beds/*"
        path "clair_output/tmp/pileup_output"
        path "clair_output/tmp/merge_output"
        path "clair_output/tmp/phase_output"
        path "clair_output/tmp/gvcf_tmp_output"
        path "clair_output/tmp/full_alignment_output"
        path "clair_output/tmp/phase_output/phase_vcf"
        path "clair_output/tmp/phase_output/phase_bam"
        path "clair_output/tmp/full_alignment_output/candidate_bed"
        path "clair_output/tmp/pileup_output/pileup_*_*.vcf"
    shell:
        '''
        mkdir -p clair_output/log
        touch clair_output/log/parallel_1_call_var_bam_pileup.log
        . env.sh !{bam} !{bed} !{ref} "clair_output" !{params.threads}
        !{params.python} $(which clair3.py) CheckEnvs \
            --bam_fn ${BAM} \
            --bed_fn !{bed} \
            --output_fn_prefix \$OUTPUT_FOLDER \
            --ref_fn !{ref} \
            --vcf_fn !{params.vcf_fn} \
            --ctg_name !{params.ctg_name} \
            --chunk_num 0 \
            --chunk_size 5000000 \
            --include_all_ctgs !{params.include_all_ctgs} \
            --threads !{params.threads}  \
            --python !{params.python} \
            --pypy !{params.pypy} \
            --samtools samtools \
            --whatshap whatshap \
            --parallel parallel \
            --qual 2 \
            --sampleName !{params.sample_name} \
            --var_pct_full !{params.var_pct_full} \
            --ref_pct_full !{params.ref_pct_full} \
            --snp_min_af !{params.snp_min_af} \
            --indel_min_af !{params.indel_min_af}

        readarray -t CHR < "${OUTPUT_FOLDER}/tmp/CONTIGS"
        cd \$OUTPUT_FOLDER
        CUDA_VISIBLE_DEVICES=""
        echo "bed, bed, ref: $BAM , $BED , $REF"
        parallel \
            --retries 4 \
            -C ' ' \
            --joblog \$OUTPUT_FOLDER/log/parallel_1_call_var_bam_pileup.log \
            -j 4 \
            --env BED \
            --env BAM \
            --env REF \
            --env OUTPUT_FOLDER \
            'echo "hello: ${BED}" && !{params.python} $(which clair3.py) CallVarBam \
                --chkpnt_fn ${CONDA_DIR}/bin/models/ont/pileup \
                --bam_fn "${BAM}" \
                --call_fn $OUTPUT_FOLDER/tmp/pileup_output/pileup_{1}_{2}.vcf \
                --sampleName !{params.sample_name} \
                --ref_fn "${REF}"\
                --extend_bed ${OUTPUT_FOLDER}/tmp/split_beds/{1} \
                --bed_fn "${BED}" \
                --vcf_fn !{params.vcf_fn} \
                --ctgName {1} \
                --chunk_id {2} \
                --chunk_num {3} \
                --snp_min_af !{params.snp_min_af} \
                --indel_min_af !{params.indel_min_af} \
                --temp_file_dir ${OUTPUT_FOLDER}/tmp/gvcf_tmp_output \
                --pileup' \
            :::: ${OUTPUT_FOLDER}/tmp/CHUNK_LIST \
            |& tee ${LOG_PATH}/1_call_var_bam_pileup.log

        !{params.pypy} $(which clair3.py) SortVcf \
            --input_dir ${OUTPUT_FOLDER}/tmp/pileup_output \
            --vcf_fn_prefix pileup \
            --output_fn ${OUTPUT_FOLDER}/pileup.vcf \
            --sampleName !{params.sample_name} \
            --ref_fn ${REF}

        if [ "$( gzip -fdc ${OUTPUT_FOLDER}/pileup.vcf.gz | grep -v '#' | wc -l )" -eq 0 ]; then echo "[INFO] Exit in pileup variant calling"; exit 0; fi
        '''
}

process run_clair3_stage_2 {
    label "clair3"
    input:
        path bam
        path bai
        path ref
        path fai
        path bed
        path "clair_output"
    output:
        path "clair_output"
        path "clair_output/tmp/phase_output/phase_vcf/phase_qual"
        path "clair_output/log/parallel_2_select_hetero_snp.log"
        path "clair_output/log/2_select_hetero_snp.log"
    shell:
        '''
        . env.sh !{bam} !{bed} !{ref} "clair_output"
        echo $''
        echo "[INFO] 2/7 Select heterozygous SNP variants for Whatshap phasing and haplotagging"
        readarray -t CHR < "${OUTPUT_FOLDER}/tmp/CONTIGS"
        gzip \
            -fdc ${OUTPUT_FOLDER}/pileup.vcf.gz | \
            !{params.pypy} $(which clair3.py) SelectQual \
                --phase \
                --output_fn ${OUTPUT_FOLDER}/tmp/phase_output/phase_vcf

        time \
            parallel \
                --retries 4 \
                --joblog ${LOG_PATH}/parallel_2_select_hetero_snp.log \
                --env OUTPUT_FOLDER \
                -j !{params.threads} \
                "echo '{1}'; !{params.pypy} $(which clair3.py) SelectHetSnp \
                    --vcf_fn ${OUTPUT_FOLDER}/pileup.vcf.gz \
                    --split_folder ${OUTPUT_FOLDER}/tmp/phase_output/phase_vcf \
                    --ctgName {1}" \
                ::: ${CHR[@]} ::: ${ALL_SAMPLE[@]} \
                |& tee ${LOG_PATH}/2_select_hetero_snp.log
        '''
}


process run_clair3_stage_3{
    label "clair3"
    input:
        path bam
        path bai
        path ref
        path fai
        path bed
        path "clair_output"
    output:
        path "clair_output"
        path "clair_output/log/parallel_3_phase.log"
        path "clair_output/log/parallel_3_phase.log"
        path "clair_output/tmp/phase_output/phase_vcf/phased_*.vcf.gz"
    shell:
        '''
        echo $''
        echo "[INFO] 3/7 Phase VCF file using Whatshap"
        . env.sh !{bam} !{bed} !{ref} "clair_output"
        readarray -t CHR < "${OUTPUT_FOLDER}/tmp/CONTIGS"
        parallel \
            --retries 4 \
            --joblog ${LOG_PATH}/parallel_3_phase.log \
            --env OUTPUT_FOLDER \
            --env BAM \
            --env REF \
            -j !{params.threads} \
            "whatshap phase \
                --output ${OUTPUT_FOLDER}/tmp/phase_output/phase_vcf/phased_{1}.vcf.gz \
                --reference ${REF} \
                --chromosome {1} \
                --distrust-genotypes \
                --ignore-read-groups \
                ${OUTPUT_FOLDER}/tmp/phase_output/phase_vcf/{1}.vcf \
                ${BAM}" \
            ::: ${CHR[@]} |& tee ${LOG_PATH}/3_phase.log

        parallel \
            --j !{params.threads} \
            "tabix \
                -f \
                -p vcf \
                ${OUTPUT_FOLDER}/tmp/phase_output/phase_vcf/phased_{1}.vcf.gz" \
            ::: ${CHR[@]}

        '''
}

process run_clair3_stage_4{
    label "clair3"
    input:
        path bam
        path bai
        path ref
        path fai
        path bed
        path "clair_output"
    output:
        path "clair_output"
        path "clair_output/log/parallel_4_haplotag.log"
        path "clair_output/log/4_haplotag.log"
        path "clair_output/tmp/phase_output/phase_bam/*.bam"
        path "clair_output/tmp/phase_output/phase_vcf/phased_*.vcf.gz"
    shell:
        '''
        echo $''
        echo "[INFO] 4/7 Haplotag input BAM file using Whatshap"
        . env.sh !{bam} !{bed} !{ref} "clair_output"
        readarray -t CHR < "${OUTPUT_FOLDER}/tmp/CONTIGS"
        time \
            parallel \
                --retries 4 \
                --joblog ${LOG_PATH}/parallel_4_haplotag.log \
                --env REF \
                --env OUTPUT_FOLDER \
                -j !{params.threads} \
            "whatshap haplotag \
                --output ${OUTPUT_FOLDER}/tmp/phase_output/phase_bam/{1}.bam  \
                --reference ${REF} \
                --ignore-read-groups \
                --regions {1} \
                ${OUTPUT_FOLDER}/tmp/phase_output/phase_vcf/phased_{1}.vcf.gz \
                ${BAM}" \
            ::: ${CHR[@]} |& tee ${LOG_PATH}/4_haplotag.log

        parallel \
            -j !{params.threads} \
            --env OUTPUT_FOLDER \
            "samtools index \
                -@12 \
                ${OUTPUT_FOLDER}/tmp/phase_output/phase_bam/{1}.bam" \
        ::: ${CHR[@]}
        '''
}

process run_clair3_stage_5{
    label "clair3"
    input:
        path bam
        path bai
        path ref
        path fai
        path bed
        path "clair_output"
    output:
        path "clair_output"
        path "clair_output/log/parallel_5_select_candidate.log"
        path "clair_output/log/5_select_candidate.log"
        path "clair_output/tmp/full_alignment_output/candidate_bed/FULL_ALN_FILE_*"
        path "clair_output/tmp/full_alignment_output/candidate_bed/qual"
    shell:
        '''
        # Full alignment calling
        echo $''
        echo "[INFO] 5/7 Select candidates for full-alignment calling"
        . env.sh !{bam} !{bed} !{ref} "clair_output"
        readarray -t CHR < "${OUTPUT_FOLDER}/tmp/CONTIGS"
        gzip -fdc ${OUTPUT_FOLDER}/pileup.vcf.gz | \
            !{params.pypy} $(which clair3.py) SelectQual \
                --output_fn ${OUTPUT_FOLDER}/tmp/full_alignment_output/candidate_bed \
                --var_pct_full !{params.var_pct_full} \
                --ref_pct_full !{params.ref_pct_full} \
                --platform ont \
                --vcf_fn !{params.vcf_fn}

        time \
            parallel \
                --retries 4 \
                --joblog ${LOG_PATH}/parallel_5_select_candidate.log \
                --env OUTPUT_FOLDER \
                --env REF \
                -j !{params.threads} \
                "!{params.pypy} $(which clair3.py) SelectCandidates \
                    --pileup_vcf_fn ${OUTPUT_FOLDER}/pileup.vcf.gz \
                    --split_folder ${OUTPUT_FOLDER}/tmp/full_alignment_output/candidate_bed \
                    --ref_fn ${REF} \
                    --var_pct_full !{params.var_pct_full} \
                    --ref_pct_full !{params.ref_pct_full} \
                    --platform ont\
                    --ctgName {1}" \
                ::: ${CHR[@]}  |& tee ${LOG_PATH}/5_select_candidate.log
        '''
}

process run_clair3_stage_6 {
    label "clair3"
    input:
        path bam
        path bai
        path ref
        path fai
        path bed
        path "clair_output"
    output:
        path "clair_output"
        path "clair_output/tmp/full_alignment_output/candidate_bed/FULL_ALN_FILES"
        path "clair_output/log/parallel_6_call_var_bam_full_alignment.log"
        path "clair_output/log/6_call_var_bam_full_alignment.log"
        path "clair_output/full_alignment.vcf.gz"
        path "clair_output/full_alignment.vcf.gz.tbi"
    shell:
        '''
        echo $''
        echo "[INFO] 6/7 Call low-quality variants using full-alignment model"
        . env.sh !{bam} !{bed} !{ref} "clair_output"
        cat ${OUTPUT_FOLDER}/tmp/full_alignment_output/candidate_bed/FULL_ALN_FILE_* > ${OUTPUT_FOLDER}/tmp/full_alignment_output/candidate_bed/FULL_ALN_FILES

        time \
            parallel \
                --retries 4 \
                --joblog ${LOG_PATH}/parallel_6_call_var_bam_full_alignment.log \
                --env OUTPUT_FOLDER \
                --env REF \
                -j !{params.threads} \
                "!{params.python} $(which clair3.py) CallVarBam \
                    --chkpnt_fn ${CONDA_DIR}/bin/models/ont/full_alignment \
                    --bam_fn ${OUTPUT_FOLDER}/tmp/phase_output/phase_bam/{1/.}.bam \
                    --call_fn ${OUTPUT_FOLDER}/tmp/full_alignment_output/full_alignment_{1/}.vcf \
                    --sampleName !{params.sample_name} \
                    --vcf_fn !{params.vcf_fn} \
                    --ref_fn ${REF} \
                    --full_aln_regions {1} \
                    --ctgName {1/.} \
                    --add_indel_length \
                    --phasing_info_in_bam \
                    --gvcf !{params.GVCF} \
                    --python !{params.python} \
                    --pypy !{params.pypy} \
                    --samtools samtools \
                    --platform ont" \
                :::: ${OUTPUT_FOLDER}/tmp/full_alignment_output/candidate_bed/FULL_ALN_FILES |& tee ${LOG_PATH}/6_call_var_bam_full_alignment.log

        !{params.pypy} $(which clair3.py) SortVcf \
            --input_dir ${OUTPUT_FOLDER}/tmp/full_alignment_output \
            --vcf_fn_prefix "full_alignment" \
            --output_fn ${OUTPUT_FOLDER}/full_alignment.vcf \
            --sampleName !{params.sample_name} \
            --ref_fn ${REF}

        if [ "$( gzip -fdc ${OUTPUT_FOLDER}/full_alignment.vcf.gz | grep -v '#' | wc -l )" -eq 0 ];
            then echo "[INFO] Exit in full-alignment variant calling";
            exit 0;
        fi

        if [ !{params.GVCF} == True ];
            then cat ${OUTPUT_FOLDER}/tmp/gvcf_tmp_output/*.tmp.g.vcf | \
                ${PYPY} ${CLAIR3} SortVcf
                    --output_fn ${OUTPUT_FOLDER}/tmp/gvcf_tmp_output/non_var.gvcf;
        fi
        '''
}

process run_clair3_stage_7{
    label "clair3"
    input:
        path bam
        path bai
        path ref
        path fai
        path bed
        path "clair_output"
    output:
        path "clair_output"
        path "clair_output/log/parallel_7_merge_vcf.log"
        path "clair_output/log/7_merge_vcf.log"
        path "clair_output/tmp/merge_output/merge_*.vcf"
        path "clair_output/merge_output.vcf.gz"
        path "clair_output/merge_output.vcf.gz.tbi"
    shell:
        '''
        echo $''
        echo "[INFO] 7/7 Merge pileup VCF and full-alignment VCF"
        . env.sh !{bam} !{bed} !{ref} "clair_output"
        readarray -t CHR < "${OUTPUT_FOLDER}/tmp/CONTIGS"
        time \
            parallel \
                --retries 4 \
                --joblog ${LOG_PATH}/parallel_7_merge_vcf.log \
                --env OUTPUT_FOLDER \
                --env REF \
                -j !{params.threads} \
                "!{params.pypy} $(which clair3.py) MergeVcf \
                    --pileup_vcf_fn ${OUTPUT_FOLDER}/pileup.vcf.gz \
                    --bed_fn_prefix ${OUTPUT_FOLDER}/tmp/full_alignment_output/candidate_bed \
                    --full_alignment_vcf_fn ${OUTPUT_FOLDER}/full_alignment.vcf.gz \
                    --output_fn ${OUTPUT_FOLDER}/tmp/merge_output/merge_{1}.vcf \
                    --platform ont \
                    --print_ref_calls False \
                    --gvcf !{params.GVCF} \
                    --haploid_precise False \
                    --haploid_sensitive False \
                    --gvcf_fn ${OUTPUT_FOLDER}/tmp/merge_output/merge_{1}.gvcf \
                    --non_var_gvcf_fn ${OUTPUT_FOLDER}/tmp/gvcf_tmp_output/non_var.gvcf \
                    --ref_fn ${REF} \
                    --ctgName {1}" \
                ::: ${CHR[@]} |& tee ${LOG_PATH}/7_merge_vcf.log

        !{params.pypy} $(which clair3.py) SortVcf \
            --input_dir ${OUTPUT_FOLDER}/tmp/merge_output \
            --vcf_fn_prefix "merge" \
            --output_fn ${OUTPUT_FOLDER}/merge_output.vcf \
            --sampleName !{params.sample_name}\
            --ref_fn ${REF}

        if [ "$( gzip -fdc ${OUTPUT_FOLDER}/merge_output.vcf.gz | grep -v '#' | wc -l )" -eq 0 ]; \
            then echo "[INFO] Exit in variant merging"; \
            exit 0; \
        fi

        if [ !{params.GVCF} == True ]; \
            then cat ${OUTPUT_FOLDER}/tmp/merge_output/merge_*.gvcf | \
                !{params.pypy} $(which clair3.py) SortVcf --output_fn ${OUTPUT_FOLDER}/merge_output.gvcf; \
        fi

        echo $''
        echo "[INFO] Finish calling, output file: ${OUTPUT_FOLDER}/merge_output.vcf.gz"
        '''
}

process hap {
    label "happy"
    input:
        path clair_output
        path ref
        path fai
        path bed
        path vcf
        path vcf_bed
    output:
        path "hap_output"
    """
    mkdir hap_output
    ls -l $vcf
    ls -l ${clair_output}/full_alignment.vcf.gz
    /opt/hap.py/bin/hap.py \
        "${vcf}" \
        "${clair_output}/full_alignment.vcf.gz" \
        -f ${vcf_bed} \
        -r ${ref} \
        -o "hap_output/happy" \
        --engine=vcfeval \
        --threads=${params.threads} \
        --pass-only
    """
}


// See https://github.com/nextflow-io/nextflow/issues/1636
// This is the only way to publish files from a workflow whilst
// decoupling the publish from the process steps.
process output {
    // publish inputs to output directory
    label "clair3"
    publishDir "${params.out_dir}", mode: 'copy', pattern: "*"
    input:
        file fname
    output:
        file fname
    """
    echo "Writing output files"
    """
}


// workflow module
workflow pipeline {
    take:
        bam
        bai
        bed
        ref
        fai
    main:
        clair_output = run_clair3_stage_1(bam, bai, ref, fai, bed)
        clair_output = run_clair3_stage_2(bam, bai, ref, fai, bed, clair_output[0])
        clair_output = run_clair3_stage_3(bam, bai, ref, fai, bed, clair_output[0])
        clair_output = run_clair3_stage_4(bam, bai, ref, fai, bed, clair_output[0])
        clair_output = run_clair3_stage_5(bam, bai, ref, fai, bed, clair_output[0])
        clair_output = run_clair3_stage_6(bam, bai, ref, fai, bed, clair_output[0])
        clair_output = run_clair3_stage_7(bam, bai, ref, fai, bed, clair_output[0])
    emit:
        clair_output[0]
}


workflow pipeline_with_hap {
    take:
        bam
        bai
        bed
        ref
        fai
        vcf
        vcf_bed
    main:
        clair_output = pipeline(bam, bai, bed, ref, fai)
        hap_output = hap(clair_output, ref, fai, bed, vcf, vcf_bed)
    emit:
        clair_output
        hap_output

}

process download_demo_files {
    label "clair3"
    output:
        path "demo_data"
    """
    mkdir demo_data
    wget -O demo_data.tar.gz  https://ont-exd-int-s3-euwst1-epi2me-labs.s3.amazonaws.com/clair_tutorial/demo_data.tar.gz
    tar -xzvf demo_data.tar.gz -C demo_data
    """
}

workflow download_demo {
    main:
        output_files = download_demo_files()
    emit:
        output_files
}

// entrypoint workflow
workflow {

if (params.help) {
        helpMessage()
        exit 1
    } else if (params.download) {
        demo_files = download_demo()
        output(demo_files)
    } else {
        if (!params.bam) {
            helpMessage()
            println("")
            println("`--bam` is required")
            exit 1
        }
        bai = new File(params.bam + ".bai")
        if (!bai.exists()){
            helpMessage()
            println("")
            println("Bam file missing index: ($bai)")
            exit 1
        }
        if (!params.bed) {
            helpMessage()
            println("")
            println("`--bed` is required")
            exit 1
        }
        if (!params.ref) {
            helpMessage()
            println("")
            println("`--ref` is required")
            exit 1
        }
        fai = new File(params.ref + ".fai")
        if (!fai.exists()) {
            helpMessage()
            println("")
            println("`ref file missing index: (${fai})")
            exit 1
        }
        ref = channel.fromPath(params.ref, checkIfExists:true)
        fai = channel.fromPath(fai, checkIfExists:true)
        bam = channel.fromPath(params.bam, checkIfExists:true)
        bai = channel.fromPath(bai, checkIfExists:true)
        bed = channel.fromPath(params.bed, checkIfExists:true)
        if (params.vcf){
            vcf = channel.fromPath(params.vcf, checkIfExists:true)
            tbi = new File(params.vcf + ".tbi")
            if (!tbi.exists()){
                helpMessage()
                println("")
                println("`vcf file missing index: (${tbi})")
                exit 1
            }
            vcf_bed = channel.fromPath(params.vcf_bed, checkIfExists:true)
            pipeline_output = pipeline_with_hap(bam, bai, bed, ref, fai, vcf, vcf_bed)
            output(pipeline_output[0].concat(pipeline_output[1]))
        } else {
            pipeline_output = pipeline(bam, bai, bed, ref, fai)
            output(pipeline_output)
        }
    }
}
