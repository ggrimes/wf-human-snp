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
params.threads = 4
params.download = ""
params.vcf = ""
params.vcf_bed = ""
params.python="python3"
params.pypy="pypy3"
params.sample_name = "SAMPLE"
params.vcf_fn = "EMPTY"
params.ctg_name = "EMPTY"
params.include_all_ctgs = "False"
params.CUDA_VISIBLE_DEVICES = ""
params.ref_pct_full = 0.1
params.var_pct_full = 0.3
params.GVCF = "False"
params.snp_min_af = 0.08
params.indel_min_af = 0.15


def helpMessage(){
    log.info """
Workflow template'

Usage:
    nextflow run epi2melabs/wf-clair [options]

Script Options:
    --bam       FILE    Path to bam file
    --bed       FILE    Path to bed file
    --ref       FILE    Path to ref file
    --model     DIR     Path to model
    --output    DIR     Output path
"""
}

process run_clair3_stage_1a {
    label "clair3"
    cpus params.threads
    input:
        path bam
        path bai
        path ref
        path fai
        path bed
    output:
        path "clair_output/tmp/CONTIGS", emit: contigs_file
        path "clair_output/tmp/CHUNK_LIST", emit: chunks_file
        path "clair_output/tmp/split_beds", emit: split_beds
        path "clair_output/tmp/gvcf_tmp_output", emit: gvcf_tmp_output
    shell:
        '''
        mkdir -p clair_output/log
        touch clair_output/log/parallel_1_call_var_bam_pileup.log
        . env.sh !{bam} !{bed} !{ref} "clair_output" !{params.threads}
        !{params.python} $(which clair3.py) CheckEnvs \
            --bam_fn !{bam} \
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
        '''
}

process generate_chunks {
    // generate: chr20 1 13
    label "clair3"
    cpus params.threads
    input:
        path "clair_output/tmp/CHUNK_LIST"
    output:
        stdout
    shell:
        '''
        #!/usr/bin/env python
        with open("clair_output/tmp/CHUNK_LIST") as f:
            for chunk in f:
                print(chunk.strip())
        '''

}

process generate_contigs {
    // generate: chr20
    label "clair3"
    cpus params.threads
    input:
        path "clair_output/tmp/CONTIGS"
    output:
        stdout
    shell:
        '''
        #!/usr/bin/env python
        with open("clair_output/tmp/CONTIGS") as f:
            for contig in f:
                print(contig.strip())
        '''

}

process generate_regions {
    input:
        val region
    output:
        tuple val(contig), val(chunk_id), val(chunk_num), emit: regions
    exec: (_, contig, _, chunk_id, _, chunk_num) = (region =~ /(.+)(\s)(.+)(\s)(.+)/)[0]    
}


process run_clair3_stage_1b {
    label "clair3"
    cpus params.threads
    input:
        each region
        path bam
        path bai
        path ref
        path fai
        path bed
        path "clair_output/tmp/split_beds"
        path "clair_output/tmp/gvcf_tmp_output"
        path model
    output:
        path "clair_output/log/1_call_var_bam_pileup.log"
        path "clair_output/tmp/pileup_output/*", optional: true, emit: "tmp_pileup_output"
    script:
        contig=region[0]
        chunk_id=region[1]
        chunk_num=region[2]
        
    """
        . env.sh $bam $bed $ref "clair_output" ${params.threads}
        export OMP_NUM_THREADS=${params.threads}
        export CUDA_VISIBLE_DEVICES=""
        mkdir -p clair_output/tmp/pileup_output
        ${params.python} \$(which clair3.py) CallVarBam \
            --chkpnt_fn $model/pileup \
            --bam_fn $bam \
            --call_fn clair_output/tmp/pileup_output/pileup_"$chunk_id"_chr20.vcf \
            --ref_fn $ref \
            --ctgName $contig \
            --chunk_id $chunk_id \
            --chunk_num $chunk_num \
            --tensorflow_threads ${params.threads} \
            --pileup \
        |& tee clair_output/log/1_call_var_bam_pileup.log
    """
}
process run_clair3_stage_1c {
    label "clair3"
    cpus params.threads
    input:
        path bam
        path bai
        path ref
        path fai
        path bed
        path "clair_output/tmp/pileup_output/"
        path contigs
    output:
        path "clair_output/pileup.vcf.gz", emit: pileup_gz
        path "clair_output/pileup.vcf.gz.tbi", emit: pileup_gz_tbi
    shell:
        '''
        . env.sh !{bam} !{bed} !{ref} "clair_output" !{params.threads}
        !{params.pypy} $(which clair3.py) SortVcf \
            --input_dir clair_output/tmp/pileup_output/ \
            --vcf_fn_prefix pileup \
            --output_fn clair_output/pileup.vcf \
            --sampleName !{params.sample_name} \
            --ref_fn ${REF} \
            --contigs_fn !{contigs}

        if [ "$( gzip -fdc clair_output/pileup.vcf.gz | grep -v '#' | wc -l )" -eq 0 ]; \
        then echo "[INFO] Exit in pileup variant calling"; exit 0; fi
        '''
}

process run_clair3_stage_2a {
    label "clair3"
    cpus params.threads
    input:
        path bam
        path bai
        path ref
        path fai
        path bed
        path "clair_output/pileup.vcf.gz"
        path "clair_output/pileup.vcf.gz.tbi"
    output:
        path "clair_output/tmp/phase_output/phase_vcf/phase_qual", emit: phase_qual
    shell:
        '''
        . env.sh !{bam} !{bed} !{ref} "clair_output"
        mkdir -p "clair_output/tmp/phase_output/phase_vcf"
        echo $''
        echo "[INFO] 2/7 Select heterozygous SNP variants for Whatshap phasing and haplotagging"
        gzip \
            -fdc clair_output/pileup.vcf.gz | \
            !{params.pypy} $(which clair3.py) SelectQual \
                --phase \
                --output_fn clair_output/tmp/phase_output/phase_vcf
        '''
}

process run_clair3_stage_2b {
    label "clair3"
    cpus params.threads
    input:
        each contig
        path bam
        path bai
        path ref
        path fai
        path bed
        path "clair_output/pileup.vcf.gz"
        path "clair_output/pileup.vcf.gz.tbi"
        path "clair_output/tmp/phase_output/phase_vcf/phase_qual"
    output:
        path "clair_output/log/2_select_hetero_snp.log"
        path "clair_output/tmp/phase_output/phase_vcf", emit: phase_vcf_dir
    shell:
        '''
        . env.sh !{bam} !{bed} !{ref} "clair_output"
        mkdir -p "/tmp/phase_output/phase_vcf"
        !{params.pypy} $(which clair3.py) SelectHetSnp \
            --vcf_fn clair_output/pileup.vcf.gz \
            --split_folder clair_output/tmp/phase_output/phase_vcf \
            --ctgName !{contig} \
        |& tee ${LOG_PATH}/2_select_hetero_snp.log
        '''
}


process run_clair3_stage_3{
    label "clair3"
    cpus params.threads
    input:
        each contig
        path bam
        path bai
        path ref
        path fai
        path bed
        path "clair_output/tmp/phase_output/phase_vcf"
    output:
        path "clair_output/tmp/phase_output/phase_vcf", emit: "phase_vcf_dir"
    shell:
        '''
        echo $''
        echo "[INFO] 3/7 Phase VCF file using Whatshap"
        export OMP_NUM_THREADS=!{params.threads}
        . env.sh !{bam} !{bed} !{ref} "clair_output"
        whatshap phase \
            --output clair_output/tmp/phase_output/phase_vcf/phased_!{contig}.vcf.gz \
            --reference !{ref} \
            --chromosome !{contig} \
            --distrust-genotypes \
            --ignore-read-groups \
            clair_output/tmp/phase_output/phase_vcf/!{contig}.vcf \
            !{bam} \
        |& tee clair_output/tmp/3_phase.log

        tabix \
            -f \
            -p vcf \
            clair_output/tmp/phase_output/phase_vcf/phased_!{contig}.vcf.gz

        '''
}

process run_clair3_stage_4{
    label "clair3"
    cpus params.threads
    input:
        each contig
        path bam
        path bai
        path ref
        path fai
        path bed
        path "clair_output/tmp/phase_output/phase_vcf"
    output:
        path "clair_output/log/4_haplotag.log"
        path "clair_output/tmp/phase_output/phase_bam/*.bam", emit: phased_bam
        path "clair_output/tmp/phase_output/phase_bam/*.bam.bai", emit: phased_bam_index
        path "clair_output/tmp/phase_output/phase_vcf/*.vcf.gz", emit: phased_vcf

    shell:
        '''
        echo $''
        echo "[INFO] 4/7 Haplotag input BAM file using Whatshap"
        . env.sh !{bam} !{bed} !{ref} "clair_output"
        mkdir -p clair_output/tmp/phase_output/phase_bam
        whatshap haplotag \
            --output clair_output/tmp/phase_output/phase_bam/!{contig}.bam  \
            --reference !{ref} \
            --ignore-read-groups \
            --regions !{contig} \
            clair_output/tmp/phase_output/phase_vcf/phased_!{contig}.vcf.gz \
            !{bam} \
        |& tee ${LOG_PATH}/4_haplotag.log

        samtools index \
            -@12 \
            clair_output/tmp/phase_output/phase_bam/!{contig}.bam
        '''
}

process run_clair3_stage_5a{
    label "clair3"
    cpus params.threads
    input:
        path bam
        path bai
        path ref
        path fai
        path bed
        path "clair_output/pileup.vcf.gz"
        path "clair_output/pileup.vcf.gz.tbi"
    output:
        path "clair_output/tmp/full_alignment_output/candidate_bed/qual", emit: full_qual
    shell:
        '''
        # Full alignment calling
        echo $''
        echo "[INFO] 5/7 Select candidates for full-alignment calling"
        . env.sh !{bam} !{bed} !{ref} "clair_output"
        mkdir -p clair_output/tmp/full_alignment_output
        mkdir -p clair_output/tmp/full_alignment_output/candidate_bed
        gzip -fdc clair_output/pileup.vcf.gz | \
        !{params.pypy} $(which clair3.py) SelectQual \
                --output_fn clair_output/tmp/full_alignment_output/candidate_bed \
                --var_pct_full !{params.var_pct_full} \
                --ref_pct_full !{params.ref_pct_full} \
                --platform ont 
        '''
}

process run_clair3_stage_5b{
    label "clair3"
    cpus params.threads
    input:
        each contig
        path bam
        path bai
        path ref
        path fai
        path bed
        path "clair_output/pileup.vcf.gz"
        path "clair_output/pileup.vcf.gz.tbi"
        path "clair_output/tmp/full_alignment_output/candidate_bed/qual"
    output:
        path "clair_output", emit: clair
        path "clair_output/log/5_select_candidate.log"
        path "clair_output/tmp/full_alignment_output/candidate_bed/FULL_ALN_FILE_*", emit: aln_bed
    shell:
        '''
        . env.sh !{bam} !{bed} !{ref} "clair_output"
        
        !{params.pypy} $(which clair3.py) SelectCandidates \
            --pileup_vcf_fn clair_output/pileup.vcf.gz \
            --split_folder clair_output/tmp/full_alignment_output/candidate_bed \
            --ref_fn !{ref} \
            --var_pct_full !{params.var_pct_full} \
            --ref_pct_full !{params.ref_pct_full} \
            --platform ont\
            --ctgName !{contig} \
        |& tee clair_output/log/5_select_candidate.log
        '''
}

process run_clair3_stage_6a {
    label "clair3"
    cpus params.threads
    input:
        path bam
        path bai
        path ref
        path fai
        path bed
        path "clair_output"
        
    output:
        path "clair_output", emit: clair
        path "clair_output/tmp/full_alignment_output/candidate_bed/chr*", emit: full_alignment
        stdout
    shell:
        '''
        echo $''
        . env.sh !{bam} !{bed} !{ref} "clair_output"
        cat clair_output/tmp/full_alignment_output/candidate_bed/FULL_ALN_FILE_* > clair_output/tmp/full_alignment_output/candidate_bed/FULL_ALN_FILES
        cat clair_output/tmp/full_alignment_output/candidate_bed/FULL_ALN_FILES
        '''
}


process run_clair3_stage_6b {
    label "clair3"
    cpus params.threads
    input:
        path bam
        path bai
        path ref
        path fai
        path bed
        path model
        path phased_bam
        path phased_vcf
        each alignment_path
        path phased_bam_index
        path "alignment_sections/*"
    output:
        path "clair_output/full_alignment_*.vcf", emit: full_alignment
        path "clair_output/log/6_call_var_bam_full_alignment.log"
    script:
        File f = new File(alignment_path)
        filename = f.getName()
        """
        echo "[INFO] 6/7 Call low-quality variants using full-alignment model"
        . env.sh $bam $bed $ref "clair_output" ${params.threads}
        export OMP_NUM_THREADS=${params.threads}
        ${params.python} \$(which clair3.py) CallVarBam \
            --chkpnt_fn $model/full_alignment \
            --bam_fn $phased_bam \
            --call_fn clair_output/full_alignment_"$filename".vcf \
            --ref_fn $ref \
            --full_aln_regions alignment_sections/$filename \
            --ctgName chr20 \
            --add_indel_length \
            --phasing_info_in_bam \
            --platform ont \
            --tensorflow_threads ${params.threads} \
        |& tee clair_output/log/6_call_var_bam_full_alignment.log
        """
}

process run_clair3_stage_6c {
    label "clair3"
    cpus params.threads
    input:
        path bam
        path bai
        path ref
        path fai
        path bed
        path "full_alignment/*"
        path contigs
    output:
        path "clair_output", emit: clair
        path "full_alignment.vcf.gz", emit: merged_alignment
        path "full_alignment.vcf.gz.tbi"
    shell:
        '''
        . env.sh !{bam} !{bed} !{ref} "clair_output"

        !{params.pypy} $(which clair3.py) SortVcf \
            --input_dir full_alignment \
            --vcf_fn_prefix "full_alignment" \
            --output_fn full_alignment.vcf \
            --sampleName !{params.sample_name} \
            --ref_fn !{ref} \
            --contigs_fn !{contigs}

        if [ "$( gzip -fdc full_alignment.vcf.gz | grep -v '#' | wc -l )" -eq 0 ];
            then echo "[INFO] Exit in full-alignment variant calling";
            exit 0;
        fi

        if [ !{params.GVCF} == True ];
            then cat clair_output/tmp/gvcf_tmp_output/*.tmp.g.vcf | \
                ${PYPY} ${CLAIR3} SortVcf
                    --output_fn clair_output/tmp/gvcf_tmp_output/non_var.gvcf;
        fi
        '''
}

process bed_files_to_folder{
    label "clair3"
    cpus params.threads
    input:
        
        path full_qual
        path aln_bed
        path full_alignment
    output:
        path "candidate_bed", emit: bed_files
        
    """
     mkdir candidate_bed
     mv $full_qual candidate_bed
     mv $aln_bed candidate_bed
     mv $full_alignment candidate_bed
    """
}


process run_clair3_stage_7a{
    label "clair3"
    cpus params.threads
    input:
        each contig
        path bam
        path bai
        path ref
        path fai
        path bed
        path pile_up_vcf
        path full_aln_vcf
        path bed_files
    output:
        path "clair_output", emit: clair
        path "clair_output/log/7_merge_vcf.log"
        path "clair_output/tmp/merge_output/merge_*.vcf", emit: merge_output
    shell:
        '''
        mkdir -p clair_output/tmp/merge_output
        echo $''
        echo "[INFO] 7/7 Merge pileup VCF and full-alignment VCF"
        . env.sh !{bam} !{bed} !{ref} "clair_output"
        !{params.pypy} $(which clair3.py) MergeVcf \
            --pileup_vcf_fn !{pile_up_vcf} \
            --bed_fn_prefix !{bed_files} \
            --full_alignment_vcf_fn !{full_aln_vcf} \
            --output_fn clair_output/tmp/merge_output/merge_!{contig}.vcf \
            --platform ont \
            --print_ref_calls False \
            --gvcf !{params.GVCF} \
            --haploid_precise False \
            --haploid_sensitive False \
            --gvcf_fn clair_output/tmp/merge_output/merge_!{contig}.gvcf \
            --non_var_gvcf_fn clair_output/tmp/gvcf_tmp_output/non_var.gvcf \
            --ref_fn ${REF} \
            --ctgName !{contig} \
        |& tee ${LOG_PATH}/7_merge_vcf.log
        '''
}

process run_clair3_stage_7b{
    label "clair3"
    cpus params.threads
    input:
        path bam
        path bai
        path ref
        path fai
        path bed
        path "merge_output/*"
        path contigs
        path pileup
    output:
        path "merge_output.vcf.gz", emit: merge_output
        path "merge_output.vcf.gz.tbi"
    shell:
        '''
        . env.sh !{bam} !{bed} !{ref} "clair_output"
        
        !{params.pypy} $(which clair3.py) SortVcf \
            --input_dir merge_output \
            --output_fn merge_output.vcf \
            --sampleName !{params.sample_name} \
            --ref_fn !{ref} \
            --contigs_fn !{contigs}

        if [ "$( gzip -fdc merge_output.vcf.gz | grep -v '#' | wc -l )" -eq 0 ]; \
            then echo "[INFO] Exit in variant merging"; \
            exit 0; \
        fi

        if [ !{params.GVCF} == True ]; \
            then cat clair_output/tmp/merge_output/merge_*.gvcf | \
                !{params.pypy} $(which clair3.py) SortVcf --output_fn clair_output/merge_output.gvcf; \
        fi

        echo $''
        echo "[INFO] Finish calling, output file: clair_putput/merge_output.vcf.gz"
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
        model
    main:
        clair_output = run_clair3_stage_1a(bam, bai, ref, fai, bed)
        regions =  generate_chunks(run_clair3_stage_1a.out.chunks_file).splitText().map(it -> it.trim())
        contigs = generate_contigs(run_clair3_stage_1a.out.contigs_file).splitText().map(it -> it.trim())
        regions2 = generate_regions(regions)
        clair_output =  run_clair3_stage_1b(
            regions2.regions, 
            bam, bai, ref, fai, bed, 
            run_clair3_stage_1a.out.split_beds, 
            run_clair3_stage_1a.out.gvcf_tmp_output,
            model
            )
        clair_output = run_clair3_stage_1c(
            bam, bai, ref, fai, bed, 
            run_clair3_stage_1b.out.tmp_pileup_output.collect(),
            run_clair3_stage_1a.out.contigs_file
            )
        clair_output = run_clair3_stage_2a(
            bam, bai, ref, fai, bed,
            run_clair3_stage_1c.out.pileup_gz,
            run_clair3_stage_1c.out.pileup_gz_tbi
            )
        clair_output = run_clair3_stage_2b(
            contigs, bam, bai, ref, fai, bed,
            run_clair3_stage_1c.out.pileup_gz,
            run_clair3_stage_1c.out.pileup_gz_tbi,
            run_clair3_stage_2a.out.phase_qual 
            )
        clair_output = run_clair3_stage_3(
            contigs, bam, bai, ref, fai, bed,
            run_clair3_stage_2b.out.phase_vcf_dir 
            )
        clair_output = run_clair3_stage_4(
            contigs, bam, bai, ref, fai, bed,
            run_clair3_stage_3.out.phase_vcf_dir,
            )
        clair_output = run_clair3_stage_5a(
            bam, bai, ref, fai, bed, 
            run_clair3_stage_1c.out.pileup_gz,
            run_clair3_stage_1c.out.pileup_gz_tbi,
            )
        clair_output = run_clair3_stage_5b(
            contigs, bam, bai, ref, fai, bed, 
            run_clair3_stage_1c.out.pileup_gz,
            run_clair3_stage_1c.out.pileup_gz_tbi,
            run_clair3_stage_5a.out.full_qual)
        clair_output_test = run_clair3_stage_6a(
            bam, bai, ref, fai, bed,
            run_clair3_stage_5b.out.clair)
        full_aln_paths = clair_output[2].splitText().map(it -> it.trim())
        clair_output_6b = run_clair3_stage_6b(
            bam, bai, ref, fai, bed,
            model,
            run_clair3_stage_4.out.phased_bam,
            run_clair3_stage_4.out.phased_vcf,
            full_aln_paths,
            run_clair3_stage_4.out.phased_bam_index,
            run_clair3_stage_6a.out.full_alignment)
        clair_output = run_clair3_stage_6c(
            bam, bai, ref, fai, bed,
            run_clair3_stage_6b.out.full_alignment.collect(),
            run_clair3_stage_1a.out.contigs_file)
        candidate_bed = bed_files_to_folder(
            run_clair3_stage_5a.out.full_qual,
            run_clair3_stage_5b.out.aln_bed,
            run_clair3_stage_6a.out.full_alignment
            )
        clair_output = run_clair3_stage_7a(contigs, bam, bai, ref, fai, bed,
            run_clair3_stage_1c.out.pileup_gz,  run_clair3_stage_6c.out.merged_alignment,
            candidate_bed.bed_files)
        clair_final = run_clair3_stage_7b(bam, bai, ref, fai, bed,
            run_clair3_stage_7a.out.merge_output.collect(),
            run_clair3_stage_1a.out.contigs_file,
            run_clair3_stage_1c.out.pileup_gz
            )
    emit:
        clair_final[0]
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
        ref = file(params.ref, type: "file")
        fai = channel.fromPath(fai, checkIfExists:true)
        bam = file(params.bam, type: "file")
        bai = channel.fromPath(bai, checkIfExists:true)
        bed = file(params.bed, type: "file")
        model = file(params.model, type: "dir")
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
            pipeline_output = pipeline(bam, bai, bed, ref, fai, model)
            output(pipeline_output)
        }
    }
}
