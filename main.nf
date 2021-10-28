#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

params.help = ""
params.bam = ""
params.ref = ""
params.bed = ""
params.threads = 4
params.download = ""
params.vcf = ""
params.vcf_bed = ""
params.sample_name = "SAMPLE"
params.vcf_fn = "EMPTY"
params.ctg_name = "EMPTY"
params.include_all_ctgs = "False"
params.ref_pct_full = 0.1
params.var_pct_full = 0.7
params.GVCF = "False"
params.snp_min_af = 0.0
params.indel_min_af = 0.0


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
    --out_dir   DIR     Output path
"""
}

process run_clair3_stage_1a {
    label "clair3"
    cpus params.threads
    input:
        tuple path(bam), path(bai)
        tuple path(ref), path(fai)
        path bed
    output:
        path "clair_output/tmp/CONTIGS", emit: contigs_file
        path "clair_output/tmp/CHUNK_LIST", emit: chunks_file
        path "clair_output/tmp/gvcf_tmp_output", emit: gvcf_tmp_output
    shell:
        '''
        mkdir -p clair_output
        python $(which clair3.py) CheckEnvs \
            --bam_fn !{bam} \
            --bed_fn !{bed} \
            --output_fn_prefix clair_output \
            --ref_fn !{ref} \
            --vcf_fn !{params.vcf_fn} \
            --ctg_name !{params.ctg_name} \
            --chunk_num 0 \
            --chunk_size 5000000 \
            --include_all_ctgs !{params.include_all_ctgs} \
            --threads !{params.threads}  \
            --qual 2 \
            --sampleName !{params.sample_name} \
            --var_pct_full !{params.var_pct_full} \
            --ref_pct_full !{params.ref_pct_full} \
            --snp_min_af !{params.snp_min_af} \
            --indel_min_af !{params.indel_min_af}
        '''
}


// TODO: this doesn't need to be a process
process generate_chunks {
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


// TODO: this doesn't need to be a process
process generate_contigs {
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


// TODO: this doesn't need to be a process
process generate_regions {
    input:
        val region
    output:
        tuple val(contig), val(chunk_id), val(chunk_num), emit: regions
    exec: (_, contig, _, chunk_id, _, chunk_num) = (region =~ /(.+)(\s)(.+)(\s)(.+)/)[0]    
}


process run_clair3_stage_1b {
    // Calls variants per region ("chunk") using pileup network.
    label "clair3"
    cpus params.threads
    input:
        each region
        tuple path(bam), path(bai)
        tuple path(ref), path(fai)
        path model
    output:
        path "clair_output/tmp/pileup_output/*", optional: true, emit: "tmp_pileup_output"
    script:
        contig=region[0]
        chunk_id=region[1]
        chunk_num=region[2]
        """
        mkdir -p clair_output/tmp/pileup_output
        python \$(which clair3.py) CallVarBam \
            --chkpnt_fn $model/pileup \
            --bam_fn $bam \
            --call_fn clair_output/tmp/pileup_output/pileup_${contig}_${chunk_id}.vcf \
            --ref_fn $ref \
            --ctgName $contig \
            --chunk_id $chunk_id \
            --chunk_num $chunk_num \
            --platform ont \
            --call_snp_only False \
            --fast_mode False \
            --snp_min_af ${params.snp_min_af} \
            --indel_min_af ${params.indel_min_af} \
            --pileup
    """
}


process run_clair3_stage_1c {
    // Aggregates and sorts all variants (across all chunks of all contigs)
    // from pileup network.
    label "clair3"
    cpus params.threads
    input:
        tuple path(ref), path(fai)
        path "clair_output/tmp/pileup_output/"
        path contigs
    output:
        tuple path("clair_output/pileup.vcf.gz"), path("clair_output/pileup.vcf.gz.tbi"), emit: pileup_vcf
    shell:
        '''
        pypy $(which clair3.py) SortVcf \
            --input_dir clair_output/tmp/pileup_output/ \
            --vcf_fn_prefix "pileup" \
            --output_fn clair_output/pileup.vcf \
            --sampleName !{params.sample_name} \
            --ref_fn !{ref} \
            --contigs_fn !{contigs}

        if [ "$( gzip -fdc clair_output/pileup.vcf.gz | grep -v '#' | wc -l )" -eq 0 ]; \
        then echo "[INFO] Exit in pileup variant calling"; exit 1; fi
        '''
}

// TODO: combine this with the previous to avoid a step
process run_clair3_stage_2a {
    // Filters a VCF (by quality?) for the purposes of input to whatshap for phasing.
    label "clair3"
    cpus params.threads
    input:
        tuple path("pileup.vcf.gz"), path("pileup.vcf.gz.tbi")
    output:
        path "phase_vcf/phase_qual", emit: phase_qual
    shell:
        '''
        mkdir -p phase_vcf
        echo "[INFO] 2/7 Select heterozygous SNP variants for Whatshap phasing and haplotagging"
        gzip -fdc pileup.vcf.gz | \
            pypy $(which clair3.py) SelectQual --phase --output_fn phase_vcf
        '''
}


process run_clair3_stage_2b {
    // Filters a VCF by contig, selecting only het SNPs.
    label "clair3"
    cpus params.threads
    input:
        each contig
        tuple path("pileup.vcf.gz"), path("pileup.vcf.gz.tbi")
        path "phase_vcf/phase_qual"
    output:
        // TODO: this looks nasty: we're outputting a superset of the input.
        //       Don't we want to simply the outputs?
        path "phase_vcf", emit: phase_vcf_dir
    shell:
        '''
        pypy $(which clair3.py) SelectHetSnp \
            --vcf_fn pileup.vcf.gz \
            --split_folder phase_vcf \
            --ctgName !{contig}
        '''
}


process run_clair3_stage_3 {
    // Phases a VCF using whatshap.
    label "clair3"
    cpus params.threads
    input:
        each contig
        tuple path(bam), path(bai)
        tuple path(ref), path(fai)
        path "clair_output/tmp/phase_output/phase_vcf"
    output:
        // TODO: should we be just returning the phased vcf?
        path "clair_output/tmp/phase_output/phase_vcf", emit: "phase_vcf_dir"
    shell:
        '''
        echo "[INFO] 3/7 Phase VCF file using Whatshap"
        whatshap phase \
            --output clair_output/tmp/phase_output/phase_vcf/phased_!{contig}.vcf.gz \
            --reference !{ref} \
            --chromosome !{contig} \
            --distrust-genotypes \
            --ignore-read-groups \
            clair_output/tmp/phase_output/phase_vcf/!{contig}.vcf \
            !{bam}

        tabix -f -p vcf \
            clair_output/tmp/phase_output/phase_vcf/phased_!{contig}.vcf.gz
        '''
}


process run_clair3_stage_4 {
    // Tags reads in a BAM with haplotype information using a VCF.
    label "clair3"
    cpus params.threads
    input:
        // TODO: THIS IS INCORRECT
        //       the third input here is a channel containing phased vcfs, one item
        //       per contig. We shouldn't be using "each contig" here, the op is
        //       a no-op in the case of a single contig.
        each contig
        tuple path(bam), path(bai)
        tuple path(ref), path(fai)
        path "clair_output/tmp/phase_output/phase_vcf"
    output:
        tuple path("clair_output/tmp/phase_output/phase_bam/*.bam"), path("clair_output/tmp/phase_output/phase_bam/*.bam.bai"), emit: phased_bam
        path "clair_output/tmp/phase_output/phase_vcf/*.vcf.gz", emit: phased_vcf
    shell:
        '''
        echo "[INFO] 4/7 Haplotag input BAM file using Whatshap"
        mkdir -p clair_output/tmp/phase_output/phase_bam
        whatshap haplotag \
            --output clair_output/tmp/phase_output/phase_bam/!{contig}.bam  \
            --reference !{ref} \
            --ignore-read-groups \
            --regions !{contig} \
            clair_output/tmp/phase_output/phase_vcf/phased_!{contig}.vcf.gz \
            !{bam}

        samtools index \
            -@!{task.cpus} \
            clair_output/tmp/phase_output/phase_bam/!{contig}.bam
        '''
}


process run_clair3_stage_5a {
    // Determines quality filter for selecting candidate variants for
    // second stage "full alignment" calling
    label "clair3"
    cpus params.threads
    input:
        tuple path("pileup.vcf.gz"), path("pileup.vcf.gz.tbi")
    output:
        path "qual", emit: full_qual
    shell:
        '''
        # Full alignment calling
        echo "[INFO] 5/7 Select candidates for full-alignment calling"
        mkdir -p clair_output/tmp/full_alignment_output
        mkdir -p clair_output/tmp/full_alignment_output/candidate_bed
        gzip -fdc pileup.vcf.gz | \
        pypy $(which clair3.py) SelectQual \
                --output_fn . \
                --var_pct_full !{params.var_pct_full} \
                --ref_pct_full !{params.ref_pct_full} \
                --platform ont 
        '''
}


process run_clair3_stage_5b {
    // Create BED files for candidate variants for "full alignment" network.
    label "clair3"
    cpus params.threads
    input:
        each contig
        tuple path(ref), path(fai)
        tuple path("pileup.vcf.gz"), path("pileup.vcf.gz.tbi")
        path "candidate_bed/qual"
    output:
        path "candidate_bed/${contig}.*", emit: candidate_bed
    shell:
        // this creates BED files as candidate_bed/<ctg>.0_14 with candidates
        // along with a file the FULL_ALN_FILE_<ctg> listing all of the BED files.
        // All we really want are the BEDs, the file of filenames is used for the
        // purposes of parallel in the original workflow.
        // https://github.com/HKU-BAL/Clair3/blob/329d09b39c12b6d8d9097aeb1fe9ec740b9334f6/scripts/clair3.sh#L218
        '''
        mkdir -p clair_output/tmp/full_alignment_output/candidate_bed 
        pypy $(which clair3.py) SelectCandidates \
            --pileup_vcf_fn pileup.vcf.gz \
            --split_folder candidate_bed \
            --ref_fn !{ref} \
            --var_pct_full !{params.var_pct_full} \
            --ref_pct_full !{params.ref_pct_full} \
            --platform ont \
            --ctgName !{contig}
        '''
}



process run_clair3_stage_6b {
    label "clair3"
    cpus params.threads
    input:
        tuple path(ref), path(fai)
        path model
        tuple path(phased_bam), path(phased_bam_index)
        each path(candidate_bed)
    output:
        path "clair_output/full_alignment_*.vcf", emit: full_alignment
    script:
        // TODO: clean this up. 
        f = file(candidate_bed)
        filename = f.getName()
        contig = "chr20"
        """
        mkdir clair_output
        echo "[INFO] 6/7 Call low-quality variants using full-alignment model"
        python \$(which clair3.py) CallVarBam \
            --chkpnt_fn $model/full_alignment \
            --bam_fn $phased_bam \
            --call_fn clair_output/full_alignment_${filename}.vcf \
            --ref_fn $ref \
            --full_aln_regions $candidate_bed \
            --ctgName $contig \
            --phasing_info_in_bam \
            --snp_min_af ${params.snp_min_af} \
            --indel_min_af ${params.indel_min_af} \
            --platform ont \
            --add_indel_length \
            --gvcf False \
        """
}


process run_clair3_stage_6c {
    label "clair3"
    cpus params.threads
    input:
        tuple path(ref), path(fai)
        path "full_alignment/*"
        path contigs
    output:
        tuple path("full_alignment.vcf.gz"), path("full_alignment.vcf.gz.tbi"), emit: full_aln_vcf
    shell:
        '''
        pypy $(which clair3.py) SortVcf \
            --input_dir full_alignment \
            --output_fn full_alignment.vcf \
            --sampleName !{params.sample_name} \
            --ref_fn !{ref} \
            --contigs_fn !{contigs}

        if [ "$( gzip -fdc full_alignment.vcf.gz | grep -v '#' | wc -l )" -eq 0 ];
            then echo "[INFO] Exit in full-alignment variant calling";
            exit 0;
        fi

        # TODO: WHAT IS THIS?
        if [ !{params.GVCF} == True ];
            then cat clair_output/tmp/gvcf_tmp_output/*.tmp.g.vcf | \
                ${PYPY} ${CLAIR3} SortVcf
                    --output_fn clair_output/tmp/gvcf_tmp_output/non_var.gvcf;
        fi
        '''
}


//TODO: this is pointless
process bed_files_to_folder{
    label "clair3"
    cpus params.threads
    input:
        path full_qual
        path aln_bed
        path full_alignment
    output:
        path "candidate_bed", emit: bed_files
    script:
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
        tuple path(ref), path(fai)
        tuple path(pile_up_vcf), path(pile_up_vcf_tbi)
        tuple path(full_aln_vcf), path(full_aln_vcf_tbi)
        path "candidate_beds/*"
    output:
        path "clair_output", emit: clair
        path "clair_output/tmp/merge_output/merge_*.vcf", emit: merge_output
    shell:
        '''
        mkdir -p clair_output/tmp/merge_output
        echo "[INFO] 7/7 Merge pileup VCF and full-alignment VCF"
        pypy $(which clair3.py) MergeVcf \
            --pileup_vcf_fn !{pile_up_vcf} \
            --bed_fn_prefix candidate_beds \
            --full_alignment_vcf_fn !{full_aln_vcf} \
            --output_fn clair_output/tmp/merge_output/merge_!{contig}.vcf \
            --platform ont \
            --print_ref_calls False \
            --gvcf !{params.GVCF} \
            --haploid_precise False \
            --haploid_sensitive False \
            --gvcf_fn clair_output/tmp/merge_output/merge_!{contig}.gvcf \
            --non_var_gvcf_fn clair_output/tmp/gvcf_tmp_output/non_var.gvcf \
            --ref_fn !{ref} \
            --ctgName !{contig}
        '''
}


process run_clair3_stage_7b{
    label "clair3"
    cpus params.threads
    input:
        tuple path(ref), path(fai)
        path "merge_output/*"
        path contigs
    output:
        tuple path("merge_output.vcf.gz"), path("merge_output.vcf.gz.tbi"), emit: final_vcf
    shell:
        '''
        pypy $(which clair3.py) SortVcf \
            --input_dir merge_output \
            --output_fn merge_output.vcf \
            --sampleName !{params.sample_name} \
            --ref_fn !{ref} \
            --contigs_fn !{contigs}

        if [ "$( gzip -fdc merge_output.vcf.gz | grep -v '#' | wc -l )" -eq 0 ]; \
            then echo "[INFO] Exit in variant merging"; \
            exit 0; \
        fi

        # TODO: this looks broken
        if [ !{params.GVCF} == True ]; \
            then cat clair_output/tmp/merge_output/merge_*.gvcf | \
                pypy $(which clair3.py) SortVcf --output_fn clair_output/merge_output.gvcf; \
        fi

        echo "[INFO] Finish calling, output file: merge_output.vcf.gz"
        '''
}


// TODO: test this
process hap {
    label "happy"
    input:
        path "clair.vcf.gz"
        tuple path("ref.fasta"), path("ref.fasta.fai")
        path "truth.vcf"
        path "truth.bed"
    output:
        path "happy"
    """
    mkdir hap_output
    ls -l $vcf
    ls -l ${clair_output}/full_alignment.vcf.gz
    /opt/hap.py/bin/hap.py \
        truth.vcf \
        clair.vcf.gz \
        -f truth.bed \
        -r ref.fasta \
        -o happy \
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
        bed
        ref
        model
    main:
        // TODO: what are regions, contigs, regions2?
        run_clair3_stage_1a(bam, ref, bed)
        regions =  generate_chunks(run_clair3_stage_1a.out.chunks_file).splitText().map(it -> it.trim())
        contigs = generate_contigs(run_clair3_stage_1a.out.contigs_file).splitText().map(it -> it.trim())
        regions2 = generate_regions(regions)
        regions2.regions.view()

        run_clair3_stage_1b(regions2.regions, bam, ref, model)
        run_clair3_stage_1c(ref, run_clair3_stage_1b.out.tmp_pileup_output.collect(), run_clair3_stage_1a.out.contigs_file)

        // TODO: these should be combined
        run_clair3_stage_2a(run_clair3_stage_1c.out.pileup_vcf)
        run_clair3_stage_2b(
            contigs,
            run_clair3_stage_1c.out.pileup_vcf,
            run_clair3_stage_2a.out.phase_qual)

        run_clair3_stage_3(
            contigs, bam, ref,
            run_clair3_stage_2b.out.phase_vcf_dir)

        run_clair3_stage_4(
            contigs, bam, ref,
            run_clair3_stage_3.out.phase_vcf_dir)

        run_clair3_stage_5a(run_clair3_stage_1c.out.pileup_vcf)
        run_clair3_stage_5b(
            contigs, ref, 
            run_clair3_stage_1c.out.pileup_vcf,
            run_clair3_stage_5a.out.full_qual)

        candidate_beds = run_clair3_stage_5b.out.candidate_bed.flatten()
        // TODO: need to match candidate beds (multiple per contig) with
        //       the correct phased bam
        run_clair3_stage_6b(
            ref, model,
            run_clair3_stage_4.out.phased_bam,
            candidate_beds)
        run_clair3_stage_6c(
            ref,
            run_clair3_stage_6b.out.full_alignment.collect(),
            run_clair3_stage_1a.out.contigs_file)

        //bed_files_to_folder(
        //    run_clair3_stage_5a.out.full_qual,
        //    run_clair3_stage_5b.out.aln_bed,
        //    run_clair3_stage_6a.out.full_alignment
        //    )

        // Merge VCFs, could this been done at the contig level
        // for added parallelism?
        run_clair3_stage_7a(
            contigs, ref,
            run_clair3_stage_1c.out.pileup_vcf,
            run_clair3_stage_6c.out.full_aln_vcf,
            candidate_beds.collect())

        clair_final = run_clair3_stage_7b(
            ref,
            run_clair3_stage_7a.out.merge_output.collect(),
            run_clair3_stage_1a.out.contigs_file)
    emit:
        clair_final[0]
}
    
   
workflow pipeline_with_hap {
    take:
        bam
        bed
        ref
        vcf
        vcf_bed
    main:
        clair_vcf = pipeline(bam, bed, ref)
        hap_output = hap(clair_vcf, ref, vcf, vcf_bed)
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
        // TODO: why do we need a fai? is this a race condition on
        //       on multiple processes trying to create it?
        //       We should likely use https://www.nextflow.io/docs/latest/process.html#storedir
        //       for these things
        ref = Channel.fromPath(params.ref, checkIfExists: true)
        fai = Channel.fromPath(params.ref + ".fai", checkIfExists: true)
        ref = ref.concat(fai).buffer(size: 2)

        bam = Channel.fromPath(params.bam, checkIfExists: true)
        bai = Channel.fromPath(params.bam + ".bai", checkIfExists: true)
        bam = bam.concat(bai).buffer(size: 2)

        // TODO: do we really need a bed? Whats it used for?
        bed = Channel.fromPath(params.bed, checkIfExists: true)
        model = Channel.fromPath(params.model, type: "dir", checkIfExists: true)
 
        if (params.vcf){  // TODO: rename vcf -- its truth VCF for benchmarking
            vcf = Channel.fromPath(params.vcf, checkIfExists:true)
            vcf_bed = Channel.fromPath(params.vcf_bed, checkIfExists:true)
            outputs = pipeline_with_hap(bam, bed, ref, vcf, vcf_bed)
            output(
                outputs[0].concat(outputs[1]).flatten())
        } else {
            pipeline_output = pipeline(bam, bed, ref, model)
            output(pipeline_output)
        }
    }
}
