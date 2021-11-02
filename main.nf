#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

params.help = ""
params.bam = ""
params.ref = ""
params.bed = ""
params.threads = 4
params.download = ""

params.sample_name = "SAMPLE"
params.vcf_fn = "EMPTY"  // use as external candidates
params.ctg_name = "EMPTY"
params.include_all_ctgs = "False"
params.ref_pct_full = 0.1
params.var_pct_full = 0.7
params.GVCF = "False"
params.snp_min_af = 0.0
params.indel_min_af = 0.0

params.truth_vcf = ""
params.truth_bed = ""

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
        path "CHUNK_LIST"
    output:
        stdout
    shell:
        '''
        #!/usr/bin/env python
        with open("CHUNK_LIST") as f:
            for chunk in f:
                print(chunk.strip())
        '''
}


// TODO: this doesn't need to be a process
process generate_contigs {
    label "clair3"
    cpus params.threads
    input:
        path "CONTIGS"
    output:
        stdout
    shell:
        '''
        #!/usr/bin/env python
        with open("CONTIGS") as f:
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
        // TODO: make this explicit, why is pileup VCF optional?
        path "pileup_*.vcf", optional: true, emit: pileup_vcf_chunks
        path "gvcf_tmp_path/*", optional: true, emit: pileup_gvcf_chunks
    script:
        // note: the VCF output here is required to use the contig
        //       name since that's parsed in the SortVcf step
        contig=region[0]
        chunk_id=region[1]
        chunk_num=region[2]
        """
        # TODO: plumb in option gvcf output path
        python \$(which clair3.py) CallVarBam \
            --chkpnt_fn ${model}/pileup \
            --bam_fn ${bam} \
            --call_fn pileup_${contig}_${chunk_id}.vcf \
            --ref_fn ${ref} \
            --ctgName ${contig} \
            --chunk_id ${chunk_id} \
            --chunk_num ${chunk_num} \
            --platform ont \
            --fast_mode False \
            --snp_min_af ${params.snp_min_af} \
            --indel_min_af ${params.indel_min_af} \
            --call_snp_only False \
            --gvcf ${params.GVCF} \
            --temp_file_dir gvcf_tmp_path \
            --pileup
        """ 
}


process run_clair3_stage_1c {
    // Aggregates and sorts all variants (across all chunks of all contigs)
    // from pileup network. Determines quality filter for selecting variants
    // to use for phasing
    label "clair3"
    cpus params.threads
    input:
        tuple path(ref), path(fai)
        // these need to be named as original, as program uses info from
        // contigs file to filter
        path "input_vcfs/*"
        path contigs
    output:
        tuple path("pileup.vcf.gz"), path("pileup.vcf.gz.tbi"), emit: pileup_vcf
        path "phase_qual", emit: phase_qual
    shell:
        '''
        pypy $(which clair3.py) SortVcf \
            --input_dir input_vcfs/ \
            --vcf_fn_prefix pileup \
            --output_fn pileup.vcf \
            --sampleName !{params.sample_name} \
            --ref_fn !{ref} \
            --contigs_fn !{contigs}

        if [ "$( gzip -fdc pileup.vcf.gz | grep -v '#' | wc -l )" -eq 0 ]; \
        then echo "[INFO] Exit in pileup variant calling"; exit 1; fi

        gzip -fdc pileup.vcf.gz | \
            pypy $(which clair3.py) SelectQual --phase --output_fn .
        '''
}


process run_clair3_stage_2b {
    // Filters a VCF by contig, selecting only het SNPs.
    label "clair3"
    cpus params.threads
    input:
        each contig
        tuple path("pileup.vcf.gz"), path("pileup.vcf.gz.tbi")
        // this is used implicitely by the program
        // https://github.com/HKU-BAL/Clair3/blob/329d09b39c12b6d8d9097aeb1fe9ec740b9334f6/preprocess/SelectHetSnp.py#L29
        path "split_folder/phase_qual"
    output:
        tuple val("chr20"), path("split_folder/chr20.vcf"), emit: het_snps_vcf
    shell:
        '''
        pypy $(which clair3.py) SelectHetSnp \
            --vcf_fn pileup.vcf.gz \
            --split_folder split_folder \
            --ctgName !{contig}
        '''
}


process run_clair3_stage_3 {
    // Phases a VCF using whatshap.
    label "clair3"
    cpus params.threads
    input:
        tuple val(contig), path(het_snps), path(bam), path(bai), path(ref), path(fai)
    output:
        tuple val(contig), path("${contig}.bam"), path("${contig}.bam.bai"), emit: phased_bam
    shell:
        '''
        echo "[INFO] 3/7 Phase VCF file using Whatshap"
        whatshap phase \
            --output phased_!{contig}.vcf.gz \
            --reference !{ref} \
            --chromosome !{contig} \
            --distrust-genotypes \
            --ignore-read-groups \
            !{het_snps} \
            !{bam}

        tabix -f -p vcf phased_!{contig}.vcf.gz

        whatshap haplotag \
            --output !{contig}.bam  \
            --reference !{ref} \
            --ignore-read-groups \
            --regions !{contig} \
            phased_!{contig}.vcf.gz \
            !{bam}

        samtools index -@!{task.cpus} !{contig}.bam
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
        path "output/qual", emit: full_qual
    shell:
        '''
        echo "[INFO] 5/7 Select candidates for full-alignment calling"
        mkdir output
        gzip -fdc pileup.vcf.gz | \
        pypy $(which clair3.py) SelectQual \
                --output_fn output \
                --var_pct_full !{params.var_pct_full} \
                --ref_pct_full !{params.ref_pct_full} \
                --platform ont 
        '''
}


process run_clair3_stage_5b {
    // Create BED files for candidate variants for "full alignment" network
    // from the previous full "pileup" variants across all chunks of all chroms
    //
    // Performed per chromosome; output a list of bed files one for each chunk.
    label "clair3"
    cpus params.threads
    input:
        each contig
        tuple path(ref), path(fai)
        tuple path("pileup.vcf.gz"), path("pileup.vcf.gz.tbi")
        // this is used implicitely by the program
        // https://github.com/HKU-BAL/Clair3/blob/329d09b39c12b6d8d9097aeb1fe9ec740b9334f6/preprocess/SelectCandidates.py#L146
        path "candidate_bed/qual"
    output:
        tuple val(contig), path("candidate_bed/${contig}.*"), emit: candidate_bed
    shell:
        // this creates BED files as candidate_bed/<ctg>.0_14 with candidates
        // along with a file the FULL_ALN_FILE_<ctg> listing all of the BED files.
        // All we really want are the BEDs, the file of filenames is used for the
        // purposes of parallel in the original workflow.
        // https://github.com/HKU-BAL/Clair3/blob/329d09b39c12b6d8d9097aeb1fe9ec740b9334f6/scripts/clair3.sh#L218
        '''
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
    // Run "full alignment" network for variants in candidate bed files.
    label "clair3"
    cpus params.threads
    input:
        tuple val(contig), path(phased_bam), path(phased_bam_index)
        tuple val(contig), path(candidate_bed)
        tuple path(ref), path(fai)
        path(model)
    output:
        path "output/full_alignment_*.vcf", emit: full_alignment
    script:
        // TODO: clean this up. 
        f = file(candidate_bed)
        filename = f.getName()
        """
        mkdir output
        echo "[INFO] 6/7 Call low-quality variants using full-alignment model"
        python \$(which clair3.py) CallVarBam \
            --chkpnt_fn $model/full_alignment \
            --bam_fn $phased_bam \
            --call_fn output/full_alignment_${filename}.vcf \
            --sampleName ${params.sample_name} \
            --ref_fn ${ref} \
            --full_aln_regions ${candidate_bed} \
            --ctgName ${contig} \
            --add_indel_length \
            --phasing_info_in_bam \
            --gvcf ${params.GVCF} \
            --snp_min_af ${params.snp_min_af} \
            --indel_min_af ${params.indel_min_af} \
            --platform ont
        """
}


process run_clair3_stage_6c {
    // Sort and merge all "full alignment" variants
    label "clair3"
    cpus params.threads
    input:
        tuple path(ref), path(fai)
        path "full_alignment/*"
        path contigs
        path "gvcf_tmp_path/*"
    output:
        tuple path("full_alignment.vcf.gz"), path("full_alignment.vcf.gz.tbi"), emit: full_aln_vcf
        path "non_var.gvcf", optional: true, emit: non_var_gvcf
    shell:
        '''
        pypy $(which clair3.py) SortVcf \
            --input_dir full_alignment \
            --output_fn full_alignment.vcf \
            --sampleName !{params.sample_name} \
            --ref_fn !{ref} \
            --contigs_fn !{contigs}

        if [ "$( gzip -fdc full_alignment.vcf.gz | grep -v '#' | wc -l )" -eq 0 ]; then
            echo "[INFO] Exit in full-alignment variant calling"
            exit 0
        fi

        # TODO: this could be a separate process
        if [ !{params.GVCF} == True ]; then
            pypy $(which clair3.py) SortVcf \
                --input_dir gvcf_tmp_path \
                --vcf_fn_suffix .tmp.gvcf \
                --output_fn non_var.gvcf \
                --sampleName !{params.sample_name} \
                --ref_fn !{ref} \
                --contigs_fn !{contigs}
        fi
        '''
}


process run_clair3_stage_7a{
    label "clair3"
    cpus params.threads
    input:
        each contig
        tuple path(ref), path(fai)
        tuple path(pile_up_vcf), path(pile_up_vcf_tbi)
        tuple path(full_aln_vcf), path(full_aln_vcf_tbi)
        path "non_var.gvcf"
        path "candidate_beds/*"
    output:
        path "output/merge_${contig}.vcf", emit: merged_vcf
        path "output/merge_${contig}.gvcf", optional: true, emit: merged_gvcf
    shell:
        '''
        # TODO: plumb in non_var.gvcf file channel
        mkdir output
        echo "[INFO] 7/7 Merge pileup VCF and full-alignment VCF"
        pypy $(which clair3.py) MergeVcf \
            --pileup_vcf_fn !{pile_up_vcf} \
            --bed_fn_prefix candidate_beds \
            --full_alignment_vcf_fn !{full_aln_vcf} \
            --output_fn output/merge_!{contig}.vcf \
            --platform ont \
            --print_ref_calls False \
            --gvcf !{params.GVCF} \
            --haploid_precise False \
            --haploid_sensitive False \
            --gvcf_fn output/merge_!{contig}.gvcf \
            --non_var_gvcf_fn non_var.gvcf \
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
        path "merge_outputs_gvcf/*"
        path contigs
    output:
        tuple path("all_contigs.vcf.gz"), path("all_contigs.vcf.gz.tbi"), emit: final_vcf
    shell:
        '''
        pypy $(which clair3.py) SortVcf \
            --input_dir merge_output \
            --vcf_fn_prefix merge \
            --output_fn all_contigs.vcf \
            --sampleName !{params.sample_name} \
            --ref_fn !{ref} \
            --contigs_fn !{contigs}

        if [ "$( gzip -fdc all_contigs.vcf.gz | grep -v '#' | wc -l )" -eq 0 ]; then
            echo "[INFO] Exit in all contigs variant merging"
            exit 0
        fi

        # TODO: this could be a separate process
        if [ !{params.GVCF} == True ]; then
            pypy $(which clair3.py) SortVcf \
                --input_dir merge_outputs_gvcf \
                --vcf_fn_prefix merge \
                --vcf_fn_suffix .gvcf \
                --output_fn all_contigs.gvcf \
                --sampleName !{params.sample_name} \
                --ref_fn !{ref} \
                --contigs_fn !{contigs}
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
        regions2.regions.dump(tag: 'regions')

        run_clair3_stage_1b(regions2.regions, bam, ref, model)
        run_clair3_stage_1c(ref, run_clair3_stage_1b.out.pileup_vcf_chunks.collect(), run_clair3_stage_1a.out.contigs_file)

        run_clair3_stage_2b(
            contigs,
            run_clair3_stage_1c.out.pileup_vcf,
            run_clair3_stage_1c.out.phase_qual)

        // `each` doesn't work with tuples, so we have to make the product ourselves
        phase_inputs = run_clair3_stage_2b.out.het_snps_vcf.combine(bam).combine(ref)
        run_clair3_stage_3(phase_inputs)

        run_clair3_stage_5a(run_clair3_stage_1c.out.pileup_vcf)
        run_clair3_stage_5b(
            contigs, ref, 
            run_clair3_stage_1c.out.pileup_vcf,
            run_clair3_stage_5a.out.full_qual)

        // Have to go through a bit of a song and dance here, stay with me.
        candidate_beds = run_clair3_stage_5b.out.candidate_bed.flatMap {
            x ->
                y = x[1]; if(! (y instanceof java.util.ArrayList)){y = [y]}
                y.collect { [x[0], it] } }  // [chr, bed]
        bams_beds_and_stuff = run_clair3_stage_3.out.phased_bam
            .cross(candidate_beds)
            .combine(ref.map {it->[it]})
            .combine(model) // [[chr, bam, bai], [chr20, bed], [ref, fai], model] 

        bams_beds_and_stuff.multiMap { it ->
            bams: it[0]
            candidates: it[1]
            ref: it[2]
            model: it[3]
        }.set { mangled }
        run_clair3_stage_6b(
            mangled.bams, mangled.candidates, mangled.ref, mangled.model)
        // merge and sort all files for all chunks for all contigs
        // gvcf is optional, stuff an empty file in, so we have at least one
        // item to flatten/collect and tthis stage can run.
        gvcfs = run_clair3_stage_1b.out.pileup_gvcf_chunks
            .flatten()
            .ifEmpty(file("$projectDir/data/OPTIONAL_FILE"))
            .collect()
        run_clair3_stage_1b.out.pileup_gvcf_chunks.flatten().collect()
        run_clair3_stage_6c(
            ref,
            run_clair3_stage_6b.out.full_alignment.collect(),
            run_clair3_stage_1a.out.contigs_file,
            gvcfs)

        // merge "pileup" and "full alignment" variants, per contig
        // note: the candidate beds aren't actually used by the program for ONT
        non_var_gvcf = run_clair3_stage_6c.out.non_var_gvcf
            .ifEmpty(file("$projectDir/data/OPTIONAL_FILE"))
        run_clair3_stage_7a(
            contigs, ref,
            run_clair3_stage_1c.out.pileup_vcf,
            run_clair3_stage_6c.out.full_aln_vcf,
            non_var_gvcf,
            candidate_beds.map {it->it[1] }.collect())

        // merge all files across contigs
        gvcfs = run_clair3_stage_7a.out.merged_gvcf
            .ifEmpty(file("$projectDir/data/OPTIONAL_FILE"))
        clair_final = run_clair3_stage_7b(
            ref,
            run_clair3_stage_7a.out.merged_vcf.collect(),
            gvcfs.collect(),
            run_clair3_stage_1a.out.contigs_file)
    emit:
        clair_final
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
 
        if (params.truth_vcf){
            vcf = Channel.fromPath(params.truth_vcf, checkIfExists:true)
            vcf_bed = Channel.fromPath(params.truth_bed, checkIfExists:true)
            outputs = pipeline_with_hap(bam, bed, ref, vcf, vcf_bed)
            output(outputs[0].concat(outputs[1]).flatten())
        } else {
            pipeline_output = pipeline(bam, bed, ref, model)
            output(pipeline_output)
        }
    }
}
