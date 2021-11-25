#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

import groovy.json.JsonBuilder


process make_chunks {
    // Do some preliminaries. Ordinarily this would setup a working directory
    // that all other commands would make use off, but all we need here are the
    // list of contigs and chunks.
    label "clair3"
    cpus 1
    input:
        tuple path(bam), path(bai)
        tuple path(ref), path(fai)
        path bed
    output:
        path "clair_output/tmp/CONTIGS", emit: contigs_file
        path "clair_output/tmp/CHUNK_LIST", emit: chunks_file
    script:
        def bedargs = bed.name != 'OPTIONAL_FILE' ? "--bed_fn ${bed}" : ''
        """
        mkdir -p clair_output
        python \$(which clair3.py) CheckEnvs \
            --bam_fn ${bam} \
            ${bedargs} \
            --output_fn_prefix clair_output \
            --ref_fn ${ref} \
            --vcf_fn ${params.vcf_fn} \
            --ctg_name ${params.ctg_name} \
            --chunk_num 0 \
            --chunk_size 5000000 \
            --include_all_ctgs ${params.include_all_ctgs} \
            --threads 1  \
            --qual 2 \
            --sampleName ${params.sample_name} \
            --var_pct_full ${params.var_pct_full} \
            --ref_pct_full ${params.ref_pct_full} \
            --snp_min_af ${params.snp_min_af} \
            --indel_min_af ${params.indel_min_af}
        """
}


process pileup_variants {
    // Calls variants per region ("chunk") using pileup network.
    label "clair3"
    cpus 2
    errorStrategy 'retry'
    maxRetries params.call_retries
    input:
        each region
        tuple path(bam), path(bai)
        tuple path(ref), path(fai)
        path model
    output:
        // TODO: make this explicit, why is pileup VCF optional?
        path "pileup_*.vcf", optional: true, emit: pileup_vcf_chunks
        path "gvcf_tmp_path/*", optional: true, emit: pileup_gvcf_chunks
    shell:
        // note: the VCF output here is required to use the contig
        //       name since that's parsed in the SortVcf step
        // Set TF threads to 1, default is 4 but less than 100%
        // observed in tests.
        '''
        python $(which clair3.py) CallVarBam \
            --tensorflow_threads 1 \
            --chkpnt_fn !{model}/pileup \
            --bam_fn !{bam} \
            --call_fn pileup_!{region.contig}_!{region.chunk_id}.vcf \
            --ref_fn !{ref} \
            --ctgName !{region.contig} \
            --chunk_id !{region.chunk_id} \
            --chunk_num !{region.total_chunks} \
            --platform ont \
            --fast_mode False \
            --snp_min_af !{params.snp_min_af} \
            --indel_min_af !{params.indel_min_af} \
            --call_snp_only False \
            --gvcf !{params.GVCF} \
            --temp_file_dir gvcf_tmp_path \
            --pileup
        ''' 
}


process aggregate_pileup_variants {
    // Aggregates and sorts all variants (across all chunks of all contigs)
    // from pileup network. Determines quality filter for selecting variants
    // to use for phasing.
    label "clair3"
    cpus 2
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

        if [ "$( bgzip -fdc pileup.vcf.gz | grep -v '#' | wc -l )" -eq 0 ]; \
        then echo "[INFO] Exit in pileup variant calling"; exit 1; fi

        bgzip -fdc pileup.vcf.gz | \
            pypy $(which clair3.py) SelectQual --phase --output_fn .
        '''
}


process select_het_snps {
    // Filters a VCF by contig, selecting only het SNPs.
    label "clair3"
    cpus 2
    input:
        each contig
        tuple path("pileup.vcf.gz"), path("pileup.vcf.gz.tbi")
        // this is used implicitely by the program
        // https://github.com/HKU-BAL/Clair3/blob/329d09b39c12b6d8d9097aeb1fe9ec740b9334f6/preprocess/SelectHetSnp.py#L29
        path "split_folder/phase_qual"
    output:
        tuple val(contig), path("split_folder/${contig}.vcf.gz"), path("split_folder/${contig}.vcf.gz.tbi"), emit: het_snps_vcf
    shell:
        '''
        pypy $(which clair3.py) SelectHetSnp \
            --vcf_fn pileup.vcf.gz \
            --split_folder split_folder \
            --ctgName !{contig}

        bgzip -c split_folder/!{contig}.vcf > split_folder/!{contig}.vcf.gz
        tabix split_folder/!{contig}.vcf.gz
        '''
}


process phase_contig {
    // Phases a VCF using whatshap and selected het SNPS. Tags reads in the
    // input BAM.
    label "clair3"
    cpus 2
    input:
        tuple val(contig), path(het_snps), path(het_snps_tbi), path(bam), path(bai), path(ref), path(fai)
    output:
        tuple val(contig), path("${contig}.bam"), path("${contig}.bam.bai"), emit: phased_bam
    shell:
        '''
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


process phase_region {
    label "clair3"
    cpus 1
    input:
        tuple val(contig), path(het_snps), path(het_snps_tbi), path(bam), path(bai), path(ref), path(fai)
        val(region)
    output:
        tuple val(contig), path("phased_${region.replace(':','-')}.vcf.gz"), emit: region_vcf
    shell:
        '''
        bcftools view !{het_snps} -r !{region} | bgzip -c > region_hets.vcf.gz
        tabix region_hets.vcf.gz

        # if VCF is empty, whatshap doesn't write a file
        whatshap phase \
            --output phased_!{region.replace(":","-")}.vcf.gz \
            --reference !{ref} \
            --chromosome !{contig} \
            --distrust-genotypes \
            --ignore-read-groups \
            region_hets.vcf.gz \
            !{bam}
        '''
}


process haplotag_phase_regions {
    label "clair3"
    cpus 1
    input:
        tuple val(contig), path(bam), path(bai), path(ref), path(fai), path("vcfs/*")
    output:
        tuple val(contig), path("${contig}.bam"), path("${contig}.bam.bai"), emit: phased_bam
    shell:
        '''
        # Run a python program to join VCFs
        bcftools concat vcfs/* | bcftools sort | bgzip -c > phased_!{contig}.vcf.gz
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


process get_contig_chunks {
    label "clair3"
    cpus 1
    input:
        val(contig)
        path("fasta_index.fai")
        val(chunk_size)
        val(min_chunk_size)
    output:
        path("regions.bed")
    shell:
        '''
        #!/usr/bin/env python
        import math

        contig = "!{contig}"
        chunk_size = !{chunk_size}
        min_chunk_size = !{min_chunk_size}

        def make_chunks(size):
            # don't let last chunk become too small
            chunks = math.ceil(size / chunk_size)
            remain = size % chunk_size
            if remain < min_chunk_size:
                chunks -= 1
            chsize = math.ceil(size / chunks)
            start = 0
            while start < size:
                end = min(start + chsize, size)
                yield start, end
                start = end

        chr_length = None
        with open("fasta_index.fai") as fh:
            for line in fh.readlines():
                chr, length, *_ = line.split()
                if chr == contig:
                    chr_length = int(length)
                    break

        with open("regions.bed", "w") as out:
            for start, end in make_chunks(chr_length): 
                out.write(f"{chr}:{start}-{end}\\n")
        '''
}


workflow sharded_phase {
    take:
        phase_inputs
    main:
        //phase_outputs = phase_contig(phase_inputs)
        phase_inputs.multiMap {
            it ->
                contig: it[0]
                fai: it[6]
        }.set { chunks }
        res = get_contig_chunks(
			chunks.contig, chunks.fai,
			params.phase_chunk, params.phase_chunk_min)
        regions = res.splitText() {
            chr = (it =~ /(.+):/)[0]
            [chr[1], it.trim()] }
        // make an auxiliary channel counting the number of regions per-chromosome
        reg_counts = regions
            .groupTuple()
            .map {chr, regs -> tuple(chr, groupKey(chr, regs.size())) }
       
        // create inputs to phase sub-contig regions
        phase_inputs
            .map { [it[0], it] }  // add contig key
            .cross(regions)       // duplicate inputs across regions
            .multiMap { it -> 
                phase_inputs: it[0][1]
                region: it[1][1]
            }.set { region_inputs }
        region_vcfs = phase_region(region_inputs.phase_inputs, region_inputs.region)
        
        // collect region VCFs by chrom and join to phase_inputs
        // we first join on chr, and then group by the groupKey
        // with knowledge of the sizes so big onward tasks for
        // different chroms don't wait for all region_vcf tasks.
        reg_counts
            .cross(phase_region.out.region_vcf)  // cross is performed using chr
            .map { i, j -> [i[1], i[0], j[1]] }  // i[1] is *chr
            .groupTuple()                        // uses *chr as key
            .map { i, j, k -> [j[0], k]}         // dispense with *chr, only need one chr!
            .set { vcfs_by_chrom }
        stuff = phase_inputs
            .map { [it[0], it] }  // add join key
            .join(vcfs_by_chrom)  // join
            .map {                // pull out the bits we need
                // chrom, bam, bai, ref, fai, vcfs
                [it[0], it[1][3], it[1][4], it[1][5], it[1][6], it[2]]
            }
        output = haplotag_phase_regions(stuff)
    emit:
        output
}


process get_qual_filter {
    // Determines quality filter for selecting candidate variants for second
    // stage "full alignment" calling.
    label "clair3"
    cpus 2
    input:
        tuple path("pileup.vcf.gz"), path("pileup.vcf.gz.tbi")
    output:
        path "output/qual", emit: full_qual
    shell:
        '''
        echo "[INFO] 5/7 Select candidates for full-alignment calling"
        mkdir output
        bgzip -fdc pileup.vcf.gz | \
        pypy $(which clair3.py) SelectQual \
                --output_fn output \
                --var_pct_full !{params.var_pct_full} \
                --ref_pct_full !{params.ref_pct_full} \
                --platform ont 
        '''
}


process create_candidates {
    // Create BED files for candidate variants for "full alignment" network
    // from the previous full "pileup" variants across all chunks of all chroms
    //
    // Performed per chromosome; output a list of bed files one for each chunk.
    label "clair3"
    cpus 2
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
        // This creates BED files as candidate_bed/<ctg>.0_14 with candidates
        // along with a file the FULL_ALN_FILE_<ctg> listing all of the BED
        // files.  All we really want are the BEDs, the file of filenames is
        // used for the purposes of parallel in the original workflow.
        // https://github.com/HKU-BAL/Clair3/blob/329d09b39c12b6d8d9097aeb1fe9ec740b9334f6/scripts/clair3.sh#L218

        // TODO: would be nice to control the number of BEDs produced to enable
        // better parallelism.
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


process evaluate_candidates {
    // Run "full alignment" network for variants in a candidate bed file.
    label "clair3"
    cpus 2
    errorStrategy 'retry'
    maxRetries params.call_retries
    input:
        tuple val(contig), path(phased_bam), path(phased_bam_index)
        tuple val(contig), path(candidate_bed)
        tuple path(ref), path(fai)
        path(model)
    output:
        path "output/full_alignment_*.vcf", emit: full_alignment
    script:
        filename = candidate_bed.name
        """
        mkdir output
        echo "[INFO] 6/7 Call low-quality variants using full-alignment model"
        python \$(which clair3.py) CallVarBam \
            --tensorflow_threads 1 \
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


process aggregate_full_align_variants {
    // Sort and merge all "full alignment" variants
    label "clair3"
    cpus 2
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

        if [ "$( bgzip -fdc full_alignment.vcf.gz | grep -v '#' | wc -l )" -eq 0 ]; then
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


process merge_pileup_and_full_vars{
    // Merge VCFs
    label "clair3"
    cpus 2
    input:
        each contig
        tuple path(ref), path(fai)
        tuple path(pile_up_vcf), path(pile_up_vcf_tbi)
        tuple path(full_aln_vcf), path(full_aln_vcf_tbi)
        path "non_var.gvcf"
        path "candidate_beds/*"
    output:
        tuple val(contig), path("output/merge_${contig}.vcf.gz"), path("output/merge_${contig}.vcf.gz.tbi"), emit: merged_vcf
        path "output/merge_${contig}.gvcf", optional: true, emit: merged_gvcf
    shell:
        '''
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

        bgzip -c output/merge_!{contig}.vcf > output/merge_!{contig}.vcf.gz
        tabix output/merge_!{contig}.vcf.gz
        '''
}


process aggregate_all_variants{
    label "clair3"
    cpus 4
    input:
        tuple path(ref), path(fai)
        path "merge_output/*"
        path "merge_outputs_gvcf/*"
        val(phase_vcf)
        path contigs
    output:
        tuple path("all_contigs.vcf.gz"), path("all_contigs.vcf.gz.tbi"), emit: final_vcf
    shell:
        '''
        prefix="merge"
        phase_vcf=!{params.phase_vcf}
        if [[ $phase_vcf == "true" ]]; then
            prefix="phased"
        fi
        ls merge_output/*.vcf.gz | parallel --jobs 4 "bgzip -d {}"

        pypy $(which clair3.py) SortVcf \
            --input_dir merge_output \
            --vcf_fn_prefix $prefix \
            --output_fn all_contigs.vcf \
            --sampleName !{params.sample_name} \
            --ref_fn !{ref} \
            --contigs_fn !{contigs}

        if [ "$( bgzip -fdc all_contigs.vcf.gz | grep -v '#' | wc -l )" -eq 0 ]; then
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


process hap {
    label "happy"
    input:
        tuple path("clair.vcf.gz"), path("clair.vcf.gz.tbi")
        tuple path("ref.fasta"), path("ref.fasta.fai")
        path "truth.vcf"
        path "truth.bed"
    output:
        path "happy"
    shell:
        '''
        mkdir happy
        /opt/hap.py/bin/hap.py \
            truth.vcf \
            clair.vcf.gz \
            -f truth.bed \
            -r ref.fasta \
            -o happy \
            --engine=vcfeval \
            --threads=4 \
            --pass-only
        '''
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


process readStats {
    label "clair3"
    cpus 1
    input:
        tuple path("alignments.bam"), path("alignments.bam.bai")
    output:
        path "readstats.txt", emit: stats
    """
    # using multithreading inside docker container does strange things
    bamstats --threads 1 alignments.bam > readstats.txt
    """
}


process getVersions {
    label "clair3"
    cpus 1
    output:
        path "versions.txt"
    script:
        """
        run_clair3.sh --version | sed 's/ /,/' >> versions.txt
        """
}

process getParams {
    label "clair3"
    cpus 1
    output:
        path "params.json"
    script:
        def paramsJSON = new JsonBuilder(params).toPrettyString()
        """
        # Output nextflow params object to JSON
        echo '$paramsJSON' > params.json
        """
}


process vcfStats {
    label "clair3"
    cpus 1
    input:
        tuple path(vcf), path(index)
    output:
        file "variants.stats"
    """
    bcftools stats $vcf > variants.stats
    """
}


process makeReport {
    label "clair3"
    input:
        file read_summary
        file vcfstats
        path versions
        path "params.json"
    output:
        path "wf-human-snp-*.html"
    script:
        report_name = "wf-human-snp-" + params.report_name + '.html'
        """
        report.py $report_name --versions $versions --params params.json \
        --read_stats $read_summary \
        --vcf_stats $vcfstats
        """
}



// workflow module
workflow clair3 {
    take:
        bam
        bed
        ref
        model
    main:

        // Run preliminaries to find contigs and generate regions to process in
        // parallel.
        // > Step 0
        make_chunks(bam, ref, bed)
        chunks = make_chunks.out.chunks_file
            .splitText(){ 
                cols = (it =~ /(.+)\s(.+)\s(.+)/)[0]
                ["contig": cols[1], "chunk_id":cols[2], "total_chunks":cols[3]]}
        contigs = make_chunks.out.contigs_file.splitText() { it.trim() }
        
        // Run the "pileup" caller on all chunks and collate results
        // > Step 1 
        pileup_variants(chunks, bam, ref, model)
        aggregate_pileup_variants(
            ref, pileup_variants.out.pileup_vcf_chunks.collect(),
            make_chunks.out.contigs_file)

        // Filter collated results to produce per-contig SNPs for phasing.
        // > Step 2
        select_het_snps(
            contigs,
            aggregate_pileup_variants.out.pileup_vcf,
            aggregate_pileup_variants.out.phase_qual)

        // Perform whatshap phasing for each contig.  
        // `each` doesn't work with tuples, so we have to make the product ourselves
        phase_inputs = select_het_snps.out.het_snps_vcf
            .combine(bam).combine(ref)
        // > Step 3 + Step 4
        if (params.parallel_phase) {
            sharded_phase(phase_inputs).set { phased_bam }
        } else {
            phase_contig(phase_inputs)
            phase_contig.out.phased_bam.set { phased_bam }
        }

        // Find quality filter to select variants for "full alignment"
        // processing, then generate bed files containing the candidates.
        // > Step 5
        get_qual_filter(aggregate_pileup_variants.out.pileup_vcf)
        create_candidates(
            contigs, ref, 
            aggregate_pileup_variants.out.pileup_vcf,
            get_qual_filter.out.full_qual)

        // Run the "full alignment" network on candidates. Have to go through a
        // bit of a song and dance here to generate our input channels here
        // with various things duplicated (again because of limitations on 
        // `each` and tuples).
        // > Step 6
        candidate_beds = create_candidates.out.candidate_bed.flatMap {
            x ->
                // output globs can return a list or single item
                y = x[1]; if(! (y instanceof java.util.ArrayList)){y = [y]}
                // effectively duplicate chr for all beds - [chr, bed]
                y.collect { [x[0], it] } }
        // produce something emitting: [[chr, bam, bai], [chr20, bed], [ref, fai], model]
        bams_beds_and_stuff = phased_bam
            .cross(candidate_beds)
            .combine(ref.map {it->[it]})
            .combine(model)
        // take the above and destructure it for easy reading
        bams_beds_and_stuff.multiMap {
            it ->
                bams: it[0]
                candidates: it[1]
                ref: it[2]
                model: it[3]
            }.set { mangled }
        // phew! Run all-the-things
        evaluate_candidates(
            mangled.bams, mangled.candidates, mangled.ref, mangled.model)

        // merge and sort all files for all chunks for all contigs
        // gvcf is optional, stuff an empty file in, so we have at least one
        // item to flatten/collect and tthis stage can run.
        gvcfs = pileup_variants.out.pileup_gvcf_chunks
            .flatten()
            .ifEmpty(file("$projectDir/data/OPTIONAL_FILE"))
            .collect()
        pileup_variants.out.pileup_gvcf_chunks.flatten().collect()
        aggregate_full_align_variants(
            ref,
            evaluate_candidates.out.full_alignment.collect(),
            make_chunks.out.contigs_file,
            gvcfs)

        // merge "pileup" and "full alignment" variants, per contig
        // note: we never create per-contig VCFs, so this process
        //       take the whole genome VCFs and the list of contigs
        //       to produce per-contig VCFs which are then finally
        //       merge to yield the whole genome results.

        // First merge whole-genome results from pileup and full_alignment
        //   for each contig...
        // note: the candidate beds aren't actually used by the program for ONT
        // > Step 7
        non_var_gvcf = aggregate_full_align_variants.out.non_var_gvcf
            .ifEmpty(file("$projectDir/data/OPTIONAL_FILE"))
        merge_pileup_and_full_vars(
            contigs, ref,
            aggregate_pileup_variants.out.pileup_vcf,
            aggregate_full_align_variants.out.full_aln_vcf,
            non_var_gvcf,
            candidate_beds.map {it->it[1] }.collect())

        if (params.phase_vcf) {
            data = merge_pileup_and_full_vars.out.merged_vcf
                .combine(bam).combine(ref).view()
            chrom = data.map { it -> it[0] }
            phase_region(data, chrom)
                .map { it -> [it[1]] }
                .set { final_vcfs }
        } else {
            merge_pileup_and_full_vars.out.merged_vcf
                .map { it -> [it[1]] }
                .set { final_vcfs }
        }

        // ...then collate final per-contig VCFs for whole genome results
        gvcfs = merge_pileup_and_full_vars.out.merged_gvcf
            .ifEmpty(file("$projectDir/data/OPTIONAL_FILE"))
        clair_final = aggregate_all_variants(
            ref,
            final_vcfs.collect(),
            gvcfs.collect(),
            params.phase_vcf,
            make_chunks.out.contigs_file)

        // reporting
        read_stats = readStats(bam)
        software_versions = getVersions()
        workflow_params = getParams()
        vcf_stats = vcfStats(clair_final)
        report = makeReport(
            read_stats.stats, vcf_stats[0],
            software_versions.collect(), workflow_params)

    emit:
        clair_final.concat(report)
}
    
   
workflow happy_evaluation {
    take:
        clair_vcf
        ref
        vcf
        vcf_bed
    main:
        hap_output = hap(clair_vcf, ref, vcf, vcf_bed)
    emit:
        hap_output
}


// entrypoint workflow
WorkflowMain.initialise(workflow, params, log)
workflow {
    // TODO: why do we need a fai? is this a race condition on
    //       on multiple processes trying to create it?
    //       We should likely use https://www.nextflow.io/docs/latest/process.html#storedir
    //       for these things
    ref = Channel.fromPath(params.ref)
    fai = Channel.fromPath(params.ref + ".fai")
    ref = ref.concat(fai).buffer(size: 2)

    bam = Channel.fromPath(params.bam)
    bai = Channel.fromPath(params.bam + ".bai")
    bam = bam.concat(bai).buffer(size: 2)

    if(params.bed == null) {
        bed = Channel.fromPath("${projectDir}/data/OPTIONAL_FILE")
    } else {
        bed = Channel.fromPath(params.bed)
    }
    model = Channel.fromPath(params.model, type: "dir")

    if (params.truth_vcf){
        // TODO: test this
        truth_vcf = Channel.fromPath(params.truth_vcf)
        truth_bed = Channel.fromPath(params.truth_bed)
        clair_vcf = clair3(bam, bed, ref)
        happy_results = happy_evaluation(clair_vcf, ref, truth_vcf, truth_bed)
        output(clair_vcf.concat(happy_results).flatten())
    } else {
        clair_vcf = clair3(bam, bed, ref, model)
        output(clair_vcf.flatten())
    }
}
