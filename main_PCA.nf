#!/usr/bin/env nextflow
// LASER_PCA_AncestryInference — main pipeline (user-facing template)

nextflow.enable.dsl = 2

/*
 * Author: Justin Pelletier (inspired by Peyton McClelland)
 * Year: 2025
 * Pipeline: LASER PCA Projection with QC
 *
 * REQUIREMENTS
 * - Input VCFs are bgzipped + tabix-indexed, one per chromosome named "chrN.vcf.gz" (+ .tbi)
 * - Sample include-lists are text files (one ID per line) that match VCF sample names
 * - `params.lowcomplexity_bed` is hg38; set `params.header_bed` to "TRUE" or "FALSE"
 * - LASER directory contains `laser`, `trace`, and `vcf2geno` binaries (see `params.path_to_laser`)
 * - `bin/PCA.R` exists (Nextflow adds `bin/` to PATH automatically)
 *
 * NOTES
 * - This script uses generic `module load` examples; remove or adapt to your env/containers.
 */

// ---- User-tunable defaults (can be overridden on CLI) -----------------------
params.nPCs     = params.nPCs     ?: 20
params.min_prob = params.min_prob ?: 0
params.seed     = params.seed     ?: 11
params.outdir   = params.outdir   ?: "results/LASER_PCA"

// ----------------------------------------------------------------------------
// 1) Normalize & QC: Reference
// ----------------------------------------------------------------------------
process qc_norm_ref {
    tag "$chr"
    input:
        tuple val(chr), path(vcf), path(vcf_tbi)
    output:
        tuple val(chr), path("ref_${chr}.qc.vcf.gz"), path("ref_${chr}.qc.vcf.gz.tbi")
    script:
    """
    # Adapt/remove module lines if using containers or a local toolchain
    module load bcftools

    bcftools norm -m -both -f "${params.ref_fasta}" "$vcf" | \
      bcftools view -f PASS -q 0.05 -Q 0.95 | \
      bcftools annotate -x INFO,^GT | \
      bcftools view -S "${params.qc_ref_list}" --force-samples | \
      bcftools annotate --set-id '%CHROM:%POS:%REF:%ALT' -Oz -o ref_${chr}.qc.vcf.gz

    tabix -f ref_${chr}.qc.vcf.gz
    """
}

// ----------------------------------------------------------------------------
// 2) Normalize & QC: Study
// ----------------------------------------------------------------------------
process qc_norm_study {
    tag "$chr"
    input:
        tuple val(chr), path(vcf), path(vcf_tbi)
    output:
        tuple val(chr), path("study_${chr}.qc.vcf.gz"), path("study_${chr}.qc.vcf.gz.tbi")
    script:
    """
    module load bcftools

    bcftools norm -m -both -f "${params.ref_fasta}" "$vcf" | \
      bcftools view -q 0.05 -Q 0.95 | \
      bcftools annotate -x INFO,^GT | \
      bcftools view -S "${params.qc_study_list}" --force-samples | \
      bcftools annotate --set-id '%CHROM:%POS:%REF:%ALT' -Oz -o study_${chr}.qc.vcf.gz

    tabix -f study_${chr}.qc.vcf.gz
    """
}

// ----------------------------------------------------------------------------
// 3) Intersect shared variants (reference ∩ study)
// ----------------------------------------------------------------------------
process intersect {
  tag "$chr"
  input:
    tuple val(chr), path(ref), path(ref_tbi), path(study), path(study_tbi)
  output:
    tuple val(chr), path("${chr}.merged.vcf.gz"), path("${chr}.merged.vcf.gz.tbi")
  script:
  """
  module load bcftools
  set -euo pipefail

  tmp=\$(mktemp -d)
  trap 'rm -rf "\$tmp"' EXIT

  bcftools isec -n=2 -p "\$tmp" -Oz "$ref" "$study"
  bcftools merge -m all -Oz -o "${chr}.merged.vcf.gz" "\$tmp/0000.vcf.gz" "\$tmp/0001.vcf.gz"
  tabix -f "${chr}.merged.vcf.gz"
  """
}

// ----------------------------------------------------------------------------
// 4) Final QC + LD pruning (PLINK), then convert back to VCF
// ----------------------------------------------------------------------------
process final_qc_and_prune {
  tag "$chr"
  input:
    tuple val(chr), path(vcf), path(vcf_tbi)
  output:
    tuple val(chr), path("${chr}.pruned.vcf.gz"), path("${chr}.pruned.vcf.gz.tbi")
  script:
  """
  module load plink/1.9 bcftools

  # Prepare low-complexity BED (drop header if present; keep first 3 cols)
  if [ "${params.header_bed}" = "TRUE" ]; then
    tail -n +2 ${params.lowcomplexity_bed} | cut -f1-3 > lowcomplexity.clean.bed
  else
    cut -f1-3 ${params.lowcomplexity_bed} > lowcomplexity.clean.bed
  fi

  # Exclude low-complexity/MHC regions
  bcftools view -T ^lowcomplexity.clean.bed "$vcf" -Oz -o ${chr}.filtered.vcf.gz

  # PLINK QC
  plink --vcf ${chr}.filtered.vcf.gz \
        --double-id \
        --real-ref-alleles \
        --snps-only \
        --maf 0.01 \
        --geno 0.001 \
        --biallelic-only strict \
        --make-bed \
        --out ${chr}.final_qc

  # LD pruning
  plink --bfile ${chr}.final_qc \
        --indep-pairwise 200 50 0.2 \
        --out ${chr}.prune

  plink --bfile ${chr}.final_qc \
        --extract ${chr}.prune.prune.in \
        --make-bed \
        --out ${chr}.pruned

  # Convert back to VCF
  plink --bfile ${chr}.pruned \
        --recode vcf bgz \
        --out ${chr}.pruned

  tabix -f ${chr}.pruned.vcf.gz
  """
}

// ----------------------------------------------------------------------------
// 5) Concatenate per-chromosome pruned VCFs
// ----------------------------------------------------------------------------
process concat_pruned_vcfs {
    input:
      path vcf_files
    output:
      path("allchr.pruned.vcf.gz")
    script:
    """
    module load bcftools
    bcftools concat -Oz -o allchr.pruned.vcf.gz ${vcf_files.join(' ')}
    tabix -f allchr.pruned.vcf.gz
    """
}

// ----------------------------------------------------------------------------
// 6) Double the sample IDs (FID IID) for PLINK/LASER conventions
// ----------------------------------------------------------------------------
process double_ref_ids {
  input:
    path qc_ref_list
  output:
    path "ref_list.id"
  script:
  """
  awk '{print \$1\"_\"\$1}' ${qc_ref_list} > ref_list.id
  """
}

process double_study_ids {
  input:
    path qc_study_list
  output:
    path "study_list.id"
  script:
  """
  awk '{print \$1\"_\"\$1}' ${qc_study_list} > study_list.id
  """
}

// ----------------------------------------------------------------------------
// 7) Split the doubled study list into chunks (default: 1000 per batch)
// ----------------------------------------------------------------------------
process split_study_list {
  input:
    path study_list
  output:
    path "study_batch_*.id"
  script:
  """
  split -d -l ${params.batch_size} --additional-suffix .id ${study_list} study_batch_
  """
}

// ----------------------------------------------------------------------------
// 8) Prepare reference LASER inputs + compute reference PCA (once)
// ----------------------------------------------------------------------------
process prepare_reference {
  tag "prepare_reference"
  publishDir "${params.outdir}", mode: "copy"

  input:
    path vcf
    path ref_list

  output:
    path "reference_pruned.vcf.gz"
    path "reference_pruned.vcf.gz.tbi"
    path "reference_pruned.geno"
    path "reference_pruned.site"
    path "reference.RefPC.coord"

  script:
  """
  module load bcftools

  # Extract only reference samples
  bcftools view -S ${ref_list} --force-samples -Oz -o reference_pruned.vcf.gz ${vcf}
  tabix -f reference_pruned.vcf.gz

  # Convert to LASER geno/site
  ${params.path_to_laser}/vcf2geno/vcf2geno --inVcf reference_pruned.vcf.gz --out reference_pruned

  # Run PCA on reference
  ${params.path_to_laser}/laser \
    -g reference_pruned.geno \
    -k ${params.nPCs} \
    -pca 1 \
    -o reference
  """
}

// ----------------------------------------------------------------------------
// 9) Project each N-sample batch onto the reference PCs (parallel)
// ----------------------------------------------------------------------------
process laser_pca_projection_batch {
  publishDir "${params.outdir}/batches", mode: "copy"
  input:
    tuple path(batch_file),
          path(vcf),
          path(ref_geno),
          path(ref_site),
          path(ref_pca)
  output:
    path "trace_*.ProPC.coord"
  script:
  """
  module load bcftools

  batch_id=\$(basename ${batch_file} .id)

  # Extract this batch
  bcftools view -S ${batch_file} --force-samples -Oz -o study_\${batch_id}.vcf.gz ${vcf}
  tabix -f study_\${batch_id}.vcf.gz

  # Convert to LASER geno/site
  ${params.path_to_laser}/vcf2geno/vcf2geno --inVcf study_\${batch_id}.vcf.gz --out study_\${batch_id}

  # LASER trace config for projection
  cat <<EOF > trace_\${batch_id}.conf
STUDY_FILE    study_\${batch_id}.geno
GENO_FILE     ${ref_geno}
COORD_FILE    ${ref_pca}
FIRST_IND     1
LAST_IND      \$(wc -l < study_\${batch_id}.geno)
OUT_PREFIX    trace_\${batch_id}
DIM           ${params.nPCs}
EOF

  ${params.path_to_laser}/trace -p trace_\${batch_id}.conf
  """
}

// ----------------------------------------------------------------------------
// 10) Merge all batches with the reference PCA
//     (writes a trimmed coord file + a full concat)
// ----------------------------------------------------------------------------
process merge_proj_coords {
  publishDir "${params.outdir}", mode: "copy"

  input:
    path ref_pca
    path proj_files
  output:
    path "final.ProPC.coord"

  script:
  """
  # 1) Start with cleaned reference PCA (keep IDs + PC columns)
  awk 'NR==1 { print; next }
       {
         # Clean FID and IID: ref1_ref1 → ref1
         for (i=1; i<=2; i++) { split(\$i, parts, "_"); \$i = parts[1]; }
         printf "%s\\t%s", \$1, \$2;
         for (i=3; i<=NF; i++) printf "\\t%s", \$i;
         printf "\\n";
       }' ${ref_pca} > final.ProPC.coord

  # 2) Append cleaned batch projections (keep IDs + PC columns)
  for f in ${proj_files.join(' ')}; do
    tail -n +2 \$f | awk '{
      # Clean FID and IID: sample1_sample1_sample1_sample1 → sample1_sample1
      for (i=1; i<=2; i++) { split(\$i, parts, "_"); \$i = parts[1] "_" parts[2]; }
      printf "%s\\t%s", \$1, \$2;
      for (i=7; i<=NF; i++) printf "\\t%s", \$i;  # PCs start at col 7 in LASER outputs
      printf "\\n";
    }' >> final.ProPC.coord
  done

  # Optional: also produce a full merged file of all projected batches
  head -n 1 ${proj_files[0]} > final.ProPC.full.coord
  for f in ${proj_files.join(' ')}; do
    tail -n +2 \$f | awk '{
      for (i=1; i<=2; i++) { split(\$i, parts, "_"); \$i = parts[1] "_" parts[2]; }
      printf "%s\\t%s", \$1, \$2;
      for (i=3; i<=NF; i++) printf "\\t%s", \$i;
      printf "\\n";
    }' >> final.ProPC.full.coord
  done
  """
}

// ----------------------------------------------------------------------------
// 11) Plot PCA using repo's bin/PCA.R (no absolute paths)
// ----------------------------------------------------------------------------
process run_pca_analysis {
  tag "run_pca_analysis"
  publishDir "results/PCA_plots", mode: "copy"

  input:
    path pca_file
    val mode

  output:
    path "*"


  script:
  """
  module load r
  mkdir -p plot_PCA

  Rscript PCA.R ${pca_file} ${params.meta_file} ${params.qc_study_list} ${mode} ${params.k} ${params.n_pcs} ${params.threshold_N}
  """
}

// ============================================================================
// WORKFLOW
// ============================================================================
workflow {

    // Build per-chromosome tuples for autosomes 1..22 from input globs
    def ref_ch = Channel.fromPath(params.input_reference, checkIfExists: true)
                        .map { f ->
                          def m = (f.name =~ /chr(\d+)/)
                          if (m) {
                            def c = m[0][1] as Integer
                            if (c >= 1 && c <= 22) tuple(c.toString(), f, file(f.toString() + '.tbi'))
                          }
                        }
                        .filter { it != null }

    def study_ch = Channel.fromPath(params.input_study, checkIfExists: true)
                          .map { f ->
                            def m = (f.name =~ /chr(\d+)/)
                            if (m) {
                              def c = m[0][1] as Integer
                              if (c >= 1 && c <= 22) tuple(c.toString(), f, file(f.toString() + '.tbi'))
                            }
                          }
                          .filter { it != null }

    // Parameter files as channels
    def qc_ref_list_ch   = Channel.fromPath(params.qc_ref_list)
    def qc_study_list_ch = Channel.fromPath(params.qc_study_list)

    // Normalize/QC both panels
    def ref_qc   = qc_norm_ref(ref_ch)
    def study_qc = qc_norm_study(study_ch)

    // Join by chr → intersect → QC+prune
    def ref_study_joined = ref_qc.join(study_qc)
    def intersect_results = intersect(ref_study_joined)
    def pruned_all = final_qc_and_prune(intersect_results)

    // Concatenate all pruned VCFs (collect into a list of paths)
    def concat_results = concat_pruned_vcfs( pruned_all.map { it[1] }.collect() )

    // Double reference IDs
    def ref_ids = double_ref_ids(qc_ref_list_ch)
    
    // Prepare reference (VCF → geno/site → PCA)
    def ( ref_vcf, ref_tbi, ref_geno, ref_site, ref_pca ) = prepare_reference(concat_results, ref_ids)
    
    if (params.run_projection) {
    
        // Double and split study IDs only when projection is requested
        def study_ids = double_study_ids(qc_study_list_ch)
        def study_batches = split_study_list(study_ids).flatMap { it }
    
        // Build combined input tuples for projection
        def proj_input = study_batches
                          .combine(concat_results)
                          .combine(ref_geno)
                          .combine(ref_site)
                          .combine(ref_pca)
                          .map { it -> tuple(it[0], it[1], it[2], it[3], it[4]) }
    
        // Project each batch
        def projected = laser_pca_projection_batch(proj_input)
    
        // Merge reference PCA + projected samples
        def final_pca = merge_proj_coords(ref_pca, projected.collect())
    
        // Plot merged PCA
        run_pca_analysis(final_pca, "projection")
    
    } else {
    
        // Plot reference PCA only
        run_pca_analysis(ref_pca, "reference_only")
    }
}
