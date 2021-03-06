# snakemake -n
# snakemake --cores 1

rule all:
  input:
    'data/expression_data',
    'data/metadata',
    'data/circadian_human_blood.qs',
    'data/circadian_human_blood_emat.qs'

rule process_data:
  input:
    'code/process_data.R',
    'data/expression_data',
    'data/metadata'
  output:
    'data/subj_norm_esetList.qs',
    'data/circadian_human_blood.qs',
    'data/circadian_human_blood_emat.qs'
  shell:
    'Rscript {input[0]}'

rule load_perturbations:
  input:
    'code/load_perturbations.R',
    'data/expression_data',
    'data/metadata',
  output:
    'data/subj_norm_pert_esetList.qs',
    'data/perturb_esetList.qs',
    'data/perturb_emat.qs'
  shell:
    'Rscript {input[0]}'

rule analyze_cv:
  input:
    'code/analyze_cv.R',
    'data/metadata',
    'data/circadian_human_blood_emat.qs',
    'data/genes2017.qs'
  output:
    'output/zeitzeiger_cv.pdf',
    'data/zeitzeiger_coefs.qs',
    'output/gene_zeitzeiger_coefs.pdf',
    'output/glmnet_cv.pdf',
    'output/bloodCCD_cv.pdf',
    'output/fig1.pdf',
    'data/glmnet_coefs.qs',
    'output/gene_glmnet_coefs.pdf',
    'output/suppFig1.pdf',
    'output/bloodCCD_coefs.pdf'
  shell:
    'Rscript {input[0]}'

rule analyze_cor:
  input:
    'code/analyze_cor.R',
    'data/metadata',
    'data/circadian_human_blood_emat.qs',
    'data/circadian_human_blood.qs',
    'data/subj_norm_pert_esetList.qs',
    'data/zeitzeiger_coefs.qs',
    'data/glmnet_coefs.qs',
    'data/genes2017.qs',
  output:
    'data/glmnet_cor_dt.qs',
    'output/overall_correlations.pdf',
    'output/study_cors_glmnet.pdf',
    'output/fig3.pdf',
    'output/glmnet_study_cond_corr.pdf',
    'output/study_correlations.pdf',
    'output/gene_venn.pdf',
    'output/suppFig2.pdf',
    'output/ccd_plot.pdf',
    'output/fig5.pdf',
    'output/suppFig3.pdf',
    'output/fig2.pdf',
    'data/result_blood_ref.qs',
    'data/genes_blood.csv'
  shell:
    'Rscript {input[0]}'

rule analyze_cell:
  input:
    'code/analyze_cell.R',
    'data/glmnet_cor_dt.qs',
    'data/rna_blood_cell_monaco.tsv.gz',
    'data/rna_blood_cell_schmiedel.tsv.gz'
  output:
    'output/cell_heatmap.pdf',
    'output/fig4.pdf'
  shell:
    'Rscript {input[0]}'
