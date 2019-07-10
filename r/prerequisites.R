
# Required packages -----------------------------------------------------------------------------------------------
suppressPackageStartupMessages({
    library(plyr)
    library(tidyverse)
    library(patchwork)
    library(ggsignif)
    library(binom)
    library(stringi)
    library(ggrepel)
    library(grid)
    library(readxl)
    library(ggbeeswarm)
    library(here)
    library(survival)
    library(survminer)
})

# Variables, etc. -------------------------------------------------------------------------------------------------
brca_associated = c('Breast Cancer', 'Breast cancer', 'Breast',
                    'Prostate Cancer', 'Prostate cancer', 'Prostate',
                    'Pancreatic Cancer', 'Pancreatic cancer', 'Pancreatic',
                    'Ovarian Cancer', 'Ovarian cancer', 'Ovarian')

# Convenience functions -------------------------------------------------------------------------------------------
as_perc = function(x) { 100*x }
as_pts = function(x, input_unit = 'mm') { as.numeric(convertUnit(unit(x, 'pt'), input_unit)) }
ci_upper = function(n, N) { binom.confint(n, N, methods = 'wilson')[['upper']] }
ci_lower = function(n, N) { binom.confint(n, N, methods = 'wilson')[['lower']] }
mcn_test = function(a, b, c, total) {
    mcnemar.test(matrix(c(a, b, c, total - (a + b + c)), ncol = 2))[['p.value']]
}
f_test = function(a, b, c, d) {
    fisher.test(matrix(c(a, b, c, d), ncol = 2))[['p.value']]
}

# Load data -------------------------------------------------------------------------------------------------------
si_tables = here('data/Taylor_SItables.xlsx')

brca_tbl = read_xlsx(si_tables, sheet = 'Table S3. BRCA mutation rates', skip = 1, n_max = 56)
names(brca_tbl) = c('cancer_type', 'total',
                    paste0('brca1_', c('germline', 'somatic_lof', 'somatic_del', 'somatic_lof_hm', 'somatic_vus')),
                    paste0('brca2_', c('germline', 'somatic_lof', 'somatic_del', 'somatic_lof_hm', 'somatic_vus')),
                    paste0('germline_', c('total_evaluable', 'loh', 'second_somatic')),
                    paste0('somatic_lof_', c('total_evaluable', 'loh')),
                    paste0('somatic_lof_hm_', c('total_evaluable', 'loh'))) 

brca_tbl = mutate_at(brca_tbl, vars(matches('total|^brca1|^brca2|^germline|^somatic')),
                     funs(as.integer(str_extract(., '^[0-9]+')))) %>% 
    mutate(total_brca_germline = brca1_germline + brca2_germline,
           total_brca_somatic = brca1_somatic_lof + brca2_somatic_lof +
                                brca1_somatic_del + brca2_somatic_del,
           total_brca = total_brca_germline + total_brca_somatic)

clinical_tbl = read_xlsx(si_tables, sheet = 'Table S2. Clinical') 
names(clinical_tbl) = c('patient', 'sample', 'sex', 'cancer_type', 'detailed_cancer_type',
                        'dx_age_first', 'dx_age_breast', 'dx_age_ovarian', 'dx_age_prostate', 'dx_age_pancreatic',
                        'dx_age_other', 'germline_hrd', 'germline_other',
                        'immuno_tx', 'immuno_tx_efs_months', 'immuno_tx_efs_status',
                        'parpi_tx_months', 'parpi_tx_status')

clinical_tbl = mutate(clinical_tbl, parpi_tx_months = as.numeric(parpi_tx_months),
                                    parpi_tx_status = as.integer(parpi_tx_status),
                                    immuno_tx_efs_months = as.numeric(immuno_tx_efs_months),
                                    immuno_tx_efs_status = as.integer(immuno_tx_efs_status))

wes_tbl = read_xlsx(si_tables, sheet = 'Table S6. WES data') 
names(wes_tbl) = c('patient', 'sample', 'sex', 'cancer_type', 'detailed_cancer_type',
                   'purity', 'ploidy', 'genome_doubled', 'hypermutated',
                   'germline_brca1', 'germline_brca2', 'somatic_brca1', 'somatic_brca2', 'somatic_brca_zygosity',
                   'germline_hrd', 'somatic_hrd',
                   'lst', 'hrd_loh', 'ntai', 'signature3', 'hrd_score', 'other_signatures')

wes_tbl = mutate(wes_tbl, germline = germline_brca1 == 1 | germline_brca2 == 1,
                          somatic = somatic_brca1 == 1 | somatic_brca2 == 1,
                          hrd_wt = germline == FALSE & somatic == FALSE &
                          germline_hrd == FALSE & somatic_hrd == FALSE & cancer_type != 'Ovarian cancer',
                          brca_associated = cancer_type %in% brca_associated,
                          hrd_score = as.numeric(hrd_score))

file_path1 = here('data/germline_mutations.maf') # manually enter path here
germline_brca_mutations = read_tsv(file_path1) %>% 
    filter(Hugo_Symbol %in% c('BRCA1', 'BRCA2')) %>% 
    select(patient, Hugo_Symbol)

file_path2 = here('data/germline_cnvs.txt') # manually enter path here

germline_brca_cnvs = read_tsv(file_path2)
names(germline_brca_cnvs) = c('patient', 'gene', 'exons', 'penetrance')
germline_brca_cnvs = filter(germline_brca_cnvs, gene %in% c('BRCA1', 'BRCA2'))

homdels = c("A-d2dd72ef5d9d-T01-IM5", "A-3f040d05cbee165-T02-IM6", "A-98cfcf2f0474e48-T01-IM3", "A-918cc83af322-T01-IM5", "A-36940743ba06-T01-IM6", "A-67def83769cf-T01-IM6", "A-eb2a875f1112-T01-IM3", "A-627c419aaf84914-T01-IM5", "A-ca63c5835392-T01-IM5", "A-9d92ae09fe14-T01-IM3", "A-54389a7ee98c2b9-T01-IM5", "A-a11a737243ac-T01-IM5", "A-2498b941f5ee-T01-IM5", "A-5146ec45e314eb5-T01-IM5", "A-2a263920e463-T01-IM5", "A-ad022433deb8-T01-IM3", "A-8e180655ec59-T01-IM6", "A-4249fde6d47b-T01-IM5", "A-1138840cf921-T01-IM5", "A-68b2356782a1-T01-IM5", "A-6d60606c4dd2-T01-IM6", "A-cab735cc920a570-T01-IM5", "A-d5fc353383b7564-T01-IM5", "A-2bbcb2d9b4ab-T01-IM3", "A-962aa5aa2b6901c-T01-IM6", "A-999b8d8fd837-T01-IM3", "A-afcfbbb78924fa4-T01-IM6", "A-ce906800a9e3-T01-IM5", "A-a1ef3971e9104dd-T01-IM6", "A-2948c162d428-T02-IM5", "A-79ec2b5e9dbb-T01-IM6", "A-b4b7cc57e2bc5cf-T01-IM3", "A-68cd95b74d4f-T01-IM3", "A-fc680f2ac226-T01-IM5", "A-acb632a71886-T01-IM5", "A-dc704de4813b69c-T01-IM5", "A-fae80f8f1373ccb-T01-IM6", "A-d401517f1107f16-T01-IM6", "A-2836c9de828a06b-T01-IM6", "A-3408a07474e1-T01-IM3", "A-65294968558e-T01-IM6", "A-1b95c0df7d81-T01-IM5", "A-9d8d0d601fc0f26-T01-IM5", "A-f225b6e98841-T01-IM3", "A-894eeafc547c0ff-T01-IM5", "A-e7108324b7e6-T01-IM3", "A-5f89a59b3477-T01-IM5", "A-dca4d80e9eec-T01-IM6", "A-1c88f3d1a5a2-T01-IM3", "A-156f6abe3f84-T02-IM6", "A-85bb17dd78bab87-T01-IM5", "A-c0bb079831609a3-T01-IM5", "A-3fa4c8979a99-T01-IM5", "A-ed4255ed3b11fef-T01-IM5", "A-40289a6b84f4144-T01-IM5", "A-8b7e97d889fa-T01-IM3", "A-934a9e9e39c1-T01-IM5", "A-7ca07d381032692-T01-IM6", "A-9d80442f2428-T01-IM5", "A-d84a1b7c8380-T01-IM5", "A-d51eca1d4e66-T01-IM5", "A-c8eef46b7e3d-T01-IM5", "A-54ff9f48d591-T01-IM6", "A-43c4bff754c7374-T01-IM3", "A-851f0c19958f-T02-IM6", "A-db0d57b3650a4b2-T02-IM6", "A-dfafa618071c-T02-IM5", "A-1be65363f71f1b5-T02-IM3", "A-cc6fd81d89bf-T01-IM6", "A-40118574ec2d-T01-IM6", "A-713cad65f5fa-T01-IM3", "A-c101ffe2320b-T01-IM3", "A-0b8bc8007ea7-T01-IM6", "A-830066eeac1a-T01-IM5", "A-5bcb0b2b3bdf-T03-IM6")

somatic_brca_mutations = read_xlsx(si_tables, sheet = 'Table S4. Somatic BRCA') 
names(somatic_brca_mutations) = c('patient', 'sample', 'gene', 'mutation', 'category', 'sample_mutation_class')

parpi_samples = readRDS(here('data/parpi-samples.rds'))

parpi_tx = filter(clinical_tbl, parpi_tx_status %in% c(1, 0)) %>% 
    mutate(brca_associated = cancer_type %in% brca_associated) %>% 
    select(-sample) %>% 
    left_join(., parpi_samples) %>% 
    mutate(germline_somatic_class = case_when(
        patient %in% c(germline_brca_mutations$patient, germline_brca_cnvs$patient) ~ 'Germline',
        sample %in% somatic_brca_mutations$sample[which(somatic_brca_mutations$category == 'LoF')] ~ 'Somatic',
        sample %in% somatic_brca_mutations$sample[which(somatic_brca_mutations$mutation == 'Deletion')] ~ 'Somatic',
        sample %in% homdels ~ 'Somatic',
        TRUE ~ 'WT'),
        brca_class = case_when(
            germline_somatic_class %in% c('Germline', 'Somatic') ~ 'Germline/Somatic',
            TRUE ~ 'WT')
    )

icb_patients = readRDS(here('data/icb-patients.rds'))

icb_tx = filter(clinical_tbl,
                immuno_tx %in% c('CTLA4', 'PD(L)1', 'PD(L)1/CTLA4'),
                patient %in% icb_patients$patient) %>%
    mutate(brca_class= case_when(
        patient %in% c(germline_brca_mutations$patient, germline_brca_cnvs$patient) ~ 'BRCA1/2',
        sample %in% somatic_brca_mutations$sample[which(somatic_brca_mutations$category == 'LoF' &
                                                            somatic_brca_mutations$sample_mutation_class == 'MSS')] ~ 'BRCA1/2',
        sample %in% somatic_brca_mutations$sample[which(somatic_brca_mutations$mutation == 'Deletion')] ~ 'BRCA1/2',
        sample %in% homdels ~ 'BRCA1/2',
        TRUE ~ 'WT')
    )

# Set plotting theme ----------------------------------------------------------------------------------------------
theme_set(theme_bw())
theme_update(
    text = element_text(family = 'ArialMT', color = 'black', size = 12),
    axis.text = element_text(family = 'ArialMT', color = 'black', size = 12),
    panel.grid.minor = element_blank(),
    panel.grid.major = element_blank(),
    legend.background = element_blank(),
    legend.key = element_blank(),
)

sign_colors = c('BRCA' = 'violetred1',
                'APOBEC' = 'dodgerblue3',
                'MMR/MSI' = 'darkseagreen3',
                'Smoking' = 'royalblue4',
                'UV' = 'gold2',
                '18' = 'brown4',
                'Other' = 'grey75')
