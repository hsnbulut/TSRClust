source('R/library_environment.R')
initialize_revision_library()
source('config/paths.R')
legacy <- new.env(parent = globalenv())
for (f in c('core/utils.R','core/core_tsr.R','core/data_generation.R',
            'evaluation/evaluation_metrics.R','competitors/competitors.R',
            'competitors/competitors_v2.R','config/config.R','config_v2.R')) {
  source(file.path(SOURCE_ROOT,'reproducibility',f), local=legacy)
}
# Immutable source functions remain in legacy; revision functions use separate names.
for (f in c('common.R','generation.R','validation.R','metrics.R','fitting.R')) source(file.path('R',f))
