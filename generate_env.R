library(rix)

my_token <- gitcreds::gitcreds_get()$password
Sys.setenv(GITHUB_PAT = my_token)

rix(
  date = "2025-07-07",
  r_pkgs = c(
    "languageserver",
    "tidyverse",
    "Seurat",
    "readxl",
    "writexl",
    "qs2",
    "targets",
    "tarchetypes",
    "visNetwork",
    "hdf5r",
    "scDblFinder"
  ),
  system_pkgs = NULL,
  git_pkgs = list(
    list(
      package_name = "httpgd",
      repo_url = "https://github.com/nx10/httpgd",
      commit = "dd6ed3a687a2d7327bb28ca46725a0a203eb2a19"
    ),
    list(
      package_name = "scMisc",
      repo_url = "https://github.com/mihem/scMisc",
      commit = "e2ebddcb779b935551f14216514c0429616fc91d"
    ),
    list(
      package_name = "presto",
      repo_url = "https://github.com/immunogenomics/presto",
      commit = "7636b3d0465c468c35853f82f1717d3a64b3c8f6"
    )
  ),
  ide = "none",
  project_path = ".",
  overwrite = TRUE,
  print = TRUE
)
