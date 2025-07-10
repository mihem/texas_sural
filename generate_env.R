library(rix)

my_token <- gitcreds::gitcreds_get()$password
Sys.setenv(GITHUB_PAT = my_token)

rix(
  date = "2025-07-07",
  r_pkgs = c("languageserver", "dplyr", "Seurat"),
  system_pkgs = NULL,
  git_pkgs = list(
    package_name = "httpgd",
    repo_url = "https://github.com/nx10/httpgd",
    commit = "dd6ed3a687a2d7327bb28ca46725a0a203eb2a19"
  ),
  ide = "none",
  project_path = ".",
  overwrite = TRUE,
  print = TRUE
)
