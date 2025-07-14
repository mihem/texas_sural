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
    "qs",
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
    ),
    list(
      package_name = "BPCells",
      repo_url = "https://github.com/bnprks/BPCells/r",
      commit = "9d2a036af9128d34936ad08a43e60a4e4916049c"
    ),
    list(
      package_name = "STACAS",
      repo_url = "https://github.com/carmonalab/STACAS",
      commit = "1fde4ef460643406683bfcd40909a33264ce3a88"
    )
  ),
  ide = "none",
  project_path = ".",
  overwrite = TRUE,
  print = TRUE
)

# --- BPCells block fix ---
nix_file <- "default.nix"
nix_lines <- readLines(nix_file)

bp_start <- grep("BPCells = \\(pkgs.rPackages.buildRPackage \\{", nix_lines)
bp_end <- grep("^\\s*\\}\\);", nix_lines)
bp_end <- bp_end[bp_end > bp_start][1]

if (length(bp_start) == 1 && !is.na(bp_end)) {
  # Extract repo URL, rev, sha256 from BPCells block
  url_line <- grep("url = ", nix_lines[bp_start:bp_end], fixed = TRUE) +
    bp_start -
    1
  rev_line <- grep("rev = ", nix_lines[bp_start:bp_end], fixed = TRUE) +
    bp_start -
    1
  sha_line <- grep("sha256 = ", nix_lines[bp_start:bp_end], fixed = TRUE) +
    bp_start -
    1

  url <- sub(".*url = \"(.*)\";.*", "\\1", nix_lines[url_line])
  url_src <- sub("/r$", "", url) # Remove trailing /r for BPCells-src
  rev <- sub(".*rev = \"(.*)\";.*", "\\1", nix_lines[rev_line])
  sha <- sub(".*sha256 = \"(.*)\";.*", "\\1", nix_lines[sha_line])

  # Create complete replacement block
  replacement_block <- c(
    "    BPCells-src = pkgs.fetchgit {",
    sprintf("      url = \"%s\";", url_src),
    sprintf("      rev = \"%s\";", rev),
    sprintf("      sha256 = \"%s\";", sha),
    "    };",
    "",
    "    BPCells = (pkgs.rPackages.buildRPackage {",
    "      name = \"BPCells\";",
    "      src = \"${BPCells-src}/r\";",
    "      postPatch = \"patchShebangs configure\";",
    "      nativeBuildInputs = [ pkgs.hdf5.dev ];"
  )

  propagatedBuildInputs <- grep(
    "propagatedBuildInputs = ",
    nix_lines[bp_start:bp_end],
    fixed = TRUE
  ) +
    bp_start -
    1

  # Replace BPCells block with corrected version
  nix_lines <- c(
    nix_lines[1:(bp_start - 1)],
    replacement_block,
    nix_lines[propagatedBuildInputs:bp_end],
    nix_lines[(bp_end + 1):length(nix_lines)]
  )

  writeLines(nix_lines, nix_file)
}

cat(
  "###################################################################################################\n"
)
cat(
  "###################################################################################################\n"
)
cat(
  "###################################################################################################\n"
)
cat("Updated nix files successfully:\n")
print(nix_lines)
