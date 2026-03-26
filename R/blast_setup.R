#' Download and Install BLAST+ Binaries
#' @export
install_blast <- function(packagge_name) {
  # 1. Define storage directory
  pkg_name <- "blasting" # Change this to your actual package name
  dest_dir <- tools::R_user_dir(pkg_name, which = "data")
  if (!dir.exists(dest_dir)) dir.create(dest_dir, recursive = TRUE)

  # 2. Identify the Operating System
  os <- Sys.info()["sysname"]

  # 3. Use the exact filename found in the LATEST folder (2.17.0)
  # Note the 'x64' in the name; this is the portable version
  file_name <- switch(os,
                      "Windows" = "ncbi-blast-2.17.0+-x64-win64.tar.gz",
                      "Darwin"  = "ncbi-blast-2.17.0+-x64-macosx.tar.gz",
                      "Linux"   = "ncbi-blast-2.17.0+-x64-linux.tar.gz",
                      stop("Operating system not supported.")
  )

  # 4. Construct the URL (ensure the slash is there!)
  full_url <- paste0("https://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/LATEST/", file_name)
  file_path <- file.path(dest_dir, file_name)

  message("Downloading BLAST+ from: ", full_url)

  # 5. Download and Extract
  tryCatch({
    # 'wb' is critical for Windows to avoid file corruption
    download.file(url = full_url, destfile = file_path, mode = "wb")

    message("Extracting files...")
    untar(file_path, exdir = dest_dir)

    # Clean up the downloaded archive
    unlink(file_path)

    message("\nSUCCESS: BLAST+ 2.17.0 installed in: ", dest_dir)
  }, error = function(e) {
    stop("\nInstallation failed. Error: ", e$message)
  })

  return(dest_dir)
}


#' Run BLASTp to find protein function
#' @param protein_seq A string of amino acids
#' @param db_path Path to a BLAST database (e.g., Swiss-Prot)
#' @export
run_protein_blast <- function(protein_seq, db_path) {
  # 1. Find the path where BLAST was installed
  bin_dir <- tools::R_user_dir("blasting", which = "data")

  # 2. Find the actual blastp.exe file
  # We assign it to the variable 'blastp_exe' here
  blastp_exe <- list.files(bin_dir, pattern = "blastp\\.exe$", recursive = TRUE, full.names = TRUE)

  if (length(blastp_exe) == 0) {
    stop("blastp.exe not found. Please run install_blast() first.")
  }

  # 3. Create temporary files for the query and the output
  query_file <- tempfile(fileext = ".fasta")
  writeLines(c(">query", protein_seq), query_file)
  out_file <- tempfile(fileext = ".txt")

  # 4. Run the command using shQuote to prevent space/path errors
  # Note: we use blastp_exe[1] just in case it finds multiple versions
  system2(blastp_exe[1],
          args = c("-query", shQuote(query_file),
                   "-db", shQuote(db_path),
                   "-out", shQuote(out_file),
                   "-outfmt", "6"))

  # 5. Read the results
  # Use tryCatch in case there are no matches (empty file)
  results <- tryCatch({
    read.table(out_file, sep = "\t", stringsAsFactors = FALSE)
  }, error = function(e) {
    message("No matches found in the database.")
    return(data.frame())
  })

  if (nrow(results) > 0) {
    colnames(results) <- c("query", "subject", "identity", "length", "mismatch",
                           "gap", "qstart", "qend", "sstart", "send", "eval", "bitscore")
  }

  return(results)
}


#' Download and Index Swiss-Prot Database
#' @export
setup_swissprot <- function() {
  db_dir <- tools::R_user_dir("blasting", which = "data")
  if (!dir.exists(db_dir)) dir.create(db_dir, recursive = TRUE)

  # Full FTP path to the Swiss-Prot FASTA file
  fasta_url <- "https://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/complete/uniprot_sprot.fasta.gz"
  dest_file <- file.path(db_dir, "swiss.fasta.gz")

  # Ensure a long timeout for the ~90MB download
  old_opts <- options(timeout = 1200)
  on.exit(options(old_opts))

  message("Attempting to download Swiss-Prot from UniProt FTP...")

  # Try downloading using the ftp protocol
  tryCatch({
    download.file(url = fasta_url, destfile = dest_file, mode = "wb")

    # Check if we got the actual file (approx 80-90MB)
    if (file.size(dest_file) < 50000000) {
      stop("Downloaded file is too small. It might be a connection error.")
    }

    message("Extracting...")
    # Requires library(R.utils)
    R.utils::gunzip(dest_file, overwrite = TRUE, remove = TRUE)

    message("Indexing database (this takes a few minutes)...")
    makeblastdb_exe <- list.files(db_dir, pattern = "makeblastdb\\.exe$", recursive = TRUE, full.names = TRUE)

    system2(makeblastdb_exe,
            args = c("-in", shQuote(file.path(db_dir, "swiss.fasta")),
                     "-dbtype", "prot",
                     "-out", shQuote(file.path(db_dir, "swissprot"))))

    message("SUCCESS: Swiss-Prot database is ready for use.")

  }, error = function(e) {
    stop("Failed to download or index Swiss-Prot: ", e$message)
  })
}


#' Identify Protein Function from Sequence
#' @export
identify_protein <- function(protein_seq) {

  # 1. Check for BLAST+
  if (!is_blast_installed()) {
    ans <- readline("BLAST+ is required but not found. Download it now? (y/n): ")
    if (tolower(ans) != "y") stop("Cannot proceed without BLAST+.")
    install_blast()
  }

  # 2. Check for Swiss-Prot Database
  db_dir <- tools::R_user_dir("blasting", which = "data")
  db_path <- file.path(db_dir, "swissprot")

  if (!file.exists(paste0(db_path, ".pin"))) {
    ans <- readline("Swiss-Prot database (90MB) is missing. Download it now? (y/n): ")
    if (tolower(ans) != "y") stop("Cannot proceed without the database.")
    setup_swissprot()
  }

  # 3. Run the search
  message("Searching...")
  return(run_protein_blast(protein_seq, db_path = db_path))
}

#' Check if BLAST+ is installed
#' @export
is_blast_installed <- function() {
  # Get the path where your package stores data
  bin_dir <- tools::R_user_dir("blasting", which = "data")

  # Look for the executable (blastp.exe on Windows, blastp on Mac/Linux)
  exe_name <- if (.Platform$OS.type == "windows") "blastp.exe" else "blastp"
  exe_path <- list.files(bin_dir, pattern = exe_name, recursive = TRUE, full.names = TRUE)

  return(length(exe_path) > 0)
}
