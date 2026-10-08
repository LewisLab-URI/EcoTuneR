# Builds the browser version of EcoTuneR into docs/ (served by GitHub Pages).
# Run from the repository folder:  Rscript build_web.R
app_dir <- file.path(tempdir(), "ecotuner_app")
unlink(app_dir, recursive = TRUE)
dir.create(file.path(app_dir, "www"), recursive = TRUE)
file.copy(c("app.R", "EcoTuneR.R", file.path("Example Datasets", "Example_combined_output.csv")), app_dir)
keep <- c("ecotuner0.png", "LewisLab.png", "LSU.png", "USM.png", "UCSC.png", "NASEM.png")
file.copy(file.path("www", keep), file.path(app_dir, "www"))
unlink("docs", recursive = TRUE)
shinylive::export(app_dir, "docs")
