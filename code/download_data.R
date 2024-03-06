#### download data from Zenodo ####
# download latest version of the database from Zenodo, https://zenodo.org/doi/10.5281/zenodo.10782428

cat("Downloading raw data\n")
zen4R::download_zenodo("10.5281/zenodo.10782428", path = "data/")

cat("Data successfully downloaded! \n")