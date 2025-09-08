
packages <- c(
  "tidyr", "gee", "dplyr", "data.table", "utils",
  "caTools", "caret", "GMDH2", "correlation",
  "fastDummies", "pROC"
)

installed_packages <- rownames(installed.packages())
for (p in packages) {
  if (!(p %in% installed_packages)) {
    install.packages(p)
  }
}

lapply(packages, library, character.only = TRUE)
