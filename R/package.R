
# Add . and .data to the list of global variables so that
# R CMD check --as-cran does not WARN that . and .data do not exist.
utils::globalVariables(".")
utils::globalVariables(".data")
