# Runs the whole analysis (Tasks 1-3) in order.
# Run from the project folder:  Rscript run_all.R
# Regenerates every figure in figures/ and writes the console output to outputs/console_output.txt

dir.create("outputs", showWarnings = FALSE)
log_file <- file("outputs/console_output.txt", open = "wt")
sink(log_file)
sink(log_file, type = "message")

options(width = 150)
for (f in c("scripts/task1_data_exploration.R",
            "scripts/task2_regression.R",
            "scripts/task3_abc.R")) {
  cat("\n#########", f, "\n")
  source(f, echo = FALSE, print.eval = TRUE)
}

sink(type = "message")
sink()
close(log_file)
cat("Done. Figures in figures/, console output in outputs/console_output.txt\n")
