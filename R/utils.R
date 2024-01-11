.file_ext <- function (x) {
  pos <- regexpr("\\.([[:alnum:]]+)$", x)
  ifelse(pos > -1L, substring(x, pos + 1L), "")
}

.generate_random_string <- function (
  n = 5
) {
  paste(sample(c(0:9, letters, LETTERS), n, replace = TRUE), collapse = "")
}

.tidy <- function (
  test_result
) {
  stopifnot(inherits(test_result, what = "htest"))
  tibble(estimate = test_result$estimate,
         statistic = test_result$statistic,
         p.value = test_result$p.value,
         conf.low = test_result$conf.int[1],
         conf.high = test_result$conf.int[2],
         method = test_result$method,
         alternative = test_result$alternative)
}
