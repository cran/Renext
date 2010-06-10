
## To be used later in a summary or coef

coefftable <- function(fit) {
  s <- sqrt(diag(fit$cov))
  res <- cbind(fit$estimate, s, fit$estimate/s)
  colnames(res) <- c("est.", "sd", "t")
  res
}
