contrast.ci <- function(model.lm, level = .05){
  ci_margin <- function(x){
    denom <- var(x) * (N-1)
    err <- sqrt(fcrit * msw / denom)
    err
  }
  rss <- resid(model.lm)
  rssDF <- model.lm$df.residual
  msw <- t(rss) %*% rss  / rssDF
  fcrit <- qf(level, 1, rssDF, lower.tail = FALSE)
  model.coef <- coef(model.lm)[-1]
  cont_mat <- model.matrix(model.lm)[,-1]
  N <- dim(cont_mat)[1]
  ci_mat <- matrix(NA, length(model.coef), 2L,
                   dimnames = list(names(model.coef), c("lower", "upper")))
  cont_se <- apply(cont_mat, 2, ci_margin)
  ci_mat[,1] <- model.coef - cont_se
  ci_mat[,2] <- model.coef + cont_se
  ci_mat
}
