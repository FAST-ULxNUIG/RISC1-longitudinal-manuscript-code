split_train_test <- function(N_total, test_prop) {
  if(!(test_prop > 0 & test_prop < 1)) stop("test proportion (test_prop) must be between 0 and 1")
  total_seq <- seq_len(N_total)
  test_size <- round(test_prop * N_total)
  test_inds <- sample(total_seq, size = test_size)
  train_inds <- total_seq[-test_inds]
  list(test_inds = test_inds, train_inds = train_inds)
}
