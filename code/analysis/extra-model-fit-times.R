model_fit_results <- readRDS(
  here::here("outputs", "model-fit-results.rds"))

model_fit_results$spline_ri_model$lme_fit_list
model_times_used <- model_fit_results$times[c("spline_ri_time", "fpca_time", "naive_time")]

model_times_extracted <- lapply(model_times_used, function(x) {
  round(x[["elapsed"]] / 60, 2)
})

times_df <- data.frame(
  model = names(model_times_extracted),
  time = unlist(model_times_extracted)
)

saveRDS(object = times_df,
        "outputs/tables/model-fit-times.rds")
