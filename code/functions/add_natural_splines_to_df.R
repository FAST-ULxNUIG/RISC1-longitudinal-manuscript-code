add_natural_splines_to_df <- function(df, spline_object) {
  # add natural spline basis to df, based on existing spline object.
  # check inputs:
  if(!(is.data.frame(df))) stop("df must be a data frame.")
  if(!("time" %in% names(df))) stop("df requires a 'time' variable.")
  if(!(is.numeric(df$time))) stop("time must be numeric.")
  if(!("ns" %in% class(spline_object))) stop("spline_object must be of class ns, created by the function ns()")
  # extract inputs
  spline_degree <- max(attr(spline_object, "degree"))
  knots <- attr(spline_object, "knots")
  Boundary.knots <- attr(spline_object, "Boundary.knots")
  intercept <- attr(spline_object, "intercept")
  time <- df$time
  spline_matrix <- ns(x = time, df = spline_degree, knots = knots, intercept = intercept, Boundary.knots = Boundary.knots)
  name_seq <- attr(spline_object, "dimnames")[[2]]
  colnames(spline_matrix) <- paste0("spline_", name_seq)
  # and return original df plus added spline columns:
  cbind(df, spline_matrix)
}