Code for the paper ‘A Multivariate Longitudinal Functional Model for
Repeatedly-Observed Human-Movement Data’
================

## Repository Structure:

------------------------------------------------------------------------

- :open_file_folder: **code**
  - :open_file_folder: **analysis** – scripts used to perform the data
    analysis.
    - :page_facing_up: [01 - Create Plot for
      Introduction](code/analysis/01-introduction-plot.R)
    - :page_facing_up: [02 - Create Plot of Stride
      Timings](code/analysis/02-strides-timing-plot.R)
    - :page_facing_up: [03 - Prepare Data for
      Analysis](code/analysis/03-data-preparation.R)
    - :page_facing_up: [04 - Split Data into Test and Training
      Samples](code/analysis/04-test-train-split.R)
    - :page_facing_up: [05 - mv-FPCA basis
      Representation](code/analysis/05-basis-representation.R)
  - :open_file_folder: **functions** – custom functions for modelling
    data analysis.
    - :page_facing_up: [Centering a (multivariate) object around a
      different mean](code/functions/center_fd_around_new_mean.R)
    - :page_facing_up: [Decentering a (multivariate) object around a
      different mean](code/functions/decenter_fd_around_new_mean.R)
    - :page_facing_up: [Helper functions for manipulating `fd`
      objects](code/functions/functions-helper-smoothing.R)
    - :page_facing_up: [Project (multivariate) functional data (`fd`
      object) onto (multivariate) FPCs (`pca.fd`
      object)](code/functions/project_data_onto_fpcs.R)  
    - :page_facing_up: [Computing the (%) of Variance Explained by a
      mv-FPCA
      reconstruction](code/functions/variance_explained_reconstruction.R)
    - :page_facing_up: [Custom theme for
      figures](code/functions/theme_gunning.R)
    - :page_facing_up: [Functions to generate data from a multilevel
      longitudinal design for the
      simulation](code/functions/generate_design.R)
    - :page_facing_up: [Generate data from a polynomial scalar
      longitudinal model for the
      simulation](code/functions/generate_polynomial_model_basis_coefficient.R)
    - :page_facing_up: [Generate multiple basis coefficients (i.e.,
      mv-FPC scores) from a polynomial scalar longitudinal model for the
      simulation](code/functions/generate-basis-coefficient-matrix.R)
    - :open_file_folder: **tests** – some basic tests for the custom
      functions.
      - :page_facing_up: [Test for
        `center_fd_around_new_mean()`](code/functions/tests/test-center_fd_around_new_mean.R)
      - :page_facing_up: [Test for
        `decenter_fd_around_new_mean()`](code/functions/tests/test-decenter_fd_around_new_mean.R)
      - :page_facing_up: [Test for
        `variance_explained_reconstruction()`](code/functions/tests/test-variance-explained-reconstruction.R)
      - :page_facing_up: [Test for
        `generate_design_multiple_subjects()`](code/functions/tests/test-generate-design.R)
      - :page_facing_up: [Test for
        `generate_polynomial_model_basis_coefficient()`](code/functions/tests/test-generate_polynomial_model_basis_coefficient.R)
      - :page_facing_up: [Test 1 for
        `generate_basis_coefficient_matrix()`](code/functions/tests/test-generate-basis-coefficient-matrix-01.R)
- :open_file_folder: **outputs**
  - :open_file_folder: **tables** – tables containing data-analysis
    results. Some are stored as `.csv` files while others have been
    exported to $TeX$ using `{xtable}` for inclusion in the paper.
  - :open_file_folder: **figures** – figures for the manuscript, all
    created in $TeX$ using `{tikzDevice}` and linked to overleaf, where
    the file [figures.tex](outputs/figures/figures.tex) compiles all the
    individual $TeX$ files.
  - 💾 also contains `.rds` objects saved at various stages of the data
    analysis
- :open_file_folder: **data** – contains the main dataset used in
  analysis and small data sets used to create the introduction and
  stride-timing plots.
