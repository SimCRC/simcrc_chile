


l_ANN_settings <- list(
  scale_type      = 1,
  verbose         = 0,
  n_batch_size    = 2000,
  n_epochs        = 10000,
  patience        = 2000,
  n_hidden_nodes  = 360,
  n_hidden_layers = 4,
  activation_fun  = 'tanh',
  init_W          = NULL,
  n_iter          = 100000,
  n_thin          = 100,
  n_chains        = 4,
  Normalize_inputs = FALSE,
  Normalize_outputs= FALSE,
  Scale_inputs    = TRUE,
  Scale_outputs   = TRUE,
  Remove_outliers = FALSE,
  Standardize_targets= FALSE,
  Saved_data      = FALSE,
  Selected_targets= TRUE
)

# Settings for SSP  model
l_ANN_settings_ssp <- list(
  scale_type      = 1,
  verbose         = 0,
  n_batch_size    = 2000,
  n_epochs        = 10000,
  patience        = 2000,
  n_hidden_nodes  = 360,
  n_hidden_layers = 4,
  activation_fun  = 'tanh',
  init_W          = NULL,
  n_iter          = 100000,
  n_thin          = 100,
  n_chains        = 4,
  Normalize_inputs = FALSE,
  Normalize_outputs= FALSE,
  Scale_inputs    = TRUE,
  Scale_outputs   = TRUE,
  Remove_outliers = FALSE,
  Standardize_targets= FALSE,
  Saved_data      = FALSE,
  Selected_targets= TRUE
)
