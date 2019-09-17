linreg <- setRefClass ("linreg",
  fields = c (
    regression_coefficients = "numeric", 
    fitted_values = "numeric",
    residuals = "numeric",
    degrees_of_freedom = "numeric",
    residual_variance = "numeric",
    variance_of_coefficients = "numeric",
    t_values = "numeric",
    p_values = "numeric"),
  methods = c (
    print = function(...) {
      print(...)
    }
  )
)