# # Make an example matrix of compositional data
# Each row is an individual. Rows sum to 1.
population_A <- matrix(
  c(
    .5, .3, .2,
    .4, .2, .4,
    .5, .4, .1,
    .6, .1, .3,
    .2, 0, .8
  ),
  nrow = 5,
  byrow = TRUE
)

population_B <- matrix(
  c(
    .9, 0, .1,
    .6, .4, 0,
    .7, 0, .3,
    .3, .4, .3,
    .5, .3, .2
  ),
  nrow = 5,
  byrow = TRUE
)


populations_AB <- cbind(
  data.frame(c(
    "A", "A", "A", "A", "A",
    "B", "B", "B", "B", "B"
  )),
  rbind(population_A, population_B)
)
colnames(populations_AB) <- c("population", "category_1", "category_2", "category_3")

test_that("plot_signature_prop works", {
  expect_no_warning( plot_signature_prop(
    relab_matrix = population_A,
    K = 3, # How many categories per vector?
    arrange = FALSE
  ))
  expect_no_warning(plot_signature_prop(
    relab_matrix = population_A,
    K = 3, # How many categories per vector?
    arrange = "horizontal"
  ))
  expect_no_warning(plot_signature_prop(
    relab_matrix = population_A,
    K = 3, # How many categories per vector?
    arrange = "vertical"
  ))
  expect_no_warning(plot_signature_prop(
    relab_matrix = population_A,
    K = 3, # How many categories per vector?
    arrange = TRUE # could also be "both"
  ))
  
  # You can modify the plot as you would any ggplot2 object
  expect_no_warning(plot_signature_prop(
    relab_matrix = population_A,
    K = 3, # How many categories per vector?
    arrange = TRUE) +
      # Below are example, optional modifications to the default plot
      ggplot2::ggtitle("Population A") +
      ggplot2::scale_fill_brewer("Blues") +
      ggplot2::scale_color_brewer("Blues") +
      ggplot2::xlab("Individuals"))
  
  
  expect_no_warning(plot_signature_prop(relab_matrix = populations_AB, group = "population"))
  expect_no_warning(plot_signature_prop(relab_matrix = populations_AB, group = "population", arrange = "vertical"))
  expect_no_warning(plot_signature_prop(relab_matrix = populations_AB, group = "population", arrange = "horizontal"))
  expect_no_warning(plot_signature_prop(relab_matrix = populations_AB, group = "population", arrange = "both"))
  
  
  expect_error(plot_signature_prop(relab_matrix = populations_AB, group = "population", arrange = "2"))
  expect_error(plot_signature_prop(relab_matrix = populations_AB, group = 1, arrange = "both"))
})
