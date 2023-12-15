
<!-- README.md is generated from README.Rmd. Please edit that file -->

# sigvar

<!-- badges: start -->
<!-- badges: end -->

The R package *sigvar* implements **sig**nature **var**iability
analysis, a framework for the analysis of mutational signature
activities within and across cancer samples. This R package accompanies
the paper [“Variability of mutational signatures is a footprint of
carcinogens’’ by Morrison et
al.](https://doi.org/10.1101/2023.11.23.23298821); please refer to the
paper for more details on the methods presented in this package.

The *sigvar* package contains two core functions to perform signature
variability analysis:

- `sigvar`: Compute the within-sample diversity and across-sample
  heterogeneity of mutational signature activity in one or multiple
  populations of samples

- `sigboot`: Use bootstrapping to statistically compare the
  within-sample diversity and across-sample heterogeneity of the
  mutational signature activity between two or more groups of samples

*sigvar* also includes accessory functions for the visualization of
mutational signature data, such as:

- `plot_SBS_spectrum`: Plot the SBS mutational spectrum of one or more
  samples of mutational signatures

- `plot_signature_prop`: Plot the relative activities of mutational
  signatures in each sample as a stacked bar plot

- `plot_dots`: Plot the mean mutational signature contributions of one
  or more groups of samples

## Installation

You can install the development version of *sigvar* from
[GitHub](https://github.com/MaikeMorrison/sigvar) with:

``` r
install.packages("devtools") # run only if devtools not already installed
devtools::install_github("MaikeMorrison/sigvar", dependencies = TRUE, build_vignettes = TRUE)
```

The package requires packages `dplyr`, `ggplot2`, `rlang`, `tidyr`, and
`readr`. They will all be installed automatically by the command above.
Installation time ranges from 1 to 5 minutes depending on whether
dependencies also need to be installed. Run time is expected to be a few
minutes on a typical desktop computer.

<!-- See vignettes for detailed examples of using sigvar on real datasets. -->

The package has been tested on R version 4.1.2 on a redhat linux
platform and a windows 10 pro platform. The package is available under
the MIT license.

## Tutorial

A tutorial on the usage of *sigvar* with a focus on the analysis of
compositional data representing microbiome samples is available in the
`tutorial` vignette, which is available at [this
link](https://maikemorrison.com/files/sigvar_tutorial.html) or via the
following code after package installation:

``` r
vignette("tutorial", package = "sigvar")
```

<!-- ## Example -->
<!-- This is a basic example which shows you how to import results from SigProfiler and plot the signature attributions: -->
<!-- ```{r example} -->
<!-- library(sigvar) -->
<!-- SPfolder = system.file("extdata", "example_SigProfiler_results", package = "sigvar") -->
<!-- Qlist = import_SigProfiler(SPfolder) -->
<!-- plot_dots(Qlist[[1]]) -->
<!-- ``` -->
<!-- ## Example analysis -->
<!-- In this example, we apply signature variability analysis (SVA) to data from [Moody et al. (2021)](https://www.nature.com/articles/s41588-021-00928-6). This data set contains 552 esophageal squamous cell carcinoma (ESCC) samples collected across eight countries which vary dramatically in their incidence of ESCC. Moody et al. (2021) reported the activities of 43 SBS, DBS, and ID mutational signatures for each sample. -->
<!-- ```{r} -->
<!-- library(sigvar) -->
<!-- library(dplyr) -->
<!-- library(ggplot2) -->
<!-- head(ESCC_sig_activity) -->
<!-- ``` -->
<!-- Each row of this data set reports the relative activity of each mutational signature for one sample. Because we are analyzing relative signature activities, each row sums to 1. -->
<!-- ### Signature variability analysis can be conducted in one line of code: -->
<!-- ```{r} -->
<!-- sva = sigvar(sig_activity = ESCC_sig_activity, K = 43, group = "Country") -->
<!-- knitr::kable(sva) -->
<!-- ``` -->
<!-- Signature variability analysis quantifies both the mean signature diversity within each sample, as well as the heterogeneity in signature activities across samples. Both variability statistics range between 0 and 1.  -->
<!-- If we add to our SVA results a column corresponding to the incidence of each country, we are able to determine if signature heterogeneity or diversity are associated with ESCC incidence.  -->
<!-- ```{r} -->
<!-- sva_incidence = ESCC_sig_activity %>%  -->
<!--   transmute(Country = as.character(Country),  -->
<!--             Incidence_Level) %>% -->
<!--   distinct %>% -->
<!--   right_join(sva) -->
<!-- knitr::kable(sva_incidence) -->
<!-- ggplot(sva_incidence, -->
<!--        aes(x = mean_within_sample_diversity,  -->
<!--            y = across_sample_heterogeneity,  -->
<!--            color = Incidence_Level)) +  -->
<!--   geom_point(size = 4) + -->
<!--   ggrepel::geom_text_repel(aes(label = Country)) + -->
<!--   theme_bw() -->
<!-- ``` -->
<!-- We see that high-ESCC-incidence countries have more within-sample signature diversity and less across-sample heterogeneity than low-incidence countries. -->
<!-- <!-- ```{r} -->

–\> <!-- <!-- sva_incidence %>%  --> –\>
<!-- <!--   tidyr::pivot_longer(cols = c(across_sample_heterogeneity, mean_within_sample_diversity), -->
–\>
<!-- <!--                       names_to = "Variability_statistic") %>% -->
–\>
<!-- <!-- ggplot(aes(x = Incidence_Level, y = value, fill = Incidence_Level)) +  -->
–\> <!-- <!--   geom_violin(color = NA, alpha = 0.75) + --> –\>
<!-- <!--   ggbeeswarm::geom_beeswarm() + --> –\>
<!-- <!--   ggpubr::stat_compare_means(method = "t.test") +  --> –\>
<!-- <!--   facet_wrap(~ Variability_statistic, scales = "free") + -->
–\> <!-- <!--   theme_bw() --> –\> <!-- <!-- ``` --> –\>

<!-- Note that this difference in diversity is missed when we analyze only country-level mean signature activities: -->
<!-- ```{r, fig.width=8} -->
<!-- plot_dots(sig_activity = ESCC_sig_activity,  -->
<!--           K = 43, group = "Country",  -->
<!--           facet = "Incidence_Level",  -->
<!--           pivot = TRUE) -->
<!-- ``` -->
<!-- We use bootstrapping (`sigboot`) to statistically compare the signature diversity or heterogeneity of cancer samples. For example, the below code compares SVA results among the low-incidence countries. -->
<!-- ```{r} -->
<!-- low_inc_boot = sigboot(sig_activity = ESCC_sig_activity %>%  -->
<!--                          filter(Incidence_Level == "Low"), -->
<!--                        K = 43,  -->
<!--                        group = "Country", -->
<!--                        n_replicates = 100,  -->
<!--                        seed = 1) -->
<!-- ``` -->
<!-- In short, the bootstrapping algorithm involves taking a pair of populations (for example, Japan and Brazil), scrambling their samples up many times (with some samples randomly duplicated and others omitted), conducting SVA on the scrambled populations, and then comparing the difference between SVA results for each scrambled population (black dots) to the true difference (red dots). We can see in the below plot that the red dots fall in the middle of the black dots, suggesting that there is not a significant difference in SVA results between Japan and Brazil.  -->
<!-- ```{r} -->
<!-- low_inc_boot$bootstrap_distribution_plot$`Japan--Brazil` -->
<!-- ``` -->
<!-- This lack of significance is quantified by the two-sided P-value comparing each pair of countries:  -->
<!-- ```{r} -->
<!-- low_inc_boot$P_values %>%  data.frame() %>% select(-pooled_diversity) -->
<!-- ``` -->
<!-- However, suppose we compare a high-incidence country like Kenya to a low-incidence country like Japan: -->
<!-- ```{r} -->
<!-- comparison_boot = sigboot(sig_activity = ESCC_sig_activity %>%  -->
<!--                          filter(Country %in% c("Japan", "Kenya")), -->
<!--                        K = 43,  -->
<!--                        group = "Country", -->
<!--                        n_replicates = 100,  -->
<!--                        seed = 1) -->
<!-- comparison_boot$bootstrap_distribution_plot -->
<!-- ``` -->
<!-- We see that, for the SVA statistics `across_sample_heterogeneity` and `mean_within_sample_diversity`, the red dots fall far outside of the range of most of the black dots, suggesting that there is a significant difference in SVA results between Kenya and Japan. This difference is reflected by low two-sided P-values: -->
<!-- ```{r} -->
<!-- comparison_boot$P_values %>%  data.frame() %>% select(-pooled_diversity) -->
<!-- ``` -->
