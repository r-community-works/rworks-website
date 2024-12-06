
# Note:
# On OSX may need to add following to ~/.Renviron 
# PKG_CXXFLAGS=-I /Applications/Xcode.app/Contents/Developer/Platforms/MacOSX.platform/Developer/SDKs/MacOSX.sdk/usr/include/c++/v1



run_cached <- function(f, arg_list, prefix='') {
  assumed_fn_name <- as.character(substitute(f))  # really name of variable in calling scope
  key = digest(
    paste(c(assumed_fn_name, length(arg_list), as.character(names(arg_list)), as.character(arg_list)), collapse="\n"),
    algo="md5")
  path = paste0('cache_', prefix, '_', key, '.RDS')
  if(file.exists(path)) {
    return(readRDS(path))
  }
  result = do.call(f, arg_list)
  saveRDS(result, path)
  return(result)
}


# code to replicate code per-IDX
idx_blocks <- function(txt, n_studies) {
  paste0(vapply(
    seq(n_studies),
    function(i) {
      gsub("{IDX}",
           as.character(i),
           txt
           ,
           fixed = TRUE)
    },
    character(1)
  ), collapse="\n")
}


define_Stan_model <- function(n_studies, c_replacement = "") {
# the Stan source code for combining multiple studies as a meta-analysis
src <- paste0("
data {
  int<lower=1> n_studies;  // number of studies
  array[n_studies] int<lower=1> nE;  // number treated examples
  vector[n_studies] meanE;  // mean observed treatment effect
  vector<lower=0>[n_studies] varE;  // observed treatment effect variance
  array[n_studies] int<lower=1> nC;  // number control examples
  vector[n_studies] meanC;  // mean observed control effect
  vector<lower=0>[n_studies] varC;  // observed control effect variance
}
parameters {
  {//} real inferred_grand_treatment_mean;  // unobserved expected treatment
  {//} real inferred_grand_control_mean;  // unobserved expected control effect
  {//} real<lower=0> inferred_between_group_stddev; // standard distance between groups
  vector[n_studies] inferred_group_treatment_mean;  // unobserved per-group ideal treatment effect
  vector[n_studies] inferred_group_control_mean;  // unobserved per-group ideal treatment effect
  vector<lower=0>[n_studies] inferred_in_group_stddev;  // unobserved per-group treatment variance
",
idx_blocks("
  vector[nE[{IDX}]] treatment_subject_{IDX};  // unobserved per-group and subject treatment effects
  vector[nC[{IDX}]] control_subject_{IDX};  // unobserved per-group and subject control effects
", n_studies = n_studies),
"
}
transformed parameters {
  vector[n_studies] sampled_meanE;
  vector<lower=0>[n_studies] sampled_varE;
  vector[n_studies] sampled_meanC;
  vector<lower=0>[n_studies] sampled_varC;
",
idx_blocks("
  sampled_meanE[{IDX}] = mean(treatment_subject_{IDX});
  sampled_varE[{IDX}] = variance(treatment_subject_{IDX});
  sampled_meanC[{IDX}] = mean(control_subject_{IDX});
  sampled_varC[{IDX}] = variance(control_subject_{IDX});
", n_studies = n_studies),
"
}
model {
  // priors
  {//} inferred_grand_treatment_mean ~ normal(0, 1);
  {//} inferred_grand_control_mean ~ normal(0, 1);
  {//} inferred_between_group_stddev ~ lognormal(0, 1);
  inferred_group_treatment_mean ~ normal(0, 1);
  inferred_group_control_mean ~ normal(0, 1);
  inferred_in_group_stddev ~ lognormal(0, 1);
  // more peaked/informative stuff
",
idx_blocks("
  // each group generates an unobserved treatment response 
  {//} inferred_group_treatment_mean[{IDX}] ~ normal(inferred_grand_treatment_mean, inferred_between_group_stddev);
  // treatment subjects experience effects a function of unobserved group response
  // in the normal distribution case could avoid forming individual observations as we know
  // summary mean should be normally distributed and variance chi-square (with propoper parameters and scaling)
  treatment_subject_{IDX} ~ normal(inferred_group_treatment_mean[{IDX}], inferred_in_group_stddev[{IDX}]);
  // match observed summaries
  sampled_meanE[{IDX}] ~ normal(meanE[{IDX}], 0.01);
  sampled_varE[{IDX}] ~ normal(varE[{IDX}], 0.01);
  // each group generates an unobserved control response
  {//} inferred_group_control_mean[{IDX}] ~ normal(inferred_grand_control_mean, inferred_between_group_stddev);
  // control subjects experience effects a function of unobserved group response
  // in the normal distribution case could avoid forming individual observations as we know
  // summary mean should be normally distributed and variance chi-square (with propoper parameters and scaling)
  control_subject_{IDX} ~ normal(inferred_group_control_mean[{IDX}], inferred_in_group_stddev[{IDX}]);
  // match observed summaries
  sampled_meanC[{IDX}] ~ normal(meanC[{IDX}], 0.01);
  sampled_varC[{IDX}] ~ normal(varC[{IDX}], 0.01);
", n_studies = n_studies),
"
}
")
gsub("{//}", c_replacement, src, fixed = TRUE)
}


extract_dual_density_plot_data <- function(fit, c1, c2) {
  treatment_result <- data.frame(
    effect = as.data.frame(fit)[ , c1, drop=TRUE],
    estimate = c1)
  treatment_result['mean_effect'] <- mean(treatment_result[['effect']])
  treatment_result['sd_effect'] <- sd(treatment_result[['effect']])
  control_result <- data.frame(
    effect = as.data.frame(fit)[ , c2, drop=TRUE],
    estimate = c2)
  control_result['mean_effect'] <- mean(control_result[['effect']])
  control_result['sd_effect'] <- sd(control_result[['effect']])
  results <- rbind(treatment_result, control_result)
  results['z'] <- (treatment_result[1, 'mean_effect'] - control_result[1, 'mean_effect']) / 
    ( (treatment_result[1, 'sd_effect'] + control_result[1, 'sd_effect']) / 2 )
  return(results)
}


dual_density_plot <- function(fit, c1, c2, title, vlines=c()) {
  # show the inferred distribution of plausible results
  results <- extract_dual_density_plot_data(fit, c1, c2)
  z <- results[1, 'z']
  plt <- (
    ggplot(
      data = results,
      mapping = aes(x = effect, color = estimate, fill = estimate, linetype = estimate))
    + geom_density(alpha = 0.5)
    + geom_vline(
      mapping = aes(xintercept = mean_effect, color = estimate, linetype = estimate))
    + theme(legend.position="bottom")
    + ggtitle(
        title,
        subtitle=paste0('(z ~ ', sprintf('%5.2f', z), ')')
      )
  )
  if(length(vlines) > 0) {
    for(x in vlines) {
      plt <- (
        plt
          + geom_vline(xintercept = x,linetype = 3, color = "darkgray")
      )
    }
  }
  return(plt)
}


double_dual_density_plot <- function(fitA, fitB, c1, c2, title, vlines=c()) {
  # show the inferred distribution of plausible results
  resultsA <- extract_dual_density_plot_data(fitA, c1, c2)
  zA <- resultsA[1, 'z']
  resultsA['method'] <- paste0('independent estimate (z ~ ', sprintf('%5.2f', zA), ')')
  resultsB <- extract_dual_density_plot_data(fitB, c1, c2)
  zB <- resultsB[1, 'z']
  resultsB['method'] <- paste0('joint estimate (z ~ ', sprintf('%5.2f', zB), ')')
  results <- rbind(resultsA, resultsB)
  plt <- (
    ggplot(
      data = results,
      mapping = aes(x = effect, color = estimate, fill = estimate, linetype = estimate))
    + geom_density(alpha = 0.5)
    + geom_vline(
      mapping = aes(xintercept = mean_effect, color = estimate, linetype = estimate))
    + facet_wrap(~method, ncol=1)
    + theme(legend.position="bottom")
    + ggtitle(title)
  )
  if(length(vlines) > 0) {
    for(x in vlines) {
      plt <- (
        plt
        + geom_vline(xintercept = x, linetype = 3, color = "darkgray")
      )
    }
  }
  return(plt)
}



define_Stan_individuals_model <- function(n_studies, c_replacement = "") {
  # the Stan source code for combining multiple studies as a meta-analysis
  src <- paste0("
data {
  int<lower=1> n_studies;  // number of studies
  array[n_studies] int<lower=1> nE;  // number treated examples
  array[n_studies] int<lower=1> nC;  // number control examples
",
                idx_blocks("
  vector[nE[{IDX}]] treatment_subject_{IDX};  // observed per-group and subject treatment effects
  vector[nC[{IDX}]] control_subject_{IDX};  // observed per-group and subject control effects
", n_studies = n_studies),
                "
}
parameters {
  {//} real inferred_grand_treatment_mean;  // unobserved expected treatment
  {//} real inferred_grand_control_mean;  // unobserved expected control effect
  {//} real<lower=0> inferred_between_group_stddev; // standard distance between groups
  vector[n_studies] inferred_group_treatment_mean;  // unobserved per-group ideal treatment effect
  vector[n_studies] inferred_group_control_mean;  // unobserved per-group ideal treatment effect
  vector<lower=0>[n_studies] inferred_in_group_stddev;  // unobserved per-group treatment variance
}
model {
  // priors
  {//} inferred_grand_treatment_mean ~ normal(0, 1);
  {//} inferred_grand_control_mean ~ normal(0, 1);
  {//} inferred_between_group_stddev ~ lognormal(0, 1);
  inferred_group_treatment_mean ~ normal(0, 1);
  inferred_group_control_mean ~ normal(0, 1);
  inferred_in_group_stddev ~ lognormal(0, 1);
  // more peaked/informative stuff
",
idx_blocks("
  // each group generates an unobserved treatment response 
  {//} inferred_group_treatment_mean[{IDX}] ~ normal(inferred_grand_treatment_mean, inferred_between_group_stddev);
  // treatment subjects experience effects a function of unobserved group response
  treatment_subject_{IDX} ~ normal(inferred_group_treatment_mean[{IDX}], inferred_in_group_stddev[{IDX}]);
  // each group generates an unobserved control response
  {//} inferred_group_control_mean[{IDX}] ~ normal(inferred_grand_control_mean, inferred_between_group_stddev);
  // control subjects experience effects a function of unobserved group response
  control_subject_{IDX} ~ normal(inferred_group_control_mean[{IDX}], inferred_in_group_stddev[{IDX}]);
", n_studies = n_studies),
"
}
")
  gsub("{//}", c_replacement, src, fixed = TRUE)
}
