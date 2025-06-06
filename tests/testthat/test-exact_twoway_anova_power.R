test_that("input checks", {
  refmean <- 1
  treatgroups <- 4
  timepoints <- 5
  treateff <- 1.25
  timeeff <- 0.85
  cellswithinteraction <- matrix(c(rep(2,3), 3:5), 3,2) #second level of factor A interacts with 3rd, 4th and 5th level of factor B
  effects_treat_time_interact <- calculate_mean_matrix(refmean = refmean, nlfA = treatgroups, nlfB = timepoints,
                                                       fAeffect = treateff, fBeffect = timeeff, label_list = list(treatment=letters[1:treatgroups], time=1:timepoints),
                                                       groupswinteraction = cellswithinteraction, interact=1.3)
  fxs <- effsize(effects_treat_time_interact)
  expect_error(exact_twoway_anova_power(a= treateff, b=timepoints, effect_sizes=fxs, n=seq(6,16,2)))
  expect_error(exact_twoway_anova_power(a= treatgroups, b=timepoints, effect_sizes=fxs[-1], n=seq(6,16,2)))
  expect_error(exact_twoway_anova_power(a= treatgroups, b=timepoints, effect_sizes=fxs, n=treateff))
  expect_error(exact_twoway_anova_power(a= treatgroups, b=timepoints, effect_sizes=fxs, n=seq(6,16,2),
                                        factor_names = c("treament")))
})

test_that("factor name assignment", {
  refmean <- 1
  treatgroups <- 4
  timepoints <- 5
  treateff <- 1.25
  timeeff <- 0.85
  cellswithinteraction <- matrix(c(rep(2,3), 3:5), 3,2) #second level of factor A interacts with 3rd, 4th and 5th level of factor B
  effects_treat_time_interact <- calculate_mean_matrix(refmean = refmean, nlfA = treatgroups, nlfB = timepoints,
                                                       fAeffect = treateff, fBeffect = timeeff, label_list = list(treatment=letters[1:treatgroups], time=1:timepoints),
                                                       groupswinteraction = cellswithinteraction, interact=1.3)
  fxs <- effsize(effects_treat_time_interact)
  res <- exact_twoway_anova_power(a= treatgroups, b=timepoints, effect_sizes=fxs, n=seq(6,16,2),
                                        factor_names = c("meds", "follow-up"))
  expect_named(res$exactpower$powercurve, c("n", "meds", "follow-up", "meds:follow-up"))
  res <- exact_twoway_anova_power(a= treatgroups, b=timepoints, effect_sizes=as.vector(fxs), n=seq(6,16,2))
  expect_named(res$exactpower$powercurve, c("n", "FactorA", "FactorB", "Interaction"))
})



test_that("alternative hypothesis", {
  refmean <- 1
  treatgroups <- 4
  timepoints <- 5
  treateff <- 1.25
  timeeff <- 0.85
  cellswithinteraction <- matrix(c(rep(2,3), 3:5), 3,2) #second level of factor A interacts with 3rd, 4th and 5th level of factor B
  effects_treat_time_interact <- calculate_mean_matrix(refmean = refmean, nlfA = treatgroups, nlfB = timepoints,
                                                       fAeffect = treateff, fBeffect = timeeff, label_list = list(treatment=letters[1:treatgroups], time=1:timepoints),
                                                       groupswinteraction = cellswithinteraction, interact=1.3)
  fxs <- effsize(effects_treat_time_interact)
  res <- exact_twoway_anova_power(a= treatgroups, b=timepoints, effect_sizes=fxs, n=seq(6,16,2))
  expect_true(all(res$exactpower$powercurve$treatment>0.9))
  expect_true(all(res$exactpower$powercurve$`treatment:time`[4:6]>0.8))
})


test_that("null hypothesis", {
  refmean <- 1
  treatgroups <- 4
  timepoints <- 5
  treateff <- 1
  timeeff <- 1
  nulleffects_treat_time <- calculate_mean_matrix(refmean = refmean, nlfA = treatgroups, nlfB = timepoints,
                                                       fAeffect = treateff, fBeffect = timeeff, label_list = list(treatment=letters[1:treatgroups], time=1:timepoints))
  fxs <- effsize(nulleffects_treat_time)
  res <- exact_twoway_anova_power(a= treatgroups, b=timepoints, effect_sizes=fxs, n=seq(6,16,2))
  expect_true(all(res$exactpower$powercurve$treatment<0.06))
  expect_true(all(res$exactpower$powercurve$`treatment:time`<0.06))
})

test_that("plot option", {
  refmean <- 1
  treatgroups <- 4
  timepoints <- 5
  treateff <- 1
  timeeff <- 1
  nulleffects_treat_time <- calculate_mean_matrix(refmean = refmean, nlfA = treatgroups, nlfB = timepoints,
                                                  fAeffect = treateff, fBeffect = timeeff, label_list = list(treatment=letters[1:treatgroups], time=1:timepoints))
  fxs <- effsize(nulleffects_treat_time)
  res <- exact_twoway_anova_power(a= treatgroups, b=timepoints, effect_sizes=fxs, n=seq(6,16,2))
  expect_s3_class(res$exactpower$powercurve, "data.frame")
  expect_s3_class(res[[2]], "gg")
  res <- exact_twoway_anova_power(a= treatgroups, b=timepoints, effect_sizes=fxs, n=seq(6,16,2), plot = FALSE)
  expect_equal(length(res), 6)
})
