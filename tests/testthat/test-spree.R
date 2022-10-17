# Test if totals after spree match with survey totals
#library(testthat)

# Load needed data
data("census")
data("survey_17")
data("P2017c")
data("survey_12")




names(census)[1]<-"Canton"

test_that("Totals after spree match with survey totals?", {
  totals_c<- colSums(spree(population_domains= "Canton", sample_domains = "Canton", population_data=census, sample_data=survey_12, type="SPREE")$updated_point[,-1])
  totals_r<- rowSums(spree(population_domains= "Canton", sample_domains = "Canton", population_data=census, sample_data=survey_12, type="SPREE")$updated_point[,-1])
  expect_equal(totals_c, colSums(survey_12[,-1]))
  expect_equal(totals_r, rowSums(survey_12[,-1]))
})


test_that("Totals after gspree match with survey totals?", {
  totals_c<- colSums(spree(population_domains= "Canton",  sample_domains = "Canton", population_data=census, sample_data=survey_12, type="GSPREE")$updated_point[,-1])
  totals_r<- rowSums(spree(population_domains= "Canton", sample_domains = "Canton", population_data=census, sample_data=survey_12, type="GSPREE")$updated_point[,-1])
  expect_equal(totals_c, colSums(survey_12[,-1]))
  expect_equal(totals_r, rowSums(survey_12[,-1]))

})

test_that("Totals after mspree with ML match with survey totals?", {
  totals_c<- colSums(spree(population_domains= "Canton",  sample_domains = "Canton", population_data=census, sample_data=survey_12, type="MSPREE", method = "ML")$updated_point[,-1])
  totals_r<- rowSums(spree(population_domains= "Canton",  sample_domains = "Canton", population_data=census, sample_data=survey_12, type="MSPREE", method = "ML")$updated_point[,-1])
  expect_equal(totals_c, colSums(survey_12[,-1]))
  expect_equal(totals_r, rowSums(survey_12[,-1]))

})




