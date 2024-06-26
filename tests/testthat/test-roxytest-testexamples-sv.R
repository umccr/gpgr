# Generated by roxytest: do not edit by hand!

# File R/sv.R: @testexamples

test_that("Function split_double_col() @ L27", {
  
  x <- tibble::tibble(
    a = letters[1:11],
    b = c("0.4,0.8", paste0(round(runif(10), 2), ",", round(runif(10), 2))),
    nacol = rep(NA, 11),
    namix = sample(c(NA, "0.4,0.6"), 11, replace = TRUE)
  )
  (b <- gpgr:::split_double_col(x, "b"))
  (nacol <- gpgr:::split_double_col(x, "nacol"))
  (namix <- gpgr:::split_double_col(x, "namix"))
  expect_equal(colnames(b), "b")
  expect_equal(nrow(x), 11)
  expect_error(gpgr:::split_double_col(x, "c"))
  expect_equal(b$b[1], "0.6 (0.4, 0.8)")
})


test_that("Function count_pieces() @ L77", {
  
  (a <- gpgr:::count_pieces("foo,bar,baz", sep = ","))
  (b <- gpgr:::count_pieces("foo", sep = ","))
  (k <- gpgr:::count_pieces("", sep = ","))
  (m <- gpgr:::count_pieces(",", sep = ","))
  expect_equal(a, 3)
  expect_equal(b, 1)
  expect_equal(k, 0)
  expect_equal(m, 2)
  expect_error(gpgr:::count_pieces("foo", NA))
})


test_that("Function abbreviate_effect() @ L112", {
  
  (e1 <- gpgr:::abbreviate_effect("3_prime_UTR_truncation&start_lost&splice_region_variant"))
  (e2 <- gpgr:::abbreviate_effect("duplication&foo&gene_fusion&BOOM&intron_variant"))
  (e3 <- gpgr:::abbreviate_effect("TF_binding_site_variant&TFBS_ablation"))
  (e4 <- gpgr:::abbreviate_effect("foo&bar&stop_gained&badaboom"))
  expect_equal(e1, "3UTRtrunc, SpliceRegV, StartLoss")
  expect_equal(e2, "BOOM, Dup, foo, FusG, IntronV")
  expect_equal(e3, "TFBSDel, TFBSVar")
  expect_equal(e4, "badaboom, bar, foo, StopGain")
})

