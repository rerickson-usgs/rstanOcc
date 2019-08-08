context("Test sample detection one probability.")
library(rstanOcc)

sdo_0    <- sample_detection_one(J = 0, K = 8, theta = 0.5, p_detection = 0.75)
sdo_1    <- sample_detection_one(J = 1, K = 8, theta = 0.5, p_detection = 0.5)
sdo_1000 <- sample_detection_one(J = 1000, K = 8, theta = 0.5, p_detection = 0.5)

sdo_1000_k0 <-
    sample_detection_one(J = 1000, K = 0, theta = 0.5, p_detection = 0.5)

sdo_1000_t0 <-
    sample_detection_one(J = 1000, K = 8, theta = 0.0, p_detection = 0.5)

sdo_1000_p0 <-
    sample_detection_one(J = 1000, K = 8, theta = 0.5, p_detection = 0.0)



test_that("sample detection one function works as expected",
{
    expect_equal( sdo_0, 0.0)
    expect_equal(sdo_1, 0.498046875)
    expect_equal(sdo_1000, 1.0)
    expect_equal(sdo_1000_k0, 0.0)
    expect_equal(sdo_1000_t0, 0.0)
    expect_equal(sdo_1000_p0, 0.0)
}
)
