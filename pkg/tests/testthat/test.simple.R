context("Simple test")

test_that("A simple example works", { 
    expect_true(file.exists("../../data/example.txt"))
    ret <- PhyDOSE("../../data/example.txt")
    expect_equal(ret$k_star, 37)
})

test_that("A non existent file returns -1", { 
    expect_false(file.exists("../../data/nonexistent.txt"))
    ret <- PhyDOSE("../../data/nonexistent.txt")
    expect_equal(ret$k_star, -1)
})