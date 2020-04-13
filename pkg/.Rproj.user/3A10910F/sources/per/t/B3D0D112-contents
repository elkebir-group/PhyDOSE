context("Simple test")

test_that("A simple example works", { 
    expect_true(file.exists("../../data/example.txt"))
    expect_equal(PhyDOSE("../../data/example.txt"), 37)
})

test_that("A non existent file returns -1", { 
    expect_false(file.exists("../../data/nonexistent.txt"))
    expect_equal(PhyDOSE("../../data/nonexistent.txt"), -1)
})