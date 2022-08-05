test_that("Data partitioning.", {
  
  # Case 1: Complete data.
  data <- data.frame(
    t = c(0, 1, 2),
    s = c(0, 1, 2),
    x = c(0, 1, 2)
  )
  
  part <- PartitionData(data$t, data$s, data$x)
  expect_equal(c(part$Dims$n0, part$Dims$n1, part$Dims$n2), c(3, 0, 0))
  
  # Case 2: Target missing.
  data <- data.frame(
    t = c(NA, NA, 2),
    s = c(0, 1, 2),
    x = c(0, 1, 2)
  )
  
  part <- PartitionData(data$t, data$s, data$x)
  expect_equal(c(part$Dims$n0, part$Dims$n1, part$Dims$n2), c(1, 2, 0))  

  # Case 3: Surrogate missing.
  data <- data.frame(
    t = c(0, 1, 2),
    s = c(NA, NA, 2),
    x = c(0, 1, 2)
  )  
  
  part <- PartitionData(data$t, data$s, data$x)
  expect_equal(c(part$Dims$n0, part$Dims$n1, part$Dims$n2), c(1, 0, 2))  
  
  # Case 4: Bilateral missingness.
  data <- data.frame(
    t = c(NA, NA, 2),
    s = c(0, NA, 2),
    x = c(0, 1, 2)
  )  
  
  part <- PartitionData(data$t, data$s, data$x)
  expect_equal(c(part$Dims$n0, part$Dims$n1, part$Dims$n2), c(1, 1, 0)) 
    
})