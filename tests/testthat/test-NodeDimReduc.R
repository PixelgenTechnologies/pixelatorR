test_that("CreateNodeDimReducObject works as expected", {
  embeddings <- matrix(
    1:6,
    nrow = 3,
    dimnames = list(c("n1", "n2", "n3"), NULL)
  )
  loadings <- matrix(
    1:4,
    nrow = 2,
    dimnames = list(c("f1", "f2"), NULL)
  )
  dr <- CreateNodeDimReducObject(
    embeddings = embeddings,
    loadings = loadings,
    stdev = c(3, 2),
    key = "PC",
    method = "pca"
  )
  expect_s4_class(dr, "NodeDimReduc")
  expect_equal(colnames(SeuratObject::Embeddings(dr)), c("PC_1", "PC_2"))
  expect_equal(SeuratObject::Stdev(dr), c(3, 2))
  expect_equal(dr@key, "PC_")
  expect_equal(SeuratObject::Cells(dr), c("n1", "n2", "n3"))
  expect_equal(nrow(SeuratObject::Loadings(dr)), 2)
})

test_that("CreateNodeDimReducObject fails when invalid input is provided", {
  expect_error(CreateNodeDimReducObject(embeddings = 1))
  expect_error(CreateNodeDimReducObject(embeddings = matrix(1:4, nrow = 2)))
  embeddings <- matrix(1:4, nrow = 2, dimnames = list(c("a", "a"), NULL))
  expect_error(CreateNodeDimReducObject(embeddings = embeddings))
})
