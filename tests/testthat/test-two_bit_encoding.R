test_that("pack_2bits and unpack_2bits work as expected", {
  expect_no_error(pack_2bits("ACTG"))
  expect_no_error(pack_2bits(c("ACTG", "ACTT")))
  expect_no_error(unpack_2bits(bit64::as.integer64(180), k = 4L))
  expect_no_error(unpack_2bits(bit64::as.integer64(c(180, 187)), k = 4L))
  expect_no_error(unpack_2bits(c("180", "187"), k = 4L))


  test_seqs <-
    c(
      "ACTG",
      "TTTGT",
      "GCGCGC",
      "GCGCGCGC",
      "TCGATCGA",
      "TCGATCGATCGATCCGA",
      "TCGTGGACCCGAGGCGGGGGGACAGGGT"
    )

  for (seq in test_seqs) {
    packed <- pack_2bits(seq)
    unpacked <- unpack_2bits(packed, k = nchar(seq))
    expect_identical(unpacked, seq)
  }

  expect_identical(
    unpack_2bits("18012321353423242", k = 27),
    "GGAGCGCTCGCTTTCCATCAGTTTTTT"
  )
})

test_that("pack_2bits and unpack_2bits fails with invalid input", {
  expect_error(pack_2bits("Invalid"))
  expect_error(unpack_2bits("Invalid"))
  expect_error(unpack_2bits(bit64::as.integer64(180), k = "Invalid"))
  expect_error(pack_2bits(c("ACTG", "ACTTA")))
  expect_error(pack_2bits(c("ACTG", "ACTN")))
})
