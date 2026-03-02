library(Seurat)
seur <- ReadMPX_Seurat(system.file("extdata/five_cells", "five_cells.pxl", package = "pixelatorR"))

test_that("CalculateDispersion works as expected", {
  expect_no_error(seur <- CalculateDispersion(seur, layer = "counts"))
  expect_equal(
    seur[[]]$dispersion_gini,
    c(
      0.874683039462636,
      0.885185848863404,
      0.897175910527191,
      0.87788441903206,
      0.89908581068758
    )
  )

  expect_no_error(seur <- CalculateDispersion(seur, method = "tau", layer = "counts", margin = 2))
  expect_equal(
    seur[[]]$dispersion_tau,
    c(
      0.983626936942866,
      0.977802893309222,
      0.980362089500532,
      0.978009070904472,
      0.974876484085913
    )
  )

  expect_no_error(seur <- CalculateDispersion(seur, layer = "counts", metadata_name = "test"))
  expect_equal(
    seur[[]]$test,
    c(
      0.874683039462636,
      0.885185848863404,
      0.897175910527191,
      0.87788441903206,
      0.89908581068758
    )
  )

  expect_no_error(disp <- CalculateDispersion(LayerData(seur, layer = "counts"), margin = 2))
  expect_equal(seur[[]]$dispersion_gini, unname(disp))

  expect_named(
    disp,
    c(
      "RCVCMP0000217",
      "RCVCMP0000118",
      "RCVCMP0000487",
      "RCVCMP0000655",
      "RCVCMP0000263"
    )
  )

  expect_no_error(disp <- CalculateDispersion(as.matrix(LayerData(seur, layer = "counts")), margin = 2))
  expect_equal(seur[[]]$dispersion_gini, unname(disp))

  expect_no_error(disp <- CalculateDispersion(as.data.frame(LayerData(seur, layer = "counts")), margin = 2))
  expect_equal(seur[[]]$dispersion_gini, unname(disp))

  expect_no_error(disp <- CalculateDispersion(seur, layer = "counts", margin = 1))
  expect_equal(
    disp,
    c(
      CD274 = 0.241025641025641, CD44 = 0.381586679725759, CD25 = 0.336,
      CD279 = 0.575438596491228, CD41 = 0.391836734693878, `HLA-ABC` = 0.221447596988998,
      CD54 = 0.233850702143385, CD26 = 0.656728778467909, CD27 = 0.470229007633588,
      CD38 = 0.705684626791893, CD16 = 0.572727272727273, CD52 = 0.403088803088803,
      CD53 = 0.246153846153846, CD11c = 0.189473684210526, CD11a = 0.4,
      CD127 = 0.345945945945946, CD29 = 0.38328530259366, CD82 = 0.326629834254144,
      CD45RB = 0.338010471204188, CD40 = 0.544329896907217, CD19 = 0.681818181818182,
      CD8 = 0.629558541266795, CD59 = 0.626916890080429, TCRb = 0.116455696202532,
      mIgG2a = 0.238805970149254, CD11b = 0.333333333333333, CD86 = 0.698637602179836,
      CD197 = 0.257931034482759, `HLA-DR` = 0.705603335939536, CD3E = 0.424242424242424,
      CD2 = 0.407458292443572, CD20 = 0.75125, CD45RA = 0.662576687116564,
      CD14 = 0.171428571428571, CD4 = 0.695182636315511, mIgG2b = 0.2125,
      mIgG1 = 0.32, CD9 = 0.375, CD69 = 0.237837837837838, B2M = 0.259154109339739,
      CD36 = 0.451948051948052, CD45 = 0.210753197536712, CD152 = 0.363636363636364,
      CD337 = 0.64, CD1d = 0.507692307692308, CD84 = 0.363636363636364,
      CD161 = 0.539130434782609, CD163 = 0.35, CD200 = 0.45, CD137 = 0.48,
      CD229 = 0.320567375886525, CD244 = 0.665306122448979, CD154 = 0.4,
      CD18 = 0.489879931389365, CD71 = 0.572307692307692, ACTB = 0.3,
      CD48 = 0.20983606557377, CD43 = 0.525301204819277, CD150 = 0.375,
      CD22 = 0.754562383612663, CD62P = 0.363636363636364, CD50 = 0.476683937823834,
      CD33 = 0.228571428571429, CD37 = 0.705882352941177, CD162 = 0.595918367346939,
      CD328 = 0.422222222222222, CD7 = 0.452705882352941, CD102 = 0.292857142857143,
      CD47 = 0.311287128712871, CD72 = 0.779679144385027, CD5 = 0.470404172099087,
      CD55 = 0.374736842105263, CD278 = 0.601980198019802, CD32 = 0.741279069767442,
      CD268 = 0.417910447761194, CD64 = 0.133333333333333, CD49D = 0.339047619047619,
      CD158 = 0.457142857142857, CD314 = 0.664705882352941, CD35 = 0.774421768707483
    )
  )

  expect_no_error(disp2 <- CalculateDispersion(LayerData(seur, layer = "counts"), margin = 1))
  expect_equal(disp, disp2)

  expect_no_error(disp <- CalculateDispersion(seur, method = "tau", layer = "counts", margin = 1))
  expect_equal(
    disp,
    c(
      CD274 = 0.402173913043478, CD44 = 0.637889688249401, CD25 = 0.698529411764706,
      CD279 = 0.884615384615385, CD41 = 0.569444444444444, `HLA-ABC` = 0.466566866267465,
      CD54 = 0.513071895424837, CD26 = 0.934725848563969, CD27 = 0.73828125,
      CD38 = 0.9525, CD16 = 0.883333333333333, CD52 = 0.682017543859649,
      CD53 = 0.476190476190476, CD11c = 0.458333333333333, CD11a = 0.734375,
      CD127 = 0.633333333333333, CD29 = 0.736686390532544, CD82 = 0.626721763085399,
      CD45RB = 0.677458033573141, CD40 = 0.876923076923077, CD19 = 0.960526315789474,
      CD8 = 0.826652221018418, CD59 = 0.903603268945022, TCRb = 0.352272727272727,
      mIgG2a = 0.452380952380952, CD11b = 0.5, CD86 = 0.967257318952234,
      CD197 = 0.640756302521008, `HLA-DR` = 0.946727790072716, CD3E = 0.613696808510638,
      CD2 = 0.613920099875156, CD20 = 0.982441471571906, CD45RA = 0.938931297709924,
      CD14 = 0.375, CD4 = 0.914595170454545, mIgG2b = 0.45, mIgG1 = 0.625,
      CD9 = 0.75, CD69 = 0.462765957446809, B2M = 0.593230063109581,
      CD36 = 0.780487804878049, CD45 = 0.452191987906274, CD152 = 0.638888888888889,
      CD337 = 0.833333333333333, CD1d = 0.708333333333333, CD84 = 0.708333333333333,
      CD161 = 0.866666666666667, CD163 = 0.583333333333333, CD200 = 0.75,
      CD137 = 0.625, CD229 = 0.530612244897959, CD244 = 0.947530864197531,
      CD154 = 0.625, CD18 = 0.835348506401138, CD71 = 0.863095238095238,
      ACTB = 0.5, CD48 = 0.556818181818182, CD43 = 0.689189189189189,
      CD150 = 0.75, CD22 = 0.984158415841584, CD62P = 0.7, CD50 = 0.815315315315315,
      CD33 = 0.666666666666667, CD37 = 0.968009478672986, CD162 = 0.889705882352941,
      CD328 = 0.607142857142857, CD7 = 0.690789473684211, CD102 = 0.654255319148936,
      CD47 = 0.670871559633028, CD72 = 0.991355463347165, CD5 = 0.715877437325905,
      CD55 = 0.710227272727273, CD278 = 0.814655172413793, CD32 = 0.977416798732171,
      CD268 = 0.742424242424242, CD64 = 0.5, CD49D = 0.656108597285068,
      CD158 = 0.75, CD314 = 0.956896551724138, CD35 = 0.989730878186969
    )
  )


  expect_error(CalculateDispersion(seur, layer = "counts", margin = 3))
  expect_error(CalculateDispersion(seur, layer = "counts", margin = "2"))
  expect_error(CalculateDispersion(seur, layer = "nolayer"))
  expect_error(CalculateDispersion(seur, metadata_name = 37))
})
