test_that("PixelgenPalette works as expected", {
  expect_no_error(PixelgenPalette(name = "Tint", n = 5))
  expect_no_error(PixelgenPalette(name = "Pastel", n = 5))
  expect_no_error(PixelgenPalette(name = "Semi-saturated", n = 5))
  expect_no_error(PixelgenPalette(name = "Saturated", n = 5))
  expect_no_error(PixelgenPalette(name = "Cells1", n = 5))
  expect_no_error(PixelgenPalette(name = "Cells2", n = 5))
  expect_no_error(PixelgenPalette(name = "Cells3", n = 5))

  expect_error(PixelgenPalette(name = "Cells2", n = 10))
  expect_error(PixelgenPalette(name = "not a palette", n = 5))
})

test_that("PixelgenGradient works as expected", {
  expect_no_error(PixelgenGradient(name = "BluesCherry", n = 5))
  expect_no_error(PixelgenGradient(name = "Cherry", n = 5))
  expect_no_error(PixelgenGradient(name = "Blues", n = 5))
  expect_no_error(PixelgenGradient(name = "Magenta", n = 5))
  expect_no_error(PixelgenGradient(name = "Cyan", n = 5))
  expect_no_error(PixelgenGradient(name = "Cyan", n = 50))
  expect_no_error(PixelgenGradient(name = "Cyan", n = 500))
  expect_no_error(PixelgenGradient(name = "Cyan", n = 5000))

  expect_error(PixelgenGradient(name = "Cyan", n = 0))
  expect_error(PixelgenGradient(name = "not a palette", n = 5))
})

test_that("PixelgenLegacyPalette works as expected", {
  expect_no_error(PixelgenLegacyPalette(name = "PNA product sheet 2025"))

  expect_error(PixelgenLegacyPalette(name = "not a palette"))
})



test_that("theme_pixelgen works as expected", {
  expect_s3_class(theme_pixelgen(), "theme")
})


test_that("PixelgenAccentColors works as expected", {
  expect_identical(
    PixelgenAccentColors(hue = "greens"),
    c(
      greens1 = "#F4FBF3", greens2 = "#D2ECD0", greens3 = "#B0DCAD",
      greens4 = "#90CA8C", greens5 = "#71B96C", greens6 = "#53A54D",
      greens7 = "#369030", greens8 = "#2A8025", greens9 = "#206F1C",
      greens10 = "#175C14", greens11 = "#0F480D", greens12 = "#093207"
    )
  )
  expect_identical(
    PixelgenAccentColors(hue = "greens", level = 1),
    c(greens1 = "#F4FBF3")
  )
  expect_identical(
    PixelgenAccentColors(hue = c("greens", "reds"), level = 1:3),
    c(
      greens1 = "#F4FBF3", greens2 = "#D2ECD0", greens3 = "#B0DCAD",
      reds1 = "#FDF5F6", reds2 = "#FDF0F2", reds3 = "#FBE2E3"
    )
  )
  expect_identical(
    PixelgenAccentColors(level = 1),
    c(
      purples1 = "#F7F5FD", blues1 = "#F3F6FD", cyans1 = "#F2FBFA",
      greens1 = "#F4FBF3", pinks1 = "#FDF7FB", reds1 = "#FDF5F6", oranges1 = "#FDFAF6",
      yellows1 = "#FEFDF2", greys1 = "#F9F9F9", beiges1 = "#FFF6EE",
      standardblues1 = "#F8F9FD"
    )
  )

  expect_error(PixelgenAccentColors())
  expect_error(PixelgenAccentColors(hue = "not a hue"))
  expect_error(PixelgenAccentColors(level = 0))
  expect_error(PixelgenAccentColors(hue = "purples", level = 0))
  expect_error(PixelgenAccentColors(hue = "purples", level = 13))
})


library(ggplot2)
p_discrete <- ggplot(mtcars, aes(cyl, wt, fill = factor(gear))) +
  geom_col()
p_seq <- ggplot(mtcars, aes(cyl, wt, color = disp)) +
  geom_point()

test_that("color_discrete_pixelgen, color_sequential_pixelgen and color_divergent_pixelgen works as expected", {
  expect_no_error(p_discrete + color_discrete_pixelgen())
  expect_no_error(p_seq + color_sequential_pixelgen())
  expect_no_error(p_seq + color_divergent_pixelgen())
})
