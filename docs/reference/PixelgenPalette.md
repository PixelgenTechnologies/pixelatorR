# Get Pixelgen palette colors

This function returns a vector of Pixelgen branded colors.

## Usage

``` r
PixelgenPalette(n, name)
```

## Arguments

- n:

  The number of colors to return

- name:

  The name of the palette to return. Options are "Tint", "Pastel",
  "Semi-saturated", "Saturated", "Cells1", "Cells2", and "Cells3"

## Value

A vector of colors

## Examples

``` r
PixelgenPalette(5, "Tint")
#> [1] "#E0E6EF" "#BECCE0" "#D1C6BB" "#DAD6D7" "#C4C4C4"
PixelgenPalette(5, "Pastel")
#> [1] "#8197BD" "#637EA5" "#D887A0" "#C86584" "#D8BA98"
PixelgenPalette(5, "Semi-saturated")
#> [1] "#4D988D" "#496389" "#1F395F" "#E05573" "#BF9871"
PixelgenPalette(5, "Saturated")
#> [1] "#1B9E8A" "#25C6F2" "#E24B7E" "#AA498D" "#FFC950"
```
