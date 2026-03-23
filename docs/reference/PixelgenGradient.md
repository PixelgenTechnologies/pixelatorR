# Get Pixelgen gradient colors

This function returns a vector of Pixelgen branded colors that can be
used to create a gradient in a plot.

## Usage

``` r
PixelgenGradient(n, name)
```

## Arguments

- n:

  The number of colors to return

- name:

  The name of the gradient to return. Options are "BluesCherry",
  "BluesGrayCherry", "GrayblueRose", "Cherry", "Blues", "Magenta", and
  "Cyan".

## Value

A vector of colors

## Examples

``` r
PixelgenGradient(5, "BluesCherry")
#> [1] "#1F395F" "#8DA2C1" "#FFFFFF" "#DB8CA6" "#781534"
PixelgenGradient(5, "Cherry")
#> [1] "#F2F2F2" "#F9D2DF" "#DB8CA6" "#AC4B69" "#781534"
PixelgenGradient(5, "Blues")
#> [1] "#F2F2F2" "#D1DAE6" "#8EA2C1" "#536D93" "#1F395F"
```
