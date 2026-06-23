# Pixelgen theme

This function returns a ggplot2 theme that is styled with Pixelgen
branding.

## Usage

``` r
theme_pixelgen()
```

## Value

A ggplot2 theme

## Examples

``` r
library(ggplot2)
ggplot(mtcars, aes(x = wt, y = mpg)) +
  geom_point() +
  theme_pixelgen()

```
