# Get Pixelgen accent colors

This function returns a vector of Pixelgen branded colors that can be
used as accent colors in a plot, given one or more hues and levels. If
only a hue is provided, all levels of that hue are returned. If only a
level is provided, all hues at that level are returned. If both a hue
and a level are provided, the colors at that hue and level are returned.

## Usage

``` r
PixelgenAccentColors(hue = NULL, level = NULL)
```

## Arguments

- hue:

  A character vector of hues to return. Options are "purples", "blues",
  "cyans", "greens", "pinks", "reds", "oranges", and "yellows". If NULL,
  all hues are returned.

- level:

  A numeric vector of levels to return. If NULL, all levels are
  returned.

## Value

A vector of colors

## Examples

``` r
PixelgenAccentColors(c("purples", "blues"), c(1, 3))
#>  purples1  purples3    blues1    blues3 
#> "#F7F5FD" "#C6B6EE" "#F3F6FD" "#A8BEE9" 
PixelgenAccentColors(c("cyans", "greens"))
#>    cyans1    cyans2    cyans3    cyans4    cyans5    cyans6    cyans7    cyans8 
#> "#F2FBFA" "#CBEFE9" "#A6E0D7" "#82D0C4" "#60BFB0" "#3FAC9B" "#209785" "#168777" 
#>    cyans9   cyans10   cyans11   cyans12   greens1   greens2   greens3   greens4 
#> "#0F7567" "#086256" "#044D44" "#013630" "#F4FBF3" "#D2ECD0" "#B0DCAD" "#90CA8C" 
#>   greens5   greens6   greens7   greens8   greens9  greens10  greens11  greens12 
#> "#71B96C" "#53A54D" "#369030" "#2A8025" "#206F1C" "#175C14" "#0F480D" "#093207" 
PixelgenAccentColors(level = 2)
#>       purples2         blues2         cyans2        greens2         pinks2 
#>      "#DED4F6"      "#CDDAF4"      "#CBEFE9"      "#D2ECD0"      "#FCF3F9" 
#>          reds2       oranges2       yellows2         greys2        beiges2 
#>      "#FDF0F2"      "#FDF6F0"      "#FDFBEC"      "#E7E7E7"      "#FCE9D6" 
#> standardblues2 
#>      "#EEF1F8" 
PixelgenAccentColors(hue = "reds")
#>     reds1     reds2     reds3     reds4     reds5     reds6     reds7     reds8 
#> "#FDF5F6" "#FDF0F2" "#FBE2E3" "#F7CACD" "#F19DA7" "#E96978" "#DD4154" "#CB2539" 
#>     reds9    reds10    reds11    reds12 
#> "#A72030" "#9B1828" "#7E0E1C" "#4F0E16" 
```
