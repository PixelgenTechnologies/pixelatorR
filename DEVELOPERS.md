## Linting

To run the linter, you need to install `lintr`. Then you can use one of the the following commands:

```r
# Lint entire pacakge
lintr::lint_package()

# Lint single file
lintr::lint("path/to/file.R")
```

Alternatively, you can run the linter from RStudio through Addins -> Lint current file or Addins -> Lint current package.

The configuration file `.lintr` is used to specify the rules that the linter should follow. For compatibility with styler, some linting rules have been disabled.

## Styler

To style the code, you need to install `styler`. Then you can use one of the the following commands:

```r
# Style entire package
styler::style_pkg(transformers = pixelatorR::pixelatorR_style())

# Style single file
styler::style_file("path/to/file.R", transformers = pixelatorR::pixelatorR_style())
```

Alternatively, you can run the styler from RStudio. Here, you need to configure `styler` to use the style guide provided in pixelatorR. Go to Addins -> Styler -> Set style and set `pixelatorR::pixelatorR_style()` as the style guide. Then you can use Addins -> Styler -> Style active file or Addins -> Styler -> Style active package.
