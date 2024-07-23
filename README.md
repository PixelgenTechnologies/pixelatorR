<!-- badges: start -->
[![DOI](https://img.shields.io/badge/DOI%3A-10.1038%2Fs41592--024--02268--9-red
)](https://www.nature.com/articles/s41592-024-02268-9)
[![codecov](https://codecov.io/gh/PixelgenTechnologies/pixelatorR/graph/badge.svg?token=ClGH1zHvuD)](https://codecov.io/gh/PixelgenTechnologies/pixelatorR)
[![R-CMD-check](https://github.com/PixelgenTechnologies/pixelatorR/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/PixelgenTechnologies/pixelatorR/actions/workflows/R-CMD-check.yaml)
![Static Badge](https://img.shields.io/badge/beta_release-v0.10.1-orange)
<!-- badges: end -->

[**Installation**](#installation) |
[**Usage**](#usage) |
[**Notes**](#notes) |
[**Contact**](#contact) |
[**Credits**](#credits) |
[**Documentation**](https://pixelgentechnologies.github.io/pixelatorR/reference/index.html) |
[**Tutorials**](https://software.pixelgen.com/mpx-analysis/introduction)

# pixelatorR

pixelatorR provides the infrastructure to process, analyze and visualize MPX data for R users.

<br/>

> [!NOTE] 
> Visit our [package website](https://pixelgentechnologies.github.io/pixelatorR/) for more details about installation. Here, you can also access the function documentation and tutorials for advanced users.

<br/>

## Installation

pixelatorR can be installed from GitHub by running the following code from R:

``` r
install.packages("remotes")
remotes::install_github("PixelgenTechnologies/pixelatorR")
```

### Installation from source

You can also install pixelatorR from source by cloning the repository.

```shell
git clone https://github.com/pixelgentechnologies/pixelatorR.git
```

Then, in an R session, run:

```r
install.packages("<path to pixelatorR directory>", repos = NULL, type = "source")
```

## Usage

Visit our [tutorials](https://software.pixelgen.com/mpx-analysis/introduction) for a step-by step guide on MPX data analysis with pixelatorR. For advanced users, we also provide additional tutorials on our [package website](https://pixelgentechnologies.github.io/pixelatorR/), which are more directed towards users who are interested in the details of how pixelatorR stores and handles MPX data in R.

## Notes

pixelatorR is designed to work with objects types from the [SeuratObject](https://github.com/satijalab/seurat-object) R package. For the best user experience, we recommend installing [Seurat v5](https://satijalab.org/seurat/).

## Contact

For feature requests or bug reports, please use the GitHub [issues](https://github.com/PixelgenTechnologies/pixelatorR/issues).

You can also email the development team at [developers@pixelgen.com](mailto:developers@pixelgen.com).

## Credits

pixelatorR is developed and maintained by the [developers](https://github.com/PixelgenTechnologies) at [Pixelgen Technologies](https://pixelgen.com).

When using pixelatorR in your research, please cite the following publication:

> Karlsson, F., Kallas, T., Thiagarajan, D. et al. Molecular pixelation: spatial proteomics of single cells by sequencing. Nat Methods (2024). https://doi.org/10.1038/s41592-024-02268-9

Main developers:

- Ludvig Larsson ([@ludvigla](https://github.com/ludvigla))
- Max Karlsson ([@maxkarlsson](https://github.com/maxkarlsson))
- Vincent van Hoef ([@vincent-van-hoef](https://github.com/vincent-van-hoef))
- Alvaro Martinez Barrio ([@ambarrio](https://github.com/ambarrio))
- Johan Dahlberg ([@johandahlberg](https://github.com/johandahlberg))

A huge thank you to all [code contributors](https://github.com/PixelgenTechnologies/pixelatorR/graphs/contributors)!
