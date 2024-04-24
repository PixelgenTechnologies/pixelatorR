# pixelatorR

pixelatorR provides the infrastructure to process, analyze and visualize MPX data for R users.

<p align="center">
    <img src="https://www.pixelgen.com/wp-content/uploads/2022/12/share-image-pixelgen.png" height=200
     alt="Pixelgen Technologies" />
</p>
<div align="center">© 2023 - Pixelgen Technologies AB</div>

## Installation

pixelatorR can be installed from GitHub by running the following code from R:

``` r
install.packages("remotes")
remotes::install_github("PixelgenTechnologies/pixelatorR")
```

### Installation from source

You can also install pixelatorR from source by cloning the repository.

```console
git clone https://github.com/pixelgentechnologies/pixelatorR.git
```

Then, in an R session, run:

```r
install.packages("<path to pixelatorR directory>", repos = NULL, type = "source")
```

## Usage

Visit our [tutorials](https://software.pixelgen.com/mpx-analysis/introduction) for a step-by step guide on MPX data analysis with pixelatorR. On this site, we also provide additional additional tutorials that are more directed towards users who are interested in the details of how pixelatorR stores and handles MPX data in R.

Function documentation can be accessed from the [reference](https://pixelgentechnologies.github.io/pixelatorR/reference/) tab above. Alternatively, you can access the documentation from an R session by running `?function_name` once the package is installed. 

## Notes

pixelatorR is designed to work with classes from the [SeuratObject](https://github.com/satijalab/seurat-object) R package. For the best analysis experience, we recommend installing [Seurat v5](https://satijalab.org/seurat/).

## Contact

For feature requests or bug reports, please use the GitHub [issues](https://github.com/PixelgenTechnologies/pixelatorR/issues).

You can also email the development team at [developers@pixelgen.com](mailto:developers@pixelgen.com).

## Credits

pixelatorR is developed and maintained by the [developers](https://github.com/PixelgenTechnologies) at [Pixelgen Technologies](https://pixelgen.com).

When using pixelator in your research, please cite the following publication:

> Karlsson, Filip, Tomasz Kallas, Divya Thiagarajan, Max Karlsson, Maud Schweitzer, Jose Fernandez Navarro, Louise Leijonancker, _et al._ “Molecular Pixelation: Single Cell Spatial Proteomics by Sequencing.” bioRxiv, June 8, 2023. https://doi.org/10.1101/2023.06.05.543770.

Main developers:

- Ludvig Larsson ([@ludvigla](https://github.com/ludvigla))
- Max Karlsson ([@maxkarlsson](https://github.com/maxkarlsson))

A huge thank you to all [code contributors](https://github.com/PixelgenTechnologies/pixelatorR/graphs/contributors)!
