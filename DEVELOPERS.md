## Installation

pixelatorR has more than 20 dependencies and depending on the system you are working on, some of these dependencies can create problems. For most systems, the easiest way to install pixelatorR is to set up a conda environment with the necessary dependencies and install it within that environment. In our GitHub repo, we provide some tasks which can be used to set up such an environment and install pixelatorR.

Below is a list of requirements to use these tasks:

- [micromamba](https://mamba.readthedocs.io/en/latest/installation/micromamba-installation.html)
- [task](https://taskfile.dev/installation/)

First you need to clone the repo and navigate to pixelatorR:

````
git clone PixelgenTechnologies/pixelatorR
cd pixelatorR
````

### Windows, Linux and Mac on intel

The following task can be used to setup a conda environment named `r-mpx` with pixelatorR installed:

````
task setup-env-pixelatorR ENV_NAME=r-mpx
````

The task will install all dependencies listed in the `environment.yml` file and then the pixelatorR package from GitHub.

Activate the environment with:

````
micromamba activate r-mpx
````

### For Mac OS on ARM

On machines with ARM processors (e.g. M1 or M2), the solution provided above will not work. The reason is that at the time of writing, conda-forge doesn't provide all required pixelatorR dependencies for ARM processors.

Instead we can setup an empty environment with the following task:

````
task setup-env-bare ENV_NAME=r-mpx
````

Then we can install the dependencies using the [pak](https://pak.r-lib.org/) instead. `pak` is a package manager for R that can install packages from CRAN, Bioconductor, GitHub, and other sources. It can also also fetch pre-built binaries from Posit Package Manager (PPM) if these are available. This is useful for systems where conda-forge doesn't provide all required dependencies. 

````
micromamba activate r-mpx
task install-pixelatorR-with-pak
````

Note that compared to the micromamba solution, pak is significantly slower.

### Rosetta 2

An alternative option is to use [Rosetta 2](https://support.apple.com/en-us/102527) to emulate an intel processor. With this emulation, you can use the first solution and add and additional flag to the `micromamba` command:

````
task setup-env-pixelatorR ENV_NAME=r-mpx MICROMAMBA_ARGS="--platform osx-64"
````

Now micromamba should be able to install all dependencies as these are all available for intel processors. 

This is not the recommended solution for ARM processors for two reasons:

- performance penalty : the emulation will slow down the execution of the code
- arrow : if R is running R under emulation, arrow might segfault without error. See [this](https://github.com/apache/arrow/pull/37777) issue on GitHub for more information. You can still use this configuration, but you will see a warnings message on package load when running under emulation on macOS (i.e., use of x86 installation of R on M1/aarch64).

## package dev tasks

Create [roxygen2](https://roxygen2.r-lib.org/) documentation for the pacakge:

````
task document
````

***

Run all package tests (`devtools::test()`):

````
task test-all
````

***

Run all examples (`devtools::run_examples()`):

````
task test-examples
````

***

Test staged R test scripts. This task will look for any test R scripts that are staged for commit and run these. The test scripts are located in the `tests/testthat` directory.

````
task test-staged-files
````

***

Style staged R scripts using [styler](https://styler.r-lib.org/). This task will look for any R scripts that are staged for commit and style these based on the pixelatorR style guide.

````
task style-staged-files
````


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
