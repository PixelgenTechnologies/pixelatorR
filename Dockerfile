# Install system dependencies
FROM rocker/r-ver:4

# Install system dependencies for HDF5 and R packages
RUN apt-get update && apt-get install -y \
    libhdf5-dev \
    hdf5-tools \
    libhdf5-serial-dev \
    libcurl4-openssl-dev \
    libssl-dev \
    libtk8.6 \
    libxml2-dev \
    patch \
    && rm -rf /var/lib/apt/lists/*

# Install R package dependencies
RUN R -e "install.packages('hdf5r', repos='https://cloud.r-project.org/')"
RUN R -e "install.packages(c('pak', 'devtools'))"

# Copy package source and install
WORKDIR /app
COPY . /app
# Install pixelatorR dependencies
RUN R -e "pak::pak()"
# Install pixelatorR
RUN R -e "devtools::install()"
