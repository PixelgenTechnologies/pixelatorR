library(Seurat)

for (assay_version in c("v3", "v5")) {
  options(Seurat.object.assay.version = assay_version)

  # Create a dummy Seurat object
  set.seed(37)
  object <-
    CreateSeuratObject(
      counts = matrix(
        c(
          rpois(100000, 40),
          rpois(100000, 5)
        )[sample(1:200000, 200000)],
        nrow = 100,
        ncol = 2000,
        dimnames = list(
          paste0("Feature", 1:100),
          paste0("Cell", 1:2000)
        )
      ) %>%
        as("dgCMatrix")
    ) %>%
      AddMetaData(
        metadata = data.frame(
          sample = rep(c("A", "B"), each = 1000),
          sample_type = rep(
            c("Unstimulated", "Stimulated"),
            each = 500,
            times = 2
          ),
          row.names = paste0("Cell", 1:2000)
        )
      )

  test_that("DensityScatterPlot works as expected", {
    # No facetting, no gating
    expect_no_error(
      DensityScatterPlot(
        object,
        marker1 = "Feature1",
        marker2 = "Feature2",
        layer = "counts"
      )
    )

    # Other colors
    expect_no_error(
      DensityScatterPlot(
        object,
        marker1 = "Feature1",
        marker2 = "Feature2",
        layer = "counts",
        colors = c("blue", "red")
      )
    )

    # Marginal density
    expect_no_error(
      DensityScatterPlot(
        object,
        marker1 = "Feature1",
        marker2 = "Feature2",
        layer = "counts",
        margin_density = T,
        coord_fixed = F
      )
    )

    # Facetting by two variables
    expect_no_error(
      DensityScatterPlot(
        object,
        marker1 = "Feature1",
        marker2 = "Feature2",
        layer = "counts",
        facet_vars = c("sample", "sample_type"),
        scale_density = F
      )
    )

    # Single common gate
    plot_gate <-
      tibble(
        xmin = 20,
        xmax = 70,
        ymin = 20,
        ymax = 50
      )

    # No facetting
    expect_no_error(
      DensityScatterPlot(
        object,
        marker1 = "Feature1",
        marker2 = "Feature2",
        layer = "counts",
        plot_gate = list(plot_gate, "rectangle")
      )
    )

    # Facetting by one variable
    expect_no_error(
      DensityScatterPlot(
        object,
        marker1 = "Feature1",
        marker2 = "Feature2",
        facet_vars = "sample",
        layer = "counts",
        plot_gate = list(plot_gate, "rectangle")
      )
    )

    # Gating by one variable
    plot_gate <-
      tibble(
        xmin = c(20, 25),
        xmax = c(70, 70),
        ymin = c(20, 25),
        ymax = c(70, 50),
        sample = c("A", "B")
      )

    # Facetting by one variable
    expect_no_error(
      DensityScatterPlot(
        object,
        marker1 = "Feature1",
        marker2 = "Feature2",
        facet_vars = "sample",
        layer = "counts",
        plot_gate = list(plot_gate, "rectangle")
      )
    )

    # Facetting by two variables
    expect_no_error(
      DensityScatterPlot(
        object,
        marker1 = "Feature1",
        marker2 = "Feature2",
        facet_vars = c("sample", "sample_type"),
        layer = "counts",
        plot_gate = list(plot_gate, "rectangle")
      )
    )

    # Multiple gates
    plot_gate <-
      tibble(
        xmin = c(-1, 21, -1),
        xmax = c(20, 71, 18),
        ymin = c(-1, 21, -1),
        ymax = c(70, 60, 18)
      )

    expect_no_error(
      DensityScatterPlot(
        object,
        marker1 = "Feature1",
        marker2 = "Feature2",
        facet_vars = c("sample", "sample_type"),
        layer = "counts",
        plot_gate = list(plot_gate, "rectangle")
      )
    )

    # Test plot_gate list format and quadrant gates
    test_that("DensityScatterPlot handles new gate formats correctly", {
      rect_gate <- tibble(
        xmin = 20,
        xmax = 70,
        ymin = 20,
        ymax = 50
      )

      expect_no_error(
        DensityScatterPlot(
          object,
          marker1 = "Feature1",
          marker2 = "Feature2",
          plot_gate = list(rect_gate, "rectangle")
        )
      )

      # Quadrant gate
      quad_gate <- tibble(
        x = 30,
        y = 40
      )

      expect_no_error(
        DensityScatterPlot(
          object,
          marker1 = "Feature1",
          marker2 = "Feature2",
          plot_gate = list(quad_gate, "quadrant")
        )
      )

      # Quadrant gate with faceting
      quad_gate_faceted <- tibble(
        x = c(30, 35),
        y = c(40, 45),
        sample = c("A", "B")
      )

      expect_no_error(
        DensityScatterPlot(
          object,
          marker1 = "Feature1",
          marker2 = "Feature2",
          facet_vars = "sample",
          plot_gate = list(quad_gate_faceted, "quadrant")
        )
      )

      # Test annotation_params
      expect_no_error(
        DensityScatterPlot(
          object,
          marker1 = "Feature1",
          marker2 = "Feature2",
          plot_gate = list(quad_gate, "quadrant"),
          annotation_params = list(
            color = "red",
            size = 4,
            fontface = "bold"
          )
        )
      )
    })

    # Test error cases for new functionality
    test_that("DensityScatterPlot errors correctly for invalid inputs", {
      # Invalid gate type
      expect_error(
        DensityScatterPlot(
          object,
          marker1 = "Feature1",
          marker2 = "Feature2",
          plot_gate = list(tibble(x = 30, y = 40), "invalid_type")
        ),
        "Gate type must be either 'rectangle' or 'quadrant'"
      )

      # Wrong format for plot_gate (not a list)
      expect_error(
        DensityScatterPlot(
          object,
          marker1 = "Feature1",
          marker2 = "Feature2",
          plot_gate = tibble(x = 30, y = 40)
        )
      )

      # Missing required columns for quadrant gate
      expect_error(
        DensityScatterPlot(
          object,
          marker1 = "Feature1",
          marker2 = "Feature2",
          plot_gate = list(tibble(x = 30), "quadrant")
        )
      )

      # Missing required columns for rectangle gate
      expect_error(
        DensityScatterPlot(
          object,
          marker1 = "Feature1",
          marker2 = "Feature2",
          plot_gate = list(tibble(xmin = 20, xmax = 70), "rectangle")
        )
      )

      # Invalid annotation_params (not a list)
      expect_error(
        DensityScatterPlot(
          object,
          marker1 = "Feature1",
          marker2 = "Feature2",
          plot_gate = list(quad_gate, "quadrant"),
          annotation_params = "invalid"
        )
      )
    })

    # Expected errors

    plot_gate <-
      tibble(
        xmin = c(20, 25),
        xmax = c(70, 70),
        ymin = c(20, 25),
        ymax = c(70, 50),
        sample = c("A", "B")
      )

    expect_error(
      DensityScatterPlot(
        object,
        marker1 = "FeatureNotHere",
        marker2 = "Feature2",
        layer = "counts"
      )
    )
    expect_error(
      DensityScatterPlot(
        object,
        marker1 = "Feature1",
        marker2 = "Feature2",
        layer = "counts",
        facet_vars = c("sample", "sample_type", "sample_type")
      )
    )
    expect_error(
      DensityScatterPlot(
        object,
        marker1 = "Feature1",
        marker2 = "Feature2",
        layer = "counts",
        facet_vars = c("sample", "columnNotHere")
      )
    )
    expect_error(
      DensityScatterPlot(
        object,
        marker1 = "Feature1",
        marker2 = "Feature2",
        layer = "counts",
        plot_gate = list(tibble(xmin = 20), "rectangle")
      )
    )
    expect_error(
      DensityScatterPlot(
        object,
        marker1 = "Feature1",
        marker2 = "Feature2",
        facet_vars = NULL,
        layer = "counts",
        plot_gate = plot_gate
      )
    )
  })
}
