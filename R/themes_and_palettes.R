#' Get Pixelgen gradient colors
#'
#' This function returns a vector of Pixelgen branded colors that can be used
#' to create a gradient in a plot.
#'
#' @param n The number of colors to return
#' @param name The name of the gradient to return. Options are "BluesCherry",
#' "Cherry", "Blues", "Magenta", and "Cyan"
#'
#' @return A vector of colors
#'
#' @examples
#' PixelgenGradient(5, "BluesCherry")
#' PixelgenGradient(5, "Cherry")
#' PixelgenGradient(5, "Blues")
#'
#' @rdname PixelgenGradient
#'
#' @export
#'
PixelgenGradient <- function(
  n,
  name
) {
  # Input validation
  pixelatorR:::assert_single_value(name, type = "string")
  pixelatorR:::assert_single_value(n, type = "integer")
  pixelatorR:::assert_within_limits(n, c(1, Inf))

  colors <-
    switch(name,
           "BluesCherry" = c(
             "#1F395F", "#496389", "#728BB1", "#AABAD1", "#DFE5EE", "#FFFFFF",
             "#FFE0EA", "#E9AABF", "#CD6F8D", "#A23F5E", "#781534"
           ),
           "GrayblueRose" = c("#798AAC", "#93A1BD", "#C4CBDB", "#FFFFFF", "#E8BFCD", "#D190A4", "#C1728B"),
           "Cherry" = c("#F2F2F2", "#FFE0EA", "#E9AABF", "#CD6F8D", "#A23F5E", "#781534"),
           "Blues" = c("#F2F2F2", "#DFE5EE", "#AABAD1", "#728BB1", "#496389", "#1F395F"),
           "Magenta" = c("#F2F2F2", "#FDE0EF", "#F1B6DA", "#DE77AE", "#C51C7D", "#8E0152"),
           "Cyan" = c("#F2F2F2", "#C2E5E1", "#9FE5DD", "#7CD5D0", "#59C5C3", "#36B5B6")
    )

  if (all(is.null(colors))) abort("Invalid palette name.")

  colors <-
    colorRampPalette(colors)(n)


  return(colors)
}

#' Get Pixelgen palette colors
#'
#' This function returns a vector of Pixelgen branded colors.
#'
#' @param n The number of colors to return
#' @param name The name of the palette to return. Options are "Tint", "Pastel",
#' "Semi-saturated", "Saturated", "Cells1", "Cells2", and "Cells3"
#'
#' @return A vector of colors
#'
#' @examples
#' PixelgenPalette(5, "Tint")
#' PixelgenPalette(5, "Pastel")
#' PixelgenPalette(5, "Semi-saturated")
#' PixelgenPalette(5, "Saturated")
#'
#' @rdname PixelgenPalette
#'
#' @export
#'
PixelgenPalette <- function(
  n,
  name
) {
  # Input validation
  pixelatorR:::assert_single_value(name, type = "string")
  pixelatorR:::assert_single_value(n, type = "integer")
  pixelatorR:::assert_within_limits(n, c(1, Inf))

  colors <- switch(name,
           "Tint" = c(
             "#E0E6EF", "#BECCE0", "#D1C6BB",
             "#DAD6D7", "#C4C4C4"
           ),
           "Pastel" = c(
             "#8197BD", "#637EA5", "#D887A0",
             "#C86584", "#D8BA98", "#E2A489",
             "#B4ADAF", "#978D89"
           ),
           "Semi-saturated" = c(
             "#4D988D", "#496389", "#1F395F",
             "#E05573", "#BF9871", "#918F8F"
           ),
           "Saturated" = c(
             "#1B9E8A", "#25C6F2", "#E24B7E",
             "#AA498D", "#FFC950", "#231F20"
           ),
           "Cells1" = c(
             "#E9CD98", "#E19DB0", "#526C92",
             "#BECCE0", "#9E9188", "#7E9EA3",
             "#DEBA95", "#C6C6C6", "#46616E",
             "#1F395F"
           ),
           "Cells2" = c(
             "#C89433", "#809EA2", "#85756C",
             "#48696E", "#BA9D80", "#556C92",
             "#EBCD97", "#21395E"
           ),
           "Cells3" = c(
             "#A3A9CC", "#8DBFB3", "#F2EBC0",
             "#F3B462", "#F06060", "#44D593",
             "#B99095", "#E5E6E6", "#5D4C52",
             "#07475A"
           )
    )

  if (all(is.null(colors))) abort("Invalid palette name.")
  if (n > length(colors)) abort(glue("Palette '{name}' only has {length(colors)} colors."))

  colors <- colors[seq_len(n)]

  return(colors)
}

#' Get Pixelgen palette colors
#'
#' This function returns palettes that have previously been used at Pixelgen.
#'
#' @param name The name of the palette to return. Options are "PNA product sheet 2025"
#' and "PNA preprint 2025".
#'
#' @return A vector of colors
#'
#' @examples
#' PixelgenLegacyPalette("PNA product sheet 2025")
#'
#' @rdname PixelgenLegacyPalette
#'
#' @export
#'
PixelgenLegacyPalette <- function(
  name
) {
  # Input validation
  pixelatorR:::assert_single_value(name, type = "string")

  colors <-
    switch(name,
           "PNA product sheet 2025" =
             c(
               "Naive CD4 T" = "#B9CDED",
               "TCM CD4 T" = "#4A73C0",
               "TEM CD4 T" = "#224792",
               "Tregs" = "#1C3A76",
               "Naive CD8 T" = "#D0EDE6",
               "TCM CD8 T" = "#4AAF9D",
               "TEM CD8 T" = "#1B7E6F",
               "MAIT" = "#156559",
               "CD56dim NK" = "#A28EDB",
               "CD56bright NK" = "#866CCD",
               "Naive B" = "#F0DAC4",
               "Intermediate B" = "#C7A989",
               "Memory B" = "#917557",
               "mDC" = "#DA94C1",
               "pDC" = "#BB5391",
               "CD14 Mono" = "#E6BB43",
               "CD16 Mono" = "#AE8A1B",
               "Neutrophils" = "#CDCDCD",
               "Basophils" = "#797979",
               "Platelets" = "#7C2628",
               "gdT" = "#5A3E9E"
             ),

           "PNA preprint 2025" =
             c(
               "Naive CD4 T" = "#B9CDED",
               "CD4 TSCM" = "#92B0E0",
               "CD4 TCM" = "#6D92D1",
               "CD4 TEM" = "#4A73C0",
               "CD4 TEFF" = "#224792",
               "CD4 TEMRA" = "#1C3A76",
               "Treg" = "#2955AE",
               "Naive CD8 T" = "#D0EDE6",
               "CD8 TSCM" = "#A2DACE",
               "CD8 TCM" = "#75C5B5",
               "CD8 TEM" = "#4AAF9D",
               "CD8 TEFF" = "#209785",
               "CD8 TEMRA" = "#1B7E6F",
               "MAIT" = "#156559",
               "DPT" = "#1B7E9F",
               "DNT" = "#7B787F",
               "CD56dim NK" = "#A28EDB",
               "CD56bright NK" = "#866CCD",
               "NKT" = "#6C4ABD",
               "Naive B" = "#F0DAC4",
               "Memory B" = "#917557",
               "Plasma cells" = "#DE9982",
               "cDC1" = "#DA94C1",
               "cDC2" = "#BB5391",
               "pDC" = "#9C4579",
               "Classical Mono" = "#F0C966",
               "Intermediate Mono" = "#DAAC22",
               "Non-classical Mono" = "#836714",
               "Neutrophils" = "#CDCDCD",
               "Basophils" = "#797979",
               "Platelets" = "#7C2628",
               "gdT" = "#5A3E9E"
             )

    )

  if (all(is.null(colors))) abort("Invalid palette name.")

  return(colors)
}

#' Pixelgen theme
#'
#' This function returns a ggplot2 theme that is styled with Pixelgen branding.
#'
#' @return A ggplot2 theme
#'
#' @examples
#' library(ggplot2)
#' ggplot(mtcars, aes(x = wt, y = mpg)) +
#'   geom_point() +
#'   theme_pixelgen()
#'
#' @rdname theme_pixelgen
#'
#' @importFrom ggplot2 theme_bw element_blank element_rect
#'
#' @export
#'
theme_pixelgen <- function() {
  theme_bw() +
    theme(
      panel.grid = element_blank(),
      strip.background = element_rect(fill = "gray95")
    )
}


#' Pixelgen ggplot2 color themes
#'
#' These functions return ggplot2 color scales that uses Pixelgen accent colors.
#' Works both for discrete and sequential data and sets colors both for fill and color
#' aesthetics.
#'
#' @param hue The hue to use for the colors. See `Pixelgen_accent_colors` for
#' available hues.
#' @param level The level of the hue to use. See `Pixelgen_accent_colors` for
#' available levels.
#' @param direction The direction of the color gradient. 1 for normal, -1 for reversed.
#' @param hue_low,hue_high The low/high hue for the color gradient. For example,
#' `hue_low = 'blues'` and `hue_high = 'reds'` will create a gradient from blue to red
#' with white in the middle.
#' @param limits The limits of the color scale.
#' @param aes_type The aesthetics type to color. Either "fill" or "color" or NULL to
#' color both (default).
#' @param shuffle Whether to shuffle the colors in the resulting palette. Default is FALSE.
#' @param indices A vector of indices to control the selection of colors.
#'
#' @return A list of ggplot2 scales that can be added to a ggplot object.
#'
#' @examples
#' library(ggplot2)
#' library(dplyr)
#'
#' # Discrete colors
#' ggplot(mtcars, aes(cyl, wt, fill = factor(gear))) +
#'   geom_col() +
#'   color_discrete_pixelgen()
#'
#' # Sequential scale
#' ggplot(mtcars, aes(cyl, wt, color = disp)) +
#'   geom_point() +
#'   color_sequential_pixelgen()
#'
#' # Divergent scale
#' ggplot(mtcars %>% mutate(disp = scale(disp)), aes(cyl, gear, fill = disp)) +
#'   geom_tile() +
#'   color_divergent_pixelgen(limits = c(-2, 2), hue_low = "cyans", hue_high = "beiges")
#'
#' @name color-themes
#' @rdname color-themes
#'
#' @export
#'
color_discrete_pixelgen <- function(
  hue = NULL,
  level = NULL,
  aes_type = NULL,
  shuffle = FALSE,
  indices = NULL
) {
  pixelatorR:::assert_single_value(aes_type, type = "string", allow_null = TRUE)
  pixelatorR:::assert_vector(indices, type = "integer", allow_null = TRUE)
  if (!is.null(hue) || !is.null(level)) {
    cols <- PixelgenAccentColors(hue, level) %>% unname()
    if (!is.null(indices)) {
      cols <- cols[indices]
    }
    if (shuffle && is.null(indices)) {
      cols <- sample(cols)
    }
  } else {
    rows <- c(7, 6, 8, 5, 9, 4, 10, 3, 11, 2, 12, 1)
    columns <- c(1, 3, 5, 7, 2, 4, 6, 8, 9, 10, 11)
    if (!is.null(indices)) {
      columns <- columns[indices]
    }
    if (shuffle && is.null(indices)) {
      columns <- sample(columns)
    }
    cols <- Pixelgen_accent_colors[rows, columns] %>%
      as.matrix() %>%
      t() %>%
      as.vector()
  }
  if (is.null(aes_type)) {
    return(list(scale_fill_manual(values = cols), scale_color_manual(values = cols)))
  }
  if (aes_type == "fill") {
    return(scale_fill_manual(values = cols))
  }
  if (aes_type == "color") {
    return(scale_color_manual(values = cols))
  }
}

#' @rdname color-themes
#'
#' @export
#'
color_sequential_pixelgen <- function(
  hue = "purples",
  direction = 1,
  aes_type = NULL
) {
  pixelatorR:::assert_single_value(direction, type = "integer")
  pixelatorR:::assert_single_value(aes_type, type = "string", allow_null = TRUE)
  if (!direction %in% c(-1, 1)) {
    cli::cli_abort(
      c(
        "x" = "{.var direction} must be either -1 or 1."
      )
    )
  }
  cols <- PixelgenAccentColors(hue, level = NULL)
  if (direction == -1) {
    cols <- rev(cols)
  }
  if (is.null(aes_type)) {
    return(list(scale_fill_gradientn(colors = cols), scale_color_gradientn(colors = cols)))
  }
  if (aes_type == "fill") {
    return(scale_fill_gradientn(colors = cols))
  }
  if (aes_type == "color") {
    return(scale_color_gradientn(colors = cols))
  }
}

#' @rdname color-themes
#'
#' @export
#'
color_divergent_pixelgen <- function(
  hue_low = "blues",
  hue_high = "reds",
  limits = NULL,
  aes_type = NULL
) {
  pixelatorR:::assert_single_value(aes_type, type = "string", allow_null = TRUE)
  cols_low <- PixelgenAccentColors(hue_low, level = NULL) %>% rev()
  cols_high <- PixelgenAccentColors(hue_high, level = NULL)
  list(
    scale_fill_gradientn(colors = c(cols_low, "#FFFFFF", cols_high), limits = limits),
    scale_color_gradientn(colors = c(cols_low, "#FFFFFF", cols_high), limits = limits)
  )
  if (is.null(aes_type)) {
    return(list(
      scale_fill_gradientn(colors = c(cols_low, "#FFFFFF", cols_high), limits = limits),
      scale_color_gradientn(colors = c(cols_low, "#FFFFFF", cols_high), limits = limits)
    ))
  }
  if (aes_type == "fill") {
    return(scale_fill_gradientn(colors = c(cols_low, "#FFFFFF", cols_high), limits = limits))
  }
  if (aes_type == "color") {
    return(scale_color_gradientn(colors = c(cols_low, "#FFFFFF", cols_high), limits = limits))
  }
}

#' Pixelgen accent colors
#'
#' With a soft base we allow accent colors to pop and reveal important data and
#' insights without complicating the picture. Accent colors are meant to be used
#' for data visualisation or as small touches in a composition. The accent colors
#' can also be used sparingly to emphasize clickable content and to create visual
#' hierarchy.
#'
#' @format A tibble with 12 rows and 11 columns
#'
#' @rdname Pixelgen_accent_colors
#'
#' @export
Pixelgen_accent_colors <-
  tibble(
    purples = c(
      "#F7F5FD", "#DED4F6", "#C6B6EE", "#AF98E4", "#987DD8", "#8263CC",
      "#7E54E4", "#6D39DA", "#5F29C9", "#3C2178", "#2D165D", "#1F0E41"
    ),
    blues = c(
      "#F3F6FD", "#CDDAF4", "#A8BEE9", "#86A3DD", "#6588CF", "#466EC0",
      "#2955AE", "#1E469B", "#143887", "#0C2B71", "#061F58", "#02143E"
    ),
    cyans = c(
      "#F2FBFA", "#CBEFE9", "#A6E0D7", "#82D0C4", "#60BFB0", "#3FAC9B",
      "#209785", "#168777", "#0F7567", "#086256", "#044D44", "#013630"
    ),
    greens = c(
      "#F4FBF3", "#D2ECD0", "#B0DCAD", "#90CA8C", "#71B96C", "#53A54D",
      "#369030", "#2A8025", "#206F1C", "#175C14", "#0F480D", "#093207"
    ),
    pinks = c(
      "#FDF7FB", "#FCF3F9", "#F9E9F5", "#F6D1EB", "#F0ADDB", "#E57AC0",
      "#D953B1", "#C63288", "#AB216D", "#8D205A", "#771E4D", "#46122E"
    ),
    reds = c(
      "#FDF5F6", "#FDF0F2", "#FBE2E3", "#F7CACD", "#F19DA7", "#E96978",
      "#DD4154", "#CB2539", "#A72030", "#9B1828", "#7E0E1C", "#4F0E16"
    ),
    oranges = c(
      "#FDFAF6", "#FDF6F0", "#FAEBDC", "#F5D6BA", "#EFC097", "#E69C6A",
      "#DD7C45", "#CE5E2E", "#AE4525", "#8D3823", "#6C2D1F", "#491F15"
    ),
    yellows = c(
      "#FEFDF2", "#FDFBEC", "#FCF6CD", "#FAEDA6", "#F7E188", "#EFC438",
      "#DDAA29", "#BE861F", "#955E1A", "#7B4C1C", "#693D1D", "#4E2E15"
    ),
    greys = c(
      "#F9F9F9", "#E7E7E7", "#D5D5D5", "#C2C2C2", "#B0B0B0", "#9D9D9D",
      "#8B8B8B", "#737373", "#5C5C5C", "#444444", "#2D2D2D", "#151515"
    ),
    beiges = c(
      "#FFF6EE", "#FCE9D6", "#F6DCC1", "#EDCFAF", "#E2C29F", "#D5B592",
      "#C4A788", "#AC9175", "#947B62", "#7B6550", "#624F3E", "#48392C"
    ),
    standardblues = c(
      "#F8F9FD", "#EEF1F8", "#E0E6EF", "#CBD5E5", "#B2C1D9", "#98ABCA",
      "#7C90B1", "#607291", "#465671", "#344157", "#242D3E", "#161D2B"
    )
  )


#' Show accent colors
#'
#' This function displays the Pixelgen accent colors in a grid format.
#'
#' @param label_discrete_order Logical, the color order used for
#' \code{\link{color_discrete_pixelgen}} are shown in the plot.
#'
#' @return A ggplot object displaying the accent colors.
#'
#' @export
#'
show_accent_colors <- function(
  label_discrete_order = FALSE
) {
  accent_colors <- Pixelgen_accent_colors %>%
    mutate(level = factor(paste0(seq_len(n())), levels = seq_len(n()))) %>%
    pivot_longer(cols = colnames(Pixelgen_accent_colors), names_to = "hue", values_to = "color") %>%
    mutate(hue = factor(hue, levels = names(Pixelgen_accent_colors)))

  if (label_discrete_order) {
    rows <- paste0(c(7, 6, 8, 5, 9, 4, 10, 3, 11, 2, 12, 1))
    columns <- names(Pixelgen_accent_colors)[c(1, 3, 5, 7, 2, 4, 6, 8, 9, 10, 11)]
    ord <- tibble(
      level = rep(rows, each = length(columns)) %>% factor(levels = levels(accent_colors$level)),
      hue = rep(columns, times = length(rows)) %>% factor(levels = levels(accent_colors$hue)),
    ) %>%
      mutate(label = seq_len(n())) %>%
      mutate(text_color = if_else(
        level %in% paste0(1:6),
        "#000000", "#FFFFFF"
      ))
    accent_colors <- accent_colors %>%
      left_join(ord, by = c("level", "hue"))
  }

  p <- ggplot(accent_colors, aes(x = level, y = hue)) +
    geom_tile(aes(fill = color), color = "white") +
    scale_fill_identity() +
    labs(title = "Pixelgen Accent Colors", x = "Level", y = "Hue") +
    theme_minimal() +
    scale_y_discrete(limits = rev) +
    theme(plot.background = element_rect(fill = "#F2F2F2"),
          axis.text = element_text(size = 10))

  if (label_discrete_order) {
    p <- p +
      geom_text(aes(label = label, color = text_color)) +
      scale_color_identity()
  }

  return(p)
}

#' Get Pixelgen accent colors
#'
#' This function returns a vector of Pixelgen branded colors that can be used
#' as accent colors in a plot, given one or more hues and levels. If only a
#' hue is provided, all levels of that hue are returned. If only a level is
#' provided, all hues at that level are returned. If both a hue and a level
#' are provided, the colors at that hue and level are returned.
#'
#' @param hue A character vector of hues to return. Options are "purples",
#' "blues", "cyans", "greens", "pinks", "reds", "oranges", and "yellows".
#' If NULL, all hues are returned.
#' @param level A numeric vector of levels to return. If NULL, all levels are
#' returned.
#'
#' @return A vector of colors
#'
#' @examples
#' PixelgenAccentColors(c("purples", "blues"), c(1, 3))
#' PixelgenAccentColors(c("cyans", "greens"))
#' PixelgenAccentColors(level = 2)
#' PixelgenAccentColors(hue = "reds")
#'
#' @rdname PixelgenAccentColors
#'
#' @export
#'
PixelgenAccentColors <- function(
  hue = NULL,
  level = NULL
) {
  # Input validation
  pixelatorR:::assert_vector(hue, type = "character", n = 1, allow_null = TRUE)
  pixelatorR:::assert_vector(level, type = "numeric", n = 1, allow_null = TRUE)

  stopifnot(
    "'hue' and 'level' cannot both be NULL." = !(is.null(hue) && is.null(level))
  )

  if (is.null(level)) level <- seq_len(nrow(Pixelgen_accent_colors))
  if (is.null(hue)) hue <- colnames(Pixelgen_accent_colors)

  pixelatorR:::assert_x_in_y(hue, colnames(Pixelgen_accent_colors))
  pixelatorR:::assert_within_limits(level, c(1, nrow(Pixelgen_accent_colors)))

  # Get colors
  colors <-
    Pixelgen_accent_colors[level, hue] %>%
    unlist() %>%
    set_names(paste0(
      rep(hue, each = length(level)),
      rep(level, length(hue))
    ))

  return(colors)
}
