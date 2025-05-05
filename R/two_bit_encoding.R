#' Encode a DNA sequence using 2 bits per base
#'
#' Converts a DNA string (containing only A, C, G, T) into a 64-bit integer using 2 bits per base.
#' This allows compact storage and efficient comparison of short DNA sequences (kmers).
#'
#' @param seq A character string of DNA bases (e.g., "ACGT"). Must only contain A, C, G, or T.
#'
#' @return An object of class `integer64` representing the 2-bit packed sequence.
#'
#' @examples
#' library(bit64)
#' pack_2bits("TCGT") # returns integer64 encoding
#'
#' @export
pack_2bits <- function(seq) {
  expect_bit64()

  assert_class(seq, "character")
  assert_single_value(seq)

  lookup <- c(A = 0L, C = 1L, G = 2L, T = 3L)

  # Convert sequence string into individual bases
  seq_chars <- rev(strsplit(seq, split = "")[[1]])

  # Initialize packed as 64-bit integer
  packed <- bit64::as.integer64(0)

  for (base in seq_chars) {
    val <- lookup[[base]]
    if (is.null(val)) {
      stop("Invalid base: ", base)
    }

    # Shift packed 2 bits to the left (multiply by 4)
    packed <- packed * 4L
    # Add the value corresponding to the current base
    packed <- packed + val
  }

  return(packed)
}

#' Decode a 2-bit encoded DNA sequence from a 64-bit integer
#'
#' Given a 64-bit integer encoding a DNA sequence using 2 bits per base (as from \code{pack_2bits}),
#' this function reconstructs the original DNA string.
#'
#' @param packed A 64-bit integer (class `integer64`) or character representing the packed DNA sequence.
#' @param k Integer. The number of bases in the original DNA sequence.
#'
#' @return A character string of DNA bases (e.g., "ACGT").
#'
#' @examples
#' library(bit64)
#' packed <- pack_2bits("TCGT")
#' unpack_2bits(packed, 4) # returns "TCGT"
#'
#' @export
unpack_2bits <- function(packed, k) {
  expect_bit64()

  assert_class(packed, c("integer64", "character"))

  if (is.character(packed)) {
    packed <- bit64::as.integer64(packed)
  }

  assert_single_value(packed, type = "numeric")
  assert_class(k, c("numeric", "integer"))
  assert_single_value(k, type = "numeric")

  lookup <- c("A", "C", "G", "T")

  seq_chars <- character(k)

  for (i in 1:k) {
    # Get last 2 bits
    val <- as.integer(packed %% 4)
    seq_chars[i] <- lookup[val + 1]

    # Divide by 4 to right shift by 2 bits
    packed <- packed %/% 4
  }

  dna_seq <- paste(seq_chars, collapse = "")

  return(dna_seq)
}
