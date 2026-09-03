# =============================================================================
# Shared helpers.
# =============================================================================

#' Locate an interaction term among a set of parameter names
#'
#' Model matrices spell interactions in whichever order the formula produced
#' them, so every ordering of the components is tried.
#'
#' @param nms available parameter names.
#' @param ... two or three component names.
#' @returns The matching parameter name.
#' @noRd
.find_int <- function(nms, ...) {
  parts <- c(...)
  found <- .find_int_quiet(nms, parts)
  if (is.null(found))
    stop("Interaction term not found for: ", paste(parts, collapse = " * "),
         "\nAvailable: ", paste(nms, collapse = ", "), call. = FALSE)
  found
}


#' Locate an interaction term, or return NULL
#'
#' The body of \code{.find_int()}, for callers that would rather handle a
#' missing term themselves than stop. lavaan input needs this: it accepts
#' two-way products in its own syntax but not three-way ones, so a missing
#' term there calls for a different explanation.
#'
#' @param nms available parameter names.
#' @param parts character vector of two or three component names.
#' @returns The matching parameter name, or \code{NULL}.
#' @noRd
.find_int_quiet <- function(nms, parts) {
  candidates <- if (length(parts) == 2L) {
    c(paste(parts, collapse = ":"), paste(rev(parts), collapse = ":"))
  } else {
    p <- parts
    c(paste(p[1], p[2], p[3], sep = ":"),
      paste(p[1], p[3], p[2], sep = ":"),
      paste(p[2], p[1], p[3], sep = ":"),
      paste(p[2], p[3], p[1], sep = ":"),
      paste(p[3], p[1], p[2], sep = ":"),
      paste(p[3], p[2], p[1], sep = ":"))
  }
  found <- candidates[candidates %in% nms]
  if (length(found)) found[1] else NULL
}
