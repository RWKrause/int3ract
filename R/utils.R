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
  if (!length(found))
    stop("Interaction term not found for: ", paste(parts, collapse = " * "),
         "\nAvailable: ", paste(nms, collapse = ", "), call. = FALSE)
  found[1]
}
