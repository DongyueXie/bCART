#' @title Mute any message from a function
#' @description source: \url{https://r.789695.n4.nabble.com/Suppressing-output-e-g-from-cat-td859876.html}
#' @export

quiet <- function(x) {
  sink(tempfile())
  on.exit(sink())
  invisible(force(x))
}
