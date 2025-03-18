#' Title
#'
#' @return
#' @export
#'
#' @examples
start_timer <- function() {
  start <- Sys.time()
  print(paste0("Start time: ", round(start)))
  return(start)
}

#' Title
#'
#' @return
#' @export
#'
#' @examples
stop_timer <- function() {
  stop <- Sys.time()
  print(paste0("Stop time: ", round(stop)))
  return(stop)
}

#' Title
#'
#' @param start
#' @param stop
#'
#' @return
#' @export
#'
#' @examples
print_run_time <- function(start, stop) {
  run_time <- stop - start
  units <- attributes(run_time)[[2]]
  return(paste0("Run time: ", round(run_time, 2), " ", units))
}
