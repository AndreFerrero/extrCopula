#' No adaptation scheme
#'
#' @return Adaptation object with update() method
adapt_none <- function() {

  list(
    type = "none",
    #' Update proposal state (identity)
    #'
    #' @param state Current proposal state
    #' @param param Current parameter
    #' @param accept Logical indicating acceptance
    #' @param iter Current iteration number
    #' @return Unchanged state
    update = function(state, param, accept, iter) {
      state
    }
  )
}
