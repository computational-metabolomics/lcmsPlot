setClassUnion("characterOrNULL", c("character", "NULL"))

setClassUnion("numericOrNULL", c("numeric", "NULL"))

create_lcmsPlot_function <- function(fn_name, default_args = list()) {
  function(...) {
    function(obj) {
      args <- list(...)
      args <- modifyList(default_args, args)
      args <- c(list(obj = obj), args)

      do.call(fn_name, args)
    }
  }
}
