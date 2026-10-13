# runs the expression and returns its value together with all messages and
# warnings emitted by toolStatusMessage (which uses warnings for "warn" status)
captureConditions <- function(expr) {
  conditions <- character()
  value <- withCallingHandlers(
    expr,
    message = function(m) {
      conditions <<- unique(c(conditions, conditionMessage(m)))
      invokeRestart("muffleMessage")
    },
    warning = function(w) {
      conditions <<- unique(c(conditions, conditionMessage(w)))
      invokeRestart("muffleWarning")
    }
  )
  return(list(value = value, conditions = conditions))
}
