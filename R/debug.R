DebugTrace <- setRefClass("DebugTrace",
  methods = list(
    dbg_trace = function(...) {
      class_ref <- getRefClass()      
      method_name <- strsplit(as.character(sys.call(sys.parent())), "\\$")[[1]][2]
      message(paste0(class_ref@className[1], "::", method_name, " args = ", paste(..., sep = ", ", collapse = ", ")))
    }
  )
)


