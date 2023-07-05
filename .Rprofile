dllunload <- function(){
  dyn.unload(
    system.file("libs", "x64", "Apollonius.dll", package = "Apollonius")
  )
}

myinstall <- function() {
  if(rstudioapi::isAvailable()) {
    rstudioapi::restartSession(
      "devtools::install(quick = TRUE, keep_source = TRUE)"
    )
  } else {
    try(dllunload())
    try(detach("package:Apollonius", character.only = TRUE))
    devtools::install(quick = TRUE, keep_source = TRUE)
  }
}

mydocument <- function() {
  if(rstudioapi::isAvailable()) {
    rstudioapi::restartSession(
      "roxygen2::roxygenise(load_code = roxygen2::load_installed)"
    )
  } else {
    roxygen2::roxygenise(load_code = roxygen2::load_installed)
  }
}
