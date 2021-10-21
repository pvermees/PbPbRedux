loadTable <- function(fn) {
    out <- list()
    out$tab <- read.csv(system.file(fn,package='PbPbRedux'))
    out$foo <- 'hello'
    out
}

freeformServer <- function(port=NULL) {
  appDir <- R.utils::getAbsolutePath("../inst/www")
  shinylight::slServer(host='0.0.0.0', port=port, appDir=appDir, daemonize=TRUE,
    interface=list(
      loadTable=loadTable
    )
  )
}

freeformServer(8000)
