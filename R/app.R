loadTable <- function(fn) {
    read.csv(system.file(fn,package='PbPbRedux'))
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
