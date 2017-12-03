#'  get operation system info
#'
getOS <- function() {
    os = toupper(.Platform$OS.type)
    sysinf = Sys.info()
    if ( !is.null(sysinf) ) {
        os = sysinf['sysname']
        if ( os == 'Darwin' ) {
            os = "OSX"
        } else if ( os == 'Linux' ) {
            os = "LINUX"
        }
    } else {
        if ( grepl("^darwin", R.version$os, perl=T) ) {
            os = "OSX"
        } else if ( grepl("linux-gnu", R.version$os, perl=T) ) {
            os = "LINUX"
        }
    }

    return(os)
}
