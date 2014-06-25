#!/usr/bin/env Rscript

## Aim: does this and that
## choose between:
## Author: Hajar CHOUIKI
## License: 
#  
#  Copyright 2014 INRA-CIRAD
#  
#  This program is free software; you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation; either version 2 of the License, or
#  (at your option) any later version.
#  
#  This program is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.
#  
#  You should have received a copy of the GNU General Public License
#  along with this program; if not, see <http://www.gnu.org/licenses/> or 
#  write to the Free Software Foundation, Inc., 
#  51 Franklin Street, Fifth Floor, Boston,
#  MA 02110-1301, USA.
#  
#  
rm(list=ls())
prog.name <- "GenotypeConcordance2.R"
prog.version <- "1.0"

R.v.maj <- as.numeric(R.version$major)
R.v.min.1 <- as.numeric(strsplit(R.version$minor, "\\.")[[1]][1])
if(R.v.maj < 2 || (R.v.maj == 2 && R.v.min.1 < 15))
    stop("require R >= 2.15 (for paste0)", call.=FALSE)

##' Display the help on stdout
##'
##' The format complies with help2man (http://www.gnu.org/s/help2man) but use --no-discard-stderr
##' @title Help
help <- function(){
  txt <- paste0("`", prog.name, "' does this and that.\n")
  txt <- paste0(txt, "\n")
  txt <- paste0(txt, "Usage: ", prog.name, " [OPTIONS] ...\n")
  txt <- paste0(txt, "\n")
  txt <- paste0(txt, "Options:\n")
  txt <- paste0(txt, " -h, --help\tdisplay the help and exit\n")
  txt <- paste0(txt, " -V, --version\toutput version information and exit\n")
  txt <- paste0(txt, " -v, --verbose\tverbosity level (0/default=1/2/3)\n")
  txt <- paste0(txt, " -n, --input\tpath to the input missing data file\n")
  txt <- paste0(txt, " -jpg, --input\tpath for the plot jpg file\n")
  txt <- paste0(txt, "\n")
  txt <- paste0(txt, "Examples:\n")
  txt <- paste0(txt, " ", prog.name, " -n <NA>\n")
  txt <- paste0(txt, " ", prog.name, " -jpg <jpg>\n")
  txt <- paste0(txt, "\n")
  txt <- paste0(txt, "Report bugs to <>.")
  message(txt)
}

##' Display version and license information on stdout
##'
##' To comply with help2man (http://www.gnu.org/s/help2man) but use --no-discard-stderr
##' @title Version
version <- function(){
  txt <- paste0(prog.name, " ", prog.version, "\n")
  txt <- paste0(txt, "\n")
  ## choose between:
  ## txt <- paste0(txt, "Not copyrighted -- provided to the public domain\n")
  ## or:
  ##txt <- paste0(txt, "Copyright (C) \n")
  txt <- paste0(txt, "License GPLv3+: GNU GPL version 3 or later <http://gnu.org/licenses/gpl.html>\n")
  txt <- paste0(txt, "This is free software; see the source for copying conditions. There is NO\n")
  txt <- paste0(txt, "warranty; not even for MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.\n")
  txt <- paste0(txt, "\n")
  txt <- paste0(txt, "Written by Hajar CHOUIKI.")
  message(txt)
}

##' Parse the command-line arguments
##'
##' Allow short and long options
##' @title Command-line
##' @param params list of parameters initialized with default values
##' @return List of parameters
parseCmdLine <- function(params){
  args <- commandArgs(trailingOnly=TRUE)
  ## print(args)
  
  i <- 0
  while(i < length(args)){ # use "while" loop for options with no argument
    i <- i + 1
    if(args[i] == "-h" || args[i] == "--help"){
      help()
      quit("no", status=0)
    }
    else if(args[i] == "-V" || args[i] == "--version"){
      version()
      quit("no", status=0)
    }
    else if(args[i] == "-v" || args[i] == "--verbose"){
      params$verbose <- as.numeric(args[i+1])
      i <- i + 1
    }
    
    else if(args[i] == "-n" || args[i] == "--NA"){
      params$na.file <- args[i+1]
      i <- i + 1
    }
    else if(args[i] == "-jpg" || args[i] == "--jpg"){
      params$jpg.file <- args[i+1]
      i <- i + 1
    }
    
    else{
      write(paste0(prog.name, ": invalid option -- ", args[i], "\n"), stderr())
      help()
      quit("no", status=1)
    }
    
  }
  
  return(params)
}

##' Check the values of the command-line parameters
##'
##' @param params list of parameters
checkParams <- function(params){
  if(is.null(params$na.file)){
    write("ERROR: missing compulsory option --NA\n", stderr())
    help()
    quit("no", status=1)
  }
  if(! file.exists(params$na.file)){
    write(paste0("ERROR: can't find file ", params$na.file, "\n"), stderr())
    help()
    quit("no", status=1)
  }
  
}

run <- function(params){
  
  ## specific code ...
	 ##print (params$report.file)
	 ##print (params$jpg.file)
	na = read.table(params$na.file, sep="\t", header=TRUE)
	png(params$jpg.file)
	SNPNA <- as.matrix(na, rownames.force = TRUE)
	SNPNA[SNPNA < 1] = NA
	image(x=1:nrow(SNPNA), y=1:ncol(SNPNA), z=is.na(SNPNA), col=c("white","black"), main="Missing values", xlab="Genes", ylab="Samples")
	dev.off()
}

main <- function(){
  params <- list(verbose=1,
                 in.file=NULL)
  
  params <- parseCmdLine(params)
  
  checkParams(params)
  
  if(params$verbose > 0){
    message(paste0("START ", prog.name, " ",
                   format(Sys.time(), "%Y-%m-%d %H:%M:%S")))
    message(paste0("cwd: ", getwd()))
  }
  
  system.time(run(params))
  
  if(params$verbose > 0){
    message(paste0("END ", prog.name, " ",
                   format(Sys.time(), "%Y-%m-%d %H:%M:%S")))
    ## print(object.size(x=lapply(ls(), get)), units="Kb") # return an error I don't understand
  }
}

if(! interactive())
    main()

