#!/usr/bin/env Rscript

## Aim: does this and that
## choose between:
## Author: Hajar CHOUIKI
## Not copyrighted -- provided to the public domain
## License: 

rm(list=ls())
prog.name <- "GenotypeConcordance.R"
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
  txt <- paste0(txt, " -i, --input\tpath to the input file\n")
  txt <- paste0(txt, " -r, --input\tpath to the output report file\n")
  txt <- paste0(txt, "\n")
  txt <- paste0(txt, "Examples:\n")
  txt <- paste0(txt, " ", prog.name, " -i <input>\n")
  txt <- paste0(txt, " ", prog.name, " -r <report>\n")
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
    else if(args[i] == "-i" || args[i] == "--input"){
      params$in.file <- args[i+1]
      i <- i + 1
    }
    else if(args[i] == "-r" || args[i] == "--report"){
      params$report.file <- args[i+1]
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
  if((is.null(params$in.file)) || (is.null(params$in.file))){
    write("ERROR: missing compulsory option --input\n", stderr())
    help()
    quit("no", status=1)
  }
  if(! file.exists(params$in.file)){
    write(paste0("ERROR: can't find file ", params$in.file, "\n"), stderr())
    help()
    quit("no", status=1)
  }

  if(is.null(params$report.file)){
    write("ERROR: missing compulsory option --report\n", stderr())
    help()
    quit("no", status=1)
  }
  if(! file.exists(params$report.file)){
    write(paste0("ERROR: can't find file ", params$report.file, "\n"), stderr())
    help()
    quit("no", status=1)
  }
  
}

run <- function(params){
  
  ## specific code ...
	 ##print (params$report.file)
	 ##print (params$jpg.file)
	sink (params$report.file, append = TRUE)
	library(gsalib)
	d = gsa.read.gatkreport(params$in.file)
	cat ("\n")
	dimnames(d$SiteConcordance_Summary)[[2]][5]<-"file1_ONLY"
	dimnames(d$SiteConcordance_Summary)[[2]][6]<-"file2_ONLY"
	cat ("Site Concordance Summary : \n")
	cat ("----------------------------\n\n")
	print (d$SiteConcordance_Summary[c(1,2,3,4,5,6)])
	cat ("\n\n")
	cat ("For strictly bi-allelic VCFs, only the ALLELES_MATCH, EVAL_ONLY, TRUTH_ONLY fields will be populated, but where multi-allelic sites are involved counts for EVAL_SUBSET_TRUTH and EVAL_SUPERSET_TRUTH will be generated.\n")
	cat ("then the site is tabulated as EVAL_SUBSET_TRUTH. Were the situation reversed, it would be EVAL_SUPERSET_TRUTH. \n")
	cat ("For more information : http://www.broadinstitute.org/gatk/gatkdocs/org_broadinstitute_sting_gatk_walkers_variantutils_GenotypeConcordance.html\n")
	cat ("\n\n")
	cat ("Genotype Concordance Summary : \n")
	cat ("-------------------------------\n\n")
	print (d$GenotypeConcordance_Summary[c(1,3,5,6)])
	cat ("\n\n")
	cat ("Genotype Concordance Eval Proportions :\n")
	cat ("----------------------------------------\n\n")
	print (d$GenotypeConcordance_EvalProportions[c(1,4,5,9,11,16,15,20)])
	cat ("----------------------------------------\n\n")
	sink()

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

