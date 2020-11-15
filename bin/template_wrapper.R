#!/usr/bin/env Rscript

############################## ARGUMENTS SECTION #############################
## Collect arguments
args <- commandArgs(TRUE)

## Default setting when no all arguments passed or help needed
if("--help" %in% args | "help" %in% args | (length(args) == 0) | (length(args) == 1) ) {
  cat("
      The helper R Script template_wrapper.R
      Mandatory arguments:
          --input=path              - The path to the input file

          --help                    - you are reading it
          
      Optionnal arguments:

          --a_number=num            - A numeric value to check parameters of type numeric
                                      Default: 12
          --a_name=chr              - A name for checking parameters of type character
                                      Default: 'this_name'
    
      Usage:
      
          The typical command for running the script is as follows:
    
          ./template_wrapper.R --input='example.txt' --a_number=8 

          To see help:
         ./template_wrapper.R --help
      
      WARNING : here put all the things the user has to know
      \n")
  
  q(save="no")
}

## Parse arguments (we expect the form --arg=value)
parseArgs    <- function(x) strsplit(sub("^--", "", x), "=")

argsL        <- as.list(as.character(as.data.frame(do.call("rbind", parseArgs(args)))$V2))
names(argsL) <- as.data.frame(do.call("rbind", parseArgs(args)))$V1
args         <- argsL
rm(argsL)

## Give some value to optional arguments if not provided
if(is.null(args$a_number)) {args$a_number = 12} else {args$a_number=as.numeric(args$a_number)}
if(is.null(args$a_name)) {args$a_name = "this_name"} else {args$a_name=as.character(args$a_name)}

############################## LIBRARIES SECTION #############################

suppressWarnings(suppressMessages(library(stats)))

# ######################### VARIABLES REASSIGNMENT SECTION ###############################

# Facilitates testing and protects from wh-spaces, irregular chars

# required
input    <- args$input

# optional
a_number <- args$a_number
a_name   <- args$a_name

cat("\n")
cat("ARGUMENTS SUMMARY")
cat("\n")
cat("input    : ", input, "\n",sep="")
cat("a_number : ", a_number,   "\n",sep="")
cat("a_name   : ", a_name,        "\n",sep="")

# ############################### SCRIPT SECTION ###############################

