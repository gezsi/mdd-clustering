#!/usr/bin/Rscript 

# check and install packages
list.of.packages <- c("data.table", "optparse", "dplyr")
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if( length(new.packages) > 0 ) 
  install.packages(new.packages)

# load packages
suppressWarnings(suppressMessages(library(data.table, quietly = T)))
suppressWarnings(suppressMessages(library(optparse, quietly = T)))
suppressWarnings(suppressMessages(library(dplyr, quietly = T, warn.conflicts = F)))

bnmcmc.default.path <- "./bnmcmc.exe"
bnmcmc.default.arguments <- "--model-init randomize --param-prior bdeu --parents 4 -q edge --burn-in 2000000 --steps 10000000 --chains 4 --heat-profile linear --swap-profile neighbors --temperature 0.1 --conv-diag-level mbmGeweke --conv-diag-steps 10000 --print-progress 10000 -r 50585086 --data-collapse-identical-rows true --score-cache nodes"

compute.prevalence <- function( dt, diseases )
{
  dt.freq.list <- list()
  # prevalence of diseases
  for( disease in diseases )
  {
    s <- summary( dt[,disease] )
    
    dt.freq.list[[disease]] <- as.data.frame( cbind( disease, t(s) ) ) %>% mutate_if(is.factor, as.character)
  }
  
  dt.freq2 <- bind_rows( dt.freq.list )
  for( i in 2:ncol(dt.freq2) )
    dt.freq2[,i] <- suppressWarnings(as.numeric(as.character(dt.freq2[,i])))
  
  # fix missing 'NA's
  dt.freq2$`NA's`[ is.na(dt.freq2$`NA's`) ] <- 0
  
  dt.freq2$Frequency <- 1.0 - (dt.freq2$`NA's` / nrow(dt))
  dt.freq2$Cases <- nrow(dt) - dt.freq2$`NA's`
  
  return( dt.freq2 )
}

option_list = list(
  make_option( c("-i", "--input"), action = "store", default = NA, type = 'character',
               help = "input file name (csv)" ),
  make_option( c("--age"), action="store", default="Age",
               help="name of variable 'age' [default: %default]" ),
  make_option( c("--sex"), action="store", default="Sex",
               help="name of variable 'sex' [default: %default]" ),
  make_option( c("--household-income"), action="store", default="Income",
               help="name of variable 'household income' [default: %default]" ),
  make_option( c("--depression"), action="store", default="F32",
               help="name of variable 'depression' (if it is composed of multiple variables, then the format should be: comma separated strings without spaces) [default: %default]" ),
  make_option( c("--min-freq"), action="store", default=0.01,
               help="minimum frequency of a disease/medication to include [default: %default]" ),
  make_option( c("--onset-intervals"), action="store", default="20,40,60,70",
               help="onset interval boundaries (format: comma separated numerical values [of age in years] without spaces) [default: %default]" ),
  make_option( c("--number-of-previous-intervals"), action="store", default=1,
               help="number of previous time intervals to consider in each time slice (>=1) [default: %default]" ),
  make_option( c("-r","--remove"), action="store", default=NA, type = 'character',
               help="names of variables to remove (format: comma separated strings without spaces)" ),
  make_option( c("--threads"), action="store", default=20,
               help="number of parallel threads [default: %default]" ) )

opt_parser = OptionParser(option_list=option_list)
opt = parse_args( opt_parser )

if( is.na(opt$input) )
{
  print_help( opt_parser )
  stop( "Input file name can not be empty.", call. = F )
}

age.var <- opt$age
gender.var <- opt$sex
household_income.var <- opt$`household-income`
depression.vars <- unique( strsplit( x = opt$depression, split = ",", fixed = T )[[1]] )
remove.vars <- unique( strsplit( x = opt$remove, split = ",", fixed = T )[[1]] )
min.freq <- opt$`min-freq`
number.of.previous.timeslices <- opt$`number-of-previous-intervals`
threads <- opt$threads

onset.intervals <- as.numeric( strsplit( opt$`onset-intervals`, split = ",", fixed = T )[[1]] )

# check arguments
stopifnot( min.freq >= 0 & min.freq <= 1 )
stopifnot( length(onset.intervals) > 0 )
stopifnot( number.of.previous.timeslices >= 1 )
stopifnot( threads >= 1 & threads <= 512 )
stopifnot( threads == round( threads ) )

# read input file
df <- as.data.frame( fread( opt$input ) )

#### check dataset ####

any.errors <- F

if( !( age.var %in% colnames(df) ) )
{
  warning( paste0( "The age variable ('", age.var, "') is missing from the dataset." ) )
  any.errors <- T
}

if( !( gender.var %in% colnames(df) ) )
{
  warning( paste0( "The sex variable ('", gender.var, "') is missing from the dataset." ) )
  any.errors <- T
}

if( !( household_income.var %in% colnames(df) ) )
{
  warning( paste0( "The houshold income variable ('", household_income.var, "') is missing from the dataset." ) )
  any.errors <- T
}

if( any( !( depression.vars %in% colnames(df) ) ) )
{
  warning( paste0( "The depression variable(s) ('", paste0( depression.vars, collapse = "," ), "') is/are missing from the dataset." ) )
  any.errors <- T
}

if( any.errors )
{
  warning( "Please fix the above issues, and run the script again." )
  stop()
}

any.errors <- F

## check sex
genders <- unique( df[,gender.var] )
if( length( genders ) != 2 ) {
  warning( "Sex should be a binary variable." )
  any.errors <- F
}
if( !( 1 %in% genders ) | !( 2 %in% genders ) ) {
  warning( "Sex should be coded with with values '1' and '2'." )
  any.errors <- F
}
if( any(is.na(genders)) ) {
  warning( "Sex contains NA values." )
  any.errors <- T
}

## check income
incomes <- unique( df[,household_income.var] )
if( length( incomes ) != 3 ) {
  warning( "Household income should be a ternary variable." )
  any.errors <- T
}
if( !( 1 %in% incomes ) | !( 2 %in% incomes ) | !( 3 %in% incomes ) ) {
  warning( "Household income should be coded with values '1', '2' and '3'." )
  any.errors <- T
}
if( any(is.na(incomes)) ) {
  warning( "Household income contains NA values. You might want to impute missing values (e.g. by random values based on the household income distribution stratified by sex and age group)." )
}

## check age
if( any(is.na(df[,age.var])) ) {
  warning( "Age contains NA values. You might want to set the age to the highest disease onset of each patient whose age is missing." )
  any.errors <- T
}
if( min(df[,age.var], na.rm = T) <= 0 ) {
  warning( "Age contains zero or negative values." )
  any.errors <- T
}

## check whether there are negative values
sum.negatives <- sum( df < 0, na.rm=T )
if( sum.negatives > 0 )
{
  warning( paste0( "There are ", sum.negatives, " negative values in the dataset. Please check your data." ) )
  any.errors <- T
}

# we assume that all other variable are diseases
vars <- setdiff( colnames( df ), c( age.var, gender.var, household_income.var ) )

# remove variables if needed
if( length(remove.vars) > 0 )
  vars <- setdiff( vars, remove.vars )

# clean dataset
df <- df[,c( age.var, gender.var, household_income.var, vars )]

df[,household_income.var] <- as.factor(df[,household_income.var])

# check the relation between age and disease onsets
max.onsets <- do.call( pmax, c( df[,vars], na.rm = TRUE ) )
wrong.ages <- df[,age.var] < max.onsets
sum.wrong.ages <- sum( wrong.ages, na.rm=T )
if( sum.wrong.ages > 0 )
{
  warning( paste0( "Age contains lower values than the highest onset in case of ", sum.wrong.ages, " patients. Please check your data. You might want to set the age to the highest disease onset of these patients." ) )
  any.errors <- T
}

if( any.errors )
{
  stop( "Please fix the above issues, and run the script again." )
}

#### ####

## compute disease prevalences
message( "Computing disease prevalences..." )

dt.freq <- compute.prevalence( df, 
                               diseases = vars )

dt.freq.depression <- compute.prevalence( df[ which( rowSums( !is.na(df[,depression.vars,drop=F]) ) > 0 ), ], 
                                          diseases = vars )


frequent.vars <- unique( na.omit( c( dt.freq$disease[ dt.freq$Frequency >= min.freq ], 
                                     dt.freq.depression$disease[ dt.freq.depression$Frequency >= min.freq ] ) ) )


message( paste0( "Number of variables that are more frequent than ", min.freq, " in the whole population or in the depressed subpopulation: ", length( frequent.vars ) ) )

bnmcmc.script <- "#!/bin/bash\n\n"

get.model.file.node <- function( df.ts, var )
{
  unique.values = length(unique(df.ts[,var]))
  
  node <- paste0( "<chanceNode id=\"",var, "\" type=\"", var, "\"><parents></parents><cpd><table><distribution values=\"", paste( rep(1/unique.values, times=unique.values), collapse="," ), "\"/></table></cpd></chanceNode>\n" )
  
  return( node )  
}

get.model.file.type <- function( df.ts, var )
{
  unique.values = unique(df.ts[,var])
  
  type <- paste0( "<enumType id=\"", var, "\"><values>" )
  for( i in unique.values[order(unique.values)] )
    type <- paste0( type, paste0( "<value id=\"",i,"\" lowerBound=\"",i,"\" upperBound=\"",i,"\"/>" ) )
  type <- paste0( type, "</values></enumType>\n" )
  
  return( type )  
}

# creating the time slice datasets
for( i in 1:length(onset.intervals) )
{
  age.start <- 1e-5 
  age.end <- onset.intervals[i]
  
  print( paste0( "interval ", i, " => [", age.start, "-", age.end, "[" ) )
  
  df.filtered <- df[ df[,age.var] >= age.end, ]
  
  df.ts <- if( length(unique(df.filtered[,gender.var])) == 2 ) {
    df.filtered[, c(gender.var, household_income.var) ] 
  } else {
    df.filtered[, household_income.var, drop=F ] 
  }
  
  prefix.current <- paste0("interval", i, ".")
  for( v in frequent.vars )
  {
    if( !( v %in% colnames(df.filtered) ) )
    {
      message( paste0( "Variable '", v, "' is not present in the dataset." ) )
      next
    }
    
    indices <- which( !is.na(df.filtered[,v]) & df.filtered[,v] >= age.start & df.filtered[,v] < age.end )
    
    if( length(indices) == 0 )
      next
    
    v.new <- paste0( prefix.current, v )
    df.ts[ , v.new ] <- 0
    df.ts[ indices, v.new ] <- 1 
  }
  
  previous.filename <- ""
  if( i >= 1 ) 
  {
    pre.age.start <- 1e-5
    pre.age.end <- onset.intervals[i-1]
    
    print( paste0( "pre interval => [", pre.age.start, "-", pre.age.end, "[" ) )
    
    vars.pre.decade <- c()
    
    if( length(unique(df.filtered[,gender.var])) == 2 ) {
      vars.pre.decade <- c(gender.var, household_income.var)
    } else {
      vars.pre.decade <- c(household_income.var)
    }
    
    s.index <- max( 1, i-number.of.previous.timeslices )
    e.index <- i-1
    
    for( index in c(s.index:e.index) )
    {
      prefix.pre <- paste0("interval", index, ".")
      
      for( v in frequent.vars )
      {
        if( !( v %in% colnames(df.filtered) ) )
        {
          message( paste0( "Variable '", v, "' is not present in the dataset." ) )
          next
        }
        
        indices <- which( !is.na(df.filtered[,v]) & df.filtered[,v] >= pre.age.start & df.filtered[,v] < pre.age.end )
        
        if( length(indices) == 0 )
          next
        
        v.pre <- paste0( prefix.pre, v )
        df.ts[ , v.pre ] <- 0
        df.ts[ indices, v.pre ] <- 1
        
        vars.pre.decade <- c( vars.pre.decade, v.pre )
      }
    }
    
    previous.filename <- paste0("previous.timeslices.vars.in_interval-",i,".txt")
    write.table( vars.pre.decade, previous.filename, quote = F, sep = "\t", row.names = F, col.names = F )
  }
  
  message( paste0( "Number of rows in time-slice [", age.start, ",", age.end, "[ : ", nrow(df.ts) ) )
  
  dataset.filename <- paste0( "dataset_ihDBN", "_interval-", i, "_freq", min.freq, ".csv" )
  
  write.table( df.ts, file = dataset.filename, quote = F, sep = ",", row.names = F, col.names = T )
  
  # create xml model file
  model.filename <- paste0( dataset.filename, ".xml" )
  
  sink( model.filename )
  cat( "<?xml version=\"1.0\" encoding=\"UTF-8\" standalone=\"yes\"?>\n" )
  cat( paste0( "<model id=\"", dataset.filename, "\">\n<nodes>\n" ) )
  for( v in colnames( df.ts ) )
    cat( get.model.file.node( df.ts, v ) )
  cat( "</nodes>\n<types>\n" )
  for( v in colnames( df.ts ) )
    cat( get.model.file.type( df.ts, v ) )
  cat( "</types><groups></groups></model>" )
  sink()
  
  targets <- c()
  if( length( depression.vars ) == 1 )
  {
    targets <- paste0( "-q MBG ", prefix.current, depression.vars )
  } else {
    for( depression.var in depression.vars )
    {
      targets <- c( targets, paste0( "-q MBG ", prefix.current, depression.var ) )
    }
    targets <- c( targets, paste0( "-q MBG ", paste0( prefix.current, depression.vars, collapse="," ) ) )
  }
  targets <- paste0( targets, collapse = " " )
  
  # append bnmcmc script
  bnmcmc.script <- paste0( bnmcmc.script, paste( bnmcmc.default.path,
                                                 paste0( "--data-parallel ", threads, " ", dataset.filename ),
                                                 paste0( "--model ", model.filename ),
                                                 paste0( "--output ", "bdmm_output_", i ),
                                                 ifelse( previous.filename != "", paste0( "--is-interventional-file ", previous.filename ), "" ),
                                                 targets, 
                                                 bnmcmc.default.arguments,
                                                 "&\n",
                                                 sep = " " ) )
  
  
  
  
}

sink( "run_BDMM.sh" )
cat( bnmcmc.script )
sink()

