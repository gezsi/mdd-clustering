#!/usr/bin/Rscript 

# check and install packages
list.of.packages <- c("data.table", "optparse", "dplyr", "openxlsx", "pdist")
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if( length(new.packages) > 0 ) 
  install.packages(new.packages)

# load packages
suppressWarnings(suppressMessages(library(data.table, quietly = T)))
suppressWarnings(suppressMessages(library(optparse, quietly = T)))
suppressWarnings(suppressMessages(library(openxlsx, quietly = T)))
suppressWarnings(suppressMessages(library(dplyr, quietly = T, warn.conflicts = F)))
suppressWarnings(suppressMessages(library(pdist, quietly = T)))

## functions

create.profile.matrix <- function( interval, comorbidity.profile, weighted = T )
{
  condition <- comorbidity.profile$Interval == interval & !is.na( comorbidity.profile$Disease )
  
  P <- comorbidity.profile$P[ condition ]

  if( !weighted ) { P[ P > 0 ] <- 1.0 }

  m <- matrix( P, ncol = 1, dimnames = list( comorbidity.profile$Variable[ condition ],
                                             paste( "ALL_", 
                                                    ifelse( weighted, "Weighted", "NotWeighted" ), 
                                                    sep = "_" ) ) )
  
  colnames(m) <- paste0( "interval", interval, "_", colnames(m) )
  
  return( m )
}


## read options

option_list = list(
  make_option( c("-i", "--input"), action = "store", default = NA, type = 'character',
               help = "input file name (csv)" ),
  make_option( c("-o", "--output"), action = "store", default = ".", type = 'character',
               help = "output directory name [default: %default]" ),
  make_option( c("--cluster-definitions"), action="store", default=NA, type = 'character',
               help="cluster definition file name [default: %default]" ),
  make_option( c("--id"), action = "store", default = NA, type = 'character',
               help = "name of id column (if set, it will be copied to the result files to identify subjects) [default: %default]" ),
  make_option( c("--age"), action="store", default="Age",
               help="name of variable 'age' [default: %default]" ),
  make_option( c("--onset-intervals"), action="store", default="20,40,60,70",
               help="onset interval boundaries (format: comma separated numerical values [of age in years] without spaces) [default: %default]" ) )

opt_parser = OptionParser(option_list=option_list)
opt = parse_args( opt_parser )

if( is.na(opt$input) )
{
  print_help( opt_parser )
  stop( "Input file name can not be empty.", call. = F )
}

if( is.na(opt$`cluster-definitions`) )
{
  print_help( opt_parser )
  stop( "Cluster definitions file name can not be empty.", call. = F )
}

input.file <- opt$input
output.dir <- opt$output
cluster.definitions.file <- opt$`cluster-definitions`
id.column <- opt$id
age.var <- opt$age
onset.intervals <- as.numeric( strsplit( opt$`onset-intervals`, split = ",", fixed = T )[[1]] )

# check arguments
stopifnot( length(onset.intervals) > 0 )

if( output.dir != "." )
{
  dir.create( output.dir, showWarnings = F, recursive = T )
}

#input.file <- "e:/Munka/Trajectome/2021-UKB-Trajectome/ukb42450_disease_onset_imputedIncome_fixedSex_consensual_variables.csv"
#cluster.definitions.file <- "e:/Munka/Trajectome/ConsensualClustering_2021nov/cluster-definitions-all-7.xlsx"

message( "Reading input file." )

# read input file
df <- fread( input.file, data.table = T )

comorbidity.profile <- openxlsx::read.xlsx( cluster.definitions.file, sheet = 1 )
cluster.definitions <- openxlsx::read.xlsx( cluster.definitions.file, sheet = 2 )
cluster.definitions.scaling <- openxlsx::read.xlsx( cluster.definitions.file, sheet = 3 )

# check comorbidity profile
stopifnot( "Comorbidity profile is not set." = !is.null( comorbidity.profile ) )

# check id column
if( !is.na(id.column) )
{
  stopifnot( "ID column is not found in input file." = id.column %in% colnames(df) )
}

message( "Creating transformed dataset." )

if( is.na( id.column ) )
{
  id.column = "rowId"
  df[[id.column]] <- 1:nrow(df) 
} 

select.cols <- c( id.column )

df.transformed.to.intervals <- df[,..select.cols]

# creating the time slice datasets
for( i in 1:length(onset.intervals) )
{
  age.start <- 1e-5 
  age.end <- onset.intervals[i]
  
  prefix.current <- paste0("interval", i, ".")
  for( v in colnames(df) )
  {
    variable.current <- paste0( prefix.current, v )
    if( !( variable.current %in% comorbidity.profile$Variable ) ) #### TODO: depression!
      next
    
    indices.of.proper.age <- which( df[,..age.var] >= age.end )
    # fill those values what we know, but not necessarily complete (i.e. the person may have been censored during the period)
    indices.of.onset <- which( !is.na(df[,..v]) & df[,..v] >= age.start & df[,..v] < age.end )  # df[,..age.var] >= age.end & 
    
    df.transformed.to.intervals[indices.of.proper.age, (variable.current) := 0L]
    df.transformed.to.intervals[indices.of.onset, (variable.current) := 1L]
    #df.transformed.to.intervals[, (variable.current) := lapply(.SD, as.factor), .SDcols = variable.current ]
  }
}

#### Compute patient profiles ####

message( "Computing patient profiles." )

select.cols <- c( id.column )

patient.profiles <- df.transformed.to.intervals[,..select.cols]

weighted <- T
for( i in 1:length(onset.intervals) )
{
    message( paste0( "Interval ", i ) )
    
    profile.matrix <- create.profile.matrix( i, comorbidity.profile, weighted = weighted )
    
    common.variables <- intersect( rownames(profile.matrix), colnames(df.transformed.to.intervals) )
    if( length( setdiff( rownames(profile.matrix), common.variables ) ) > 0 )
    {
      warning( paste0( "In interval ", i, ", certain variables are missing from the dataset. This should not happen normally." ) )
    }
    
    df.transformed.to.intervals.filtered.to.profile <- as.matrix( df.transformed.to.intervals[ , ..common.variables ] )
    profile.matrix <- profile.matrix[ common.variables, , drop=F ]
    
    patient.profiles <- cbind( patient.profiles, 
                               df.transformed.to.intervals.filtered.to.profile %*% profile.matrix )
}

# This will be used for clustering!
m <- matrix( data = NA, 
             nrow = nrow( patient.profiles ), 
             ncol = 4, dimnames = list( rownames( patient.profiles ),
                                        c("interval1__0W","interval2_vs_1__0W","interval3_vs_2__0W","interval4_vs_3__0W" ) ) )

m[,"interval1__0W"] <- patient.profiles$interval1_ALL__Weighted
m[,"interval2_vs_1__0W"] <- patient.profiles$interval2_ALL__Weighted - patient.profiles$interval1_ALL__Weighted
m[,"interval3_vs_2__0W"] <- patient.profiles$interval3_ALL__Weighted - patient.profiles$interval2_ALL__Weighted
m[,"interval4_vs_3__0W"] <- patient.profiles$interval4_ALL__Weighted - patient.profiles$interval3_ALL__Weighted

m.scaled <- t( ( t(m) - cluster.definitions.scaling$Centers ) / cluster.definitions.scaling$Scales )

pdists <- as.matrix( pdist( m.scaled, as.matrix( cluster.definitions ) ) )

clusters <- apply( pdists, MARGIN = 1, FUN = function(x) { which.min(x) } )

message( "Cluster sizes:" )
print( table(clusters) )

# 1. Cluster indices
output.clusters <- data.frame( id = patient.profiles[[ id.column ]],
                               cluster = clusters )
colnames(output.clusters)[1] <- id.column
output.clusters.file <- paste0( output.dir, "/", "cluster_indices.txt" )

write.table( x = output.clusters,
             file = output.clusters.file,
             quote = F, append = F, sep = "\t", row.names = F, col.names = T )

message( paste0( "Cluster indices: ", output.clusters.file ) )

# 2. Cluster probabilities
pdists <- exp(-pdists)
pdists <- pdists / rowSums(pdists)

output.cluster.probs <- cbind( data.frame( id = patient.profiles[[ id.column ]] ),
                               pdists )
colnames(output.cluster.probs)[1] <- id.column
output.cluster.probs.file <- paste0( output.dir, "/", "cluster_probabilities.txt" )

write.table( x = output.cluster.probs,
             file = output.cluster.probs.file,
             quote = F, append = F, sep = "\t", row.names = F, col.names = T )

message( paste0( "Cluster probabilities: ", output.cluster.probs.file ) )

# 2.b Cluster probabilities log-odds
pdists.prob.logodds <- log( pdists / (1.0 - pdists) )

output.cluster.probs.logodds <- cbind( data.frame( id = patient.profiles[[ id.column ]] ),
                                       pdists.prob.logodds )
colnames(output.cluster.probs.logodds)[1] <- id.column
output.cluster.probs.logodds.file <- paste0( output.dir, "/", "cluster_probabilities_logodds.txt" )

write.table( x = output.cluster.probs.logodds,
             file = output.cluster.probs.logodds.file,
             quote = F, append = F, sep = "\t", row.names = F, col.names = T )

message( paste0( "Cluster probabilities (log-odds): ", output.cluster.probs.logodds.file ) )


# 3. Cluster raw variables
output.clusters.raw <- cbind( data.frame( id = patient.profiles[[ id.column ]] ),
                          as.data.frame( m ) )
colnames(output.clusters.raw)[1] <- id.column
output.clusters.raw.file <- paste0( output.dir, "/", "cluster_variables.txt" )

write.table( x = output.clusters.raw,
             file = output.clusters.raw.file,
             quote = F, append = F, sep = "\t", row.names = F, col.names = T )

message( paste0( "Cluster variables (raw values): ", output.clusters.raw.file ) )

message( "Done." )



