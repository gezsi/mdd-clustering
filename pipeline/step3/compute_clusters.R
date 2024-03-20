#!/usr/bin/Rscript 

# check and install packages
list.of.packages <- c("data.table", "optparse", "dplyr", "openxlsx", "pdist", "fst", "ggplot2", "GGally", "gridExtra", "cluster", "viridis", "tidyr", "factoextra", "MASS")
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if( length(new.packages) > 0 ) 
  install.packages(new.packages)

# load packages
suppressWarnings(suppressMessages(library(data.table, quietly = T)))
suppressWarnings(suppressMessages(library(optparse, quietly = T)))
suppressWarnings(suppressMessages(library(openxlsx, quietly = T)))
suppressWarnings(suppressMessages(library(dplyr, quietly = T, warn.conflicts = F)))
suppressWarnings(suppressMessages(library(tidyr, quietly = T)))
suppressWarnings(suppressMessages(library(pdist, quietly = T)))
suppressWarnings(suppressMessages(library(ggplot2, quietly = T)))
suppressWarnings(suppressMessages(library(GGally, quietly = T)))
suppressWarnings(suppressMessages(library(gridExtra, quietly = T)))
suppressWarnings(suppressMessages(library(cluster, quietly = T)))
suppressWarnings(suppressMessages(library(viridis, quietly = T)))
suppressWarnings(suppressMessages(library(fst, quietly = T)))
suppressWarnings(suppressMessages(library(factoextra, quietly = T)))

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

read.bdmm.results <- function( directory, targets )
{
  dts.list <- list()
  
  for( i in 1:4 )
  {
    target <- paste0( paste0( paste0( "interval", i, "." ), targets ), collapse = "-" )
    filename <- paste0( directory, paste0( "bdmm_output_", i, "_query-",target,".mbm.csv" ) )
    
    if( !file.exists( filename ) )
    {
      message( paste0( "Looking for '", filename, "'." ) )
      stop( paste0( "The given directory '", directory, "' doesn't contain the appropriate MBM files, or they are empty. Is this the right directory?" ), call. = F )
    }
    
    dt <- read.table( filename, header = F, sep = ";" ) %>%
      rename( "P" = "V1", "Variable" = "V2" ) %>%
      tidyr::separate( col = Variable, into = c("IntervalOfDisease", "Disease"), sep = "[.]", remove = F, fill = "right" ) %>%
      filter( !( Disease %in% targets ) ) %>%
      mutate( Interval = i ) %>%
      select( Variable, IntervalOfDisease, Disease, Interval, P )
    
    dts.list[[ i ]] <- dt
  }
  
  return( bind_rows(dts.list) )
}

is.fst.file <- function( filename )
{
  exts <- strsplit(basename(filename), split="\\.")[[1]]
  return( exts[length(exts)] == "fst" )
}

my_fn <- function(data, mapping, N=100, ...){
  
  get_density <- function(x, y, n ) {
    dens <- MASS::kde2d(x = x, y = y, n = n)
    ix <- findInterval(x, dens$x)
    iy <- findInterval(y, dens$y)
    ii <- cbind(ix, iy)
    return(dens$z[ii])
  }
  
  X <- eval_data_col(data, mapping$x)
  Y <- eval_data_col(data, mapping$y)
  
  data$density <- get_density(x=X, y=Y, n=N)
  
  p <- ggplot(data, mapping) +
    geom_point(aes(colour=density), ...) +
    viridis::scale_color_viridis()      
  p
}

pca.plot <- function( m, clusters )
{
  m.pca <- as.matrix( m ) #,
  ir.pca <- prcomp(m.pca,
                   center = TRUE, 
                   scale. = TRUE) 
  
  #print( ir.pca )
  #plot(ir.pca, type = "l")
  print( summary(ir.pca) )
  # plot(ir.pca$x[, 1], ir.pca$x[, 2], main = "PCA", xlab = "PC1", ylab = "PC2")
  
  var.explained <- ir.pca$sdev*ir.pca$sdev / sum(ir.pca$sdev*ir.pca$sdev)
  
  
  # create data frame with correlations between variables and PCs
  correlations = as.data.frame(cor(m.pca, ir.pca$x))
  
  # data frame with arrows coordinates
  arrows = data.frame(x1 = rep(0,nrow(correlations)), y1 = rep(0,nrow(correlations)), x2 = correlations$PC1, 
                      y2 = correlations$PC2)
  
  
  # create data frame with scores
  scores = as.data.frame(ir.pca$x)
  scores$Cluster <- factor(clusters)
  scores <- plyr::count( scores )
  
  ggplot() + 
    #geom_path( data = scores[ res2$order, ], aes(x = PC1, y = PC2), alpha = 0.3 ) +
    #geom_path(data = scores[ scores$Group != 0,], aes(x = PC1, y = PC2, color = Group, group = Group), alpha = 0.4 ) +
    geom_point(data = scores, aes(x = PC1, y = PC2, color = Cluster, size = freq), alpha = 0.8 ) +
    geom_segment(data = arrows, aes(x = x1, y = y1, xend = x2 * 6, yend = y2 * 6), colour = "gray65") + 
    geom_text(data = correlations, aes(x = PC1 * 6, y = PC2 * 6, label = rownames(correlations)), size=3) + 
    geom_hline(yintercept = 0, colour = "gray65") + 
    geom_vline(xintercept = 0, colour = "gray65") +
    labs( x = paste0( "PC1 axis (", format(var.explained[1]*100, digits = 4), "%)" ), 
          y = paste0( "PC2 axis (", format(var.explained[2]*100, digits = 4), "%)" ) )
}


## read options

option_list = list(
  make_option( c("-i", "--input"), action = "store", default = NA, type = 'character',
               help = "input file name (csv)" ),
  make_option( c("-o", "--output"), action = "store", default = "clusters", type = 'character',
               help = "output directory name prefix (it will be appended by the number of clusters) [default: %default]" ),
  make_option( c("--cohort-profile"), action="store", default=NA, type = 'character',
               help="cohort specific comorbidity profile definition based on BDMM analysis results, name of directory [default: %default]" ),
  make_option( c("--number-of-clusters"), action="store", default=NA, type = 'integer',
               help="number of clusters [default: %default]" ),
  make_option( c("--age"), action="store", default="Age",
               help="name of variable 'age' [default: %default]" ),
  make_option( c("--sex"), action="store", default="Sex",
               help="name of variable 'sex' [default: %default]" ),
  make_option( c("--household-income"), action="store", default="Income",
               help="name of variable 'household income' [default: %default]" ),
  make_option( c("--depression"), action="store", default="F32",
               help="name of variable 'depression' (if it is composed of multiple variables, then the format should be: comma separated strings without spaces) [default: %default]" ),
  make_option( c("--codetable"), action="store", default=NA, type = 'character',
               help="Filename for code table containing disease names for each disease code (columns of the input file) [default: %default]" ),
  make_option( c("--id"), action = "store", default = NA, type = 'character',
               help = "name of id column (if set, it will be copied to the result files to identify subjects) [default: %default]" ),
  make_option( c("--onset-intervals"), action="store", default="20,40,60,70",
               help="onset interval boundaries (format: comma separated numerical values [of age in years] without spaces) [default: %default]" ) )

opt_parser = OptionParser(option_list=option_list)
opt = parse_args( opt_parser )

if( is.na(opt$input) )
{
  print_help( opt_parser )
  stop( "Input file name can not be empty.", call. = F )
}

if( is.na(opt$`cohort-profile`) )
{
  print_help( opt_parser )
  stop( "Cohort profile directory name can not be empty. This should be a directory of a previous BDMM analysis.", call. = F )
}

if( !is.na(opt$codetable) ) {
	codetable <- openxlsx::read.xlsx( opt$codetable, sheet = 1, colNames = T )
} else {
	codetable <- NA
}

input.file <- opt$input
output.dir <- opt$output
cohort.profile.dir <- paste0( opt$`cohort-profile`, "/" )
number.of.clusters <- as.integer( opt$`number-of-clusters` )
id.column <- opt$id
age.var <- opt$age
gender.var <- opt$sex
household_income.var <- opt$`household-income`
depression.vars <- unique( strsplit( x = opt$depression, split = ",", fixed = T )[[1]] )
onset.intervals <- as.numeric( strsplit( opt$`onset-intervals`, split = ",", fixed = T )[[1]] )

# check arguments
stopifnot( length(onset.intervals) > 0 )

stopifnot( "Number of clusters is not set." = !is.na(number.of.clusters) )
stopifnot( "Number of clusters is not appropriate." = number.of.clusters > 0 )

if( output.dir == "." ) {
  output.dir <- paste0( "_", number.of.clusters, "/" )
} else {
  output.dir <- paste0( output.dir, "_", number.of.clusters, "/" )
}

if( output.dir != "." )
{
  dir.create( output.dir, showWarnings = F, recursive = T )
}

## Cohort-specific comorbidity profile ##

message( "Reading cohort-specific comorbidity profile." )

comorbidity.profile <- read.bdmm.results( cohort.profile.dir, depression.vars )

stopifnot( "Comorbidity profile is not set." = !is.null( comorbidity.profile ) )

## Input file

if( is.fst.file( input.file ) ) 
{
  message( "Reading input summary file." )
  
  df.transformed.to.intervals <- fst::read.fst( input.file, as.data.table = F ) 
  
  # deaggregate dataset
  df.transformed.to.intervals <- df.transformed.to.intervals[ rep(row.names(df.transformed.to.intervals), df.transformed.to.intervals$freq), 
                                                              setdiff( colnames(df.transformed.to.intervals), "freq" ) ] %>%
    mutate( across(where(is.factor), function(x) { as.numeric(as.character(x)) } ) )
  
  setDT( df.transformed.to.intervals )
  
  if( is.na( id.column ) )
  {
    id.column = "rowId"
    df.transformed.to.intervals[[id.column]] <- 1:nrow(df.transformed.to.intervals) 
  } 
  
} else {
  
  message( "Reading input file." )
  
  df <- fread( input.file, data.table = T ) 
 
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
      if( !( variable.current %in% comorbidity.profile$Variable ) & !( v %in% depression.vars ) )
        next
      
      indices.of.proper.age <- which( df[,..age.var] >= age.end )
      # fill those values what we know, but not necessarily complete (i.e. the person may have been censored during the period)
      indices.of.onset <- which( !is.na(df[,..v]) & df[,..v] >= age.start & df[,..v] < age.end )  # df[,..age.var] >= age.end & 
      
      df.transformed.to.intervals[indices.of.proper.age, (variable.current) := 0L]
      df.transformed.to.intervals[indices.of.onset, (variable.current) := 1L]
      #df.transformed.to.intervals[, (variable.current) := lapply(.SD, as.factor), .SDcols = variable.current ]
    }
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

patient.profiles$MaxInterval <- NA
patient.profiles$MaxInterval[ which( !is.na( patient.profiles$interval1_ALL__Weighted ) ) ] <- 1
patient.profiles$MaxInterval[ which( !is.na( patient.profiles$interval2_ALL__Weighted ) ) ] <- 2
patient.profiles$MaxInterval[ which( !is.na( patient.profiles$interval3_ALL__Weighted ) ) ] <- 3
patient.profiles$MaxInterval[ which( !is.na( patient.profiles$interval4_ALL__Weighted ) ) ] <- 4

# This will be used for clustering!
patient.profiles$interval1__0W <- patient.profiles$interval1_ALL__Weighted
patient.profiles$interval2_vs_1__0W <- patient.profiles$interval2_ALL__Weighted - patient.profiles$interval1_ALL__Weighted
patient.profiles$interval3_vs_2__0W <- patient.profiles$interval3_ALL__Weighted - patient.profiles$interval2_ALL__Weighted
patient.profiles$interval4_vs_3__0W <- patient.profiles$interval4_ALL__Weighted - patient.profiles$interval3_ALL__Weighted

patient.profiles$WeightedBurden <- 5*patient.profiles$interval4_ALL__Weighted +  # from 65 till 70
  20*patient.profiles$interval3_ALL__Weighted + # from 50 till 70
  40*patient.profiles$interval2_ALL__Weighted + # from 30 till 70
  60*patient.profiles$interval1_ALL__Weighted   # from 10 till 70

patients.with.full.life.trajectories <- which( patient.profiles$MaxInterval == 4 )
patients.of.interest <- patients.with.full.life.trajectories

pp <- patient.profiles[ patients.of.interest, ]
df.patients.of.interest <- df.transformed.to.intervals[ patients.of.interest, ]

m <- as.matrix( pp[, c("interval1__0W","interval2_vs_1__0W","interval3_vs_2__0W","interval4_vs_3__0W" ) ] )

# SCALING!!!
m.scaled <- scale( m )
scaling.centers <- attr( m.scaled, "scaled:center")
scaling.scales <- attr( m.scaled, "scaled:scale")

NUMBER.OF.CLUSTERS = number.of.clusters
directory <- output.dir
resample = T

# Do the STUFF
if(T)
{
  message( "Computing clusters." )
  clustering <- kmeans( m.scaled, centers = NUMBER.OF.CLUSTERS, iter.max = 100, nstart = 25 )
  
  message( "Creating plots." )
  
  # reorder clustering index by burden
  cc <- clustering$centers
  cc <- as.data.frame( t( t(cc) * scaling.scales + scaling.centers ) )
  cc$interval1_ALL__Weighted <- cc$interval1__0W
  cc$interval2_ALL__Weighted <- cc$interval2_vs_1__0W + cc$interval1_ALL__Weighted
  cc$interval3_ALL__Weighted <- cc$interval3_vs_2__0W + cc$interval2_ALL__Weighted
  cc$interval4_ALL__Weighted <- cc$interval4_vs_3__0W + cc$interval3_ALL__Weighted
  cc$burden <- 5*cc$interval4_ALL__Weighted + 
    20*cc$interval3_ALL__Weighted +
    40*cc$interval2_ALL__Weighted +
    60*cc$interval1_ALL__Weighted
  cc$Clustering <- order(order(cc$burden))
  cc$Clustering <-  factor( cc$Clustering )
  
  write.xlsx( list( "scaled" = clustering$centers[order(cc$Clustering),],
                    "backscaled" = cc,
                    "scaling_parameters" = data.frame( Centers = scaling.centers,
                                                       Scales = scaling.scales ) ), 
              file = paste0( directory, "cluster_centers-", NUMBER.OF.CLUSTERS, ".xlsx" ),
              overwrite = T ) 
  
  REAL.CLUSTER.INDEX <- order(order(cc$burden))[ clustering$cluster ]
  
  m.sample.rows <- if( resample ) { sample(1:nrow(m), size = min(20000, nrow(m))) } else { 1:nrow(m) }
  m.sample <- m.scaled[ m.sample.rows, ]
  clustering.sample <- REAL.CLUSTER.INDEX[ m.sample.rows ]
  
  ## Creating silhouette plot
  
  s <- silhouette( clustering.sample, dist( m.sample, method = "euclidean" ) )
  summary(s)
  
  ggsave( fviz_silhouette( s ) + 
            scale_colour_viridis_d() + 
            scale_fill_viridis_d(), 
          filename = paste0( directory, "silhouette_plot_", NUMBER.OF.CLUSTERS, ".png" ), 
          device = "png", width = 297, height = 210, units = "mm", dpi = 300 )
  
  # This is really nice, we had better not showing it
  ggsave( pca.plot( m.sample, clusters = clustering.sample ) +
            scale_colour_viridis_d(),
          filename = paste0( directory, "PCA_plot_", NUMBER.OF.CLUSTERS, ".png" ), 
          device = "png", width = 260, height = 210, units = "mm", dpi = 300, scale = 0.9 )
  
  # Write cluster index back to the original tables
  pp$Clustering <- factor( REAL.CLUSTER.INDEX )
  clustering.variable <- paste0( "Clustering_", NUMBER.OF.CLUSTERS )
  df.patients.of.interest[[clustering.variable]] <- factor( REAL.CLUSTER.INDEX )
  
  message( "Cluster sizes:" )
  print( table( REAL.CLUSTER.INDEX ) )
  
  # Cluster centers for each cohort
  clustering.centers <- merge( cc,
                               data.frame( Clustering = 1:nrow(clustering$centers) ),
                               by = "Clustering" ) %>%
    tidyr::pivot_longer( cols = starts_with("interval"), 
                         names_to = c("interval","type"),
                         names_pattern = "interval(.*)__(.*)",
                         values_to = "Value" )
  clustering.centers$Interval <- clustering.centers$interval
  clustering.centers$Interval[clustering.centers$interval == "1_ALL"] <- "1"
  clustering.centers$Interval[clustering.centers$interval == "2_ALL"] <- "2"
  clustering.centers$Interval[clustering.centers$interval == "3_ALL"] <- "3"
  clustering.centers$Interval[clustering.centers$interval == "4_ALL"] <- "4"
  
  g <- GGally::ggpairs( cbind( as.data.frame(m.sample), 
                               Clustering = factor(clustering.sample) ),
                        columns = 1:4,  
                        ggplot2::aes(colour = Clustering, alpha = 0.7),
                        upper = "blank") + 
    scale_colour_viridis_d()
  
  ggsave( g, filename = paste0( directory, "scatterplot_of_clusterVariables_", NUMBER.OF.CLUSTERS, ".png" ), 
          device = "png", width = 297, height = 210, units = "mm", dpi = 150 )
  
  pdists <- as.matrix( pdist( m.scaled[, c("interval1__0W","interval2_vs_1__0W","interval3_vs_2__0W","interval4_vs_3__0W" ) ],
                              clustering$centers[order(cc$Clustering), c("interval1__0W","interval2_vs_1__0W","interval3_vs_2__0W","interval4_vs_3__0W" ) ] ) )
  if( !all(table(apply( pdists, MARGIN = 1, FUN = function(x) { which.min(x) } )) == table(pp$Clustering) ) )
  {
    stop( "Problem with clustering" )
  }
  pdists <- exp(-pdists)
  pdists <- pdists / rowSums(pdists)
  
  g <- GGally::ggpairs( as.data.frame(pdists)[ sample(1:nrow(pdists), size = min(10000,nrow(pdists)) ), ],
                        lower = list(continuous = my_fn, combo = "box_no_facet" ) )
  
  ggsave( g, filename = paste0( directory, "scatterplot_of_clusterProbabilities_", NUMBER.OF.CLUSTERS, ".png" ), 
          device = "png", width = 210, height = 210, units = "mm", dpi = 150 )
  
  
  # Create statistics in longer form
  pp.freq <- pp %>%
    dplyr::sample_n( min(10000,nrow(pp)) ) %>%
    plyr::count() %>%
    dplyr::mutate( Id = dplyr::row_number() ) %>%
    tidyr::pivot_longer( cols = starts_with("interval"), 
                         names_to = c("interval","type"),
                         names_pattern = "interval(.*)__(.*)",
                         values_to = "Value" )
  
  pp.freq$Interval <- pp.freq$interval
  pp.freq$Interval[pp.freq$interval == "1_ALL"] <- 1
  pp.freq$Interval[pp.freq$interval == "2_ALL"] <- 2
  pp.freq$Interval[pp.freq$interval == "3_ALL"] <- 3
  pp.freq$Interval[pp.freq$interval == "4_ALL"] <- 4
  
  g1 <- ggplot( pp.freq %>% filter( type == "Weighted" ),
                aes( x = as.character(Interval), y = Value ) ) +
    geom_line( aes( group = Id, color = Clustering ), alpha = 0.1 ) +
    geom_line( data = clustering.centers %>% filter( type == "Weighted" ), 
               aes( x = Interval, y = Value, group = Clustering ), 
               color="red" ) +
    facet_grid( . ~ Clustering ) +
    scale_color_viridis_d() +
    xlab( "Interval" ) +
    theme_light() +
    theme( legend.position = "none" ) 
  
  g1b <- ggplot( pp.freq %>% filter( type == "Weighted" ),
                 aes( x = as.character(Interval), y = Value ) ) +
    geom_line( aes( group = Id, color = Clustering ), alpha = 0.1 ) +
    geom_line( data = clustering.centers %>% filter( type == "Weighted" ), 
               aes( x = Interval, y = Value, group = Clustering ), 
               color="red" ) +
    facet_wrap( ~ Clustering ) +
    scale_color_viridis_d() +
    xlab( "Interval" ) +
    theme_light() +
    theme( legend.position = "none" ) 
  
  ggsave( filename = paste0( directory, "clustering_", NUMBER.OF.CLUSTERS, ".png"), device = "png", plot = g1b, width = 297, height = 210, units = "mm", scale = 1, dpi = 150 )
  
  
  g4b <- ggplot( pp  ) +
    geom_density( aes( x = WeightedBurden, fill = Clustering ), alpha = 0.5 ) +
    scale_fill_viridis_d() +
    scale_y_continuous( trans = "sqrt") +
    scale_x_continuous( trans = "sqrt" ) +
    labs( title = "Distribution of weighted burden in clusters" ) +
    theme_light()
  
  ggsave( filename = paste0( directory, "distribution_of_weightedBurden_", NUMBER.OF.CLUSTERS, ".png"), device = "png", plot = g4b, width = 297, height = 160, units = "mm", scale = 1, dpi = 150 )
  
  message( "Computing disease ORs." )
  
  ### Odds of diseases in various clusters ###
  diseases.test <- list()
  for( disease in c( colnames(df.patients.of.interest)[ grep( pattern = "interval", colnames(df.patients.of.interest) ) ] ) )
  {
      for( cluster in 1:NUMBER.OF.CLUSTERS )
      {
        id <- paste0(c(disease,cluster),collapse="_")
        
        # print(disease)
        # print(table( df.patients.of.interest[[disease]]))
        # print(clustering.variable)
        # print(table( df.patients.of.interest[[clustering.variable]] == cluster ))

        s <- strsplit( disease, split = ".", fixed = T )[[1]]
        x <- as.matrix( table( df.patients.of.interest[[disease]], 
                               df.patients.of.interest[[clustering.variable]] == cluster ) )
        if( nrow( x ) == 2 & ncol( x ) == 2 & all( x > 1 ) )
        {
          ft <- fisher.test( x )
          
          diseases.test[[ id ]] <- data.frame( Interval = if( length(s) == 2 ) { as.numeric( gsub( "interval", "", s[1] ) ) } else { 1 },
                                               Clustering = cluster,
                                               Disease = if( length(s) == 2 ) { s[2] } else { disease },
                                               Prevalence = if( length(s) == 2 ) { sum( df.patients.of.interest[[disease]] ) / nrow(df.patients.of.interest) } else { sum( df.patients.of.interest[[disease]] - 1 ) / nrow(df.patients.of.interest) },
                                               OR = ft$estimate,
                                               OR.CIL = ft$conf.int[1],
                                               OR.CIH = ft$conf.int[2],
                                               pvalue = ft$p.value )
        }
      }
  }

      
  diseases.test <- bind_rows( diseases.test )
  diseases.test$Clustering <- factor( diseases.test$Clustering, levels = 1:NUMBER.OF.CLUSTERS )
  diseases.test$AdjP <- p.adjust( diseases.test$pvalue, method = "bonferroni" )
  diseases.test$Significant <- diseases.test$AdjP < 0.05
  if( !is.na(opt$codetable) ) 
	  diseases.test <- diseases.test %>% left_join( codetable, by = c("Disease" = "ICD10" ) )
  else
  	diseases.test$DiseaseName <- diseases.test$Disease
  diseases.test$DiseaseName[ diseases.test$Disease == "Sex" ] <- "Sex"
  if( !is.na(opt$codetable) && length(intersect(depression.vars,codetable$ICD10)) == 0 )
    diseases.test$DiseaseName[ diseases.test$Disease %in% depression.vars ] <- "Depression"
  #diseases.test$DiseaseName[ diseases.test$Disease == "DEPRESSION" ] <- "Depression"
  
  write.xlsx( diseases.test, file = paste0( directory, "diseases_tests_", NUMBER.OF.CLUSTERS, ".xlsx" ), overwrite = T ) 
  
  max.prevalences <- diseases.test %>%
    group_by( Disease ) %>%
    summarize( MaxPrevalence = max(Prevalence) )
  
  diseases.test$DiseaseName <- factor( diseases.test$DiseaseName )
  #diseases.test$DiseaseName <- relevel( diseases.test$DiseaseName, ref = "Depression" )
  
  selected.diseases <- diseases.test %>% 
    filter( Disease %in% depression.vars | 
              ( Significant == T & 
                  ( Disease %in% ( max.prevalences %>% 
                                     filter( MaxPrevalence > 0.05 ) %>% 
                                     pull( Disease ) ) & 
                      abs(log2(OR)) > log2(3) ) ) )
  
  selected.diseases$DiseaseName <- factor( selected.diseases$DiseaseName )
  
  
  g9 <- ggplot( selected.diseases ) +
    geom_point( aes( size = abs(log2(OR)), y = DiseaseName, x = Interval, color = (OR < 1), alpha = abs(log2(OR)) ) ) +
    scale_color_manual( values = c("red","blue") ) +
    scale_y_discrete( limits = rev(levels(selected.diseases$DiseaseName)) ) + #, labels = function(x) stringr::str_wrap(x, width = 40)) +
    scale_x_discrete( limits = factor(1:length(onset.intervals))) +
    facet_grid( . ~ Clustering ) +
    ylab( NULL ) +
    theme_light()
  
  g10 <- ggplot( clustering.centers %>% filter( type == "Weighted" ) ) +
    geom_line( aes( x = Interval, y = Value, group = Clustering, color = Clustering ), 
               size = 3 ) +
    facet_grid( . ~ Clustering ) +
    scale_color_viridis_d() +
    xlab( "Interval" ) +
    ylab( "Direct multimorbidity\nfor depression (dMM4D)") +
    theme_light() +
    theme( legend.position = "none" )
  
  fg12 <- cowplot::plot_grid( g9, g10, align = "v", ncol = 1, axis = "lr", rel_heights = c(3,1) )
  
  ggsave( filename = paste0( directory, "ORs_", NUMBER.OF.CLUSTERS, ".png"), device = "png", plot = fg12, width = 297, height = 210, units = "mm", scale = 1.0, dpi = 300 )
}

## Compute clusters for all other samples
  
# This will be used for clustering!
m <- matrix( data = NA, 
             nrow = nrow( patient.profiles ), 
             ncol = 4, dimnames = list( rownames( patient.profiles ),
                                        c("interval1__0W","interval2_vs_1__0W","interval3_vs_2__0W","interval4_vs_3__0W" ) ) )

m[,"interval1__0W"] <- patient.profiles$interval1_ALL__Weighted
m[,"interval2_vs_1__0W"] <- patient.profiles$interval2_ALL__Weighted - patient.profiles$interval1_ALL__Weighted
m[,"interval3_vs_2__0W"] <- patient.profiles$interval3_ALL__Weighted - patient.profiles$interval2_ALL__Weighted
m[,"interval4_vs_3__0W"] <- patient.profiles$interval4_ALL__Weighted - patient.profiles$interval3_ALL__Weighted

m.scaled <- t( ( t(m) - scaling.centers ) / scaling.scales )

pdists <- as.matrix( pdist( m.scaled,
                            clustering$centers[ order(cc$Clustering), 
                                                c("interval1__0W","interval2_vs_1__0W","interval3_vs_2__0W","interval4_vs_3__0W" ) ] ) )

clusters <- apply( pdists, MARGIN = 1, FUN = function(x) { which.min(x) } )

message( "Cluster sizes:" )
print( table( clusters ) )
print( table( patient.profiles$MaxInterval,
              clusters ) )

# 1. Cluster indices
output.clusters <- data.frame( id = patient.profiles[[ id.column ]],
                               cluster = clusters )
colnames(output.clusters)[1] <- id.column
output.clusters.file <- paste0( output.dir, "/", "cluster_indices_", NUMBER.OF.CLUSTERS, ".txt" )

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
output.cluster.probs.file <- paste0( output.dir, "/", "cluster_probabilities_", NUMBER.OF.CLUSTERS, ".txt" )

write.table( x = output.cluster.probs,
             file = output.cluster.probs.file,
             quote = F, append = F, sep = "\t", row.names = F, col.names = T )

message( paste0( "Cluster probabilities: ", output.cluster.probs.file ) )

# 3. Cluster raw variables
output.clusters.raw <- cbind( data.frame( id = patient.profiles[[ id.column ]] ),
                          as.data.frame( m ) )
colnames(output.clusters.raw)[1] <- id.column
output.clusters.raw.file <- paste0( output.dir, "/", "cluster_variables_", NUMBER.OF.CLUSTERS, ".txt" )

write.table( x = output.clusters.raw,
             file = output.clusters.raw.file,
             quote = F, append = F, sep = "\t", row.names = F, col.names = T )

message( paste0( "Cluster variables (raw values): ", output.clusters.raw.file ) )

message( "Done." )



