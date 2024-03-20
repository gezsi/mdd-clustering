#!/usr/bin/Rscript 

# check and install packages
list.of.packages <- c("ggplot2", "data.table", "optparse", "dplyr", "openxlsx")
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if( length(new.packages) > 0 ) 
  install.packages(new.packages)

# load packages
suppressWarnings(suppressMessages(library(data.table, quietly = T)))
suppressWarnings(suppressMessages(library(ggplot2, quietly = T)))
suppressWarnings(suppressMessages(library(optparse, quietly = T)))
suppressWarnings(suppressMessages(library(dplyr, quietly = T, warn.conflicts = F)))
suppressWarnings(suppressMessages(library(openxlsx, quietly = T)))

plot.number.of.diseases.per.age <- function( df, vars, age.var = "Age", top.age = 8, max.age = 70 )
{
  pb = txtProgressBar(min = 1, max = max(df[,age.var],na.rm = T), initial = 1, width = 100, style = 3 ) 
  
  df.numOfDiseases <- data.frame()
  for( age in 1:max(df[,age.var],na.rm = T) )
  {
    # progress bar
    setTxtProgressBar( pb, age )
    
    condition <- df[,age.var] >= age
    x <- cbind( Age = age, 
                as.data.frame( table( apply( df[condition,vars] > 0 & df[condition,vars] <= age, 1, sum, na.rm=T ) ) ) )
    x$Percent <- x$Freq / sum( condition, na.rm = T )
    x$Var1 <-  as.numeric(as.character(x$Var1))
    x <- rbind( x[x$Var1 <= (top.age-1),],
                data.frame( Age = age,
                            Var1 = top.age, 
                            Freq = sum( x$Freq[x$Var1 >= top.age ] ),
                            Percent = sum( x$Percent[x$Var1 >= top.age ] ) ) )
    df.numOfDiseases <- rbind( df.numOfDiseases, x )
  }
  
  df.numOfDiseases <- df.numOfDiseases %>% filter( Var1 > 0 ) %>% mutate( Var1 = as.character(Var1) )
  df.numOfDiseases$Var1[ df.numOfDiseases$Var1 == as.character(top.age) ] <- paste0( top.age, "+" )
  
  p <- ggplot( df.numOfDiseases ,
               aes( x = Age, y = Percent, fill = as.factor(Var1) ) ) +
    geom_area( stat = "identity" ) + 
    scale_fill_brewer( palette = "Set3", name = "Number of\ndiseases" ) +
    xlim( 0, max.age ) +
    theme_minimal()
  
  close(pb)
  
  return( p )
}

plot.number.of.diseases.per.age.stratified.by.sexAndIncome <- function( df, vars, age.var = "Age", sex.var = "Sex", income.var = "Income", top.age = 8, max.age = 70 )
{
  pb = txtProgressBar(min = 1, max = max(df[,age.var],na.rm = T), initial = 1, width = 100, style = 3 ) 
  
  df.numOfDiseases <- data.frame()
  for( age in 1:max(df[,age.var],na.rm = T) )
    for( sex in unique(df[,sex.var] ) )
      for( income in unique(df[,income.var] ) )
      {
        # progress bar
        setTxtProgressBar( pb, age )
        
        condition <- df[,age.var] >= age & df[,sex.var] == sex & df[,income.var] == income
        
        if( sum(condition) == 0 )
          next
        
        x <- cbind( Sex = sex,
                    Income = income, 
                    Age = age, 
                    as.data.frame( table( apply( df[condition,vars] > 0 & df[condition,vars] <= age, 1, sum, na.rm=T ) ) ) )
        x$Percent <- x$Freq / sum( condition, na.rm = T )
        x$Var1 <-  as.numeric(as.character(x$Var1))
        x <- rbind( x[x$Var1 <= (top.age-1),],
                    data.frame( Sex = sex,
                                Income = income, 
                                Age = age,
                                Var1 = top.age, 
                                Freq = sum( x$Freq[x$Var1 >= top.age ] ),
                                Percent = sum( x$Percent[x$Var1 >= top.age ] ) ) )
        df.numOfDiseases <- rbind( df.numOfDiseases, x )
      }
  
  df.numOfDiseases <- df.numOfDiseases %>% filter( Var1 > 0 ) %>% mutate( Var1 = as.character(Var1) )
  df.numOfDiseases$Var1[ df.numOfDiseases$Var1 == as.character(top.age) ] <- paste0( top.age, "+" )
  
  p <- ggplot( df.numOfDiseases,
               aes( x = Age, y = Percent, fill = as.factor(Var1) ) ) +
    geom_area( stat = "identity" ) + 
    facet_grid( Sex ~ Income, labeller = label_both ) +
    scale_fill_brewer( palette = "Set3", name = "Number of\ndiseases" ) +
    xlim( 0, max.age ) +
    theme_minimal()
  
  close(pb)
  
  return( p )
}

compute.prevalence <- function( dt, diseases, codetable = NA )
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
  
  dt.freq2$Frequency <- 1.0 - (dt.freq2$`NA's` / nrow(dt))
  dt.freq2$Cases <- nrow(dt) - dt.freq2$`NA's`
  
  if( any(!is.na(codetable)) )
	dt.freq2 <- dt.freq2 %>% left_join( codetable, by = c("disease" = "ICD10") )
  
  return( dt.freq2 )
}

option_list = list(
  make_option( c("-i", "--input"), action = "store", default = NA, type = 'character',
               help = "input file name (csv)" ),
  make_option( c("--age"), action="store", default="Age",
               help="Name of variable 'age' [default: %default]" ),
  make_option( c("--sex"), action="store", default="Sex",
               help="Name of variable 'sex' [default: %default]" ),
  make_option( c("--household-income"), action="store", default="Income",
               help="Name of variable 'household income' [default: %default]" ),
  make_option( c("--depression"), action="store", default="F32",
               help="Name of variable 'depression' (if it is composed of multiple variables, then the format should be: comma separated strings without spaces) [default: %default]" ),
  make_option( c("--codetable"), action="store", default=NA, type = 'character',
               help="Filename for code table containing disease names for each disease code (columns of the input file) [default: %default]" ),
  make_option( c("-r","--remove"), action="store", default=NA, type = 'character',
               help="names of variables to remove (format: comma separated strings without spaces)" ),
  make_option(c("--plot"), action="store_true", default = T,
              help="Create plots [default: %default]." ),
  make_option(c("--no-plot"), action="store_false", dest="plot",
              help="Don't create plots." ) )

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

# read input file
message( "Reading input file..." )
df <- fread( opt$input, data.table = F, verbose = F, showProgress = T )

if( !is.na(opt$codetable) ) {
	codetable <- openxlsx::read.xlsx( opt$codetable, sheet = 1, colNames = T )
} else {
	codetable <- NA
}

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

# check dataset
## check sex
genders <- unique( df[,gender.var] )
if( length( genders ) != 2 ) {
  warning( "Sex should be a binary variable." )
  any.errors <- T
}
if( !( 1 %in% genders ) | !( 2 %in% genders ) ) {
  warning( "Sex should be coded with with values '1' and '2'." )
  any.errors <- T
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

## compute disease prevalences
message( "Computing disease prevalences..." )

dt.freq <- compute.prevalence( df, 
                               diseases = vars, 
                               codetable = codetable )

dt.freq.depression <- compute.prevalence( df[ which( rowSums( !is.na(df[,depression.vars,drop=F]) ) > 0 ), ], 
                                          diseases = vars, 
                                          codetable = codetable )

openxlsx::write.xlsx( list( "Prevalences" = dt.freq, 
                            "Prevalences in depression" = dt.freq.depression ),
                      file = "disease-prevalences.xlsx",
                      colNames = T, asTable = T, firstRow = T, colWidths = "auto", overwrite = T )

if( opt$plot )
{
  message( "Creating plots..." )
  # compute disease burden plots
  p <- plot.number.of.diseases.per.age( df, 
                                        vars = vars, 
                                        age.var = age.var, 
                                        top.age = 9, 
                                        max.age = quantile( df[,age.var], probs = 0.95, na.rm = T ) )
  
  suppressWarnings( ggsave( "number_of_diseases_per_age.pdf", plot = p, device = "pdf", width = 210, height = 210, units = "mm" ) )
  
  p2 <- plot.number.of.diseases.per.age.stratified.by.sexAndIncome( df, 
                                                                    vars = vars, 
                                                                    age.var = age.var, 
                                                                    sex.var = gender.var,
                                                                    income.var = household_income.var,
                                                                    top.age = 9, 
                                                                    max.age = quantile( df[,age.var], probs = 0.95, na.rm = T ) )
  
  suppressWarnings( ggsave( "number_of_diseases_per_age_stratified.pdf", plot = p2, device = "pdf", width = 210, height = 210, units = "mm" ) )
  
  
  # plot disease counts
  disease.counts <- rowSums( df[,vars] > 0, na.rm = T )
  
  p3 <- ggplot( cbind( df[,c(gender.var,household_income.var)], Count = disease.counts ) ) +
    geom_violin( aes( x = 1, y = Count ) ) +
    facet_grid( as.formula( paste0( gender.var, " ~ ", household_income.var ) ), 
                labeller = label_both ) +
    theme_minimal() +
    labs( x = NULL, y = "Count of diseases" ) +
    theme( axis.text.x = element_blank() )
  
  suppressWarnings( ggsave( "number_of_diseases_stratified.pdf", plot = p3, device = "pdf", width = 210, height = 210, units = "mm" ) )
  
  p4 <- ggplot( data = df[ !is.na(df[,household_income.var]), ], 
                aes_( x = as.name(age.var), fill = as.name(household_income.var) ) ) + 
    geom_density( alpha = 0.5 ) +
    scale_fill_discrete( name = "Household income" ) +
    facet_wrap( as.formula( paste0( "~ ", gender.var ) ), labeller = label_both ) +
    theme_minimal() +
    theme( legend.position = "bottom" )
  
  suppressWarnings( ggsave( "income_distribution_stratified.pdf", plot = p4, device = "pdf", width = 210, height = 210, units = "mm" ) )
  
  if( file.exists("Rplots.pdf") )
    file.remove("Rplots.pdf")
}

## compute pairwise association between variables
message( "Computing pairwise association between variables..." )
disease.comorbidities.list <- list()

filtered.diseases <- unique( c( dt.freq$disease[ dt.freq$Frequency > 0.005 ], 
                                dt.freq.depression$disease[ dt.freq.depression$Frequency > 0.005 ] ) )

filtered.vars <- c( gender.var, household_income.var, filtered.diseases )

pb = txtProgressBar(min = 1, max = length(filtered.vars)-1, initial = 1, width = 100, style = 3) 

for( i in 1:(length(filtered.vars)-1) )
{
  # progress bar
  setTxtProgressBar(pb,i)
  
  var.i <- filtered.vars[i]
  
  if( !(var.i %in% colnames(df)) )
    next
  
  for( j in (i+1):length(filtered.vars) )
  {
    var.j <- filtered.vars[j]
    
    if( !(var.j %in% colnames(df)) )
      next
    
    cond.i <- if( var.i %in% filtered.diseases ) { !is.na(df[,var.i]) } else { df[,var.i] }
    cond.j <- if( var.j %in% filtered.diseases ) { !is.na(df[,var.j]) } else { df[,var.j] }
    
    binary.i <- length(unique(cond.i)) == 2
    binary.j <- length(unique(cond.j)) == 2
    
    if( binary.i & binary.j )
    {
      test <- suppressWarnings( fisher.test( cond.i, cond.j ) )
      
      
      disease.comorbidities.list[[paste0(i,"-",j)]] <- data.frame( Var.A = as.character(var.i),
                                                                   Var.B = as.character(var.j),
                                                                   p.value = test$p.value,
                                                                   OR = test$estimate ) %>% mutate_if(is.factor, as.character)
  
    } else {
      test <- suppressWarnings( chisq.test( cond.i, cond.j ) )
      
      m <- as.matrix( table( cond.i, cond.j ) )
      if( nrow(m) > ncol(m) )
        m <- t(m)
		
	    if( nrow(m) < 2 | ncol(m) < 2 )
		    next
      
      max.or <- 1.0
      min.or <- 1.0
      for( l in 2:nrow(m) )
        for( k in 2:ncol(m) )
        {
          or <- (m[l,k]/m[1,k])/(m[l,1]/m[1,1])
          
          if( is.na(or) | is.nan(or) )
            next
          
          if( or < min.or )
            min.or <- or
          if( or > max.or )
            max.or <- or
        }
      
      or = if( (1.0 / min.or) > max.or ) { min.or } else { max.or }
      
      disease.comorbidities.list[[paste0(i,"-",j)]] <- data.frame( Var.A = as.character(var.i),
                                                                   Var.B = as.character(var.j),
                                                                   p.value = test$p.value,
                                                                   OR = or ) %>% mutate_if(is.factor, as.character)
    }
  }  
}

disease.comorbidities <- bind_rows( disease.comorbidities.list )
disease.comorbidities$adj.p.value <- p.adjust( disease.comorbidities$p.value, method = "bonferroni" )

openxlsx::write.xlsx( list( "Pairwise associations" = disease.comorbidities ),
                      file = "disease-comorbidities.xlsx",
                      colNames = T, asTable = T, firstRow = T, colWidths = "auto", overwrite = T )

close(pb)
