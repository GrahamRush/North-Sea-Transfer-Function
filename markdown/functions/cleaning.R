###    split the data
  df.split <- function (x, n)
  {
    spec <<-subset(x[, 1:(n-5)])
    env <<- subset(x[, (n-4):n])
  }


###    function to clean the data for low counts and rare soecies
clean_data <- function (data, region, min.count = 75, sites = NULL, file = TRUE, rare = TRUE) {
  df <- (data [data$Site %in% sites, ])
  ts.basic <- subset(df, select = -c(Agglut.unidentified,  Unidentified.spp.))
  ts.50 <- subset(ts.basic, Total.dead >= min.count) # remove samples with counts < c 
  ts.species <- subset(ts.50, select = -c(Site, Total.dead, Elevation, SWLI, HoF)) # split the data
  ts.enviro <- subset(ts.50, select = c(Site, Total.dead, Elevation, SWLI, HoF)) # split the data
  ts.species <- ts.species [, colSums(ts.species != 0) > 0] # remove species with no occurences
  ts.clean <- cbind(ts.species, ts.enviro) # join back with the environmental data
  if(rare == TRUE) { # remove rare species 
    m <- apply(ts.species, 2, max) # calculate maximum abundance
    m.2 <- (100/nrow(ts.species)) * (apply(ts.species, 2, function(x) sum( x > 0, na.rm=TRUE))) # calculate % of species occurences 
    ts.rare.species <- ts.species[,!(m.2 < 2 & m < 5)] # remove species that occur in < 2 % of samples and with a maximum abundance < 5 %
    ts.clean <- cbind(ts.rare.species, ts.enviro) # join back with the environmental data
  }
  if(file == TRUE){ # produce a csv of output
    file.name <- paste(region, "_clean", ".csv", sep = "")
    write.csv(ts.clean, file = file.name)
  }
  return(ts.clean)
}
  

# function to clean the core data
clean_core <- function(core, n = 75, td = 1, s = 15, rare = FALSE, f = FALSE, dom = TRUE, n2 = 50) {
  # 1. remove samples with less than n
  c.50 <- subset (core, Total.dead >= n | Total.dead < td)
  # 1a retain samples > n2 and with a dominant species (> 50 %)
  if (dom == TRUE){
    c.50 <- subset (core, Total.dead >= n2 | Total.dead < td)
    c.50 <- subset (c.50, Total.dead >= n | (Total.dead < n & apply(c.50[, 1:s], MARGIN = 1, function(x) any(x > 50))))
  } 
  # 2. subset data
  c.species <- c.50 [c(1:s)] # s is the last column of species data
  c.enviro <- c.50 [c(-1:-s)]
  # 3. remove all species with 0 counts
  c.species <- c.species [, colSums(c.species != 0) > 0]
  # 4. rejoin data
  c.clean <- cbind(c.species, c.enviro)
  # 5. remove rare species
  if(rare == TRUE) { # remove rare species 
      m <- apply(c.species, 2, max) # calculate maximum abundance
      m.2 <- (100/nrow(c.species)) * (apply(c.species, 2, function(x) sum(x>0, na.rm=TRUE))) # calculate % of species occurences 
      c.rare.species <- c.species[,!(m.2 < 2 & m < 5)] # remove species that occur in < 2 % of samples and with a maximum abundance < 5 %
      c.clean <- cbind(c.rare.species, c.enviro) # join back with the environmental data
    }
  if(f == TRUE){ # produce a csv of output
      write.csv(c.clean, "core-clean.csv")
  }
  return(c.clean)
}

# split the core
c.split <- function (x = core, s = 15) {
  spec <-subset(x[, 1:s])
  env <- subset(x[, (s + 1):ncol(x)])
  split <- list("spec" = spec, "env" = env)
  return(split)
}
