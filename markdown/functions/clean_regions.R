#' read data and set up
#'
#' @title set up training set
#'
#' @description takes the regional data and sets it up for analysis
#'
#' @param n is the minimum number of individuals. 50 or 75 normally
#' @param rare_species is logical to remove rare species. i.e. that occur in < 2 % of samples and with a maximum abundance < 5 % 
#'
#' @return list to add to the training_set data
#'
#' @keywordsset up
#' @export
#' @examples
#' clean_regions(n = 50, rare_species = FALSE, save_file = FALSE)
#' 

clean_regions <- function(n = 50, rare_species = FALSE, save_file = FALSE) {
  data <- list()
  data$northsea <- clean_data(data = mod_data, region = "northsea", min.count = n, sites = training_sets$names$northsea, file = save_file, rare = rare_species)    # run grahams function in 'cleaning' to clean the data and remove samples with counts less than 'min.count' and the option to remove rare species
  data$east <- clean_data(data = mod_data, region = "east", min.count = n, sites = training_sets$names$east, file = save_file, rare = rare_species)
  data$west <- clean_data(data = mod_data, region = "west", min.count = n, sites = training_sets$names$west, file = save_file, rare = rare_species)
  data$northwest <- clean_data(data = mod_data, region = "northwest", min.count = n, sites = training_sets$names$northwest, file = save_file, rare = rare_species)
  data$southwest <- clean_data(data = mod_data, region = "southwest", min.count = n, sites = training_sets$names$southwest, file = save_file, rare = rare_species)
  data$ythan <- clean_data(data = mod_data, region = "ythan", min.count = n, sites = training_sets$names$ythan, file = save_file, rare = rare_species)
  return(data)
}
