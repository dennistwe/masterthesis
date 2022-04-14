##############################################################################

# function to load raw fluorescence titration data 
 
###############################################################################

#function to read read raw fluorescence titration data
#data are stored in an xlsx format prespecified by the manufacturer
#of the software to handle fluorescence intensity experiments
#data need to have the same number of measurment points 

#required packages:
#readxl
#dplyr
#purrr

#function arguments:
#data.path: rel. path to the files e.g. "RawData/FluorescenceSpec"
#file.index: vector containing the position of the files in the folder
#trna.conc: starting concentration of the trna as int
###################################################################

read_fi_raw <- function(data.path, file.index, trna.conc){
# gives the paths of the files to be read (# files specified by n_rep)
  n_rep <- sapply(1:length(file.index), function(x) paste0("replica", x))
  input.files <- paste0(data.path, "/", list.files(path = data.path)[file.index])
  
  # xlsx files are stored in a list and intensity measurements
  # that are stored rowwise are first transposed and subsequently
  # extracted from the list entries and stored in a tibble
  input.list <- map(input.files, read_excel, skip = 11, col_names = FALSE)
  names(input.list) <- n_rep
  input.list1 <- input.list %>% 
    map( function(.x) pivot_longer(.x, c(2:ncol(.x)), values_to = c("fi"))) %>% 
    map(function(.x)  dplyr::select(.x, fi)) %>% 
    map(function(.x) unique(.x)) %>% 
    map(function(.x) filter(.x, fi != "NA")) %>%
    map_df(`[[`, "fi")
  
  input.list2 <- input.list1 %>%  
    mutate(trna = c(trna.conc, sapply(c(1:(nrow(input.list1)-2)), 
                                      function(x) trna.conc/2^x), 0))
  return(input.list2)
}
