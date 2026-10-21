# Claudia 2022

library(readr)
library(dplyr)
library(tidyr)

# Read in the mortality rates file from Andrea
# Jorge will need to update this using the Chilean life table
mort_rates <- readr::read_csv("df_mortrate_state_age_Mexico_City_2023_2050.csv")

# First, create/or look at the birth year grid (Jorge and Claudia will meet to work on this)
# to figure out the first and last birth year of the LONGEST FIT program

first_birth_cohort <- 1926 
last_birth_cohort <- 1999

first_year_of_program <- 2025

first_cohort_age <- 45
last_age <- 100

# Calculate the birth cohort by subtracting the year by the age
mort_rates$birth_cohort <- mort_rates$year - mort_rates$age

# Filter the data to only look at these birth cohorts you picked above based on the program years
# and also filter ages greater than or equal to the first screening age (earliest cohort age)
mort_rates <- mort_rates %>% filter(birth_cohort >= first_birth_cohort & birth_cohort <= last_birth_cohort & age >= first_cohort_age & age <= last_age & year >= first_year_of_program)

# Now select the columns of data you need to create the ann haz file, age, birth_cohort, and mortality rates
mort_rates <- mort_rates %>% select(age, birth_cohort, death_rate)

# DID NOT ADAPT WELL FROM OLD CODE
# # If there are values missing for certain ages, then fill in with data the same age, but prevous birth cohort
# #Loop over age
# for(bc in seq (first_birth_cohort,last_birth_cohort,1)){
#   bc <- 1997
#   temp <- mort_rates %>% filter(birth_cohort == bc)
#   # Find last age with data
#   last_val <- max(temp$age)
#   # Grab the same data from the previous birth cohort to use for any future ages.
#   value_at_age <- LT[age,last_val]
#   # print(c(age, last_val))
#   # print(value_at_age)
#   for(length_to_extend in seq(last_val,num_col,1)){
#     #Extend last val
#     LT[age,length_to_extend] <- value_at_age
#   }
# }

mort_rates

# If the age=100, set the death rates to 1000 (SimCRC trick)
mort_rates$death_rate <- ifelse(mort_rates$age == 100, 1000, mort_rates$death_rate)


# Reformat it as is needed for this analysis, fill missing values with zeros
reformatted_mort_rates <- mort_rates %>% spread(birth_cohort, death_rate, fill = 0)



# BE CAREFUL TO FIX MISSING VALUES AT THE OLDEST BIRTH COHORTS111111
# SEE FILLED IN VERSION!!!
write_csv(reformatted_mort_rates,"Formatted ann haz file/reformatted_mort_rates.csv")



