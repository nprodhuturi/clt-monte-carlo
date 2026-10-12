# setup.R: variables shared by each analysis

student_id <- 8120

set.seed(student_id)  #last 4 digits of my student ID - for reproducible random results

B <- 10000                 # n = sample size, B= number of repetitions (how many sample means we simulate)
n_values <- c(1:40,45,50,75,100)

