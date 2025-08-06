getwd()
setwd("/home/samuel/Desktop/AI_Omics_Internship_2025/Module_I/")



#Creating directories (raw_data, clean_data, scripts, results or Tasks, plots)

dir.create("raw_data")
dir.create("clean_data")
dir.create("scripts")
dir.create("results")
dir.create("plots")


#Loading my csv file 
data <- read.csv(file.choose())

#Examining my dataset
summary(data)


#Variables with inconsistent data types
#The variables with inconsistent data types are gender, smoker and diagnosis


#convert gender to factor
data$gender_fac <- as.factor(data$gender)
str(data)

#convert factor to numeric using ifelse statement (Female = 1 and Male = 0)
data$gender_num <- ifelse(data$gender_fac == "Female", 1, 0)
class(data$gender_num)

#convert numeric gender code to factor
data$gender_num <- as.factor(data$gender_num)
class(data$gender_num)



#convert smker to factor 
data$smoker_fac <- as.factor(data$smoker)
str(data)

#convert factor to numeric using ifelse statement (Yes = 1 and Male = 0)
data$smoker_num <- ifelse(data$smoker_fac == "Yes", 1, 0)
class(data$smoker_num)

#convert numeric gender code to factor
data$smoker_num <- as.factor(data$smoker_num)
class(data$smoker_num)



#convert diagnosis to factor 
data$diagnosis_fac <- as.factor(data$diagnosis)
str(data)

data$diagnosis_fac <- factor(data$diagnosis_fac,
                             levels = c("Normal", "Cancer"))


View(data)

#Saving my cleaned csv file
write.csv(data, file = "clean_data/patient_info_clean.csv")

#Saving my Rdata data file 
save.image(file = "scripts/clean_Ib.RData")
