# 1. Set Working Directory
getwd()
setwd("/home/samuel/Desktop/AI_Omics_Internship_2025/Module_I/")

# 2. Create required directories
dir.create("raw_data")
dir.create("clean_data")
dir.create("scripts")
dir.create("results")
dir.create("plots")

# 3. Load dataset
data <- read.csv(file.choose())

# Inspect dataset
summary(data)
str(data)

# Identify variables to fix: gender, smoker, diagnosis

# Convert gender to factor
data$gender_fac <- as.factor(data$gender)

# Create binary factor for gender (Female = 1, Male = 0)
data$gender_num <- as.factor(ifelse(data$gender_fac == "Female", 1, 0))

# Convert smoker to factor
data$smoker_fac <- as.factor(data$smoker)

# Create binary factor for smoker (Yes = 1, No = 0)
data$smoker_num <- as.factor(ifelse(data$smoker_fac == "Yes", 1, 0))

# Convert diagnosis to factor with specified levels
data$diagnosis_fac <- factor(data$diagnosis,
                             levels = c("Normal", "Cancer"))

# View dataset
View(data)

# Save cleaned CSV
write.csv(data, file = "clean_data/patient_info_clean.csv", row.names = FALSE)

# Save workspace
save.image(file = "scripts/clean_Ib.RData")
