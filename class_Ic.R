# Set working directory
getwd()
setwd("/home/samuel/Desktop/AI_Omics_Internship_2025/Module_I/")

# ----------------------------------------------------------------------------------------------------------------
# Exercise 1 - Checking Cholesterol Level (using if)
cholesterol <- 230
if (cholesterol > 240) {
  print("High Cholesterol")
}

# ----------------------------------------------------------------------------------------------------------------
# Exercise 2 - Blood Pressure Status (if...else)
Systolic_bp <- 130
if (Systolic_bp < 120) {
  print("Blood Pressure is normal")
} else {
  print("Blood Pressure is high")
}

# ----------------------------------------------------------------------------------------------------------------
# Exercise 3 - Automating Data Type Conversion with for loop

## --- Working with patient_info.csv ---
# Load dataset
patient_data <- read.csv(file.choose())

# Create a copy
copy_patient_data <- patient_data

# Identify categorical columns (excluding PatientID)
factor_cols <- c("gender", "diagnosis", "smoker") 

# Convert to factors
for (col in factor_cols) {
  copy_patient_data[[col]] <- as.factor(copy_patient_data[[col]])
}

# Check structure
str(copy_patient_data)
View(copy_patient_data)

## --- Working with metadata.csv ---
# Load dataset
metadata <- read.csv(file.choose())

# Create a copy
copy_metadata <- metadata

# Identify categorical columns
factor_cols <- c("gender", "height") 

# Convert to factors
for (col in factor_cols) {
  copy_metadata[[col]] <- as.factor(copy_metadata[[col]])
}

# Check structure
str(copy_metadata)
View(copy_metadata)

# ----------------------------------------------------------------------------------------------------------------
# Exercise 4 - Converting Factors to Numeric Codes

## For patient_info.csv
# Choose binary factor columns (Yes/No)
binary_cols <- c("smoker")

for (col in binary_cols) {
  copy_patient_data[[col]] <- ifelse(copy_patient_data[[col]] == "Yes", 1, 0)
}

# Check structure
str(copy_patient_data)
View(copy_patient_data)


#Checking to compare both original and modified csv files for patient_info
str(patient_data)
str(copy_patient_data)


#Checking to compare both original and modified csv files for Metadata
str(metadata)
str(copy_metadata)




