read.csv("~/Downloads/preclinical_rates_US_Chile (1).csv")

# *****************************************************************************
#### Compare symptom-detection rates for distal colon cancer (D) ####
# *****************************************************************************

df_D <- NULL

# Parameters that fit the US-F stage distribution
pSxDetS1_D <- 0.06280371
hr_SxDetS2S1_D <- 5.17765749
hr_SxDetS3S2_D <- 2.20033026
hr_SxDetS4S3_D <- 3.86768297

# Convert parameters to rates
rSxDetS1_D <- -log(1 - pSxDetS1_D)
rSxDetS2_D <- hr_SxDetS2S1_D * rSxDetS1_D
rSxDetS3_D <- hr_SxDetS3S2_D * rSxDetS2_D
rSxDetS4_D <- hr_SxDetS4S3_D * rSxDetS3_D

# Store US-F rates
df_D <- rbind(df_D, data.frame(
  id = "US-F",
  rSxDetS1_D = rSxDetS1_D,
  rSxDetS2_D = rSxDetS2_D,
  rSxDetS3_D = rSxDetS3_D,
  rSxDetS4_D = rSxDetS4_D
))

# Parameters that fit the Chile stage distribution
pSxDetS1_D <- 0.03076337
hr_SxDetS2S1_D <- 4.90499694
hr_SxDetS3S2_D <- 5.22166926
hr_SxDetS4S3_D <- 10.78068293

# Convert parameters to rates
rSxDetS1_D <- -log(1 - pSxDetS1_D)
rSxDetS2_D <- hr_SxDetS2S1_D * rSxDetS1_D
rSxDetS3_D <- hr_SxDetS3S2_D * rSxDetS2_D
rSxDetS4_D <- hr_SxDetS4S3_D * rSxDetS3_D

# Store Chile rates
df_D <- rbind(df_D, data.frame(
  id = "Chile",
  rSxDetS1_D = rSxDetS1_D,
  rSxDetS2_D = rSxDetS2_D,
  rSxDetS3_D = rSxDetS3_D,
  rSxDetS4_D = rSxDetS4_D
))

df_D


# *****************************************************************************
#### Repeat for proximal colon cancer (P) ####
# *****************************************************************************

df_P <- NULL

# Parameters that fit the US-F stage distribution
pSxDetS1_P <- 0.02773729
hr_SxDetS2S1_P <- 8.47015693
hr_SxDetS3S2_P <- 2.55182920
hr_SxDetS4S3_P <- 3.78244819

# Convert parameters to rates
rSxDetS1_P <- -log(1 - pSxDetS1_P)
rSxDetS2_P <- hr_SxDetS2S1_P * rSxDetS1_P
rSxDetS3_P <- hr_SxDetS3S2_P * rSxDetS2_P
rSxDetS4_P <- hr_SxDetS4S3_P * rSxDetS3_P

# Store US-F rates
df_P <- rbind(df_P, data.frame(
  id = "US-F",
  rSxDetS1_P = rSxDetS1_P,
  rSxDetS2_P = rSxDetS2_P,
  rSxDetS3_P = rSxDetS3_P,
  rSxDetS4_P = rSxDetS4_P
))

# Parameters that fit the Chile stage distribution
pSxDetS1_P <- 0.02048440
hr_SxDetS2S1_P <- 7.64661458
hr_SxDetS3S2_P <- 4.68534224
hr_SxDetS4S3_P <- 9.76350316

# Convert parameters to rates
rSxDetS1_P <- -log(1 - pSxDetS1_P)
rSxDetS2_P <- hr_SxDetS2S1_P * rSxDetS1_P
rSxDetS3_P <- hr_SxDetS3S2_P * rSxDetS2_P
rSxDetS4_P <- hr_SxDetS4S3_P * rSxDetS3_P

# Store Chile rates
df_P <- rbind(df_P, data.frame(
  id = "Chile",
  rSxDetS1_P = rSxDetS1_P,
  rSxDetS2_P = rSxDetS2_P,
  rSxDetS3_P = rSxDetS3_P,
  rSxDetS4_P = rSxDetS4_P
))

df_P


# *****************************************************************************
#### Repeat for rectal cancer (R) ####
# *****************************************************************************

df_R <- NULL

# Parameters that fit the US-F stage distribution
pSxDetS1_R <- 0.13094501
hr_SxDetS2S1_R <- 1.67106714
hr_SxDetS3S2_R <- 5.26794323
hr_SxDetS4S3_R <- 1.44296312

# Convert parameters to rates
rSxDetS1_R <- -log(1 - pSxDetS1_R)
rSxDetS2_R <- hr_SxDetS2S1_R * rSxDetS1_R
rSxDetS3_R <- hr_SxDetS3S2_R * rSxDetS2_R
rSxDetS4_R <- hr_SxDetS4S3_R * rSxDetS3_R

# Store US-F rates
df_R <- rbind(df_R, data.frame(
  id = "US-F",
  rSxDetS1_R = rSxDetS1_R,
  rSxDetS2_R = rSxDetS2_R,
  rSxDetS3_R = rSxDetS3_R,
  rSxDetS4_R = rSxDetS4_R
))

# Parameters that fit the Chile stage distribution
pSxDetS1_R <- 0.02113080
hr_SxDetS2S1_R <- 8.12764921
hr_SxDetS3S2_R <- 8.68509103
hr_SxDetS4S3_R <- 9.06540889

# Convert parameters to rates
rSxDetS1_R <- -log(1 - pSxDetS1_R)
rSxDetS2_R <- hr_SxDetS2S1_R * rSxDetS1_R
rSxDetS3_R <- hr_SxDetS3S2_R * rSxDetS2_R
rSxDetS4_R <- hr_SxDetS4S3_R * rSxDetS3_R

# Store Chile rates
df_R <- rbind(df_R, data.frame(
  id = "Chile",
  rSxDetS1_R = rSxDetS1_R,
  rSxDetS2_R = rSxDetS2_R,
  rSxDetS3_R = rSxDetS3_R,
  rSxDetS4_R = rSxDetS4_R
))

df_R