if(!require(pacman)) install.packages("pacman")
pacman::p_load("dplyr",         # Data management
               "qgraph",        # Network analysis
               "psych",         # Psychometrics
               "corpcor",       # Partial correlations
               "bootnet",       # Bootstrap on networks
               "graphicalVAR",  # N=1 time series VAR network
               "mlVAR",         # Multilevel time series VAR network
               "lm.beta")       # Standardized regression coefficients

# ==============================================================================
# DATA IMPORT
# ==============================================================================

# Load data from Borsboom et al. (2021) Nature Methods Primer
load(url("https://github.com/DennyBorsboom/NatureMethodsPrimer_NetworkAnalysis/raw/main/Psychopathology/data/clean_network.RData"))

# Variable labels
varLabs <- c("Relax","Irritable","Worry","Nervous","Future","Anhedonia",
             "Tired","Hungry","Alone","Angry","Social_offline","Social_online",
             "Music","Procrastinate","Outdoors","C19_occupied","C19_worry","Home")

Data5b <- Data2

# ==============================================================================
# DATA PREPARATION
# ==============================================================================

# Rename variables with meaningful labels
vars <- paste0("Q", 1:18)
names(Data5b)[names(Data5b) %in% vars] <- varLabs

# Select variables for network analysis (removing non-psychological items)
excluded_vars <- c("Hungry","Angry","Music","Procrastinate","Tired",
                   "Outdoors","Home","C19_occupied","C19_worry")
Data5b <- Data5b %>% select(-all_of(excluded_vars))
varLabs <- varLabs[!varLabs %in% excluded_vars]

# ==============================================================================
# DETRENDING SETUP
# ==============================================================================

# Significance threshold for trend removal
alpha <- 0.05

# Create data frames for storing fitted values
fitted_all <- expand.grid(
  beep = seq(min(Data5b$beep), max(Data5b$beep)),
  day = seq(min(Data5b$day), max(Data5b$day))
)

fitted_day <- data.frame(
  day = seq(min(Data5b$day), max(Data5b$day))
)

fitted_beep <- data.frame(
  beep = seq(min(Data5b$beep), max(Data5b$beep))
)

# Initialize storage
p_values <- data.frame(var = c("day", "beep"))
testStatistics <- list()
coefficients <- list()
stdcoefficients <- list()

# ==============================================================================
# TIME VARIABLE CREATION
# ==============================================================================

# Convert beep to factor with time labels
beep_labels <- c("09:00 - 12:00","12:00 - 15:00","15:00 - 18:00","18:00 - 21:00")
Data5b$beepFactor <- factor(Data5b$beep, levels = 0:3, labels = beep_labels)
fitted_all$beepFactor <- factor(fitted_all$beep, levels = 0:3, labels = beep_labels)
fitted_beep$beepFactor <- factor(fitted_beep$beep, levels = 0:3, labels = beep_labels)

# Create date variables (starting March 15, 2020)
Data5b$date <- as.Date("2020-03-15") + Data5b$day
fitted_all$date <- as.Date("2020-03-15") + fitted_all$day
fitted_day$date <- as.Date("2020-03-15") + fitted_day$day

# Add midpoint times for each beep
midpoint_times <- c("10:30","13:30","16:30","19:30")
Data5b$midTime <- as.POSIXct(paste(Data5b$date, 
                                   midpoint_times[Data5b$beep + 1]), 
                             format = "%Y-%m-%d %H:%M", 
                             tz = "Europe/Amsterdam")

fitted_all$midTime <- as.POSIXct(paste(fitted_all$date, 
                                       midpoint_times[fitted_all$beep + 1]), 
                                 format = "%Y-%m-%d %H:%M", 
                                 tz = "Europe/Amsterdam")

# ==============================================================================
# DETRENDING PROCEDURE
# ==============================================================================

data_detrended <- Data5b

# Detrend each variable
for (v in seq_along(varLabs)){
  # Model: variable ~ intercept + day trend + beep effects
  formula <- as.formula(paste0(varLabs[v], " ~ 1 + day + factor(beep)"))
  lmRes <- lm(formula, data = Data5b)
  
  # Extract coefficients
  fixed <- coef(lmRes)
  
  # Test significance and zero out non-significant trends
  p_values[[varLabs[v]]] <- anova(lmRes)[["Pr(>F)"]][1:2]
  
  if (p_values[[varLabs[v]]][1] > alpha){
    fixed[2] <- 0  # Remove day trend if not significant
  }
  if (p_values[[varLabs[v]]][2] > alpha){
    fixed[3:5] <- 0  # Remove beep effects if not significant
  }
  
  # Calculate fitted values
  fitted_all[,varLabs[v]] <- fixed[1] + 
    fixed[2] * fitted_all[["day"]] +  
    fixed[3] * (fitted_all[["beep"]] == 1) + 
    fixed[4] * (fitted_all[["beep"]] == 2) + 
    fixed[5] * (fitted_all[["beep"]] == 3)
  
  fitted_day[,varLabs[v]] <- fixed[1] + fixed[2] * fitted_day[["day"]]
  
  fitted_beep[,varLabs[v]] <- fixed[1] + 
    fixed[2] * median(fitted_day[["day"]]) +  
    fixed[3] * (fitted_beep[["beep"]] == 1) + 
    fixed[4] * (fitted_beep[["beep"]] == 2) + 
    fixed[5] * (fitted_beep[["beep"]] == 3)
  
  # Detrend by subtracting fitted values
  data_detrended[,varLabs[v]] <- Data5b[,varLabs[v]] - 
    (fixed[1] + 
       fixed[2] * Data5b[["day"]] +  
       fixed[3] * (Data5b[["beep"]] == 1) + 
       fixed[4] * (Data5b[["beep"]] == 2) + 
       fixed[5] * (Data5b[["beep"]] == 3))
  
  # Store test statistics
  testStatistics[[v]] <- cbind(
    data.frame(var = varLabs[v], effect = rownames(anova(lmRes))), 
    anova(lmRes)
  )
  
  # Store coefficients
  coefficients[[v]] <- data.frame(
    var = varLabs[v],
    type = names(coef(lmRes)),
    coef = coef(lmRes),
    std = coef(lm.beta(lmRes))
  )
}

# ==============================================================================
# DATA CLEANING
# ==============================================================================

# Remove rows with excessive missing values (>50% missing)
nas <- data_detrended %>%
  select(all_of(varLabs)) %>%
  is.na() %>%
  rowSums()

data_detrended <- data_detrended %>%
  filter(nas < 9) %>%
  distinct(id, day, beep, .keep_all = TRUE)