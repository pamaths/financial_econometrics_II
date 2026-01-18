suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(haven)
  library(lubridate)
  library(modelsummary)
  library(fixest)
  library(rstudioapi)
  library(tibble)
  library(ggplot2)
  library(tidyverse)
})

start_year <- 1994
end_year <- 2013
chartertype_commercial_bank <- 200  # only commercial banks
reference_year <- 1976
trim_bool <- TRUE
q1 <- 0.01 # triming quantile

# Importing data
callreports_path <- "callreports_1976_2020_WRDS.dta"
herf_path    <- "l1_herfdepcty.csv"
fed_path <- "fredgraph.csv"

call_raw <- read_dta(callreports_path)
herf_raw <- read_csv(herf_path)
fed_raw  <- read_csv(fed_path)


# Cleaning main database
df_call <- call_raw %>%
  mutate(
    reportdate = ymd(as.character(date)),
    year = year(reportdate),
    quarter = quarter(reportdate),
    reference_quarter = (year - reference_year) * 4 + quarter
  ) %>%
  filter(chartertype == chartertype_commercial_bank,
         year >= start_year,
         year <= end_year)

# Cleaning hhi
df_hhi <- herf_raw %>%
  mutate(
    year = as.numeric(substr(dateq, 1, 4)),
    quarter = as.numeric(substr(dateq, 6, 6)),
    reference_quarter = (year - reference_year) * 4 + quarter
  )

# Cleaning target
df_ff <- fed_raw %>%
  mutate(
    date = as.Date(observation_date),
    year = year(date),
    quarter = quarter(date),
    target =
      case_when(
        date < as.Date("2008-12-01") ~ DFEDTAR / 100,
        # before crisis
        date == as.Date("2008-12-01") ~ 0.25 / 100,
        # fix empty data
        date > as.Date("2008-12-01") ~ DFEDTARU / 100,
        # after crisis
      )
  ) %>%
  group_by(year, quarter) %>%
  summarise(target = mean(target), .groups = "drop") %>% #We take the mean FF on the quarter
  mutate(reference_quarter = (year - 1976) * 4 + quarter) %>%
  select(reference_quarter, target, year, quarter)

# Merge
df <- df_call %>%
  left_join(df_hhi, by = c("cert", "reference_quarter"))

df <- df %>%
  group_by(cert, reference_quarter) %>%
  slice_max(order_by = assets,
            n = 1,
            with_ties = FALSE) %>% #Keep only the duplicate with highest asset
  ungroup()

# Remove duplicates
#df <- df %>%
#  group_by(cert, reference_quarter) %>%
#  filter(n() == 1) %>%
#  ungroup()

# Add target
df <- df %>%
  left_join(df_ff, by = "reference_quarter")

# Construct regressors
log <- function(x)
  log(pmax(x, 1e-6))

# Filter data is not contingent in time, quarters that cannot be lagged by 1
df <- df %>%
  arrange(cert, reference_quarter) %>%
  group_by(cert) %>%
  mutate(rq_diff = lead(reference_quarter) - reference_quarter)
df <- df %>%
  filter(rq_diff == 1)

data <- df %>%
  select(
    cert,
    reference_quarter,
    year,
    deposits,
    totsavdep,
    timedep,
    fedfundsrepoliab,
    tradingliabilities,
    otherborrowedmoney,
    liabilities,
    assets,
    cash,
    securities,
    loans,
    reloans,
    ciloans,
    intexpalldep,
    intexpfordep,
    foreigndep,
    l1_herfdepcty,
    target
  ) %>%
  arrange(cert, reference_quarter) %>%                # sort by bank & time
  group_by(cert) %>%
  mutate(
    # ----------------------------------------------------------
    # DEPENDENT VARIABLES: FORWARD log-changes (t → t+1)
    # ----------------------------------------------------------
    d_total_deposits = log(lead(deposits))      - log(deposits),
    d_savings        = log(lead(totsavdep))        - log(totsavdep),
    d_time           = log(lead(timedep))       - log(timedep),
    wholesale_funding = otherborrowedmoney + replace_na(fedfundsrepoliab, 0) +
      replace_na(tradingliabilities, 0),
    d_wholesale      = log(lead(wholesale_funding)) - log(wholesale_funding),
    d_total_liab     = log(lead(liabilities))   - log(liabilities),
    
    # Panel B: Assets side
    d_total_assets = log(lead(assets))      - log(assets),
    d_cash         = log(lead(cash))        - log(cash),
    d_securities   = log(lead(securities))  - log(securities),
    d_total_loans  = log(lead(loans))       - log(loans),
    d_real_estate  = log(lead(reloans)) - log(reloans),
    d_CI_loans     = log(lead(ciloans))      - log(ciloans),
    
    # ----------------------------------------------------------
    # RHS REGRESSORS: Lag structure EXACTLY per DSS (t and t-1)
    # ----------------------------------------------------------
    d_ff    = lead(target) - target,
    # Δtarget_t
    hhi_lag = l1_herfdepcty,
    # HHI_{t−1}
    interaction = d_ff * l1_herfdepcty,
    # Already lagged
    
    # -------------------------
    # Deposit rate and deposit spread (DSS definition)
    # -------------------------
    deposit_rate     = 4 * ((intexpalldep - replace_na(intexpfordep, 0)) / pmax(deposits-replace_na(foreigndep, 0), 1e-6)),
    # annualized, avoid /0
    d_deposit_rate   = lead(deposit_rate) - deposit_rate,
    d_deposit_spread = d_ff - d_deposit_rate,
    
    # ZLB control
    post2008 = if_else(year >= 2009, 1, 0)
  ) %>%
  ungroup()

dv_list <- c(
  "d_total_deposits",
  "d_savings",
  "d_time",
  "d_wholesale",
  "d_total_liab",
  "d_total_assets",
  "d_cash",
  "d_securities",
  "d_total_loans",
  "d_real_estate",
  "d_CI_loans",
  "d_deposit_spread"   
)
# Trim data
trim <- function(x) {
  q <- quantile(x, probs = c(q1, 1 - q1), na.rm = TRUE)
  x[x < q[1] | x > q[2]] <- NA
  x
}
if (trim_bool){
data <- data %>%
  mutate(across(all_of(dv_list), trim))
}


m1 <- feols(
  d_total_deposits ~ hhi_lag + interaction |
    cert + cert^post2008 + reference_quarter,
  cluster = ~ cert,
  data = data
)
m2 <- feols(
  d_deposit_spread ~ hhi_lag + interaction |
    cert + cert^post2008 + reference_quarter,
  cluster = ~ cert,
  data = data
)
m3 <- feols(
  d_savings ~ hhi_lag + interaction |
    cert + cert^post2008 + reference_quarter,
  cluster = ~ cert,
  data = data
)
m4 <- feols(
  d_time ~ hhi_lag + interaction |
    cert + cert^post2008 + reference_quarter,
  cluster = ~ cert,
  data = data
)
m5 <- feols(
  d_wholesale ~ hhi_lag + interaction |
    cert + cert^post2008 + reference_quarter,
  cluster = ~ cert,
  data = data
)
m6 <- feols(
  d_total_liab ~ hhi_lag + interaction |
    cert + cert^post2008 + reference_quarter,
  cluster = ~ cert,
  data = data
)

m7  <- feols(
  d_total_assets ~ hhi_lag + interaction |
    cert + cert^post2008 + reference_quarter,
  cluster = ~ cert,
  data = data
)
m8  <- feols(
  d_cash ~ hhi_lag + interaction |
    cert + cert^post2008 + reference_quarter,
  cluster = ~ cert,
  data = data
)
m9  <- feols(
  d_securities ~ hhi_lag + interaction |
    cert + cert^post2008 + reference_quarter,
  cluster = ~ cert,
  data = data
)
m10 <- feols(
  d_total_loans ~ hhi_lag + interaction |
    cert + cert^post2008 + reference_quarter,
  cluster = ~ cert,
  data = data
)
m11 <- feols(
  d_real_estate ~ hhi_lag + interaction |
    cert + cert^post2008 + reference_quarter,
  cluster = ~ cert,
  data = data
)
m12 <- feols(
  d_CI_loans ~ hhi_lag + interaction |
    cert + cert^post2008 + reference_quarter,
  cluster = ~ cert,
  data = data
)

modelsummary(list(m1, m2, m3, m4, m5, m6, m7, m8, m9 , m10 , m11, m12), stars = TRUE)

etable(
  m1,
  m2,
  m3,
  m4,
  m5,
  m6,
  keep = "%interaction",
  # only show ΔFF × Bank HHI
  dict = c("interaction" = "$\\Delta FF_t \\times \\text{Bank HHI}_{t-1}$"),
  se.below = TRUE,
  digits = 3,
  signif.code = c("***" = 0.01, "**" = 0.05, "*" = 0.1),
  fitstat = c("n", "r2"),
  # Observations, R²
  tex = FALSE,
  #file = "tableVIII_panelA.tex",
  title = "Panel A: Liabilities",
  style.tex = "aer"                      # nice economics style
)

etable(
  m7,
  m8,
  m9,
  m10,
  m11,
  m12,
  keep = "%interaction",
  dict = c("interaction" = "$\\Delta FF_t \\times \\text{Bank HHI}_{t-1}$"),
  se.below = TRUE,
  digits = 3,
  signif.code = c("***" = 0.01, "**" = 0.05, "*" = 0.1),
  fitstat = c("n", "r2"),
  tex = FALSE,
  #file = "tableVIII_panelB.tex",
  title = "Panel B: Assets",
  style.tex = "aer"
)