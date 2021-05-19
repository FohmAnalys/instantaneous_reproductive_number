rm(list = ls())

#heading <- "
#------------------------------------------------------------------------------------
#
#
# CONTACT:      Lisa Brouwers
# EMAIL:        analysenheten@folkhalsomyndigheten.se
# PROJECT:      COVID-19 modelling
# STUDY:        Rt
#
#
#
# R version:   	4.0.2
# What:   	    Modelling of epidemic trajectory of COVID-19 following analysis of Timothy Churches.
#				
#------------------------------------------------------------------------------------
#\n"



analysis_date <- "2021-04-09"




#-----------------------------------------------------------------------------------
# Directory and paths
#-----------------------------------------------------------------------------------



# Input your own path where folder (end with "/") is located in project.path.
# E.g. project.path 	     <- "C:/Users/Modelling/"


project.path 	<- ""

data.path 		<- paste(project.path, "Data", sep="")
output.path 	<- paste(project.path, "Output",sep="")
script.path 	<- paste(project.path, "Script",sep="")








## ----label="options", echo=FALSE---------------------------------

options(replace.assign=TRUE, width=90)

options(scipen=999)



#---------------------------------------------------------------------------------------------------
# LIBRARIES
#---------------------------------------------------------------------------------------------------
library(xtable)
library(tidyr)
library(ggplot2)
library(EpiEstim)
library(lubridate)
library(incidence)
library(openxlsx)


## ----'functions'--------------------------------------------------------------

# custom results plotting function to avoid the ugly
# TableGrob messages returned by the plotting function in the
# EpiEstim package
plot_Ri <- function(estimate_R_obj) {
                p_I  <- plot(estimate_R_obj, "incid", add_imported_cases = TRUE, col=1:2)  # plots the incidence
                p_SI <- plot(estimate_R_obj, "SI")  # plots the serial interval distribution
                p_Ri <- plot(estimate_R_obj, "R")
                #return(gridExtra::grid.arrange(p_I, p_SI, p_Ri, ncol = 1))
                g <- arrangeGrob(p_I, p_SI, p_Ri, ncol = 1)
                #return(list(I=p_I, SI=p_SI, R=p_Ri))
}

plotit <- function(region, lower, upper){

    daily_tick  <- seq(1, nrow(incidence_daily), 2)
    daily_label <- incidence_daily$dates[daily_tick]
    y.lim       <- c(-1, max(incidence_daily$Total)*1.2)

    plot(gc_fit, ylab="Antal fall", xlab="Datum", axes=FALSE, ylim=y.lim)
    axis(1, at=daily_tick, labels=daily_label, pos=-0.5)
    axis(2, at=seq(0, max(incidence_daily$Total)*1.2, 1))

    lines(tt, lower, type="l", col="black", lty=2)
    lines(tt, upper, type="l", col="black", lty=2)


}




## ----'parameters'-------------------------------------------------------------

today <- analysis_date

## Report delay: 2 days.
yesterday <- ymd(today) - ddays(2)




# start 7 days before so that first estimation is on 2020-08-01 due to 7 day window.
start.date <- ymd("2020-08-01") - ddays(7)
start.curve <- ymd("2020-08-03")



#-----------------------------------------------------------------------------------
# Read data
#-----------------------------------------------------------------------------------
#


# Incidence data with no exclusion
# R_t_data contains two data frames
# 1. line_list which contains dates for each reported case a date and whether it was imported or local
# 2. incidence_daily which contains a summation for each day from line_list, 
# e.g. july 27 there were 3 imported cases, 68 local, and therefore in total 71 cases.



load(file.path(data.path, "R_t_data.Rdata")) 

#-----------------------------------------------------------------------------------
# Edit data/Create variables
#-----------------------------------------------------------------------------------



long_dates <- as.Date(incidence_daily$dates)

# format dates: abbreviated month - day

incidence_daily$dates <- format(as.Date(incidence_daily$dates), "%b %d")

# 
# inc_date   <- as.Date(long_dates)
# times_each <- incidence_daily$Total
# inc        <- NULL
# inc$Date   <- rep(inc_date, times_each)



## ----'plot_Rt_I', out.width='0.7\\linewidth'----------------------------------

    #ref:
    #Hiroshi Nishiura, Natalie M. Linton, Andrei R. Akhmetzhanov,
    #Serial interval of novel coronavirus (COVID-19) infections, Int J Infectious Diseases
    mu_si    <- 4.8  # days
    sigma_si <- 2.3  # days


    dat <- data.frame(dates=long_dates, local=incidence_daily$local, imported=incidence_daily$imported)

    incid <- incidence(line_list$Statistikdatum, group=line_list$utlandssmitta)

    res_parametric_si <- estimate_R(incid,
                            method = "parametric_si",
                            config = make_config(list(mean_si = mu_si,
                                                     std_si = sigma_si)))

    plot(res_parametric_si, "incid", add_imported_cases = TRUE, col=1:2, ylim=c(-1, 10))  # plots the incidence




## ----'plot_Rt_R', out.width='0.9\\linewidth'----------------------------------
    plot(res_parametric_si, "R",
            options_R=list(xlim=c(start.curve, yesterday)))  # plots the incidence


## ----'epi_si_Rt_7_I', results='asis'------------------------------------------

R_t <- res_parametric_si$R[c('t_start', 't_end', 'Mean(R)', 'Quantile.0.025(R)', 'Quantile.0.975(R)')]

end_7day_window <- res_parametric_si$dates[R_t$t_end]

R_t <- data.frame(end_7day_window=end_7day_window, R_t)

R_t$end_7day_window <- ymd(R_t$end_7day_window)

# whole Rt list for saving
R_t_all <- R_t
R_t_all <- R_t_all[R_t_all$end_7day_window >= start.curve,]
R_t_all <- R_t_all[, -c(which(colnames(R_t_all) %in% c("t_start", "t_end")))]

colnames(R_t_all) <- c("Slutdatum_fönster", "Medelvärdet", "2.5% kvantil", "97.5% kvantil")


# Rt last 15 days
start_date_15 <- yesterday - 15
R_t <- R_t[R_t$end_7day_window >= start_date_15, -c(which(colnames(R_t) %in% c("t_start", "t_end")))]

R_t[, -which(colnames(R_t) == "end_7day_window")] <- round(R_t[, -which(colnames(R_t) == "end_7day_window")],2)

R_t[,1] <- as.character(R_t[,1])
colnames(R_t) <- c("Slutdatum_fönster", "Medelvärdet", "2.5% kvantil", "97.5% kvantil")


print(xtable(R_t,
    caption=paste("Riket: momentant reproduktionstal med 95\\% trovärdig intervall baserat på ett glidande 7-dagarsfönster:", start_date_15, " -- ", yesterday, sep=""),
    label="tab:Rt_R"),
    sanitize.text.function = function(str){
        str <- gsub("_", " ", str, fixed=TRUE)
        str <- gsub("/", " ", str, fixed=TRUE)
        str <- gsub("%", "\\%", str, fixed=TRUE)
        return(str)
      },
	 table.placement="h",
     caption.placement="top", include.colnames=TRUE, include.rownames=FALSE)		



## ----'spara_listan'-----------------------------------------------------------

info <- "Sveriges momentant reproduktionstal baserad på 7-dagars fönster. Detaljer i rapporten: https://www.folkhalsomyndigheten.se/smittskydd-beredskap/utbrott/aktuella-utbrott/covid-19/statistik-och-analyser/analys-och-prognoser/. "

write.xlsx(list(R_t=R_t_all,
                Info=as.data.frame(info)),
            file=file.path(output.path, paste0("R_t_(", analysis_date, ").xlsx")), colNames = TRUE, colWidths = c("auto", "auto", "auto"))



