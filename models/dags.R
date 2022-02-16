# Title:        Specify DAGs for regression models
# Author:       Ewan Carr
# Started:      2021-09-21

library(tidyverse)
library(here)
library(dagitty)
library(ggdag)

# Functions to process/plot DAGs ----------------------------------------------

dag_to_data <- function(dag) {
  dag %>%
    tidy_dagitty() %>%
    node_status() %>%
    dag_adjustment_sets() %>%
    as_tibble() %>%
    mutate(fill = factor(coalesce(status, adjusted),
                         levels = c("exposure", "outcome",
                                    "unadjusted", "adjusted"),
                         ordered = TRUE,
                         labels = c("Exposure", "Outcome",
                                    "Unadjusted", "Adjusted")))
}

draw_dag <- function(d, ff = "Arial Narrow") {
  d$name <- str_wrap(d$name, 10)
  ggplot(d,
         aes(x = x, 
             y = -y, 
             xend = xend, 
             yend = -yend, 
             fill = fill)) +
     theme_dag(base_family = ff) +
     geom_point(colour = "black",
                alpha = 0.8,
                stroke = 1,
                shape = 21,
                size = 20) +
     geom_dag_edges_arc(curvature = d$curve) +
     geom_dag_text(size = 3,
                   family = ff,
                   color = "black") +
     expand_plot(expand_y = expansion(c(0.1, 0.1)),
                 expand_x = expansion(c(0.1, 0.1))) +
     scale_fill_manual(values = c("#ccdc49",
                                  "#47b8e4",
                                  "#fa9a9a",
                                  "#fafafa")) +
     guides(fill = guide_legend(override.aes = list(size = 10))) +
     theme(legend.title = element_blank())
}


###############################################################################
####                                                                      #####
####        H1: Clinical variables --> Usability and acceptability       #####
####                                                                      #####
###############################################################################


d1 <- dagitty('dag {
              "Benefit receipt" [pos="0.2,0.5"]
              "Comorbid condition" [pos="0.2,0.4"]
              "Currently working" [pos="0.35,0.6"]
              "Female gender" [pos="0.3,0.4"]
              "Previous wearable use" [pos="0.57,0.4"]
              "Switch to Android" [pos="0.63,0.4"]
              "Years of education" [adjusted,pos="0.2,0.6"]
              "Life event" [pos="0.5,0.6"]
              Age [adjusted,pos="0.4,0.3"]
              Partner [pos="0.45,0.4"]
              WSAS [exposure,pos="0.4,0.5"]
              Usability [outcome,pos="0.6,0.5"]
              "Life event" -> "WSAS"
              "Benefit receipt" -> WSAS
              "Comorbid condition" -> "Benefit receipt"
              "Comorbid condition" -> WSAS
              "Currently working" -> "Benefit receipt"
              "Currently working" -> "Comorbid condition" [pos="0.297,0.401"]
              "Currently working" -> WSAS
              "Female gender" -> WSAS
              "Previous wearable use" -> Usability
              "Switch to Android" -> Usability
              "Years of education" -> "Benefit receipt"
              "Years of education" -> "Comorbid condition" [pos="0.171,0.374"]
              "Years of education" -> "Currently working"
              "Years of education" -> WSAS
              "Years of education" -> Usability
              Age -> "Comorbid condition"
              Age -> "Currently working"
              Age -> Partner
              Age -> WSAS
              Age -> Usability
              Partner -> Usability
              WSAS -> Usability
}')

d1_data <- dag_to_data(d1)

# Adjust curves
d1_data$curve <- 0
yoe <- "Years of education"
cc <- "Comorbid condition"
d1_data[d1_data$name == yoe & d1_data$to == cc, ]$curve <- 0.4
d1_data[d1_data$name == "Age" & d1_data$to == "Currently working", ]$curve <- -0.1

p_d1 <- draw_dag(d1_data)

###############################################################################
####                                                                      #####
####                                 H2: Usability --> Usage              #####
####                                                                      #####
###############################################################################

d2 <- dagitty('dag {
bb="0,0,1,1"
"Comorbid condition" [pos="0.2,0.4"]
"Currently working" [pos="0.3,0.6"]
"Female gender" [pos="0.4,0.35"]
"Life event" [pos="0.5,0.6"]
"Partner" [pos="0.45,0.3"]
"Previous wearable use" [adjusted,pos="0.55,0.3"]
"Benefit receipt" [pos="0.3,0.35"]
"Switched to Android" [pos="0.5,0.3"]
"Usability" [exposure,pos="0.5,0.4"]
"Years of education" [adjusted,pos="0.2,0.5"]
Age [adjusted,pos="0.3,0.2"]
Usage [outcome,pos="0.6,0.5"]
"Clinical variables" [adjusted,pos="0.4,0.5"]
"Comorbid condition" -> "Benefit receipt"
"Comorbid condition" -> Usage
"Comorbid condition" -> "Clinical variables"
"Currently working" -> "Comorbid condition" [pos="0.297,0.401"]
"Currently working" -> "Benefit receipt"
"Currently working" -> Usage [pos="0.483,0.573"]
"Currently working" -> "Clinical variables"
"Female gender" -> "Clinical variables"
"Life event" -> Usage
"Life event" -> "Clinical variables"
"Partner" -> "Usability"
"Previous wearable use" -> "Usability"
"Previous wearable use" -> Usage
"Benefit receipt" -> "Clinical variables"
"Switched to Android" -> "Usability"
"Usability" -> Usage
"Years of education" -> "Comorbid condition"
"Years of education" -> "Currently working"
"Years of education" -> "Benefit receipt"
"Years of education" -> "Usability"
"Years of education" -> "Clinical variables"
Age -> "Comorbid condition"
Age -> "Currently working"
Age -> "Partner"
Age -> "Usability"
Age -> "Clinical variables"
"Clinical variables" -> "Usability"
"Clinical variables" -> Usage
}')
d2_data <- dag_to_data(d2)

# Adjust curves
d2_data$curve <- 0
d2_data[d2_data$name == "Age" & 
        d2_data$to == "Currently working",]$curve <- 0.15

p_d2 <- draw_dag(d2_data)

###############################################################################
####                                                                      #####
####                    H3: Clinical variables --> Usage                  #####
####                                                                      #####
###############################################################################

d3 <- dagitty('dag {
bb="0,0,1,1"
"Benefit receipt" [pos="0.3,0.5"]
"Comorbid condition" [adjusted,pos="0.2,0.6"]
"Currently working" [adjusted,pos="0.35,0.6"]
"Female gender" [pos="0.4,0.4"]
"Life event" [adjusted,pos="0.5,0.6"]
"Previous wearable use" [adjusted,pos="0.55,0.3"]
"Switched to Android" [pos="0.5,0.3"]
"Years of education" [adjusted,pos="0.2,0.4"]
Age [adjusted,pos="0.3,0.3"]
Partner [pos="0.45,0.3"]
Usability [pos="0.5,0.4"]
Usage [outcome,pos="0.6,0.5"]
WSAS [exposure,pos="0.4,0.5"]
"Benefit receipt" -> WSAS
"Comorbid condition" -> "Benefit receipt"
"Comorbid condition" -> Usage
"Comorbid condition" -> WSAS
"Currently working" -> "Benefit receipt"
"Currently working" -> "Comorbid condition"
"Currently working" -> Usage [pos="0.483,0.573"]
"Currently working" -> WSAS
"Female gender" -> WSAS
"Life event" -> Usage
"Life event" -> WSAS
"Previous wearable use" -> Usability
"Previous wearable use" -> Usage
"Switched to Android" -> Usability
"Years of education" -> "Benefit receipt"
"Years of education" -> "Comorbid condition" [pos="0.171,0.374"]
"Years of education" -> "Currently working"
"Years of education" -> Usability
"Years of education" -> WSAS
Age -> "Comorbid condition"
Age -> "Currently working"
Age -> Partner
Age -> Usability
Age -> WSAS
Partner -> Usability
Usability -> Usage
WSAS -> Usability
WSAS -> Usage
}')

d3_data <- dag_to_data(d3)

# Adjust curves
d3_data$curve <- 0
yoe <- "Years of education"
cw <- "Currently working"
d3_data[d3_data$name == yoe & d3_data$to == "Usability", ]$curve <- 0.15
d3_data[d3_data$name == yoe & d3_data$to == cw, ]$curve <- -0.2
d3_data[d3_data$name == cw & d3_data$to == "Usage", ]$curve <- -0.1

p_d3 <- draw_dag(d3_data)

# Save ------------------------------------------------------------------------

# Save DAGs
save(d1, d2, d3, file = here("outputs", "dags.Rdata"))

# Save images
p <- list(p_d1, p_d2, p_d3)
names(p) <- paste0("d", 1:3)

walk2(p, names(p),
      ~ ggsave(.x,
               filename = here("writing", "figures", paste0(.y, ".png")),
               dev = "png",
               width = 9,
               height = 6,
               dpi = 300))

