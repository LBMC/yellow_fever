install.packages("plotly")
install.packages("rjson")
library(plotly)
library(rjson)
library(tidyverse)

# load the data
data <- read_csv("data/2018_09_14_TCR_Component_Test_Laurent_Sept2018.csv") %>%
  filter(!is.na(pDEA_cell_type)) %>%
  filter(!is.na(TRBV)) %>%
  filter(!is.na(TRAV)) %>%
  mutate(chain_id = as.factor(paste0(TRBV, ",", TRAV)))

# get donor id
donors <- data %>% pull(DONOR_ID) %>% as.factor() %>% levels()

# build node_info that contain the name of the node their size and their color
node_TRBV <- data %>%
  filter(DONOR_ID %in% donors[1]) %>%
  group_by(TRBV) %>%
  summarise(mean = mean(pDEA_cell_type),
            n = n()) %>%
  rename(node = TRBV)
node_TRAV <- data %>%
  filter(DONOR_ID %in% donors[1]) %>%
  group_by(TRAV) %>%
  summarise(mean = mean(pDEA_cell_type),
            n = n()) %>%
  rename(node = TRAV)
node_info <- rbind(node_TRBV, node_TRAV)
my_palette <- colorRamp(c("blue", "gray", "red"))
node_info <- node_info %>%
  mutate(color = rgb(my_palette(mean)/255))
node_info %>% dim

node_id <- node_info %>% pull(node) %>% as.factor()
get_node_id <- function(node_id, node){
  x <- apply(as.matrix(node), 1, FUN = function(x, node_id){
    which(node_id %in% x)
  }, node_id = node_id)
  unlist(x)
}

# build x the link infos between node
x <- data %>%
  filter(DONOR_ID %in% donors[1]) %>%
  group_by(chain_id) %>%
  summarise(mean = mean(pDEA_cell_type),
            n = n()) %>%
  mutate(chain_id = as.vector(chain_id),
         TRBV = gsub("(.*),.*", "\\1", as.vector(chain_id), perl = T),
         TRAV = gsub(".*,(.*)", "\\1", as.vector(chain_id), perl = T))
x$TRBV <- get_node_id(node_id, x$TRBV)
x$TRAV <- get_node_id(node_id, x$TRAV)
x


p <- plot_ly(
    type = "sankey",
    domain = list(
      x =  c(0,1),
      y =  c(0,1)
    ),
    orientation = "h",
    valueformat = ".0f",
    valuesuffix = "pMEM",

    node = list(
      label = node_info$node,
      color = node_info$color,
      pad = 15,
      thickness = 15,
      line = list(
        color = "black",
        width = 0.5
      )
    ),

    link = list(
      source = x$TRBV,
      target = x$TRAV,
      value =  x$n,
      label = x$mean
    )
  ) %>%
  layout(
    title = "My nice sankey plot",
    font = list(
      size = 10
    ),
    xaxis = list(showgrid = F, zeroline = F),
    yaxis = list(showgrid = F, zeroline = F)
)
print(p)
chart_link = api_create(p, filename="sankey-basic")
chart_link
