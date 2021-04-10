.libPaths();
.libPaths("C:/Users/Guanghui He/Documents/R/win-library/4.0");

getwd();
setwd("C:/Users/Guanghui He/Documents/R/win-library/4.0");

# install.packages("gbm");
# install.packages("ggplot2", lib="C:/Users/Guanghui He/Documents/R/win-library/4.0");

library(DiagrammeR);


# plot a flow chart of multiple imputation analysis;

DiagrammeR::grViz("digraph {
  
graph[layout = dot, rankdir = LR]

a
b
c

a -> b -> c

}"
)


# data <- list(a=1000, b=800, c=600, d=400)


DiagrammeR::grViz("
digraph graph2 {

graph [layout = dot]

# node definitions with substituted label text
node [shape = rectangle, width = 4, fillcolor = Biege]
a [label = '@@1']
b [label = '@@2']
c [label = '@@3']
d [label = '@@4']

a -> b -> c -> d

}

[1]: paste0('X comes from', italic)
[2]: paste0('Remove Errors')
[3]: paste0('Identify Potential Customers')
[4]: paste0('Select Top Priorities')
")

render_graph(graph);

graph1 <-
  create_graph() %>%
    add_path(
      n = 5,
      edge_aes = edge_aes(
        arrowhead = c(
          "normal", "vee",
          "tee", "dot"
        ),
        color = c(
        "red", "blue",
        "orange", "purple"
        )
      )
    )

# graph %>% render_graph();

# graph %>% 
#  generate_dot() %>%
#  cat()

export_graph(graph1,
file_name = "diagrammertestyulei",
title = 'Simple Graph', height=100, width=100
)


# example code;

# With `create_graph()` we can
# simply create an empty graph (and
# add in nodes and edges later
# with other functions)
graph <- create_graph()
# A graph can be created with
# nodes and without having any edges;
# this can be done in 2 steps:
# 1. create a node data frame (ndf)
# using `create_node_df()`
ndf <-
create_node_df(n = 4)
# 2. create a new graph object with
# `create_graph()` and then pass
# in the ndf to `nodes_df`
graph <-
create_graph(
nodes_df = ndf)
# Get information on the graph's nodes
graph %>%
get_node_info()
# You can create a similar graph with
# just nodes but also providing a
# range of attributes for the nodes
# (e.g., types, labels, or arbitrary
# 'values')
ndf <-
create_node_df(
n = 4,
label = TRUE,
type = c("type_1", "type_1",
"type_5", "type_2"),
shape = c("circle", "circle",
"rectangle", "rectangle"),
values = c(3.5, 2.6, 9.4, 2.7))
graph <-
create_graph(nodes_df = ndf)
# Get information on the graph's
# internal node data frame (ndf)
graph %>%
get_node_df()
# A graph can also be created by
# specifying both the nodes and
# edges; create an edge data frame

# (edf) using the `create_edge_df()`
# function:
edf <-
create_edge_df(
from = c(1, 2, 3),
to = c(4, 3, 1),
rel = "leading_to",
values = c(7.3, 2.6, 8.3))
# Create the graph object with
# `create_graph()` and pass in the
# ndf and edf objects
graph <-
create_graph(
nodes_df = ndf,
edges_df = edf)
# Get information on the graph 
graph %>% get_edge_df()
# Get information on the graph's
# internal node data frame (ndf)
graph %>% get_node_df()


node_df <-
  create_node_df(
    n = 4,
    type = "a",
    label = c(2384, 3942, 8362, 2194),
    style = "filled",
    color = "aqua",
    shape = c("circle", "circle",
              "rectangle", "rectangle"),
    value = c(3.5, 2.6, 9.4, 2.7))


########### successful example for the book cover of imputation;


DiagrammeR::grViz("digraph {

graph [layout = dot, rankdir = LR]

# define the global styles of the nodes. We can override these in box if we wish
node [shape = rectangle, style = filled, fillcolor = Linen]

oridata[label= ' 1.1 3.1 ... \n 2.2  ? ...  \n ?   3.8 ... \n ... ', 
shape= folder, width=1, height=1, fillcolor=Red]

data1 [label = ' 1.1 3.1 ... \n 2.2 1.3 ... \n 2.3 3.8 ... \n ... ', 
width=1, height=1, shape = folder, fillcolor = Orange]
data2 [label = ' 1.1 3.1 ... \n 2.2 1.1 ... \n 1.9 3.8 ... \n ... ', 
width=1, height=1, shape = folder, fillcolor = Orange]
data3 [label = ' 1.1 3.1 ... \n 2.2 1.7 ... \n 3.4 3.8 ... \n ... ', 
width=1, height=1, shape = folder, fillcolor = Orange]
data4 [label = ' 1.1 3.1 ... \n 2.2 0.9 ... \n 1.6 3.8 ... \n ... ', 
width=1, heigh=1, shape = folder, fillcolor = Orange]
data5 [label = ' 1.1 3.1 ... \n 2.2 1.4 ... \n 2.7 3.8 ... \n ... ', 
width=1, height=1, shape = folder, fillcolor = Orange]

Results1 [shape=oval, label = 'Results \n set 1  ', shape = oval, 
width=1, height=1, fillcolor = Yellow]
Results2 [label = 'Results \n set 2  ', shape = oval, 
width=1, height=1, fillcolor = Yellow]
Results3 [label = 'Results \n set 3  ', shape = oval, 
width=1, height=1, fillcolor = Yellow]
Results4 [label = 'Results \n set 4  ', shape = oval, 
width=1, height=1, fillcolor = Yellow]
Results5 [label = 'Results \n set 5  ', shape = oval, 
width=1, height=1, fillcolor = Yellow]

results [label= ' Multiple imputation \n analysis results  ', shape=oval, fillcolor=Green]

# edge definitions with the node IDs

edge[label='Imputation']

oridata -> data1;

edge[label='']

oridata -> data2;
oridata -> data3;
oridata -> data4;
oridata -> data5;

edge[label='Completed-data analysis']

data1 -> Results1;

edge[label='']

data2 -> Results2;
data3 -> Results3;
data4 -> Results4;
data5 -> Results5;

edge[label='Combining']

Results1 -> results;

edge[label='']

Results2 -> results;
Results3 -> results;
Results4 -> results;
Results5 -> results;
# Results5 -> combining [style=dashed];



}")


###################  to create the figure in Chapter 3 ############################

DiagrammeR::grViz("digraph {

graph [layout = dot, rankdir = LR]

# define the global styles of the nodes. We can override these in box if we wish
node [shape = rectangle, style = filled, fillcolor = Linen]

oridata[label='Incomplete \n dataset', shape= folder, fillcolor=Beige]

impute1 [label = 'Imputation 1', shape = oval, fillcolor = Beige]
impute2 [label = 'Imputation 2', shape = oval, fillcolor = Beige]
impute3 [label = 'Imputation 3', shape = oval, fillcolor = Beige]
impute4 [label = 'Imputation 4', shape = oval, fillcolor = Beige]
impute5 [label = 'Imputation 5', shape = oval, fillcolor = Beige]

data1 [label = 'Completed \n dataset 1', shape = folder, fillcolor = Beige]
data2 [label = 'Completed \n dataset 2', shape = folder, fillcolor = Beige]
data3 [label = 'Completed \n dataset 3', shape = folder, fillcolor = Beige]
data4 [label = 'Completed \n dataset 4', shape = folder, fillcolor = Beige]
data5 [label = 'Completed \n dataset 5', shape = folder, fillcolor = Beige]

Results1 [label = 'Estimates 1', shape = folder, fillcolor = Beige]
Results2 [label = 'Estimates 2', shape = folder, fillcolor = Beige]
Results3 [label = 'Estimates 3', shape = folder, fillcolor = Beige]
Results4 [label = 'Estimates 4', shape = folder, fillcolor = Beige]
Results5 [label = 'Estimates 5', shape = folder, fillcolor = Beige]

combining [shape=oval, label =  'Combining \n estimates', fillcolor=Beige]

analysis1 [shape=oval, label='Analysis 1', fillcolor=Beige]
analysis2 [shape=oval, label='Analysis 2', fillcolor=Beige]
analysis3 [shape=oval, label='Analysis 3', fillcolor=Beige]
analysis4 [shape=oval, label='Analysis 4', fillcolor=Beige]
analysis5 [shape=oval, label='Analysis 5', fillcolor=Beige]

results [label= 'A single set \n of final inferences', shape=oval, fillcolor=Beige]

# edge definitions with the node IDs

oridata -> impute1;
oridata -> impute2;
oridata -> impute3;
oridata -> impute4;
oridata -> impute5;

impute1 -> data1;
impute2 -> data2;
impute3 -> data3;
impute4 -> data4;
impute5 -> data5;

data1 -> analysis1;
data2 -> analysis2;
data3 -> analysis3;
data4 -> analysis4;
data5 -> analysis5;

analysis1 -> Results1;
analysis2 -> Results2;
analysis3 -> Results3;
analysis4 -> Results4;
analysis5 -> Results5;

Results1 -> combining;
Results2 -> combining;
Results3 -> combining;
Results4 -> combining;
Results5 -> combining;
# Results5 -> combining [style=dashed];

combining -> results;



}")