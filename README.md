# Ring : Pipeline for the analysis of multiplex immunofluorescence stainings.

Benidovskaya Elena, Beyaert Simon and Huyghe Nicolas (22.02.2023)

# Packages to load for this script

```{r packages}
library(tidyverse)
library(ggplot2)
library(tibble)
library(plyr)
library(spatial)
library(spatstat)
library(contoureR)
library(sp)
library(concaveman)
library(maptools)
library(sf)
library(ggpubr)
library(fpc)
library(dbscan)
library(tmap)
library(ComplexHeatmap) ##install from Bioconductor
library(reshape)
library(SummarizedExperiment) ##install from Bioconductor
library(patchwork)
library(BBmisc)
library(readxl)
library(dplyr)
```

# Multiplex data transformation from HALO to R

## Setting the inputs and outputs

```{r variables setup}
output_data <- "~/project/tables/" ## set directory to output data
output_graphs <- "~/project/graphs/" ## set directory to output graphs

wd_data <- "~/project/data/"
wd_tables <- "~/project/tables/"

x <- list.files("~/project/data", pattern = ".csv") ## creates a vector with all file names that ends with .csv in your folder

y <- "df" ## creates a variable where you store the name of the final data frame summarizing the densities of cells of interest for each sample you have
```

## Creating a function to transform your data automatically

When exporting your data from an image analysis software, patient by patient, you end up with a .csv type file that you are going to analyze in R. This data contains the location of each cell and the markers (Alexa FLuors for example) for which your cell is positive. 

When giving you the location,the program gives you a X/Y min and X/Y max for the borders of the cell. To simplify the analysis, we are going to calculate the mean of each point.

Then we are going to select only the columns of the data frame that are interesting for us (for which Alexa FLuor is the cell positive?) and we are going to rename those columns for a easier usage. Here for example, we used Alexa Fluor 647 for the CD8 T cell receptor, ...

After what we can start to create loops according to the cell types we want to study. In fact, different cells are defined according to the markers we use, for example CD8+ T cells should be positive for CD8 and CD3 marker; epithelial tumor cells should be positive to hPanCK marker; but if a cell is positive to hPanCK and CD3/CD8 then it's probably an immune cell lost in the tumor tissue (so it should be counted as a immune cell). While defining the different cells, don't forget to define markers that should be positive (== 1) and negative (== 0). If you have problems with overlapping markers, you can try to manage them in your loops. Then you create a column flag where you put the result of your loops and you save the new tables in the folder you created earlier. Finally you apply those loops to all the list defined in the variable x.

```{r transformation function}
setwd(wd_data) ##set working directory (from the R project)

main_table <- function(a){
  Main_table <- read.csv(a)
  Main_table <- mutate(Main_table, X = (Main_table$XMin + Main_table$XMax)/2)
  Main_table <- mutate(Main_table, Y = (Main_table$YMin + Main_table$YMax)/2)
  Main_table <- select(Main_table, -c(XMin, XMax, YMin, YMax))
  Main_table <- select(Main_table, c("Alexa.Fluor.647.Positive.Classification","Alexa.Fluor.555.Positive.Classification","Alexa.Fluor.488.Positive.Classification", "Classifier.Label" , "X", "Y" ))
  colnames(Main_table) <- c("CD8", "CD3", "hPanCK", "Classifier Label","X","Y")
  
  output <- character(dim(Main_table)[1])
  condition <- Main_table$CD8 == 1 & Main_table$CD3 == 1
  
  for(i in (1:dim(Main_table)[1])[condition]) {
    if(condition[i]) {
      output[i] <- "CD8+ CD3+"
    } else {
      output[i] <- "0"
    }}
  
  condition <- Main_table$hPanCK == 1
  
  for(i in (1:dim(Main_table)[1])[condition]) {
    if(condition[i]) {
      output[i] <- "Tumor cells"
    }}
  
  condition <- Main_table$hPanCK == 1 & Main_table$CD8 == 1 & Main_table$CD3 == 1
  
  for(i in (1:dim(Main_table)[1])[condition]) {
    if(condition[i]) {
      output[i] <- "CD8+ CD3+"
    }}
  
  condition <- Main_table$hPanCK == 0 & Main_table$CD8 == 0 & Main_table$CD3 == 0
  
  for(i in (1:dim(Main_table)[1])[condition]) {
    if(condition[i]) {
      output[i] <- "Stromal cells"
    }}
  
  condition <- Main_table$CD8 == 0 & Main_table$CD3 == 1
  
  for(i in (1:dim(Main_table)[1])[condition]) {
    if(condition[i]) {
      output[i] <- "CD3+"
    }}
  
  condition <- Main_table$CD8 == 1 & Main_table$CD3 == 0
  
  for(i in (1:dim(Main_table)[1])[condition]) {
    if(condition[i]) {
      output[i] <- "CD8+ CD3-"
    }}
  
  Main_table[,"flag"] <- output
  
  output_table <- paste(output_data, str_replace(a, ".csv", ""),".csv", sep="")
  write_csv(Main_table, file=output_table, col_names = TRUE)
  
}

lapply(FUN = main_table, x)
```

## Transforming your data in mm 

Chech the resolution of your scanner an dtransform your data in mm. Here as an example, the resolution of the scanner is 0.325 so you need to multiply your X/Y by the resolution and transform in mm (x 0.001).

```{r transformation mm}
setwd(wd_tables)

main <- function(a){
  
  Main_table <- read.csv(a)
  Main_table <- Main_table %>%
    mutate(Xadj = X * 0.325 * 0.001) %>%
    mutate(Yadj = Y * 0.325 * 0.001)
  
  output_table <- paste(output_data, str_replace(a, ".csv", ""),".csv", sep="")
  write_csv(Main_table, file=output_table, col_names = TRUE)
}

lapply(FUN = main, x)
```

## Visualisation of your biopsies

You can then create a function that is going to take the coordinates and create a visualisation of your biopsies (the axis represent mm and you can assign a color to each cell type you defined earlier).

```{r graph mm}
setwd(wd_tables)

histo_plot <- function(a){
  
  Main_table <- read.csv(a)
  
  plotly <- ggplot(Main_table, aes(x = Xadj, y = Yadj, col = flag)) +
    geom_point(size = 0.1) + ## size of the points on the graph
    scale_color_manual(values = c('#FFFF33','#CC0033','#FF9900','#3333FF','#99FF33'),limits=c("CD3+","CD8+ CD3-", "CD8+ CD3+", "Stromal cells","Tumor cells"),drop=TRUE) + ## define the color you want to use for each cell type
    ggtitle(a) + ## add a title (a for the name of each biopsy)
    coord_fixed(ratio = 1) + ## allows you to correct the extent of the graph
    theme_bw() + 
    theme(panel.grid.major=element_blank(), panel.grid.minor=element_blank()) + ## removes the grid (= white background)
    theme(plot.title = element_text(hjust = 0.5)) + ## sets the title to the center of the graph
    xlab("X (mm)") +
    ylab("Y (mm)") + 
    theme(text = element_text(size = 15)) ## allows you to set the size of the text (legends and titles)
  
  complete_plot_path_name <- paste(output_graphs, str_replace(a, ".csv", ""), ".png", sep = "")
  ggsave(plotly, file = complete_plot_path_name, dpi = 320, units = "mm") ## save each graph
  
}

lapply(FUN = histo_plot, x)
```

## Clusterisation of the biopsy

Now, let's create a function which determines the number of tumor clusters you have on each slide. To do that, we transform each biopsy in a polygon and then calculate the surface of the later. To visualize what R is doing we can create a graph per step.

```{r graph cluster}
setwd(wd_tables)

area_plot <- function(a){
  
  Main_table <- read.csv(a)
  
  df <- Main_table[,c("Xadj","Yadj")]
  db <- dbscan(df, eps=0.5, minPts=50) ## function that clusters our cells
  plot(db, df, main = a, frame = FALSE)
  
  Main_table1 <- st_as_sf(df, coords=c("Xadj","Yadj"))
  Main_table1$cluster <- db$cluster
  Main_table1$flag <- Main_table$flag
  Main_table1$x <- Main_table$Xadj
  Main_table1$y <- Main_table$Yadj
  try(Main_table2 <- Main_table1 %>% filter(cluster!=0))
  polygons <- concaveman(Main_table1, concavity = 1, length_threshold = 0) ## creates polygons for each cluster determined by dbscan
  plot(st_geometry(Main_table1))
  plot(polygons)
  
  polygons2 <- map(unique(Main_table2$cluster), ~ concaveman(Main_table2[Main_table2$cluster %in% .,])) %>%
    map2(unique(Main_table2$cluster), ~ mutate(.x, cluster = .y)) %>%
    reduce(rbind)
  
  area <- st_area(polygons2) ## calculate the area of each cluster
  area1 <- as.data.frame(area)
  
  output_table <- paste(output_data, str_replace(a, ".csv", "_area"),".csv", sep="")
  write_csv(area1, file=output_table, col_names = TRUE)
  
  try(nicemap <- ## create a graph with the outline of each sample
        ggplot() + ## set up the framework
        geom_sf(data = polygons2, color="gray") + ## add the outline using geom_sf
        geom_point(data= Main_table2, aes(x=x, y=y, col = flag), size = 0.1) +
        scale_color_manual(values = c('#FFFF33','#CC0033','#FF9900','#3333FF','#99FF33'),limits=c("CD3+","CD8+ CD3-", "CD8+ CD3+", "Stromal cells","Tumor cells"),drop=TRUE) +
        theme_bw() +
        theme(panel.grid.major=element_blank(), panel.grid.minor=element_blank()) + 
        ggtitle(a) +
        theme(plot.title = element_text(hjust = 0.5))) + 
        theme(text = element_text(size = 15)) 
  
  try(complete_plot_path_name <- paste(output_graphs,  str_replace(a, ".csv", "_polygon"), ".png", sep = ""))
  try(ggsave(nicemap, file = complete_plot_path_name))
}

lapply(FUN = area_plot, x)
```

## Create a summary table with all your samples

Now, let's create a final table with a row per sample and each column represents the density of different cell types you are studying.

```{r final table}

setwd(wd_tables)

removePat <- ".csv"
x2 <- gsub(removePat, "", x) ## removes the .csv in the names of your files

df1 <- data.frame(x2) ## new data frame created with a sample per row
df1[,c(2:6)] <- 0 ## fill columns from 2 to 6 with zeros (create a column per cell type)
colnames(df1) <- c("biopsy", "area", "CD4","CD8", "true_CD8", "CK") ## rename those new columns
output_df1 <- paste(output_data, y, ".csv", sep="")
write_csv(df1, output_df1)

table_df <- function(a){ ## create a function that extracts every important information from your tables
  
  df1 <- read_csv(output_df1)
  
  df3 <- paste(a,"_area.csv", sep="") ## takes the area of each sample
  t <- read.csv(df3)
  area <- sum(t$area)
  
  df1[df1$biopsy==a, 2] <- area
  
  df4 <- paste(a,".csv", sep="") ## takes the data frames with the cells
  
  Main_table <- read.csv(df4)
  
  ## For each cell for which you need to calculate its density
  ## Filter the data frame according to the flag you want and create df
  ## Add this df as a column in your final df
  
  dfcd3 <-  Main_table %>% filter(flag=="CD3+") 
  CD3 <- sum(dfcd3$CD3) ##the $CD3 is the name of the column in the Main_table, not the future name of the column
  
  
  dfcd8 <-  Main_table %>% filter(flag=="CD8+ CD3-") 
  CD8 <- sum(dfcd8$CD8)
  
  
  dftruecd8 <- Main_table %>% filter(flag=="CD8+ CD3+") 
  true_CD8 <- sum(dftruecd8$CD8) ## Here we have two markers (both positive), so the length of $CD8 or $CD3 is the same
  ## We take CD8 here but it doesnt change anything
  
  
  dftumor <- Main_table %>% filter(flag=="Tumor cells")
  CK <- sum(dftumor$hPanCK)
  
  df1[df1$biopsy==a, 3] <- CD3
  df1[df1$biopsy==a, 4] <- CD8
  df1[df1$biopsy==a, 5] <- true_CD8
  df1[df1$biopsy==a, 6] <- CK
  
  write_csv(df1, output_df1)
  
}

lapply(FUN = table_df, x2)
```

In this new table, you have the information about the area of each sample and a count of your cells of interest. Now you can transform it into densities for each row (= sample).

Please note that here we defined true CD8+ T cells as cells positive for both CD3 and CD8 markers but as we do multiplex stainings, sometimes, if a cell binds an antibody, there is no place to bind a second one (which means that cells stained for CD8 could potentially not be stained for CD3). But as we still want to take every cell into accont we create a column CD8 and true CD8.

```{r final table 2}

setwd(wd_tables)

y2 <- paste(y,".csv", sep="") 
df1 <- read.csv(y2)

df1 <- df1 %>%
  mutate(CD3 = CD4 + CD8 + true_CD8)

df2 <- df1 %>%
  mutate(dCD3 = CD3/area) %>%
  mutate(dCD4 = CD4/area) %>%
  mutate(dCD8 = CD8/area) %>%
  mutate(dtrueCD8 = true_CD8/area) %>%
  mutate(dallCD8 = (true_CD8 + CD8)/area) %>%
  mutate(dCK = CK/area)


output_df2 <- paste(output_data, y, ".csv", sep="")
write_csv(df2, output_df2)
```

Please note that in this script there are several chunks which means that we need to set the working directory in each one of them or else the code is not going to run. If you are working in a basic R script, you don't need to reset it each time.

Now that you have created a data frame summarizing everything, you can start your statistical analysis.

# Statistical analysis

## Setting the inputs and outputs

```{r variables setup}
output_data <- "~/project/tables/" ## set directory to output data
output_graphs <- "~/project/graphs/" ## set directory to output graphs

wd <- "~/project/tables/"
setwd(wd)

df <- read.csv("df.csv")
c_data <- read_excel("clinical_data.xlsx") 

df2 <- full_join(df, c_data, by = "patient") ##be sure you can merge with a column in common
```

## Generating the graphs

Now, depending on the question you want to answer you may want to generate specific graphs.

In this section, we are going to generate graphs and add statistical values on them. To do this, we are going to need to load the `ggpubr` package. It works on the same principle as `ggplot2` but it helps you to add the results of statistical tests on your graphs.

### What is the difference of infiltration of CD3 cells before and after treatment?

To answer to this question, the better representation is a box plot:

#### Comparison between two groups

Let's say you have a column with the timepoint (before - after treatment) and a column with the densities of CD3 (number of CD3 in your tumor/tumor area). 
Here we make a Wilcox test which means that we compare the mean of density of CD3 in all tissues before and after treatment. Secondly, note that we do a paired test as we have a tissue before/after of the **same** patient. Pay attention, you need to have the same patients and the same number of tissues before/after,otherwise, remove the patients for whom you miss tissues.

If you want to compare the type of response between patients (complete, partial, ...), you compare different patients so don't forget to set `paired` to <span style="color:green">FALSE</span>!

In the `stat_compare_means()` function, we can also add the argument label where we define if we want to know the exact p value (`= "p.format"`) or if we want to know if it is significant or not (`= "p.signif"`).

```{r ggboxplot}
a <- ggboxplot(df2, x="timepoint", y="dCD3", add="jitter")+
 stat_compare_means(method="wilcox.test", paired = TRUE, aes(group=timepoint),label.y.npc = 0.9, label="p.format")
```

#### Comparison between three or more groups

If you have more than two groups (for example, you compare biopsies collected at timepoint 0, 2 and 15), you need to specify the tests you want to make.

First, create your graph (same as before, defining x and y and adding points on your box plot with jitter). Then create a variable with the comparisons you want to make and then add it in the function `stat_compare_means()`.

```{r ggboxplot 2}
b <- ggboxplot(df2, x = "timepoint", y = "dCD3", add="jitter")

my_comparisons <- list( c("0", "2"), c("2", "15"), c("0", "15") )

b <- b +
  stat_compare_means(method="wilcox.test", paired = TRUE, comparisons = my_comparisons, aes(group = timepoint), label.y.npc = 0.9, label="p.format")+ 
  scale_x_discrete(breaks=c("0","2","15"),labels=c("Baseline","Week 2", "Week 15"))+
  xlab("timepoint") + ylab("Density of CD3+ cells (cells/mm²)") + 
  ggtitle("Density of CD3+ cells in a biopsy (cells/mm²)", subtitle =  "at week 0, 2 and 15") +
  theme(plot.title = element_text(hjust = 0.5), plot.subtitle = element_text(hjust = 0.5)) + 
  theme(text = element_text(size = 15))
```

### Is there a difference in the proliferation of CD3+ cells between patients who respond or not to the treatment?

Now, let's imagine you have a column with the densities of cells CD3+ RORC+ and a second column with the densities of cells CD3+ RORC- and you want to put them in on graph, you need to create a single column.

```{r pivot long}
df_long <- df2 %>%
    pivot_longer(names_to = "cell_type", ##CD3+RORC+ or CD3+RORC-
                 values_to = "density", ##numeric values
                 cols = starts_with("dCD3")) ##define the columns you want to take

##OR

df_long <- df2 %>%
    pivot_longer(cols=c('dCD3_rorc', 'dCD3_nororc'), ##other way to define the columns you want
                 names_to = "cell_type", ##CD3+RORC+ or CD3+RORC-
                 values_to = "density") ##numeric values
```

Then this can be presented in a graph:

```{r ggboxplot 3}
c <- ggboxplot(df_long, x="cell_type", y="density", add="jitter")+
 stat_compare_means(method="wilcox.test", paired = FALSE, aes(group=cell_type),label.y.npc = 0.9, label="p.format") + 
  facet_wrap(~ response) ## duplicate your graph to see the difference in the density of CD3 RORC between patients who respond or not to the treatment
```

### What is the difference in the proportions of scoring before and after treatment ?

To answer this question, you will need to create a pie chart and not a box plott.

*First step*: Create a new data frames for each condition (before/after treatment) and counting down the number of different "immune scores" you have.

```{r pie chart df}
df_bt <- data.frame(
  immune_score = c("low","intermediate","high"),
  BeforeT = c(nrow(filter(df2, timepoint=="before") %>% filter(immune_score=="low")),
               nrow(filter(df2, timepoint=="before") %>%filter(immune_score=="intermediate")),
               nrow(filter(df2, timepoint=="before") %>%filter(immune_score=="high")))
)


df_at <- data.frame(
  immune_score = c("low","intermediate","high"),
  AfterT = c(nrow(filter(df2, timepoint=="after") %>% filter(immune_score=="low")),
             nrow(filter(df2, timepoint=="after") %>%filter(immune_score=="intermediate")),
             nrow(filter(df2, timepoint=="after") %>%filter(immune_score=="high")))
)

```

*Second step*: Create the labels (percentage of the different "immune scores"), select a color palette, define the layout of the graph (basic R plot!!) and then create the two pies (before/after treatment) with your 3 categories (here immune_score).

```{r pie chart}
# Calculate the percentage of  sections and put it in the label
alabels <- round(100*df_bt$BeforeT/sum(df_bt$BeforeT), 1)
alabels <- paste(alabels, "%", sep="")

# Calculate the percentage of  sections and put it in the label
blabels <- round(100*df_at$AfterT/sum(df_at$AfterT), 1)
blabels <- paste(blabels, "%", sep="")

palette <- c("#CCCCFF" , "#FFE5CC", "#C4F5D7") ##select a color per category

layout(matrix(c(1,2,3,3), ncol=2, byrow=TRUE), heights=c(1, 1)) ##two pies = two columns
## matrix : first row two pies thus 1 and 2, second row, legend : one legend so 3, 3 (two times cause on the upper level you have 1, 2).

par(mai=rep(0.2, 4))
pie(df_bt$BeforeT, main="Before treatment",col = palette, labels=alabels, clockwise = TRUE, cex = 2, cex.main = 2)
pie(df_at$AfterT, main="After treatment", col = palette, labels=blabels, clockwise = TRUE, cex = 2, cex.main = 2)

par(mai=c(0,0,0,0))
plot.new()

legend("center", df_bt$BeforeT, fill = palette, ncol=3, cex=1.5, title="Immune score proportions before and after treatment") ## add legend
```

*Third step*: calculate the p-value. 

Here we are going to do a Fischer test because we want to see if there is a difference of distribution of the different "immune scores" between the timepoints (before and after treatment). You could do a Chi-square test but in case you don't have a lot of values, you need to do a Fischer's exact test which will give a better value (and not an estimation).

```{r stat test pie chart}
df<-df2 %>% select(immune_score, timepoint) %>% filter(timepoint=="before" | timepoint=="after")
df <- mutate(df, biopsies=1)
df$immune_score <- as.factor(df$immune_score)
class(df$immune_score)
df <- ddply(df, .(timepoint, immune_score), summarize, biopsies=sum(biopsies))
df <- xtabs(biopsies~timepoint+immune_score, data=df)
df

chisq.test(df)$expected ##Check this value: if at least one point is lower than 5 : do a Fischer's exact test

fisher.result <- fisher.test(df)
print(fisher.result$p.value)
```

Repeat this test if you have more categories (timepoints 0, 2 and 15 for example).

### What is the general overview of the densities of all the cell types in my multiplex analysis ?

For this question, you can create a heatmap. If you have different multiplex panels, you can full_join the different tables with markers together with the `full_join()` function. Then join in the clinical data. Let's continue with our first panel with CD3 and CD8 staining that we joined with the clinical data into the *df2* variable.

For your heatmap:
- select only the columns with the cell densities you want to study and the name of the biopsy
- transform the biopsy column into rownames: you need to have an assay with only numeric values
- normalize your data, here we show one way of normalization but it may not be the best one (please pay attention that the results are going to be different if you normalize on the columns or on the rows).
- draw your heatmap.

**A note on Summarized Experiment objects**: When working with heatmaps, you need to have a matrix with only numeric values but sometimes, you will need to add comments (clinical data etc). Thus, you need to link the numeric data but not putting the columns with columns inside the matrix! You can do it easily with a Summarized Experiment type of object <https://bioconductor.org/packages/release/bioc/vignettes/SummarizedExperiment/inst/doc/SummarizedExperiment.html#:~:text=SummarizedExperiment%20is%20a%20matrix-like%20container%20where%20rows%20represent,of%20a%20SummarizedExperiment%20object%20represent%20features%20of%20interest.>. It is usually used with sequencing datas.

```{r heatmap}
df <- df2[,c("biopsy", "dCD3", "dCD4", "dCD8", "dCK")]

df <- column_to_rownames(df, "biopsy")
boxplot(df)

df1 <- normalize(df, method = "standardize")
boxplot(df1) ##just to see the overall distribution

ra <- as.matrix(df1) ##transform you assay into a matrix
rb <- na.omit(ra) ##remove all NA values

#Basic heatmap

ggp1 <- heatmap(rb) 
ggp1

rf <- t(rb) ##transpose the table

ggp2 <- Heatmap(rf, border = TRUE) ##heatmap the other way round
ggp2 ##best to have the cell's densities in rows and patients/biopsies in columns

# Complex heatmap

##create a data frame with the clinical data etc
colData <- df2 %>%
  select(timepoint, patient)

##Create a Summarized experiment (link a matrix and coldata)
se <- SummarizedExperiment(assays = rf, colData = colData)

ggpa <- Heatmap(assay(se), column_km = 3, border = TRUE) ##independant clustering
##column_km = 3 creates 3 clusters
ggpa 

colnames(assay(se)) ##check colnames: should be the same + same length as those in annotation df!

##create a variable with annotations (ha1 for treatment, ha2 for ...)
ha1 <- HeatmapAnnotation(treatment = df2$treatment, 
                         col = list(treatment = c("immunotherapy" = "#010000", "standard_of_care" = "#ffffff")))

ggpb <- Heatmap(assay(se), name = "density", column_split = colData(se)[1], border = TRUE, top_annotation = c(ha1))
                  
ggpb ##is there a difference according to the timepoint ? cause we split according to timepoint: colData(se)[1]

ggpc <- Heatmap(assay(se), column_split = colData(se)[2], border = TRUE, top_annotation = c(ha1))
ggpc ##split the graph according to patients: colData(se)[2]

ggpd <- Heatmap(assay(se), border = TRUE,top_annotation = c(ha1))
ggpd ##unsupervised clustering

ggpe <- Heatmap(assay(se), column_km = 3, border = TRUE, top_annotation = c(ha1))
ggpe 

```

Now you saw the different graphs you could generate with your data and you can replicate the principle to answer all of your questions regarding the variation of densities of you cells according to the timepoint or the response to the treatment or the type of treatment received, ...
