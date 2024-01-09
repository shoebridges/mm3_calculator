#
# This is a Shiny web application. You can run the application by clicking
# the 'Run App' button above.
#
# Find out more about building applications with Shiny here:
#
#    http://shiny.rstudio.com/
##

library(tidyverse)
library(data.table)
library(stats)
library(caret)
library(drc)
library(magic)
options(scipen = 999)

# Define UI for application that draws a histogram
ui <- fluidPage(
  
  # Application title
  titlePanel("Smart Lipid Assay (Verena)"),
  
  # Sidebar with a slider input for number of bins 
  sidebarLayout(
    sidebarPanel(numericInput("replicates", "Number of replicates", value = 3),
                 textAreaInput("groups", "Insert the group name in order, seperated by space (e.g. A A)"),
                 fileInput("data", "Insert the excel file of your data (It will take the first sheet)", accept = ".xlsx"),
                 h3("Once settings and the dataset have been entered, click Perform fit and calculate"),
                 actionButton("Calculate", "Perform fit and calculate"),
                 h2(),
                 downloadButton("corr_out", "Download the table of correlation to predicted fit"),
                 downloadButton("plot_out", "Download the fit to values plot in pdf"),
                 downloadButton("stats_out", "Download the Anova and pairwise (TukeyHSD) of the coefficents a and b"),
                 downloadButton("coeff_out", "Download the coefficents a and b"),
                 downloadButton("data_out", "Download dataset with predicted values based on regression")
    ),
    
    # Show tables for correlation and pairwise tests of coefficients, plot of predicted values and recorded values
    mainPanel(
      h1("Introduction and settings"),
      p("The following application is for fitting a line of best fit and calculation of coefficents for the lipid assays performed by Verena and Sonja.
             Statistical tests are performed for the coefficents based on the grouping of the samples. The application requires an excel sheet with colnames where the first column is labelled Time_m
             which represents time in minutes."),
      h1("Formula Info: MM.2/MM.3"),
      p("The model is defined by the three-parameter model function:",tags$br(),tags$br(), "f(x, (c, d, e)) = c + \frac{d-c}{1+(e/x)}",tags$br(),tags$br(), "It is an increasing as a function of the dose x, 
      attaining the lower limit c at dose 0 (x=0) and the upper limit d for infinitely large doses. 
      The parameter e corresponds to the dose yielding a response halfway between c and d.
      The common two-parameter Michaelis-Menten model (MM.2) is obtained by setting c equal to 0."),
p("
    To Verena, Since your controls do weird s##t, I would not set to 0 and instead implement the three-parameter model.
    All curves run off the same model but starting is assumed independently.
    C in a sense would be your lowest Absorbance, D would be the largest absorbance, E would be the time taken
    to reach halfway between C and D",tags$br(),tags$br(),
        ),
      dataTableOutput("correlation"),
      p("Here we display a plot of the recorded values and the trendline of predicted values"),
      plotOutput("plot"),
      h1("Coefficents"),
      p("The table below contains the coefficents predicted for each sample"),
      dataTableOutput("coeff"),
      h2("ANOVA (Tukey)"),
      p("The table below contains the ANOVA statistical comparisons for coefficents C, D and E. Ignore the first column or row names, its just a row order."),
      dataTableOutput("stats"),
      h2("t-test (Bonferroni)"),
      p("The table below contains a t-test using adjusted p value for coefficent C"),
      tableOutput("Bonferroni_C"),
      p("The table below contains a t-test using adjusted p value for coefficent D"),
      tableOutput("Bonferroni_D"),
      p("The table below contains a t-test using adjusted p value for coefficent E"),
      tableOutput("Bonferroni_E"),
      h2("Plot of Coefficents"),
      p("The graph below contains the plots of the coefficents C,D and E"),
      plotOutput("coeff_plot")
    )
  ))

# Define server logic required to draw a histogram
server <- function(input, output) {
  
  labels <- eventReactive(input$Calculate, {
    if (is.null(input$data)) return()
    groups <- str_split_1(input$groups, " ")
    labels <- rep(groups, input$replicates)
    labels <- sort(factor(labels, groups))
    labels <- paste(labels, seq(1:input$replicates), sep = "_")
    return(labels)
  })
  groups <- eventReactive(input$Calculate, {
    if (is.null(input$data)) return()
    groups <- str_split_1(input$groups, " ")
    labels <- rep(groups, input$replicates)
    labels <- sort(factor(labels, groups))
    return(labels)
  })
  file <- eventReactive(input$Calculate, {
    file <- input$data
    #turn off when not active
    if (is.null(input$data)) return()
    if (is.null(labels())) return()
    if (is.null(groups())) return()
    datapath <- file$datapath
    dataset <- readxl:: read_xlsx(datapath, sheet=1)
    dataset
  })
  melt_dataset <- eventReactive(input$Calculate, {
  #Melting dataframe
  df <- file()
  #if (length(labels()!=length(names(dataset)))) return()
  names(df) <- append("Time_m", labels())
  if (is.null(input$data)) return()
  if (is.null(labels())) return()
  if (is.null(groups())) return()
  melt_df <- as.data.frame(melt.data.table(data.table(df), id.vars = "Time_m", variable.name = "Treatment", value.name = "Absorbance", na.rm = T))
  #drm fit
  })
  model <- eventReactive(input$Calculate, {
    drm(Absorbance~Time_m, data = melt_dataset(), curveid = Treatment, fct = MM.3(), separate = T)})
  
  coeff <- eventReactive(input$Calculate, {
  coeff_test <- as.data.frame(model()$coefficients)
  coeff_test[, c("Coefficent", "Treatment")] <- str_split_fixed(rownames(coeff_test), ":", 2)
  names(coeff_test) <- c("Value", "Coefficent", "Treatment")
  labels <- data.frame(labels(), groups())
  coeff_test$group <- labels$groups[match(coeff_test$Treatment, labels$labels)]
  coeff_test
  })
  
  correlation <- eventReactive(input$Calculate, {
  #Predict Values
  #Correlate
  melt_df <- melt_dataset()
  model <- model()
  melt_df$Pred <-predict(model, CURVE=Treatment)
  corr_test <- melt_df %>%
    group_by(Treatment) %>%
    summarize(cor=cor(Absorbance, Pred));corr_test
  #NOTE: This will tell you the correlation between your recorded values and the fit.
  })
  
  #Plot values
  plot <- eventReactive(input$Calculate, {
    melt_df <- melt_dataset()
    model <- model()
    melt_df$Pred <-predict(model, CURVE=Treatment)
    plot1 <- ggplot(melt_df, aes(x=Time_m, y=Absorbance, color=Treatment))+geom_point(alpha=0.05)+geom_line(aes(y=Pred))+facet_wrap(~Treatment)
    plot1
  })
  stat_test <- eventReactive(input$Calculate, {
    #turn off when not active
    if (is.null(file)) return()
    coeff <- coeff()
    anova_c <- aov(formula=Value~group, data=subset(coeff, Coefficent=="c"))
    Pairwise_c<- as.data.frame(TukeyHSD(anova_c)$group)
    Pairwise_c$Co_Test <- "c"
    Pairwise_c$Comparison <- rownames(Pairwise_c)
    #
    anova_d <- aov(formula=Value~group, data=subset(coeff, Coefficent=="d"))
    Pairwise_d <- as.data.frame(TukeyHSD(anova_d)$group)
    Pairwise_d$Co_Test <- "d"
    Pairwise_d$Comparison <- rownames(Pairwise_d)
    
    anova_e <- aov(formula=Value~group, data=subset(coeff, Coefficent=="e"))
    Pairwise_e <- as.data.frame(TukeyHSD(anova_e)$group)
    Pairwise_e$Co_Test <- "e"
    Pairwise_e$Comparison <- rownames(Pairwise_e)
    
    pairwise <- rbind(Pairwise_c, Pairwise_d, Pairwise_e)
    rownames(pairwise) <- seq(1:nrow(pairwise))
    pairwise
  })
  
  Bonferroni_C <- eventReactive(input$Calculate, {
    if (is.null(file)) return()
    data <- as.data.frame(subset(coeff(), Coefficent=="c"))
    pairwise_c <- pairwise.t.test(data$Value, data$group, p.adjust.method = "bonferroni")$p.value
    pairwise_c
  })
  Bonferroni_D <- eventReactive(input$Calculate, {
    #turn off when not active
    if (is.null(file)) return()
    data <- as.data.frame(subset(coeff(), Coefficent=="d"))
    pairwise_d <- pairwise.t.test(data$Value, data$group, p.adjust.method = "bonferroni")$p.value
    pairwise_d
  })
  Bonferroni_E <- eventReactive(input$Calculate, {
    #turn off when not active
    if (is.null(file)) return()
    data <- as.data.frame(subset(coeff(), Coefficent=="e"))
    pairwise_e <- pairwise.t.test(data$Value, data$group, p.adjust.method = "bonferroni")$p.value
    pairwise_e
  })
  output$Bonferroni_C <- renderTable(Bonferroni_C(), rownames = TRUE)
  output$Bonferroni_D <- renderTable(Bonferroni_D(), rownames = TRUE)
  output$Bonferroni_E <- renderTable(Bonferroni_E(), rownames = TRUE)
  output$correlation <- renderDataTable(correlation())
  output$plot <- renderPlot(plot())
  output$coeff <- renderDataTable(coeff())
  output$stats <- renderDataTable(stat_test())
  output$coeff_plot <- renderPlot(ggplot(data=as.data.frame(coeff()), aes(x=group, y=Value, fill=group))+
                                     stat_summary(fun="mean", geom="bar")+geom_point()+theme_bw()+labs(x="group", y="Value", title = "Value of Coefficents", subtitle = "Absorbance=a*Time/(b+Time) | Treatment)")+facet_wrap(~Coefficent, ncol=1))
  #Save files
  output$plot_out <- downloadHandler(filename = function() {paste("corrplot_", Sys.Date(), ".pdf", sep = "")}, content = function(file) {ggsave(filename = file, plot = plot(), device = "pdf", units = "in", width = 12, height = 6)})
  output$corr_out <- downloadHandler(filename = function() {paste("corr_out", Sys.Date(), ".csv", sep = "")}, content = function(file) {write.csv(correlation(), file)})
  output$coeff_out <- downloadHandler(filename = function() {paste("coeff", Sys.Date(), ".csv", sep = "")}, content = function(file) {write.csv(coeff(), file)})
  output$stats_out <- downloadHandler(filename = function() {paste("stats_test",Sys.Date(), ".csv", sep = "")}, content = function(file) {write.csv(stat_test(), file)})
  output$data_out <- downloadHandler(filename = function() {paste("data_out",Sys.Date(), ".csv", sep = "")}, content = function(file) {write.csv(melt_dataset(), file)})
  
}

# Run the application 
shinyApp(ui = ui, server = server)