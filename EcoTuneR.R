#source('https://raw.githubusercontent.com/LewisLab-URI/EcoTuneR/main/EcoTuneR.R')
#combined_output <- read.csv("https://raw.githubusercontent.com/LewisLab-URI/EcoTuneR/main/Example%20Datasets/Example_combined_output.csv")

#ecotuneR(combined_output,F,T)

if (!require("svDialogs")) {
  cat(file=stderr(),"The package 'svDialogs' is required \n")
  install.packages("svDialogs")
  library(svDialogs)
  cat(file=stderr(),"'svDialogs' successfully loaded. \n")
} else {cat(file=stderr(),"'svDialogs' successfully loaded. \n")}


if (!require("ggplot2")) {
  cat(file=stderr(),"The package 'ggplot2' is required \n")
  if(ok_cancel_box("The package 'ggplot2' is required \n Do you want to install 'ggplot2'?")){
    install.packages("ggplot2")
    library(ggplot2)
    cat(file=stderr(),"'ggplot2' successfully loaded. \n")
  } else {
    msg_box("Required package 'ggplot2' is not installed. \nEcoTuneR is not available without 'ggplot2'. \n")
    cat("Required package 'ggplot2' is not installed. \nEcoTuneR is not available without 'ggplot2'. \n")
    stop("'ggplot2' not installed \n")
  }
} else {cat(file=stderr(),"'ggplot2' successfully loaded. \n")}
if (!require("bs4Dash")) {
  cat(file=stderr(),"The package 'bs4Dash' is required \n")
  if(ok_cancel_box("The package 'bs4Dash' is required \n Do you want to install 'bs4Dash'?")){
    install.packages("bs4Dash")
    library(bs4Dash)
    cat(file=stderr(),"'bs4Dash' successfully loaded. \n")
  } else {
    msg_box("Required package 'bs4Dash' is not installed. \nEcoTuneR is not available without 'bs4Dash'. \n")
    cat("Required package 'bs4Dash' is not installed. \nEcoTuneR is not available without 'bs4Dash'. \n")
    stop("'bs4Dash' not installed \n")
  }
} else {cat(file=stderr(),"'bs4Dash' successfully loaded. \n")}
if (!require("shiny")) {
  cat(file=stderr(),"The package 'shiny' is required \n")
  if(ok_cancel_box("The package 'shiny' is required \n Do you want to install 'shiny'?")){
    install.packages("shiny")
    library(shiny)
    cat(file=stderr(),"'shiny' successfully loaded. \n")
  } else {
    msg_box("Required package 'shiny' is not installed. \nEcoTuneR is not available without 'shiny'. \n")
    cat("Required package 'shiny' is not installed. \nEcoTuneR is not available without 'shiny'. \n")
    stop("'shiny' not installed \n")
  }
} else {cat(file=stderr(),"'shiny' successfully loaded. \n")}

ecotuneR <- function(combined_output=NA,export=FALSE, bypass_check= FALSE){
  #check if data is valid - if so, run this. if not, don't run it
  ready <- 0
  if(!isTruthy(combined_output) || !is.data.frame(combined_output)){
    try(
      cat(file=stderr(),"'",substitute(combined_output),"' is not a valid dataframe input for ecotuneR.\n", sep=""),
      stop())
  } else {
    if(bypass_check==FALSE){
      cat(file=stderr(),"ecotuneR requires specifically formatted input to work properly.\n", 
          "  Column 1:  Consumer\n  Column 2:  Lower Bounds\n  Column 3:  Upper Bounds\n",
          "  Column 4:  SIA Mode/Trophic Position\n  Column 5:  Ecopath Trophic Level\n",
          "\nPlease verify that '",substitute(combined_output),
          "' data are in the appopriate columns.\n \n", sep="")
      print(head(combined_output))
      msg_box(paste("ecotuneR requires specifically formatted input to work properly.\nPlease verify that '",substitute(combined_output),"' data are in the appropriate columns.\n \n", sep=""))
      if(ok_cancel_box(paste0("Please see the Console and verify that'",substitute(combined_output),"' data are in the appropriate columns."))) {
        ready = 1
      }else{
        #reformat
        cat("\n Please reformat your data")
      }
      
      
    }
    if(bypass_check==TRUE){ready <- 1}
  }
  
  if(ready == 1){
    #Run the app
    
    #Background calculations
    {
      pd <- combined_output[, c(1, 2, 3, 4, 5)]
      colnames(pd) <- c("Consumer","PD.L","PD.U","PD.M","TL") 
      ref <- list()
      ref$w_TL <- pd[!is.na(pd$TL),] #pd with Trophic Level data
      ref$wo_TL <- pd[is.na(pd$TL),] #pd without Trophic Level data
      
      ref$ci <- subset(ref$w_TL, TL <= PD.U & TL >= PD.L )# %>%nrow() Within Confidence Interval 
      ref$hi <- subset(ref$w_TL, TL > PD.U)# %>% nrow() Above Confidence Interval 
      ref$lo <- subset(ref$w_TL, TL < PD.L )# %>% nrow() Below Confidence Interval 
      
      
      
      
      
      
      
      if(nrow(ref$wo_TL) >0){
        content_pies <- fluidRow(column(6,plotOutput("piePlot1")),column(6,plotOutput("piePlot2")))
      } else { 
        content_pies <- plotOutput("piePlot2")}
      
      col1 <- (nrow(ref[[2]]) >0)
      col2 <- (nrow(ref[[3]]) >0)
      col3 <- (nrow(ref[[4]]) + nrow(ref[[5]]) >0)
      colnum<-sum(c(col1,col2,col3))
      if(col1){
        col1 <- {column(12/colnum,align = "center",
                        h4("Without DMTL"),
                        HTML(
                          paste(unlist(
                            lapply(seq_along(ref$wo_TL[,1]), function(j) {
                              paste0(
                                '<a class="action-button" href="#" id="without_',j,'">',
                                gsub("_"," ",paste(ref$wo_TL[j,1],collapse = "<br>")),
                                '</a><br>')
                            })
                          ),collapse = '<div class="menu-dropdown-divider"></div>')
                        ))}
      } else {col1 <- list()}
      if(col2){
        col2 <- {column(12/colnum,align = "center",
                        h4("Functional Groups where TL is within SIA 95% CI"),
                        HTML(
                          paste(unlist(
                            lapply(seq_along(ref$ci[,1]), function(j) {
                              paste0(
                                '<a class="action-button" href="#" id="within_',j,'">',
                                gsub("_"," ",paste(ref$ci[j,1],collapse = "<br>")),
                                '</a><br>')
                            }) 
                          ),collapse = '<div class="menu-dropdown-divider"></div>')
                          
                        ))}
      } else {col2 <- list()}
      if(col3){
        if(nrow(ref[[4]]) >0 && nrow(ref[[5]]) >0) {
          col3 <- {  column(12/colnum,align = "center",
                            h4("Overestimated "),
                            HTML(
                              paste(unlist(
                                lapply(seq_along(ref$hi[,1]), function(j) {
                                  paste0(
                                    '<a class="action-button" href="#" id="above_',j,'">',
                                    gsub("_"," ",paste(ref$hi[j,1],collapse = "<br>")),
                                    '</a><br>')
                                }) 
                              ),collapse = '<div class="menu-dropdown-divider"></div>')
                            ),
                            br(),
                            h4("Underestimated "),
                            HTML(
                              paste(unlist(
                                lapply(seq_along(ref$lo[,1]), function(j) {
                                  paste0(
                                    '<a class="action-button" href="#" id="below_',j,'">',
                                    gsub("_"," ",paste(ref$lo[j,1],collapse = "<br>")),
                                    '</a><br>')
                                }) 
                              ),collapse = '<div class="menu-dropdown-divider"></div>')
                              
                            ))}
        } else if(nrow(ref[[4]]) >0) {
          col3 <- {  column(12/colnum,align = "center",
                            h4("Overestimated "),
                            HTML(
                              paste(unlist(
                                lapply(seq_along(ref$hi[,1]), function(j) {
                                  paste0(
                                    '<a class="action-button" href="#" id="above_',j,'">',
                                    gsub("_"," ",paste(ref$hi[j,1],collapse = "<br>")),
                                    '</a><br>')
                                }) 
                              ),collapse = '<div class="menu-dropdown-divider"></div>')
                            ))}
        } else if(nrow(ref[[5]]) >0) {
          col3 <- {  column(12/colnum,align = "center",
                            h4("Underestimated "),
                            HTML(
                              paste(unlist(
                                lapply(seq_along(ref$lo[,1]), function(j) {
                                  paste0(
                                    '<a class="action-button" href="#" id="below_',j,'">',
                                    gsub("_"," ",paste(ref$lo[j,1],collapse = "<br>")),
                                    '</a><br>')
                                })  
                              ),collapse = '<div class="menu-dropdown-divider"></div>')
                            ))}
        }
      } else {col3 <- list()}
      
      content_lists <- {fluidRow(
        col1,
        col2,
        col3,
      )}
      rm(colnum,col1,col2,col3)
    }
    distPlot <<- function(match){
      #makes a graph from pd of the consumer specified by "match"
      {
        pd_subset <- subset(pd,Consumer == match) 
        mydf <- data.frame("consumer" = rep(pd_subset[1,1],4),
                           "name" = as.character(names(pd_subset[1,c(2:5)])),
                           "value" = as.numeric(pd_subset[1,c(2:5)]))
        
        mydf[mydf == "PD.L" | mydf == "PD.U"] <- "bounds"
        xmin <- min(mydf$value)
        xmax <- max(mydf$value)
      } # This sets up the data for the ggplot
      { 
        ggplot(data = mydf[c(1,2,3),]) +
          geom_density(aes(value), alpha = 0.4, fill = "slateblue") +
          theme_classic() +
          theme(legend.position = c(.9, .9)) +
          theme(legend.key = element_rect(fill = "grey90")) +
          geom_vline(data = mydf, aes(xintercept = value, #x is TL
                                      color = name,
                                      linetype = name,
                                      size = name)) +
          
          scale_size_manual(name= "         ",
                            labels=c("SIA 95% CI",
                                     "SIA Mode/Trophic Position",
                                     "Ecopath Trophic Level"),
                            values=c(.6,.6,1)) +
          scale_linetype_manual(name= "         ",
                                labels=c("SIA 95% CI", 
                                         "SIA Mode/Trophic Position",
                                         "Ecopath Trophic Level"),
                                values=c("dashed", "dashed", "dashed")) +
          scale_color_manual(name= "         ",
                             labels=c("SIA 95% CI",
                                      "SIA Mode/Trophic Position",
                                      "Ecopath Trophic Level"),
                             values=c("gray30", "white", "red")) +
          labs(x="Probability Distribution", y="Density", fill = "Legend") +
          
          expand_limits(x=c((xmin-(.8*(xmax-xmin))),
                            (xmax+(.8*(xmax-xmin)))
          )) +
          theme(axis.title.x = element_text(size=16), 
                axis.title.y = element_text(size=16),
                title = element_text(size=19),
                legend.title = element_blank(),
                legend.text = element_text(size = 13))+
          theme(plot.margin = margin(0,2,0,.5, "cm"))
        
      } #this actually does the graphy thingy
    }
    piePlot2 <<- function(){
      {
        fitdf <- data.frame(
          group = c("Above 95% CI", "Below 95% CI","Within 95% CI"),
          value = c(nrow(isolate(ref$hi)),
                    nrow(isolate(ref$lo)),
                    nrow(isolate(ref$ci)))
        )
        ggplot(fitdf, aes(x="", y=value, fill=group))+
          geom_bar(width = 1, stat = "identity") +
          coord_polar("y", start=0) +
          geom_text(aes(label = paste0(round(value / sum(value) * 100, 1), "%"), x = 1.1),
                    size = 5, position = position_stack(vjust = 0.5)) + 
          theme_void() +
          theme(legend.title=element_blank()) +
          scale_fill_manual(values=c("#fdcc8a","#b2abd2","#2b8cbe"))
        
      } #pie chart time
    }
    # Define UI
    ui = bs4Dash::dashboardPage(
      options = NULL,
      header = {dashboardHeader(
        title = dashboardBrand(
          title = "Eco-tuneR",
          color = "white",
          href = "https://web.uri.edu/lewis-lab/",
          image = "https://raw.githubusercontent.com/LewisLab-URI/EcoTuneR/main/www/ecotuner2.png"
        ),
        sidebarIcon  = shiny::icon("lemon",verify_fa=F),
        span(img(src = "https://raw.githubusercontent.com/LewisLab-URI/EcoTuneR/main/www/ecotuner0.png", height = 40, width = 115, align = "right",name = "wefwe")),
        uiOutput("headerdur")
      )},
      sidebar = {dashboardSidebar(
        collapsed = T,startExpanded = F,
        sidebarMenu(
          id = "sidebarMenu",
          menuItem(
            text = "Summary",
            tabName = "tab1",
            icon = icon("van-shuttle")
          ),
          menuItem(
            text = "Details",
            tabName = "tab2",
            icon = icon("shuttle-space")
          ),
          menuItem(
            text = "Next Steps",
            tabName = "tab3",
            icon = icon("sliders")
          )
        )
      )},
      footer = {dashboardFooter(
        span(
          a(href = "https://www.lsu.edu/",
            img(src = "https://raw.githubusercontent.com/LewisLab-URI/EcoTuneR/main/www/LSU.png", height = 40, width = 90, align = "right")),
          a(href = "https://www.rhodes.edu/",
            img(src = "https://raw.githubusercontent.com/LewisLab-URI/EcoTuneR/main/www/RC.png", height = 40, width = 90, align = "right")),
          a(href = "https://www.nationalacademies.org/gulf/gulf-research-program",
            img(src = "https://raw.githubusercontent.com/LewisLab-URI/EcoTuneR/main/www/NASEM.png", height = 40, width = 100, align = "right")),
          h6("This research was funded by the National Academy of Sciences Gulf Research Program", 
             align = "left") 
        )
      )},
      body = {dashboardBody(
        tabItems(
          
          {tabItem(
            tabName = "tab1",
            h1("Summary"),
            h3("In place of unpublished data, this tool currently references a randomly-generated example dataset. This is meant for demonstration only.", style = "color: red;"),
            content_pies,
            h3(strong(paste0(round((nrow(ref$ci)/nrow(ref$w_TL))*100,digits = 0),"% of all available EwE TLs fit within 95% CI of SIA TLs")), align = "center"),
            
            hr(),
            content_lists
            
          )}, # Tab 1 - Summary
          {tabItem(
            tabName = "tab2",
            h1("Details"),
            box(width = 12,align = "center",
                fluidRow(
                  column(2,actionButton("summary","Summary",icon("van-shuttle"))),
                  column(5,selectInput("consumer","Functional Groups",unique(pd$Consumer))),
                  column(2,br(),actionButton("back","Back",icon("backward"))),
                  column(2,br(),actionButton("forward","Next",icon("forward")))
                ), plotOutput(outputId = 'distPlot'))
          )}, # Tab 2 - Details
          {tabItem(
            tabName = "tab3",
            box(width = 12,status = "primary",title=h3("Considerations for interpreting results"),
                collapsible = F,
                h6("EcoTuneR is a visualization tool and certain trends within the data may not be apparent when using this tool.
Understanding the limitations of using Trophic Level data to validate an ecosystem model is important for a balanced interpretation of results.
Please take time to thoroughly analyze results, and refer to the suggestions found in the table below.
Species probability distributions are approximations due to the constraints of model outputs."),
            
                h5("Trophic level is within bounds of 95% CI"),
                tableOutput('nextSteps1'),
                br(),
                h5("Trophic level is outside bounds of 95% CI"),
                tableOutput('nextSteps2')#,
                #  br()
            ),
            
            
            HTML("\n\n\n"),
            br(),
            
            
            
            
            box(width = 12,title="Citations",collapsible = F,
                h6("Granjon D (2022). bs4Dash: A 'Bootstrap 4' Version of 'shinydashboard'. \nhttps://rinterface.github.io/bs4Dash/index.html, https://github.com/RinteRface/bs4Dash."),
                h6("Wickham H (2016). ggplot2: Elegant Graphics for Data Analysis. \nSpringer-Verlag New York. ISBN 978-3-319-24277-4, https://ggplot2.tidyverse.org."),
                h6("Chang W, Cheng J, Allaire J, Sievert C, Schloerke B, Xie Y, Allen J, McPherson J, Dipert A, Borges B (2023). shiny: Web Application Framework for R. R package version 1.7.4.9002, https://shiny.rstudio.com/."))
          )}  # Tab 3 - Next Steps    
          
        ))},
      title = "DashboardPage"
    )
    
    
    server <- function(input, output) {
      
      observeEvent(input$summary, {
        updateTabsetPanel(
          session = shiny::getDefaultReactiveDomain(),
          "sidebarMenu",
          selected = "tab1"
        )
      })
      observeEvent(input$forward, {
        updateSelectInput(getDefaultReactiveDomain(),inputId = "consumer",selected = pd[which(pd[,1] == input$consumer)+1,1])
      })
      observeEvent(input$back, {
        updateSelectInput(getDefaultReactiveDomain(),inputId = "consumer",selected = pd[which(pd[,1] == input$consumer)-1,1])
      })
      
      output$piePlot1 <- renderPlot({
        {
          fitdf <- data.frame(
            group = c("TL data","No TL data"),
            value = c(nrow(isolate(ref$w_TL)),
                      nrow(isolate(ref$wo_TL)))
          )
          ggplot(fitdf, aes(x="", y=value, fill=group))+
            geom_bar(width = 1, stat = "identity") +
            coord_polar("y", start=0) +
            geom_text(aes(label = paste0(round(value / sum(value) * 100, 1), "%"), x = 1.1),
                      size = 5, position = position_stack(vjust = 0.5)) +
            theme_void() +
            theme(legend.title=element_blank()) +
            scale_fill_manual(values=c("#999999","#2b8cbe"))
          
        } #pie chart time
      }, bg="transparent")
      output$piePlot2 <- renderPlot({
        piePlot2()
      }, bg="transparent")
      output$distPlot <- renderPlot({distPlot(input$consumer)}, bg="transparent")
      
      output$nextSteps1 <- renderTable({ 
        data.frame("Question"=c(
          "Are any of the credibility interval ranges larger than 5?",
          "Do all groups have an adequate sample size within the stable isotope data?"
        ),"Consideration"=c(
          "If yes, large credibility intervals increase the chance of a trophic level falling within the bounds. This may not be a result of agreement between the two methods.",
          "If no, small sample sizes of any relevant groups may lead to large credibility intervals (see above)."))
      })
      output$nextSteps2 <- renderTable({
        data.frame("Question"=c(
          "If the system is detritus based, do any species rely heavily on detritus as a part of their diet? (>50%)",
          "Were the stable isotope samples muscle or bone? How much diet information was integrated into the samples?",
          "Were the SIA and EwE data collected in generally the same time frame (i.e., year, decade)?",
          "What season were the SIA and EwE samples collected in? Is there any mismatch (i.e., summer vs fall)?"
        ),"Consideration"=c(
          "If yes, the detritus trophic level of 1 (EwE standard) may be offsetting species trophic levels to be lower than they are in the natural system.",
          "If muscle, these samples only integrate about 3 months of diet information and may not agree with EwE data collected over a longer period of time. Conversely, bone sample integrate up to 1 year of diet information",
          "If no, there may be lack of agreement in the data due to changes over time in the system of interest (i.e., fishing, construction, climate change)",
          "If yes, there may be lack of agreement in the data if data were collected in different seasons because of ecological changes in systems over the course of a year."))
      })
      
      #each of the consumer/age combos gets an individual button. can be optimized to be one picker input with all options. for now, this works.
      lapply(1:nrow(isolate(ref$ci)),function(j){
        within_com <- paste0("within_",j)
        observeEvent(input[[within_com]], {
          updateSelectInput(getDefaultReactiveDomain(),inputId = "consumer",selected = ref$ci[j,1])
          updateTabsetPanel(
            session = shiny::getDefaultReactiveDomain(),
            "sidebarMenu",
            selected = "tab2"
          )
        })
      })
      lapply(1:nrow(isolate(ref$hi)),function(j){
        above_com <- paste0("above_",j)
        observeEvent(input[[above_com]], {
          updateSelectInput(getDefaultReactiveDomain(),inputId = "consumer",selected = ref$hi[j,1])
          updateTabsetPanel(
            session = shiny::getDefaultReactiveDomain(),
            "sidebarMenu",
            selected = "tab2"
          )
        })
      })
      lapply(1:nrow(isolate(ref$lo)),function(j){
        below_com <- paste0("below_",j)
        observeEvent(input[[below_com]], {
          updateSelectInput(getDefaultReactiveDomain(),inputId = "consumer",selected = ref$lo[j,1])
          updateTabsetPanel(
            session = shiny::getDefaultReactiveDomain(),
            "sidebarMenu",
            selected = "tab2"
          )
        })
      })
      lapply(1:nrow(isolate(ref$wo_TL)),function(j){
        without_com <- paste0("without_",j)
        observeEvent(input[[without_com]], {
          updateSelectInput(getDefaultReactiveDomain(),inputId = "consumer",selected = ref$wo_TL[j,1])
          updateTabsetPanel(
            session = shiny::getDefaultReactiveDomain(),
            "sidebarMenu",
            selected = "tab2"
          )
        })
      })
    }
    
    # Run the application 
    if(export){
      cat(file=stderr(),"Exporting summary and plots")
      add_color <- function(x) {
        if (x["tl"] >= x["lower"] & x["tl"] <= x["upper"]) {
          paste0("<tr style='background-color: #c8e6c9;'>",
                 paste0("<td style='padding:10px;border: 1px solid #ddd;'>", x, "</td>", collapse = ""),
                 "</tr>")
        } else {
          paste0("<tr style='background-color: #ffcdd2;'>",
                 paste0("<td style='padding:10px;border: 1px solid #ddd;'>", x, "</td>", collapse = ""),
                 "</tr>")
        }
      }
      {
        # Create an HTML legend
        html_legend <- paste("<p style='text-align:center;'><strong>Legend:</strong></p><ul style='list-style-type:none;padding:0;margin:0;display:flex;justify-content:center;align-items:center;'><li style='margin-right:20px;'><span style='display:inline-block;height:20px;width:20px;background-color:#c8e6c9;border-radius:50%;'></span> TL within lower and upper bounds</li><li style='margin-right:20px;'><span style='display:inline-block;height:20px;width:20px;background-color:#ffcdd2;border-radius:50%;'></span> TL outside lower and upper bounds</li></ul>")
        # Create an HTML table 
        html_table <- paste0("<table style='width:80%;margin:0 auto;border-collapse: collapse;box-shadow: 0 0 20px rgba(0,0,0,0.15);'>",
                             "<thead style='background-color: #f2f2f2;text-align:center;'><tr><th style='padding:10px;border: 1px solid #ddd;'>Consumer</th><th style='padding:10px;border: 1px solid #ddd;'>Lower</th><th style='padding:10px;border: 1px solid #ddd;'>Upper</th><th style='padding:10px;border: 1px solid #ddd;'>Mode</th><th style='padding:10px;border: 1px solid #ddd;'>TL</th></tr></thead>",
                             "<tbody style='text-align:center;'>",
                             paste0(apply(combined_output, 1, add_color), collapse = ""),
                             "</tbody></table>")
      } #HTML nonsense
      # Summary Plot 
      {
        plottyplotplot <- piePlot2()
        ggsave(
          paste0("Exported Results of ",substitute(combined_output),"/Summary Plot.png"),
          plot = plottyplotplot,
          scale = 2,
          width = 700,
          height = 400,
          units = "px",
          dpi = 200,
          limitsize = F
        )} #Save the plots
      html_element<- {paste(
        h1("Summary of ",substitute(combined_output)),
        h3(strong(paste0(round((nrow(ref$ci)/nrow(ref$w_TL))*100,digits = 0),"% of all available EwE TLs fit within 95% CI of SIA TLs")), align = "center"),
        paste0(HTML("<img src=\"Summary Plot.png\" style=\"display: block; margin: 0 auto;\" width=\"40%\"/>\n")),
        br(),
        content_lists,
        tags$style({HTML('/* Create two equal columns that floats next to each other */
.col-sm-6 {
  float: left;
  width: 48%;
  padding: 10px;
}

/* Clear floats after the columns */
.row:after {
  content: "";
  display: table;
  clear: both;
}')}),
        hr(),
        
        # Combine the legend and table HTML code
        paste0("<div style='text-align:center;'>", html_legend, html_table, "</div>")
      )}
    
    write(html_element, paste0("Exported Results of ",substitute(combined_output),"/A Summary of ",substitute(combined_output),".html"))
    for(w in 1:nrow(pd)){
      plottyplotplot <- distPlot(pd[w,1]) 
      
      file_name <- paste0(pd[w,1],".png")
      file_name <- gsub("[^[:alnum:]._[:space:]]", ".", file_name)
      file_name <- paste0("Exported Results of ",substitute(combined_output),"/",file_name)
      {ggsave(
        file_name,
        plot = plottyplotplot,
        scale = 1,
        width = 2700,
        height = 1750,
        units = "px",
        limitsize = F
      )} #Save the plots
      cat(".")
    }

    cat(file=stderr(),"\nSummary and Plots exported \n")
    browseURL(paste0("Exported Results of ",substitute(combined_output)))
    
    } else {
      shinyApp(ui = ui, server = server)
    }
  
  }
}
#ecotuneR(pid)


