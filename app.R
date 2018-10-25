library(shiny)
library(ggplot2)
library(grid)
library(shinydashboard)
library(reshape)
library(reshape2)
library(stringr)
library(gridExtra)
library(gplots)
library(svglite)
library(d3heatmap)
library(shinyjs)

load_data <- function() {
  hide("loading_page")
  show("mainTabsetPanel")
}

########################################################################################
# Read Setaria Circadian Data In (See methods for details on generation of datatables)
# See setaria-circadian.R script for how data was cleaned up / normalized
########################################################################################

#ldhhf.llhcf<-read.csv(file='data/setaria-24-ldhhf.llhcf.txt',sep='\t',header=TRUE, stringsAsFactors = FALSE,strip.white = TRUE)
#ldhhf.llhcf.norm<-read.csv(file='data/setaria-24-ldhhf.llhcf.norm.txt',sep='\t',header=TRUE, stringsAsFactors = FALSE,strip.white = TRUE)
#entrainment<-read.csv(file='data/entrainment.information.txt',sep='\t',header=TRUE, stringsAsFactors = FALSE,strip.white = TRUE)

experiment<-readRDS('data/experiment.information.rds')
quinoa.data<-readRDS('data/quinoa-all-data.rds')

########################################################################################
# Setaria Shiny Application -ui.R
########################################################################################
ui <- dashboardPage(

  dashboardHeader(title="Quinoa Heat Data Explorer"),

  dashboardSidebar(disable = T),
  
  dashboardBody(
    useShinyjs(),
    div(
      id = "loading_page",
      h1("Loading...")
    ),
    hidden(
      div(id="mainTabsetPanel",
        tabsetPanel(
          
          ########################################################################################
          tabPanel(title="Welcome",
               fluidRow(
                 tags$head(includeScript("google-analytics.js")),
                 box(title="Quinoa Heat Data Explorer",width=12, solidHeader = T,status = 'primary',
                     p("This tool is brought to you by the",a(href="http://www.gehan-lab.org/",target='_blank',"Gehan Lab"),
                       "at the Donald Danforth Plant Science Center. For the code used to generate this app,
                       please visit ", a(href="https://github.com/maliagehan/diel-explorer/",target='_blank',"our Github"),". For more information on how the data was processed please refer
                       to", a(href="XXX",target='_blank',"doi: https://doi.org/10.1101/131185"),". To download the raw data go here:"),
                     
                     actionButton("rawdata","Go to Raw Data",icon=icon("th"),
                                  onclick ="window.open('XXX', '_blank')")
                     
                     ),
                 
                 box(title="Using this Tool",width=12, solidHeader = T,status = 'primary',
                      p("In the SAMPLE INFO tab see the available datasets and conditions"),
                      hr(),
                      p("In the SEARCH AND BROWSE DATA tab above, you can search for genes of interest using  by either
                      the search bar or by uploading a .txt file with gene ids. Alternatively, you can use the data filters
                      to browse the data.Once you have loaded or searched for your selections the data can be viewed and dowloaded. 
                      A plot of the data can be seen on the Plot Data tab once the plot button is hit."),
                      hr(),
                      p("In the PLOT DATA tab above, data selected from Select Data or Browse Data tabs is plotted")
                    )
                   )
                  ),
          ########################################################################################
          
          tabPanel(title="Sample Info",
                   fluidRow(
                     box(title="Sample Info",width=12, solidHeader = T,status = 'danger',
                         p("Notes for Reference XXX: Quinoa was grown xxx. Again, for more details sampling conditions, please refer to XXXXX(Pub/Preprint).
                         Differentially expressed genes are between corresponding control (C1 and C11) and experimental samples:"),
                         div(style = 'overflow-x: scroll', dataTableOutput("conditions.table"))
                     )
                   )
          ),
          ########################################################################################
          
          tabPanel(title="Search and Browse Data",
                   fluidRow(
                     box(title="Search Data with GENEID or GO",width=12, solidHeader = T,status = 'success',collapsible = TRUE, collapsed = TRUE,
                         h4("Search using small sets of GENEIDs"),
                         p("GENEIDs, Orthologs, or GO separated by a comma are allowed"),
                         textInput('gene_search_text', label="example: AUR62000001, AUR62000003", value = ""),
                         h4("Search using small sets of GO TERMS"),
                         textInput('go_search_text', label="example:GO:GO:0008270", value = ""),
                         h4("Search using small sets of ORTHOLOG GENEIDs"),
                         textInput('orth_search_text', label="example:AT5G19850.1,AT5G26667.3", value = ""),
                         p("refresh page to clear search")
                     )),
                   
                   fluidRow(
                     box(title="Search Data with File",width=12, solidHeader = T,status = 'success',collapsible = TRUE,collapsed=TRUE,
                         h4("Upload a file of GENEIDs, the GENEIDs should not be quoted, one line per geneid"),
                         fileInput('file.geneid', 'Choose file to upload',
                                   accept = c('text/csv','text/comma-separated-values','text/tab-separated-values','text/plain', '.csv','.tsv')),
                         checkboxInput('header', 'Header', FALSE),
                         h4("Upload a file of ORTHOLOG GENEIDs, the ORTHOLOG GENEIDs should not be quoted, one line per geneid"),
                         fileInput('file.ortholog', 'Choose file to upload',
                                   accept = c('text/csv','text/comma-separated-values','text/tab-separated-values','text/plain', '.csv','.tsv')),
                         checkboxInput('header1', 'Header', FALSE),
                         p("refresh page to clear search")
                     )),
                   fluidRow(
                     box(title="Genes, Orthologs, or GO Selected with Search",width=12, solidHeader = T,status = 'success',collapsible = TRUE,collapsed =TRUE,
                         verbatimTextOutput("selected.genes"),
                         verbatimTextOutput("selected.orthologs"),
                         verbatimTextOutput("selected.go")
                         )),
                   fluidRow(
                     box(title="Browse and Filter Data",width=12, solidHeader = T,status = 'success', collapsible = TRUE,collapsed =TRUE,
                         fluidRow(
                                column(6,
                                       radioButtons('cuttype', 'Cut-Off Type',c('Meets Cut-Off In All Treatments: Day 1'='day1', 
                                                                                'Meets Cut-Off In All Treatments: Day 11'='day11',
                                                                                'Meets Cut-Off All Treatments All Days'='allt',
                                                                                'Meets Cut-Off In At Least One Treatment: Day 1'='day1t',
                                                                                'Meets Cut-Off In At Least One Treatment Day 11'='day11t',
                                                                                'Meets Cut-Off In At Least One Treatment All Days'='onedayt'
                                                                                ),'allt')),
                                column(6,
                                       selectInput("bhq.cutoff1",
                                                   "Benjamini-Hochberg Q-Value Cut-Off:",
                                                   c("All","0.5", "0.05","0.005","0.0005","0.00005","0.000005",
                                                     "0.0000005","0.00000005","0.000000005","0.0000000005",
                                                     "0.00000000005","0.000000000005","0.0000000000005","0.00000000000005")))
                      ))),
                   fluidRow(
                     box(title="Data",width=12, solidHeader = T,status = 'success', collapsible = TRUE,collapsed =FALSE,
                           div(style = 'overflow-x: scroll', dataTableOutput("quinoa.data"))
                     )),
                     fluidRow(
                       box(title="Download Selected Data",width=12, solidHeader = T,status = 'success',
                             column(4,
                                    downloadButton('download.selected', "Download Selected Data"))
                           ))
                     ),
          ########################################################################################
          tabPanel(title="Plot Data",
                   fluidRow(
                     box(title="Plot Data",width=12, solidHeader = T,status = 'info',
                        h4(textOutput("numbergenes")),
                        h4("We restrict plotting to 1000 genes. 
                            Please download data from 'Search and Browse Data' tab or run a local installation of", 
                            a(href="https://github.com/maliagehan/diel-explorer/",target='_blank',"Quinoa Heat Explorer"), 
                            "if you want to graph more genes."),
                         actionButton("plot.heat", label="Plot Selected Data as Heatmap" ),
                         radioButtons('rowcol', 'Scale Heatmap',c(Row='TRUE',Column='FALSE'),'TRUE'),
                         downloadButton('download.heat',"Download Heatmap"),
                         d3heatmapOutput("circadian.expression.heat")
                    )
                   )
                 ),
          ########################################################################################
          tabPanel(title="Contact Us",
                   fluidRow(
                     box(title="Contact US",width=12, solidHeader = T,status = 'danger',
                         p("For questions or if you are interested in adding your own data contact ",
                           a(href="https://github.com/maliagehan/diel-explorer/issues",target='_blank',"Malia Gehan")),
                         p("For mor information on the Gehan lab visit our ",
                           a(href="http://www.gehan-lab.org/",target='_blank',"website"))
                     )
                   )
                )
        ))
    )
  )
)
    
########################################################################################
# Setaria Shiny Application -server.R
########################################################################################
server<-function(input,output,session){
  
  #output for sample info tab ########################################################
    
  output$conditions.table<-renderDataTable(experiment, options=list(paging=FALSE,searching=FALSE))
  
  #output for select data tab #########################################################
  
  #output for browse data tab #########################################################
  
  #get the search terms
  searchgenes<-reactive({
    genes<-gsub(" ","",input$gene_search_text)
    genes1<-data.frame(strsplit(genes,","))
    colnames(genes1)<-c("GENEID")
    genes1
    })
  
  #get the file contents
  searchfile<-reactive({
    geneids<-input$file.geneid
    
    if(is.null(geneids))
      return(NULL)
    
    genes1<-read.csv(geneids$datapath,header=input$header, strip.white = TRUE)
    colnames(genes1)<-c("GENEID")
    genes1
  })
  
  #join the gene search and the file contents
  joinedsearch<-reactive({
    rbind(searchgenes(),searchfile())
  })
  
  #output the list of genes so the user can see it
  output$selected.genes<-renderPrint({
    joinedsearch()
  })
  
  #get the search ortholog terms
  searchorthologs<-reactive({
    orth<-gsub(" ","",input$orth_search_text)
    orth1<-data.frame(strsplit(orth,","))
    colnames(orth1)<-c("ORTHOLOGS")
    orth1
  })
  
  #get the ortholog file terms
  searchfileorth<-reactive({
    orth<-input$file.ortholog
    
    if(is.null(orth))
      return(NULL)
    
    orth1<-read.csv(orth$datapath,header=input$header1,strip.white =TRUE)
    colnames(orth1)<-c("ORTHOLOGS")
    orth1
  })
  
  #join the ortholog search and the ortholog file contents
  joinedsearchorth<-reactive({
    rbind(searchorthologs(),searchfileorth())
  })
  
  #output the list of genes so the user can see it
  output$selected.orthologs<-renderPrint({
    joinedsearchorth()
  })
  
  #get list of go terms
  searchgo<-reactive({
    go<-gsub(" " ,"",input$go_search_text)
    go1<-data.frame(strsplit(go,","))
    colnames(go1)<-c("GO-TERM")
    go1
  })
  
  #output the list of go so the user can see it
  output$selected.go<-renderPrint({
    searchgo()
  })
  
  #filter data
  quinoa.input<-reactive({
    
    quinoa.table<-quinoa.data
    
    if(nrow(joinedsearch())!=0 | nrow(searchgo())!=0 | nrow(joinedsearchorth())!=0){
       colnum<-ncol(quinoa.table)
       quinoa.subset <- data.frame(matrix(ncol=colnum, nrow = 0))
       colnames(quinoa.subset) <- paste0(c(colnames(quinoa.table)))
    
       if(nrow(joinedsearch())!=0){
         for(x in 1:nrow(joinedsearch())){
           search<-as.character(joinedsearch()[x,1])
           row<-quinoa.table[(grep(search,quinoa.table$target_id)),]
           quinoa.subset<-rbind(row,quinoa.subset)
         }
       }
       
       if(nrow(joinedsearchorth())!=0){
         for(x in 1:nrow(joinedsearchorth())){
           search<-as.character(joinedsearchorth()[x,1])
           row<-quinoa.table[(grep(search,quinoa.table$best_hit_arabidopsis)),]
           quinoa.subset<-rbind(row,quinoa.subset)
         }
       }
    
       if(nrow(searchgo())!=0){
         for(x in 1:nrow(searchgo())){
           search<-as.character(searchgo()[x,1])
           row<-quinoa.table[(grep(search,quinoa.table$PANTHER_GO.slim_biological_process)),]
           row1<-quinoa.table[(grep(search,quinoa.table$PANTHER_GO.slim_molecular_function)),]
           row2<-quinoa.table[(grep(search,quinoa.table$PANTHER_GO.slim_cellular_component)),]
           quinoa.subset<-rbind(row,quinoa.subset)
           quinoa.subset<-rbind(row1,quinoa.subset)
           quinoa.subset<-rbind(row2,quinoa.subset)
         }
       }
       
       quinoa.table<-quinoa.subset
    }
    
    if(input$cuttype=="day1"){
      quinoa.table.sub<-subset(quinoa.table,select=-c(hr11_qval,hr11_de,hs11_qval,hs11_de,hrs11_qval,hrs11_de))
      quinoa.table<-quinoa.table.sub
      if (input$bhq.cutoff1!="All"){
        quinoa.table6<-quinoa.table[quinoa.table$hr1_qval<=as.numeric(input$bhq.cutoff1),]
        quinoa.table5<-quinoa.table6[quinoa.table6$hs1_qval<=as.numeric(input$bhq.cutoff1),]
        quinoa.table4<-quinoa.table5[quinoa.table5$hrs1_qval<=as.numeric(input$bhq.cutoff1),]
        quinoa.table3<-na.omit(quinoa.table4)
        quinoa.table<-quinoa.table3
      }
    }
    
    if(input$cuttype=="day11"){
      quinoa.table.sub<-subset(quinoa.table,select=-c(hr1_qval,hr1_de,hs1_qval,hs1_de,hrs1_qval,hrs1_de))
      quinoa.table<-quinoa.table.sub
      if (input$bhq.cutoff1!="All"){
        quinoa.table6<-quinoa.table[quinoa.table$hr11_qval<=as.numeric(input$bhq.cutoff1),]
        quinoa.table5<-quinoa.table6[quinoa.table6$hs11_qval<=as.numeric(input$bhq.cutoff1),]
        quinoa.table4<-quinoa.table5[quinoa.table5$hrs11_qval<=as.numeric(input$bhq.cutoff1),]
        quinoa.table3<-na.omit(quinoa.table4)
        quinoa.table<-quinoa.table3
      }
    }
    
    if(input$cuttype=="allt"){
      if (input$bhq.cutoff1!="All"){
        quinoa.table6<-quinoa.table[quinoa.table$hr1_qval<=as.numeric(input$bhq.cutoff1),]
        quinoa.table5<-quinoa.table6[quinoa.table6$hs1_qval<=as.numeric(input$bhq.cutoff1),]
        quinoa.table4<-quinoa.table5[quinoa.table5$hrs1_qval<=as.numeric(input$bhq.cutoff1),]
        quinoa.table3<-quinoa.table4[quinoa.table4$hr1_qval<=as.numeric(input$bhq.cutoff1),]
        quinoa.table2<-quinoa.table3[quinoa.table3$hs1_qval<=as.numeric(input$bhq.cutoff1),]
        quinoa.table1<-quinoa.table2[quinoa.table2$hrs1_qval<=as.numeric(input$bhq.cutoff1),]
        quinoa.table0<-na.omit(quinoa.table1)
        quinoa.table<-quinoa.table0
      }
    }
    
    if(input$cuttype=="day1t"){
      quinoa.table.sub<-subset(quinoa.table,select=-c(hr11_qval,hr11_de,hs11_qval,hs11_de,hrs11_qval,hrs11_de))
      quinoa.table<-quinoa.table.sub
      if (input$bhq.cutoff1!="All"){
        quinoa.table6<-quinoa.table[(quinoa.table$hr1_qval<=as.numeric(input$bhq.cutoff1) | 
                                       quinoa.table$hs1_qval<=as.numeric(input$bhq.cutoff1) |  
                                       quinoa.table$hrs1_qval<=as.numeric(input$bhq.cutoff1)),]
        quinoa.table5<-na.omit(quinoa.table6)
        quinoa.table<-quinoa.table5
      }
    }
    
    if(input$cuttype=="day11t"){
      quinoa.table.sub<-subset(quinoa.table,select=-c(hr1_qval,hr1_de,hs1_qval,hs1_de,hrs1_qval,hrs1_de))
      quinoa.table<-quinoa.table.sub
      if (input$bhq.cutoff1!="All"){
        quinoa.table6<-quinoa.table[(quinoa.table$hr11_qval<=as.numeric(input$bhq.cutoff1) | 
                                       quinoa.table$hs11_qval<=as.numeric(input$bhq.cutoff1) |  
                                       quinoa.table$hrs11_qval<=as.numeric(input$bhq.cutoff1)),]
        quinoa.table5<-na.omit(quinoa.table6)
        quinoa.table<-quinoa.table5
      }
    }
    
    if(input$cuttype=="onedayt"){
      if (input$bhq.cutoff1!="All"){
        quinoa.table6<-quinoa.table[(  quinoa.table$hr1_qval<=as.numeric(input$bhq.cutoff1) | 
                                       quinoa.table$hs1_qval<=as.numeric(input$bhq.cutoff1) |  
                                       quinoa.table$hrs1_qval<=as.numeric(input$bhq.cutoff1) |
                                       quinoa.table$hr11_qval<=as.numeric(input$bhq.cutoff1) | 
                                       quinoa.table$hs11_qval<=as.numeric(input$bhq.cutoff1) |  
                                       quinoa.table$hrs11_qval<=as.numeric(input$bhq.cutoff1)),]
        quinoa.table5<-na.omit(quinoa.table6)
        quinoa.table<-quinoa.table5
      }
    }
    
    quinoa.table
  })
    
  output$quinoa.data<-renderDataTable(({
    quinoa.input()}),options=list(searching=FALSE, na="NA"))
  
  output$download.selected <- downloadHandler(
    filename = function() { paste('quinoa_heat_tovar_',strftime(Sys.time(),"%m-%d-%y_%H%M%S"),'.csv', sep='') },
    content = function(file) {
      write.csv(quinoa.input(), file)
      }
    )

  #output for plot data tab ###########################################################  
  output$numbergenes<-reactive({
    paste(nrow(quinoa.input()), "Genes are currently selected")
  })
  
  plot.input3<-reactive({
    
    quinoa.graph<-quinoa.input()
    
    if(nrow(quinoa.graph)>1000){
      quinoa.plot1=quinoa.graph[1:1000,]
    }else{quinoa.plot1=quinoa.graph}
    
    quinoa.plot<-subset(quinoa.plot1,select=-c(target_id, best_hit_arabidopsis))
    quinoa.plot<-quinoa.plot[, grepl( "_de" , names(quinoa.plot ) )]
    rownames(quinoa.plot)<-quinoa.plot1$target_id
    print(quinoa.plot)
    str(quinoa.plot)
    return(quinoa.plot)
    })
  
  plot.input2<-reactive({
    setaria.matrix=plot.input3()
    output$circadian.expression.heat<-renderD3heatmap({
      if(input$rowcol=='FALSE'){
        color.palette = colorRampPalette(c("lightyellow","lightblue","blue"),space="rgb")
        d3heatmap(setaria.matrix, scale="column",dendrogram = "none",margins=c(40, 200), color=color.palette(256))
      }else{
        color.palette = colorRampPalette(c("lightyellow","lightblue","blue"),space="rgb")
        d3heatmap(setaria.matrix, scale="row",dendrogram = "none",margins=c(40, 200),color=color.palette(256))
      }})
  })

  observeEvent(input$plot.data,{
    output$circadian.expression.plot<-renderPlot({
      withProgress(message = 'Making plot', value =NULL, {
        incProgress()
        grid.arrange(plot.input1(),ncol=1)
        })
    })})
  
  observeEvent(input$plot.heat,{
    withProgress(message='Making Heatmap', value=NULL,{
    incProgress()
    plot.input2()
    })})
  
  output$download.plot <- downloadHandler(
    filename = function() { paste('quinoa_plot_',strftime(Sys.time(),"%m-%d-%y_%H%M%S"),".svg", sep='') },
    content=function(file){
      ggsave(file, plot = plot.input1(), device = "svg",height =8, width=10)
    }, contentType='image/svg')
  
  output$download.heat <- downloadHandler(
    filename = function() { paste('quinoa_plot_',strftime(Sys.time(),"%m-%d-%y_%H%M%S"),".pdf", sep='') },
    content=function(file){
      color.palette = colorRampPalette(c("lightyellow","lightblue","blue"),space="rgb")
      if(input$rowcol=='FALSE'){scale1="column"}else{scale1="none"}
      pdf(file, width=8, height=7, useDingbats =FALSE)
      heatmap.2(as.matrix(plot.input3()),
                Rowv=FALSE,
                Colv=FALSE,
                dendrogram='none',
                scale=scale1,
                col=color.palette(256),
                trace='none',
                margins=c(8,20),
                symbreaks=FALSE,
                symm=FALSE,
                cex.main=0.75,
                density.info="none"
                )
      dev.off()
    },contentType='image/pdf')
  
  hide(id = "loading_page", anim = TRUE, animType = "fade")    
  show("mainTabsetPanel")
}

########################################################################################
#  Setaria Shiny Application - run app.
########################################################################################

# Run the application 
shinyApp(ui = ui, server = server)

########################################################################################
