# MOVICS RShiny UI
# Created by Junkai Zhu on 2022/07/06

library(shiny)
library(shinyBS)
library(shinyjs)
library(shinydashboard)
library(rhandsontable)
library(DT)
library(ggplot2)
library(jsonlite)
library(TCGAbiolinks)
library(data.table)
library(ChAMPdata)
library(tidyverse)
library(magrittr)
library(readxl)
library(stringr)
library(forcats)
library(iClusterPlus)
library(data.table)
library(RColorBrewer)
library(ComplexHeatmap)
library(circlize)
library(survival)
library(survminer)
library(ggpp)
library(MOVICS)
library(patchwork)
library(ggalluvial)
library(cowplot)
library(GSVA)
library(HelpersMG)
library(CMScaller)
library(pamr)
library(clusterRepro)
library(limma)
#library(filesstrings) #安装的话使用github安装install_github('rorynolan/filesstrings')

options(timeout = 360000,shiny.maxRequestSize=30000*1024^2) #设置文件上传大小限制为30GB且响应时间为100h
# data.table::setDTthreads(threads = 1) #设置程序运行线程为1
#Sys.setlocale(category = "LC_ALL", locale = "zh_CN.utf8") #allow chinese
Sys.setlocale(category = "LC_ALL",locale = "chinese")
AppVersion="Version: V1.0"
UpdateTime="Last Updated: 11/10/2023"

ui <- fluidPage(
    useShinyjs(),
    withMathJax(), #math symbols
    
    #Title Section
    fixedRow(
        column(12,
               fixedRow(
                   div(style="color:#00008B",
                       align = "center",
                       headerPanel("MOVICS: Multi-Omics integration and VIsualization in Cancer Subtyping"
                       )
                   ),
                   column(12,tags$p(tags$h4(tags$em(tags$strong("State Key Laboratory of Natural Medicines, Research Center of Biostatistics and Computational Pharmacy, China Pharmaceutical University, Nanjing, China")),align="center",style = "color:#317EAC"))),
                   column(12,tags$em(tags$p(tags$h5(AppVersion,";",UpdateTime,align = "Center", style = "color:grey"))))
               )
        )
    ),
    
    tags$style(HTML("
        .tabs-above > .tabset > li[class=active] > a {
                    background-color: #000;
                    color: #FFF;
                    }")),
    
    #Main Page with parameters input box on left and results output box on right
    tabsetPanel(type="tabs",id="panels",
                
                #Data Preparation Section
                tabPanel(
                    title=tags$p(tags$h4(tags$strong("Data Preparation"),style="color:purple")),
                    value="Data Preparation",
                    
                    column(5,
                           
                           #1.Basic Settings
                           wellPanel(style="background-color:#47AEE9;",
                                     
                                     fixedRow(
                                         column(12,tags$h4(strong("Basic Settings"),align="center",style="color:#834005")),
                                         # column(12,textAreaInput(inputId = "cancerType",
                                         #                         label = "Cancer Types From TCGA (Splitted by english semicolons)",
                                         #                         value = "",
                                         #                         width = '100%', height = '35px')
                                         # ),
                                         # column(12,radioButtons("continue_analysis", h5(strong("Continue analysis on previous results or not")), inline=TRUE,
                                         #                        choices = c("Yes"=1,"No"=2),selected=2)),
                                         # conditionalPanel(condition="input.continue_analysis==1",column(12,textInput(inputId="continue_id",label="Please input the previously used user account",value=""))),
                                         column(12,textAreaInput(inputId = "user_name",
                                                                 label = "Username (E-mail Address)",
                                                                 value = "",
                                                                 width = '100%', height = '35px')
                                         ),
                                         column(12,textInput(inputId="user_account",label="User Account",value="")),
                                         column(12,radioButtons("log_in_first_time", h5(strong("Whether log in for the first time using this user account")), inline=TRUE,
                                                                choices = c("Yes"=1,"No"=2),selected=1)),
                                         column(12,numericInput(inputId="cancerType_number",label="The Number of Cancer Types",value=1,min=1,step=1)),
                                         column(12,offset=0,tags$p(h5(strong("Input Cancer Types From TCGA",align="left")))),
                                         column(9,rHandsontableOutput("table.cancerType")),
                                         column(12,offset=0,tags$p(h5(strong("",align="left")))),  #空一行
                                         column(12,textAreaInput(inputId = "projectName",
                                                                 label = "Project Name",
                                                                 value = "",
                                                                 width = '100%', height = '35px')
                                         ),
                                         column(12,offset=0,tags$p(h5("Now let's start to prepare the tcga datasets first.",align="left"))),
                                         column(12,checkboxGroupInput(inputId = "multiOmics",
                                                                      label = "Multi-omics Types (Choose at least two types)",
                                                                      choices = list("mRNA"=1,"lncRNA"=2,"DNA methylation"=3,"copy number alterations"=4,"binary somatic mutation"=5,"radiomics"=6),
                                                                      selected = 1:5
                                                                      )),
                                         column(12,radioButtons("data_sources", h5(strong("Data Sources")), inline=TRUE,
                                                                choices = c("Internal resources"=1,"External resources"=2),selected=2)),
                                         conditionalPanel(condition="input.data_sources==2",
                                                          column(12,radioButtons("use_integrated_datasets", h5(strong("Use integrated datasets (.RData) or not")), inline=TRUE,
                                                                                 choices = c("Yes"=1,"No"=2),selected=1)),
                                                          conditionalPanel(condition="input.use_integrated_datasets==1",
                                                                           column(12,radioButtons("use_integrated_datasets_approaches", h5(strong("The approach to obtain integrated datasets")), inline=TRUE,
                                                                                                  choices = c("Download from specified urls"=1,"Upload manually"=2),selected=2)),
                                                                           conditionalPanel(condition="input.use_integrated_datasets_approaches==1",
                                                                                            
                                                                                            column(12,textAreaInput(inputId = "integrated_datasets_downloadUrl",
                                                                                                                    label = "Download url for the integrated datasets (.RData)",
                                                                                                                    value = "",
                                                                                                                    width = '100%', height = '35px'))
                                                                           ),
                                                                           conditionalPanel(condition="input.use_integrated_datasets_approaches==2",
                                                                                            
                                                                                            column(12,fileInput('integrated_datasets_fileUpload',tags$p(tags$b("Upload the integrated datasets (.RData)")),
                                                                                                                accept=c('.RData') ))
                                                                           )
                                                                          
                                                          )
                                         ),
                                         column(12,radioButtons("external_validation", h5(strong("External validation or not")), inline=TRUE,
                                                                choices = c("Yes"=1,"No"=2),selected=1)),
                                         conditionalPanel(condition="input.data_sources==1",
                                                          column(6,actionButton("basic_create1","Create",width="100%",class="btn btn-primary")),
                                                          column(6,actionButton("basic_process","Process",width="100%",class="btn btn-primary"))
                                                          ),
                                         conditionalPanel(condition="input.data_sources==2 & input.use_integrated_datasets==1",
                                                          column(6,actionButton("basic_create2","Create",width="100%",class="btn btn-primary")),
                                                          column(6,actionButton("use_integrated_datasets_process","Process",width="100%",class="btn btn-primary"))
                                         ),
                                         conditionalPanel(condition="input.data_sources==2 & input.use_integrated_datasets==2",
                                                          column(12,actionButton("basic_create3","Create",width="100%",class="btn btn-primary"))
                                         )
                                     )
                                     ),
                           
                           conditionalPanel(condition="input.data_sources==1",
                                            
                                            #radiomics
                                            conditionalPanel(condition="input.multiOmics.indexOf('6') > -1",
                                                             wellPanel(style="background-color:#47AEE9;",
                                                                       fixedRow(
                                                                           column(12,tags$h4(strong("radiomics"),align="center",style="color:#834005")),
                                                                           column(12,radioButtons("internal_radiomics_approaches", h5(strong("Approaches")), inline=TRUE,
                                                                                                  choices = c("Download from specified urls"=1,"Upload manually"=2),selected=2)),
                                                                           
                                                                           conditionalPanel(condition="input.internal_radiomics_approaches==1",
                                                                                            column(12,radioButtons("internal_radiomics_url_sources", h5(strong("Url sources")), inline=TRUE,
                                                                                                                   choices = c("Combined dataset"=1,"Separated datasets"=2),selected=1)),
                                                                                            conditionalPanel(condition="input.internal_radiomics_url_sources==1",
                                                                                                             column(12,textAreaInput(inputId = "internal_radiomics_downloadUrl",
                                                                                                                                     label = "Download url for combined radiomics dataset",
                                                                                                                                     value = "",
                                                                                                                                     width = '100%', height = '35px'))
                                                                                            ),
                                                                                            conditionalPanel(condition="input.internal_radiomics_url_sources==2",
                                                                                                             column(12,offset=0,tags$p(h5(strong("Input radiomics dataset url for each cancer type",align="left")))),
                                                                                                             column(9,rHandsontableOutput("table.internalRadiomicsUrl")),
                                                                                                             column(12,offset=0,tags$p(h5(strong("",align="left"))))  #空一行
                                                                                            ),
                                                                                            column(12,actionButton("internal_radiomics_process1","Process",width="100%",class="btn btn-primary"))
                                                                           ),
                                                                           conditionalPanel(condition="input.internal_radiomics_approaches==2",
                                                                                            column(12,radioButtons("internal_radiomics_upload_sources", h5(strong("Upload sources")), inline=TRUE,
                                                                                                                   choices = c("Combined dataset"=1,"Separated datasets"=2),selected=1)),
                                                                                            conditionalPanel(condition="input.internal_radiomics_upload_sources==1",
                                                                                                             column(12,fileInput('internal_radiomics_fileUpload1',tags$p(tags$b("Upload the combined radiomics dataset (a .txt file)")),
                                                                                                                                 accept=c('text/csv', 'text/comma-separated-values,text/plain', '.csv') ))
                                                                                            ),
                                                                                            conditionalPanel(condition="input.internal_radiomics_upload_sources==2",
                                                                                                             column(12,fileInput('internal_radiomics_fileUpload2',tags$p(tags$b("Upload the radiomics dataset sequentially for each cancer type (a .txt file)")),
                                                                                                                                 accept=c('text/csv', 'text/comma-separated-values,text/plain', '.csv'), multiple=TRUE))
                                                                                            ),
                                                                                            column(12,actionButton("internal_radiomics_process2","Process",width="100%",class="btn btn-primary"))
                                                                           )
                                                                       )
                                                             )        
                                            )      
                           ),
                           
                           conditionalPanel(condition="input.data_sources==2 & input.use_integrated_datasets==2",
                           #2.Clinical & Survival
                           wellPanel(style="background-color:#47AEE9;",
                                     fixedRow(
                                         column(12,tags$h4(strong("Clinical & Survival"),align="center",style="color:#834005")),
                                         column(12,radioButtons("cli_sur_approaches", h5(strong("Approaches")), inline=TRUE,
                                                                choices = c("Download automatically"=1,"Download from specified urls"=2,"Upload manually"=3),selected=1)),
                                         conditionalPanel(condition="input.cli_sur_approaches==1",
                                                          column(12,actionButton("cli_sur_process1","Process",width="100%",class="btn btn-primary"))
                                                          ),
                                         conditionalPanel(condition="input.cli_sur_approaches==2",
                                                          column(12,radioButtons("cli_sur_url_sources", h5(strong("Url sources")), inline=TRUE,
                                                                                 choices = c("Combined dataset"=1,"Separated datasets"=2),selected=1)),
                                                          conditionalPanel(condition="input.cli_sur_url_sources==1",
                                                                           column(12,textAreaInput(inputId = "cli_sur_downloadUrl",
                                                                                                   label = "Download url for combined clinical dataset",
                                                                                                   value = "",
                                                                                                   width = '100%', height = '35px'))
                                                                           ),
                                                          conditionalPanel(condition="input.cli_sur_url_sources==2",
                                                                           column(12,offset=0,tags$p(h5(strong("Input clinical dataset url for each cancer type",align="left")))),
                                                                           column(9,rHandsontableOutput("table.cliSurUrl")),
                                                                           column(12,offset=0,tags$p(h5(strong("",align="left"))))  #空一行
                                                                           ),
                                                          column(12,actionButton("cli_sur_process2","Process",width="100%",class="btn btn-primary"))
                                                          ),
                                         conditionalPanel(condition="input.cli_sur_approaches==3",
                                                          column(12,radioButtons("cli_sur_upload_sources", h5(strong("Upload sources")), inline=TRUE,
                                                                                 choices = c("Combined dataset"=1,"Separated datasets"=2),selected=1)),
                                                          conditionalPanel(condition="input.cli_sur_upload_sources==1",
                                                                           column(12,fileInput('cli_sur_fileUpload1',tags$p(tags$b("Upload the combined clinical dataset (a .txt file)")),
                                                                                               accept=c('text/csv', 'text/comma-separated-values,text/plain', '.csv') ))
                                                          ),
                                                          conditionalPanel(condition="input.cli_sur_upload_sources==2",
                                                                           column(12,fileInput('cli_sur_fileUpload2',tags$p(tags$b("Upload the clinical dataset sequentially for each cancer type at once (a .txt file)")),
                                                                                               accept=c('text/csv', 'text/comma-separated-values,text/plain', '.csv'), multiple=TRUE))
                                                          ),
                                                          column(12,actionButton("cli_sur_process3","Process",width="100%",class="btn btn-primary"))
                                        )
                                     )
                           ),
                          
                           #3.mRNA & lncRNA
                           conditionalPanel(condition="input.multiOmics.includes('1') | input.multiOmics.includes('2')",
                               wellPanel(style="background-color:#47AEE9;",
                                         fixedRow(
                                             conditionalPanel(condition="input.multiOmics.includes('1') & input.multiOmics.includes('2')",column(12,tags$h4(strong("mRNA & lncRNA"),align="center",style="color:#834005"))),
                                             conditionalPanel(condition="input.multiOmics.includes('1') & !input.multiOmics.includes('2')",column(12,tags$h4(strong("mRNA"),align="center",style="color:#834005"))),
                                             conditionalPanel(condition="!input.multiOmics.includes('1') & input.multiOmics.includes('2')",column(12,tags$h4(strong("lncRNA"),align="center",style="color:#834005"))),
                                             
                                             column(12,radioButtons("RNA_approaches", h5(strong("Approaches")), inline=TRUE,
                                                                    choices = c("Download automatically"=1,"Download from specified urls"=2,"Upload manually"=3),selected=1)),
                                             conditionalPanel(condition="input.RNA_approaches==1",
                                                              column(12,actionButton("RNA_process1","Process",width="100%",class="btn btn-primary"))
                                             ),
                                             conditionalPanel(condition="input.RNA_approaches==2",
                                                              column(12,radioButtons("RNA_url_sources", h5(strong("Url sources")), inline=TRUE,
                                                                                     choices = c("Combined dataset"=1,"Separated datasets"=2),selected=1)),
                                                              conditionalPanel(condition="input.RNA_url_sources==1",
                                                                               column(12,textAreaInput(inputId = "RNA_downloadUrl",
                                                                                                       label = "Download url for combined RNA dataset",
                                                                                                       value = "",
                                                                                                       width = '100%', height = '35px'))
                                                              ),
                                                              conditionalPanel(condition="input.RNA_url_sources==2",
                                                                               column(12,offset=0,tags$p(h5(strong("Input RNA dataset url for each cancer type",align="left")))),
                                                                               column(9,rHandsontableOutput("table.RNAUrl")),
                                                                               column(12,offset=0,tags$p(h5(strong("",align="left"))))  #空一行
                                                              ),
                                                              column(12,actionButton("RNA_process2","Process",width="100%",class="btn btn-primary"))
                                             ),
                                             conditionalPanel(condition="input.RNA_approaches==3",
                                                              column(12,radioButtons("RNA_upload_sources", h5(strong("Upload sources")), inline=TRUE,
                                                                                     choices = c("Combined dataset"=1,"Separated datasets"=2),selected=1)),
                                                              conditionalPanel(condition="input.RNA_upload_sources==1",
                                                                               column(12,fileInput('RNA_fileUpload1',tags$p(tags$b("Upload the combined RNA dataset (a .txt file)")),
                                                                                                   accept=c('text/csv', 'text/comma-separated-values,text/plain', '.csv') ))
                                                              ),
                                                              conditionalPanel(condition="input.RNA_upload_sources==2",
                                                                               column(12,fileInput('RNA_fileUpload2',tags$p(tags$b("Upload the RNA dataset sequentially for each cancer type (a .txt file)")),
                                                                                                   accept=c('text/csv', 'text/comma-separated-values,text/plain', '.csv'), multiple=TRUE))
                                                              ),
                                                              column(12,actionButton("RNA_process3","Process",width="100%",class="btn btn-primary"))
                                             )
                                         )
                               )
                           ),
                           
                           #4.DNA methylation
                           conditionalPanel(condition="input.multiOmics.includes('3')",
                               wellPanel(style="background-color:#47AEE9;",
                                         fixedRow(
                                             column(12,tags$h4(strong("DNA methylation"),align="center",style="color:#834005")),
                                             column(12,radioButtons("DNA_methylation_approaches", h5(strong("Approaches")), inline=TRUE,
                                                                    choices = c("Download automatically"=1,"Download from specified urls"=2,"Upload manually"=3),selected=1)),
                                             conditionalPanel(condition="input.DNA_methylation_approaches==1",
                                                              column(12,actionButton("DNA_methylation_process1","Process",width="100%",class="btn btn-primary"))
                                             ),
                                             conditionalPanel(condition="input.DNA_methylation_approaches==2",
                                                              column(12,radioButtons("DNA_methylation_url_sources", h5(strong("Url sources")), inline=TRUE,
                                                                                     choices = c("Combined dataset"=1,"Separated datasets"=2),selected=1)),
                                                              conditionalPanel(condition="input.DNA_methylation_url_sources==1",
                                                                               column(12,textAreaInput(inputId = "DNA_methylation_downloadUrl",
                                                                                                       label = "Download url for combined DNA methylation dataset",
                                                                                                       value = "",
                                                                                                       width = '100%', height = '35px'))
                                                              ),
                                                              conditionalPanel(condition="input.DNA_methylation_url_sources==2",
                                                                               column(12,offset=0,tags$p(h5(strong("Input DNA methylation dataset url for each cancer type",align="left")))),
                                                                               column(9,rHandsontableOutput("table.DNAMethylationUrl")),
                                                                               column(12,offset=0,tags$p(h5(strong("",align="left"))))  #空一行
                                                              ),
                                                              column(12,actionButton("DNA_methylation_process2","Process",width="100%",class="btn btn-primary"))
                                             ),
                                             conditionalPanel(condition="input.DNA_methylation_approaches==3",
                                                              column(12,radioButtons("DNA_methylation_upload_sources", h5(strong("Upload sources")), inline=TRUE,
                                                                                     choices = c("Combined dataset"=1,"Separated datasets"=2),selected=1)),
                                                              conditionalPanel(condition="input.DNA_methylation_upload_sources==1",
                                                                               column(12,fileInput('DNA_methylation_fileUpload1',tags$p(tags$b("Upload the combined DNA methylation dataset (a .gz file)")),
                                                                                                   accept=".gz" ))
                                                              ),
                                                              conditionalPanel(condition="input.DNA_methylation_upload_sources==2",
                                                                               column(12,fileInput('DNA_methylation_fileUpload2',tags$p(tags$b("Upload the DNA methylation dataset sequentially for each cancer type (a .gz file)")),
                                                                                                   accept=".gz", multiple=TRUE))
                                                              ),
                                                              column(12,actionButton("DNA_methylation_process3","Process",width="100%",class="btn btn-primary"))

                                             )
                                         )
                                )
                           ),
                           
                           #5.copy number alterations
                           conditionalPanel(condition="input.multiOmics.indexOf('4') > -1",
                               wellPanel(style="background-color:#47AEE9;",
                                         fixedRow(
                                             column(12,tags$h4(strong("copy number alterations"),align="center",style="color:#834005")),
                                             column(12,radioButtons("copy_number_alterations_approaches", h5(strong("Approaches")), inline=TRUE,
                                                                    choices = c("Download automatically"=1,"Download from specified urls"=2,"Upload manually"=3),selected=1)),
                                             conditionalPanel(condition="input.copy_number_alterations_approaches==1",
                                                              column(12,actionButton("copy_number_alterations_process1","Process",width="100%",class="btn btn-primary"))
                                             ),
                                             conditionalPanel(condition="input.copy_number_alterations_approaches==2",
                                                              column(12,radioButtons("copy_number_alterations_url_sources", h5(strong("Url sources")), inline=TRUE,
                                                                                     choices = c("Combined dataset"=1,"Separated datasets"=2),selected=1)),
                                                              conditionalPanel(condition="input.copy_number_alterations_url_sources==1",
                                                                               column(12,textAreaInput(inputId = "copy_number_alterations_downloadUrl",
                                                                                                       label = "Download url for combined copy number alterations dataset",
                                                                                                       value = "",
                                                                                                       width = '100%', height = '35px'))
                                                              ),
                                                              conditionalPanel(condition="input.copy_number_alterations_url_sources==2",
                                                                               column(12,offset=0,tags$p(h5(strong("Input copy number alterations dataset url for each cancer type",align="left")))),
                                                                               column(9,rHandsontableOutput("table.copyNumberAlterationsUrl")),
                                                                               column(12,offset=0,tags$p(h5(strong("",align="left"))))  #空一行
                                                              ),
                                                              column(12,actionButton("copy_number_alterations_process2","Process",width="100%",class="btn btn-primary"))
                                             ),
                                             conditionalPanel(condition="input.copy_number_alterations_approaches==3",
                                                              column(12,radioButtons("copy_number_alterations_upload_sources", h5(strong("Upload sources")), inline=TRUE,
                                                                                     choices = c("Combined dataset"=1,"Separated datasets"=2),selected=1)),
                                                              conditionalPanel(condition="input.copy_number_alterations_upload_sources==1",
                                                                               column(12,fileInput('copy_number_alterations_fileUpload1',tags$p(tags$b("Upload the combined copy number alterations dataset (a .tar.gz file)")),
                                                                                                   accept=".tar.gz" ))
                                                              ),
                                                              conditionalPanel(condition="input.copy_number_alterations_upload_sources==2",
                                                                               column(12,fileInput('copy_number_alterations_fileUpload2',tags$p(tags$b("Upload the copy number alterations dataset sequentially for each cancer type (a .tar.gz file)")),
                                                                                                   accept=".tar.gz", multiple=TRUE))
                                                              ),
                                                              column(12,actionButton("copy_number_alterations_process3","Process",width="100%",class="btn btn-primary"))
                                             )
                                         )
                                )        
                           ),
                           
                           #6.binary somatic mutation
                           conditionalPanel(condition="input.multiOmics.indexOf('5') > -1",
                               wellPanel(style="background-color:#47AEE9;",
                                         fixedRow(
                                             column(12,tags$h4(strong("binary somatic mutation"),align="center",style="color:#834005")),
                                             column(12,radioButtons("binary_somatic_mutation_approaches", h5(strong("Approaches")), inline=TRUE,
                                                                    choices = c("Download automatically"=1,"Download from specified urls"=2,"Upload manually"=3),selected=1)),
                                             conditionalPanel(condition="input.binary_somatic_mutation_approaches==1",
                                                              column(12,actionButton("binary_somatic_mutation_process1","Process",width="100%",class="btn btn-primary"))
                                             ),
                                             conditionalPanel(condition="input.binary_somatic_mutation_approaches==2",
                                                              column(12,radioButtons("binary_somatic_mutation_url_sources", h5(strong("Url sources")), inline=TRUE,
                                                                                     choices = c("Combined dataset"=1,"Separated datasets"=2),selected=1)),
                                                              conditionalPanel(condition="input.binary_somatic_mutation_url_sources==1",
                                                                               column(12,textAreaInput(inputId = "binary_somatic_mutation_downloadUrl",
                                                                                                       label = "Download url for combined binary somatic mutation dataset",
                                                                                                       value = "",
                                                                                                       width = '100%', height = '35px'))
                                                              ),
                                                              conditionalPanel(condition="input.binary_somatic_mutation_url_sources==2",
                                                                               column(12,offset=0,tags$p(h5(strong("Input binary somatic mutation dataset url for each cancer type",align="left")))),
                                                                               column(9,rHandsontableOutput("table.binarySomaticMutationUrl")),
                                                                               column(12,offset=0,tags$p(h5(strong("",align="left"))))  #空一行
                                                              ),
                                                              column(12,actionButton("binary_somatic_mutation_process2","Process",width="100%",class="btn btn-primary"))
                                             ),
                                             conditionalPanel(condition="input.binary_somatic_mutation_approaches==3",
                                                              column(12,radioButtons("binary_somatic_mutation_upload_sources", h5(strong("Upload sources")), inline=TRUE,
                                                                                     choices = c("Combined dataset"=1,"Separated datasets"=2),selected=1)),
                                                              conditionalPanel(condition="input.binary_somatic_mutation_upload_sources==1",
                                                                               column(12,fileInput('binary_somatic_mutation_fileUpload1',tags$p(tags$b("Upload the combined binary somatic mutation dataset (a .tar.gz file)")),
                                                                                                   accept=".tar.gz" ))
                                                              ),
                                                              conditionalPanel(condition="input.binary_somatic_mutation_upload_sources==2",
                                                                               column(12,fileInput('binary_somatic_mutation_fileUpload2',tags$p(tags$b("Upload the binary somatic mutation dataset sequentially for each cancer type (a .tar.gz file)")),
                                                                                                   accept=".tar.gz", multiple=TRUE))
                                                              ),
                                                              column(12,actionButton("binary_somatic_mutation_process3","Process",width="100%",class="btn btn-primary"))
                                             )
                                         )
                               )
                           ),
                           
                           #7.radiomics
                           conditionalPanel(condition="input.multiOmics.indexOf('6') > -1",
                                            wellPanel(style="background-color:#47AEE9;",
                                                      fixedRow(
                                                          column(12,tags$h4(strong("radiomics"),align="center",style="color:#834005")),
                                                          column(12,radioButtons("external_radiomics_approaches", h5(strong("Approaches")), inline=TRUE,
                                                                                 choices = c("Download from specified urls"=1,"Upload manually"=2),selected=2)),
                                                          
                                                          conditionalPanel(condition="input.external_radiomics_approaches==1",
                                                                           column(12,radioButtons("external_radiomics_url_sources", h5(strong("Url sources")), inline=TRUE,
                                                                                                  choices = c("Combined dataset"=1,"Separated datasets"=2),selected=1)),
                                                                           conditionalPanel(condition="input.external_radiomics_url_sources==1",
                                                                                            column(12,textAreaInput(inputId = "external_radiomics_downloadUrl",
                                                                                                                    label = "Download url for combined radiomics dataset",
                                                                                                                    value = "",
                                                                                                                    width = '100%', height = '35px'))
                                                                           ),
                                                                           conditionalPanel(condition="input.external_radiomics_url_sources==2",
                                                                                            column(12,offset=0,tags$p(h5(strong("Input radiomics dataset url for each cancer type",align="left")))),
                                                                                            column(9,rHandsontableOutput("table.externalRadiomicsUrl")),
                                                                                            column(12,offset=0,tags$p(h5(strong("",align="left"))))  #空一行
                                                                           ),
                                                                           column(12,actionButton("external_radiomics_process1","Process",width="100%",class="btn btn-primary"))
                                                          ),
                                                          conditionalPanel(condition="input.external_radiomics_approaches==2",
                                                                           column(12,radioButtons("external_radiomics_upload_sources", h5(strong("Upload sources")), inline=TRUE,
                                                                                                  choices = c("Combined dataset"=1,"Separated datasets"=2),selected=1)),
                                                                           conditionalPanel(condition="input.external_radiomics_upload_sources==1",
                                                                                            column(12,fileInput('external_radiomics_fileUpload1',tags$p(tags$b("Upload the combined radiomics dataset (a .txt file)")),
                                                                                                                accept=c('text/csv', 'text/comma-separated-values,text/plain', '.csv') ))
                                                                           ),
                                                                           conditionalPanel(condition="input.external_radiomics_upload_sources==2",
                                                                                            column(12,fileInput('external_radiomics_fileUpload2',tags$p(tags$b("Upload the radiomics dataset sequentially for each cancer type (a .txt file)")),
                                                                                                                accept=c('text/csv', 'text/comma-separated-values,text/plain', '.csv'), multiple=TRUE))
                                                                           ),
                                                                           column(12,actionButton("external_radiomics_process2","Process",width="100%",class="btn btn-primary"))
                                                          )
                                                      )
                                            )        
                           )
                        ),
                        
                        #8.Datasets Integration
                        conditionalPanel(condition="input.data_sources==1 | (input.data_sources==2 & input.use_integrated_datasets==2)",
                                         
                                         wellPanel(style="background-color:#47AEE9;",
                                                   
                                                   fixedRow(
                                                       column(12,tags$h4(strong("TCGA Datasets Integration"),align="center",style="color:#834005")),
                                                       column(12,offset=0,tags$p(h5("Now all omics datasets you specified have been well prepared, then click the 'Integrate' button below to integrate all omics datasets.",align="left"))),
                                                       column(12,actionButton("datasets_integrate","Integrate",width="100%",class="btn btn-primary"))
                                                   )
                                         )          
                        ),
                        
                        #9.Validation datasets preparation
                        conditionalPanel(condition="input.external_validation==1",
                             
                                         wellPanel(style="background-color:#47AEE9;",
                                                   
                                                   fixedRow(
                                                       column(12,tags$h4(strong("Validation Datasets Preparation"),align="center",style="color:#834005")),
                                                       column(12,offset=0,tags$p(h5("Now let's start to prepare the validation datasets which will be used in 'RUN Module' and 'COMP Module'.",align="left"))),
                                                       column(12,checkboxGroupInput(inputId = "validation_multiOmics",
                                                                                    label = "Multi-omics Types",
                                                                                    choices = list("mRNA"=1,"lncRNA"=2,"DNA methylation"=3,"copy number alterations"=4,"binary somatic mutation"=5,"radiomics"=6),
                                                                                    selected = 1:5
                                                       )),
                                                       column(12,radioButtons("validation_use_integrated_datasets", h5(strong("Use integrated datasets (.RData) or not")), inline=TRUE,
                                                                              choices = c("Yes"=1,"No"=2),selected=1)),
                                                       conditionalPanel(condition="input.validation_use_integrated_datasets==1",
                                                                        column(12,radioButtons("validation_use_integrated_datasets_approaches", h5(strong("The approach to obtain integrated datasets")), inline=TRUE,
                                                                                               choices = c("Download from specified urls"=1,"Upload manually"=2),selected=2)),
                                                                        conditionalPanel(condition="input.validation_use_integrated_datasets_approaches==1",
                                                                                         
                                                                                         column(12,textAreaInput(inputId = "validation_integrated_datasets_downloadUrl",
                                                                                                                 label = "Download url for the integrated datasets (.RData)",
                                                                                                                 value = "",
                                                                                                                 width = '100%', height = '35px'))
                                                                        ),
                                                                        conditionalPanel(condition="input.validation_use_integrated_datasets_approaches==2",
                                                                                         
                                                                                         column(12,fileInput('validation_integrated_datasets_fileUpload',tags$p(tags$b("Upload the integrated datasets (.RData)")),
                                                                                                             accept=c('.RData') ))
                                                                        )
                                                       ),
                                                       conditionalPanel(condition="input.validation_use_integrated_datasets==1",
                                                                        column(6,actionButton("validation_basic_create","Create",width="100%",class="btn btn-primary")),
                                                                        column(6,actionButton("validation_use_integrated_datasets_process","Process",width="100%",class="btn btn-primary"))
                                                       ),
                                                       conditionalPanel(condition="input.validation_use_integrated_datasets==2",
                                                                        column(12,actionButton("validation_basic_create","Create",width="100%",class="btn btn-primary"))
                                                       )
                                                   )
                                            ),     
                                            
                                         conditionalPanel(condition="input.validation_use_integrated_datasets==2",
                                                          
                                                          #Clinical & Survival(Validation)
                                                          wellPanel(style="background-color:#47AEE9;",
                                                                    fixedRow(
                                                                        column(12,tags$h4(strong("Clinical & Survival(Validation)"),align="center",style="color:#834005")),
                                                                        column(12,radioButtons("validation_cli_sur_approaches", h5(strong("Approaches")), inline=TRUE,
                                                                                               choices = c("Download from specified urls"=1,"Upload manually"=2),selected=2)),
                                                                        conditionalPanel(condition="input.validation_cli_sur_approaches==1",
                                                                                         column(12,radioButtons("validation_cli_sur_url_sources", h5(strong("Url sources")), inline=TRUE,
                                                                                                                choices = c("Combined dataset"=1,"Separated datasets"=2),selected=1)),
                                                                                         conditionalPanel(condition="input.validation_cli_sur_url_sources==1",
                                                                                                          column(12,textAreaInput(inputId = "validation_cli_sur_downloadUrl",
                                                                                                                                  label = "Download url for combined clinical dataset",
                                                                                                                                  value = "",
                                                                                                                                  width = '100%', height = '35px'))
                                                                                         ),
                                                                                         conditionalPanel(condition="input.validation_cli_sur_url_sources==2",
                                                                                                          column(12,offset=0,tags$p(h5(strong("Input clinical dataset url for each cancer type",align="left")))),
                                                                                                          column(9,rHandsontableOutput("table.validation_cliSurUrl")),
                                                                                                          column(12,offset=0,tags$p(h5(strong("",align="left"))))  #空一行
                                                                                         ),
                                                                                         column(12,actionButton("validation_cli_sur_process1","Process",width="100%",class="btn btn-primary"))
                                                                        ),
                                                                        conditionalPanel(condition="input.validation_cli_sur_approaches==2",
                                                                                         column(12,radioButtons("validation_cli_sur_upload_sources", h5(strong("Upload sources")), inline=TRUE,
                                                                                                                choices = c("Combined dataset"=1,"Separated datasets"=2),selected=1)),
                                                                                         conditionalPanel(condition="input.validation_cli_sur_upload_sources==1",
                                                                                                          column(12,fileInput('validation_cli_sur_fileUpload1',tags$p(tags$b("Upload the downloaded combined clinical dataset (a .txt file)")),
                                                                                                                              accept=c('text/csv', 'text/comma-separated-values,text/plain', '.csv') ))
                                                                                         ),
                                                                                         conditionalPanel(condition="input.validation_cli_sur_upload_sources==2",
                                                                                                          column(12,fileInput('validation_cli_sur_fileUpload2',tags$p(tags$b("Upload the downloaded clinical dataset sequentially for each cancer type at once (a .txt file)")),
                                                                                                                              accept=c('text/csv', 'text/comma-separated-values,text/plain', '.csv'), multiple=TRUE))
                                                                                         ),
                                                                                         column(12,actionButton("validation_cli_sur_process2","Process",width="100%",class="btn btn-primary"))
                                                                        )
                                                                    )
                                                          ),
                                                          
                                                          #mRNA & lncRNA(Validation)
                                                          conditionalPanel(condition="input.validation_multiOmics.includes('1') | input.validation_multiOmics.includes('2')",
                                                                           wellPanel(style="background-color:#47AEE9;",
                                                                                     fixedRow(
                                                                                         conditionalPanel(condition="input.validation_multiOmics.includes('1') & input.validation_multiOmics.includes('2')",column(12,tags$h4(strong("mRNA & lncRNA(Validation)"),align="center",style="color:#834005"))),
                                                                                         conditionalPanel(condition="input.validation_multiOmics.includes('1') & !input.validation_multiOmics.includes('2')",column(12,tags$h4(strong("mRNA(Validation)"),align="center",style="color:#834005"))),
                                                                                         conditionalPanel(condition="!input.validation_multiOmics.includes('1') & input.validation_multiOmics.includes('2')",column(12,tags$h4(strong("lncRNA(Validation)"),align="center",style="color:#834005"))),
                                                                                         
                                                                                         column(12,radioButtons("validation_RNA_approaches", h5(strong("Approaches")), inline=TRUE,
                                                                                                                choices = c("Download from specified urls"=1,"Upload manually"=2),selected=2)),
                                                                                         conditionalPanel(condition="input.validation_RNA_approaches==1",
                                                                                                          column(12,radioButtons("validation_RNA_url_sources", h5(strong("Url sources")), inline=TRUE,
                                                                                                                                 choices = c("Combined dataset"=1,"Separated datasets"=2),selected=1)),
                                                                                                          conditionalPanel(condition="input.validation_RNA_url_sources==1",
                                                                                                                           column(12,textAreaInput(inputId = "validation_RNA_downloadUrl",
                                                                                                                                                   label = "Download url for combined RNA dataset",
                                                                                                                                                   value = "",
                                                                                                                                                   width = '100%', height = '35px'))
                                                                                                          ),
                                                                                                          conditionalPanel(condition="input.validation_RNA_url_sources==2",
                                                                                                                           column(12,offset=0,tags$p(h5(strong("Input RNA dataset url for each cancer type",align="left")))),
                                                                                                                           column(9,rHandsontableOutput("table.validation_RNAUrl")),
                                                                                                                           column(12,offset=0,tags$p(h5(strong("",align="left"))))  #空一行
                                                                                                          ),
                                                                                                          column(12,actionButton("validation_RNA_process1","Process",width="100%",class="btn btn-primary"))
                                                                                         ),
                                                                                         conditionalPanel(condition="input.validation_RNA_approaches==2",
                                                                                                          column(12,radioButtons("validation_RNA_upload_sources", h5(strong("Upload sources")), inline=TRUE,
                                                                                                                                 choices = c("Combined dataset"=1,"Separated datasets"=2),selected=1)),
                                                                                                          conditionalPanel(condition="input.validation_RNA_upload_sources==1",
                                                                                                                           column(12,fileInput('validation_RNA_fileUpload1',tags$p(tags$b("Upload the downloaded combined RNA dataset (a .txt file)")),
                                                                                                                                               accept=c('text/csv', 'text/comma-separated-values,text/plain', '.csv') ))
                                                                                                          ),
                                                                                                          conditionalPanel(condition="input.validation_RNA_upload_sources==2",
                                                                                                                           column(12,fileInput('validation_RNA_fileUpload2',tags$p(tags$b("Upload the downloaded RNA dataset sequentially for each cancer type (a .txt file)")),
                                                                                                                                               accept=c('text/csv', 'text/comma-separated-values,text/plain', '.csv'), multiple=TRUE))
                                                                                                          ),
                                                                                                          column(12,actionButton("validation_RNA_process2","Process",width="100%",class="btn btn-primary"))
                                                                                         )
                                                                                     )
                                                                           )
                                                          ),
                                                          
                                                          #DNA methylation(Validation)
                                                          conditionalPanel(condition="input.validation_multiOmics.includes('3')",
                                                                           wellPanel(style="background-color:#47AEE9;",
                                                                                     fixedRow(
                                                                                         column(12,tags$h4(strong("DNA methylation(Validation)"),align="center",style="color:#834005")),
                                                                                         column(12,radioButtons("validation_DNA_methylation_approaches", h5(strong("Approaches")), inline=TRUE,
                                                                                                                choices = c("Download from specified urls"=1,"Upload manually"=2),selected=2)),
                                                                                         conditionalPanel(condition="input.validation_DNA_methylation_approaches==1",
                                                                                                          column(12,radioButtons("validation_DNA_methylation_url_sources", h5(strong("Url sources")), inline=TRUE,
                                                                                                                                 choices = c("Combined dataset"=1,"Separated datasets"=2),selected=1)),
                                                                                                          conditionalPanel(condition="input.validation_DNA_methylation_url_sources==1",
                                                                                                                           column(12,textAreaInput(inputId = "validation_DNA_methylation_downloadUrl",
                                                                                                                                                   label = "Download url for combined DNA methylation dataset",
                                                                                                                                                   value = "",
                                                                                                                                                   width = '100%', height = '35px'))
                                                                                                          ),
                                                                                                          conditionalPanel(condition="input.validation_DNA_methylation_url_sources==2",
                                                                                                                           column(12,offset=0,tags$p(h5(strong("Input DNA methylation dataset url for each cancer type",align="left")))),
                                                                                                                           column(9,rHandsontableOutput("table.validation_DNAMethylationUrl")),
                                                                                                                           column(12,offset=0,tags$p(h5(strong("",align="left"))))  #空一行
                                                                                                          ),
                                                                                                          column(12,actionButton("validation_DNA_methylation_process1","Process",width="100%",class="btn btn-primary"))
                                                                                         ),
                                                                                         conditionalPanel(condition="input.validation_DNA_methylation_approaches==2",
                                                                                                          column(12,radioButtons("validation_DNA_methylation_upload_sources", h5(strong("Upload sources")), inline=TRUE,
                                                                                                                                 choices = c("Combined dataset"=1,"Separated datasets"=2),selected=1)),
                                                                                                          conditionalPanel(condition="input.validation_DNA_methylation_upload_sources==1",
                                                                                                                           column(12,fileInput('validation_DNA_methylation_fileUpload1',tags$p(tags$b("Upload the downloaded combined DNA methylation dataset (a .gz file)")),
                                                                                                                                               accept=".gz" ))
                                                                                                          ),
                                                                                                          conditionalPanel(condition="input.validation_DNA_methylation_upload_sources==2",
                                                                                                                           column(12,fileInput('validation_DNA_methylation_fileUpload2',tags$p(tags$b("Upload the downloaded DNA methylation dataset sequentially for each cancer type (a .gz file)")),
                                                                                                                                               accept=".gz", multiple=TRUE))
                                                                                                          ),
                                                                                                          column(12,actionButton("validation_DNA_methylation_process2","Process",width="100%",class="btn btn-primary"))
                                                                                                          
                                                                                         )
                                                                                     )
                                                                           )
                                                          ),
                                                          
                                                          #copy number alterations(Validation)
                                                          conditionalPanel(condition="input.validation_multiOmics.indexOf('4') > -1",
                                                                           wellPanel(style="background-color:#47AEE9;",
                                                                                     fixedRow(
                                                                                         column(12,tags$h4(strong("copy number alterations(Validation)"),align="center",style="color:#834005")),
                                                                                         column(12,radioButtons("validation_copy_number_alterations_approaches", h5(strong("Approaches")), inline=TRUE,
                                                                                                                choices = c("Download from specified urls"=1,"Upload manually"=2),selected=2)),
                                                                                         conditionalPanel(condition="input.validation_copy_number_alterations_approaches==1",
                                                                                                          column(12,radioButtons("validation_copy_number_alterations_url_sources", h5(strong("Url sources")), inline=TRUE,
                                                                                                                                 choices = c("Combined dataset"=1,"Separated datasets"=2),selected=1)),
                                                                                                          conditionalPanel(condition="input.validation_copy_number_alterations_url_sources==1",
                                                                                                                           column(12,textAreaInput(inputId = "validation_copy_number_alterations_downloadUrl",
                                                                                                                                                   label = "Download url for combined copy number alterations dataset",
                                                                                                                                                   value = "",
                                                                                                                                                   width = '100%', height = '35px'))
                                                                                                          ),
                                                                                                          conditionalPanel(condition="input.validation_copy_number_alterations_url_sources==2",
                                                                                                                           column(12,offset=0,tags$p(h5(strong("Input copy number alterations dataset url for each cancer type",align="left")))),
                                                                                                                           column(9,rHandsontableOutput("table.validation_copyNumberAlterationsUrl")),
                                                                                                                           column(12,offset=0,tags$p(h5(strong("",align="left"))))  #空一行
                                                                                                          ),
                                                                                                          column(12,actionButton("validation_copy_number_alterations_process1","Process",width="100%",class="btn btn-primary"))
                                                                                         ),
                                                                                         conditionalPanel(condition="input.validation_copy_number_alterations_approaches==2",
                                                                                                          column(12,radioButtons("validation_copy_number_alterations_upload_sources", h5(strong("Upload sources")), inline=TRUE,
                                                                                                                                 choices = c("Combined dataset"=1,"Separated datasets"=2),selected=1)),
                                                                                                          conditionalPanel(condition="input.validation_copy_number_alterations_upload_sources==1",
                                                                                                                           column(12,fileInput('validation_copy_number_alterations_fileUpload1',tags$p(tags$b("Upload the downloaded combined copy number alterations dataset (a .tar.gz file)")),
                                                                                                                                               accept=".tar.gz" ))
                                                                                                          ),
                                                                                                          conditionalPanel(condition="input.validation_copy_number_alterations_upload_sources==2",
                                                                                                                           column(12,fileInput('validation_copy_number_alterations_fileUpload2',tags$p(tags$b("Upload the downloaded copy number alterations dataset sequentially for each cancer type (a .tar.gz file)")),
                                                                                                                                               accept=".tar.gz", multiple=TRUE))
                                                                                                          ),
                                                                                                          column(12,actionButton("validation_copy_number_alterations_process2","Process",width="100%",class="btn btn-primary"))
                                                                                         )
                                                                                     )
                                                                           )        
                                                          ),
                                                          
                                                          #binary somatic mutation(Validation)
                                                          conditionalPanel(condition="input.validation_multiOmics.indexOf('5') > -1",
                                                                           wellPanel(style="background-color:#47AEE9;",
                                                                                     fixedRow(
                                                                                         column(12,tags$h4(strong("binary somatic mutation(Validation)"),align="center",style="color:#834005")),
                                                                                         column(12,radioButtons("validation_binary_somatic_mutation_approaches", h5(strong("Approaches")), inline=TRUE,
                                                                                                                choices = c("Download from specified urls"=1,"Upload manually"=2),selected=2)),
                                                                                         conditionalPanel(condition="input.validation_binary_somatic_mutation_approaches==1",
                                                                                                          column(12,radioButtons("validation_binary_somatic_mutation_url_sources", h5(strong("Url sources")), inline=TRUE,
                                                                                                                                 choices = c("Combined dataset"=1,"Separated datasets"=2),selected=1)),
                                                                                                          conditionalPanel(condition="input.validation_binary_somatic_mutation_url_sources==1",
                                                                                                                           column(12,textAreaInput(inputId = "validation_binary_somatic_mutation_downloadUrl",
                                                                                                                                                   label = "Download url for combined binary somatic mutation dataset",
                                                                                                                                                   value = "",
                                                                                                                                                   width = '100%', height = '35px'))
                                                                                                          ),
                                                                                                          conditionalPanel(condition="input.validation_binary_somatic_mutation_url_sources==2",
                                                                                                                           column(12,offset=0,tags$p(h5(strong("Input binary somatic mutation dataset url for each cancer type",align="left")))),
                                                                                                                           column(9,rHandsontableOutput("table.validation_binarySomaticMutationUrl")),
                                                                                                                           column(12,offset=0,tags$p(h5(strong("",align="left"))))  #空一行
                                                                                                          ),
                                                                                                          column(12,actionButton("validation_binary_somatic_mutation_process1","Process",width="100%",class="btn btn-primary"))
                                                                                         ),
                                                                                         conditionalPanel(condition="input.validation_binary_somatic_mutation_approaches==2",
                                                                                                          column(12,radioButtons("validation_binary_somatic_mutation_upload_sources", h5(strong("Upload sources")), inline=TRUE,
                                                                                                                                 choices = c("Combined dataset"=1,"Separated datasets"=2),selected=1)),
                                                                                                          conditionalPanel(condition="input.validation_binary_somatic_mutation_upload_sources==1",
                                                                                                                           column(12,fileInput('validation_binary_somatic_mutation_fileUpload1',tags$p(tags$b("Upload the downloaded combined binary somatic mutation dataset (a .tar.gz file)")),
                                                                                                                                               accept=".tar.gz" ))
                                                                                                          ),
                                                                                                          conditionalPanel(condition="input.validation_binary_somatic_mutation_upload_sources==2",
                                                                                                                           column(12,fileInput('validation_binary_somatic_mutation_fileUpload2',tags$p(tags$b("Upload the downloaded binary somatic mutation dataset sequentially for each cancer type (a .tar.gz file)")),
                                                                                                                                               accept=".tar.gz", multiple=TRUE))
                                                                                                          ),
                                                                                                          column(12,actionButton("validation_binary_somatic_mutation_process2","Process",width="100%",class="btn btn-primary"))
                                                                                         )
                                                                                     )
                                                                           )
                                                          ),
                                                          
                                                          #radiomics
                                                          conditionalPanel(condition="input.validation_multiOmics.indexOf('6') > -1",
                                                                           wellPanel(style="background-color:#47AEE9;",
                                                                                     fixedRow(
                                                                                         column(12,tags$h4(strong("radiomics"),align="center",style="color:#834005")),
                                                                                         column(12,radioButtons("validation_external_radiomics_approaches", h5(strong("Approaches")), inline=TRUE,
                                                                                                                choices = c("Download from specified urls"=1,"Upload manually"=2),selected=2)),
                                                                                         
                                                                                         conditionalPanel(condition="input.validation_external_radiomics_approaches==1",
                                                                                                          column(12,radioButtons("validation_external_radiomics_url_sources", h5(strong("Url sources")), inline=TRUE,
                                                                                                                                 choices = c("Combined dataset"=1,"Separated datasets"=2),selected=1)),
                                                                                                          conditionalPanel(condition="input.validation_external_radiomics_url_sources==1",
                                                                                                                           column(12,textAreaInput(inputId = "validation_external_radiomics_downloadUrl",
                                                                                                                                                   label = "Download url for combined radiomics dataset",
                                                                                                                                                   value = "",
                                                                                                                                                   width = '100%', height = '35px'))
                                                                                                          ),
                                                                                                          conditionalPanel(condition="input.validation_external_radiomics_url_sources==2",
                                                                                                                           column(12,offset=0,tags$p(h5(strong("Input radiomics dataset url for each cancer type",align="left")))),
                                                                                                                           column(9,rHandsontableOutput("table.validation_externalRadiomicsUrl")),
                                                                                                                           column(12,offset=0,tags$p(h5(strong("",align="left"))))  #空一行
                                                                                                          ),
                                                                                                          column(12,actionButton("validation_external_radiomics_process1","Process",width="100%",class="btn btn-primary"))
                                                                                         ),
                                                                                         conditionalPanel(condition="input.validation_external_radiomics_approaches==2",
                                                                                                          column(12,radioButtons("validation_external_radiomics_upload_sources", h5(strong("Upload sources")), inline=TRUE,
                                                                                                                                 choices = c("Combined dataset"=1,"Separated datasets"=2),selected=1)),
                                                                                                          conditionalPanel(condition="input.validation_external_radiomics_upload_sources==1",
                                                                                                                           column(12,fileInput('validation_external_radiomics_fileUpload1',tags$p(tags$b("Upload the combined radiomics dataset (a .txt file)")),
                                                                                                                                               accept=c('text/csv', 'text/comma-separated-values,text/plain', '.csv') ))
                                                                                                          ),
                                                                                                          conditionalPanel(condition="input.validation_external_radiomics_upload_sources==2",
                                                                                                                           column(12,fileInput('validation_external_radiomics_fileUpload2',tags$p(tags$b("Upload the radiomics dataset sequentially for each cancer type (a .txt file)")),
                                                                                                                                               accept=c('text/csv', 'text/comma-separated-values,text/plain', '.csv'), multiple=TRUE))
                                                                                                          ),
                                                                                                          column(12,actionButton("validation_external_radiomics_process2","Process",width="100%",class="btn btn-primary"))
                                                                                         )
                                                                                     )
                                                                           )        
                                                          ),
                                                          
                                                          #Validation Datasets Integration
                                                          wellPanel(style="background-color:#47AEE9;",
                                                                    
                                                                    fixedRow(
                                                                        column(12,tags$h4(strong("Validation Datasets Integration"),align="center",style="color:#834005")),
                                                                        column(12,offset=0,tags$p(h5("Now all omics datasets you specified have been well prepared, then click the 'Integrate' button below to integrate all omics datasets.",align="left"))),
                                                                        column(12,actionButton("validation_datasets_integrate","Integrate",width="100%",class="btn btn-primary"))
                                                                    )
                                                          )        
                                         )           
                                    )
                    ),
                    br(),
                    column(7,
                           tabsetPanel(
                               
                               #Data Preparation (TCGA)
                               tabPanel(tags$h5(tags$strong("Data Preparation (TCGA)"),align="center",style="color:#FFA500"),br(),
                                        
                                        ###创建用户
                                        column(12,uiOutput("user_create")),
                                        
                                        ###TCGA
                                        ####使用内置数据
                                        column(12,uiOutput("internal_process_finish")),
                                        ####使用整合外部数据
                                        column(12,uiOutput("use_integrated_datasets_finish")),
                                        ####临床生存数据集
                                        column(12,uiOutput("clisur_download_finish")),
                                        column(12,uiOutput("clisur_finish")),
                                        ####RNA数据集
                                        column(12,uiOutput("RNA_download_finish")),
                                        column(12,uiOutput("RNA_finish")),
                                        ####DNA methylation数据集
                                        column(12,uiOutput("DNA_methylation_download_finish")),
                                        column(12,uiOutput("DNA_methylation_finish")),
                                        ####copy number alterations数据集
                                        column(12,uiOutput("copy_number_alterations_download_finish")),
                                        column(12,uiOutput("copy_number_alterations_finish")),
                                        ####binary somatic mutation数据集
                                        column(12,uiOutput("binary_somatic_mutation_download_finish")),
                                        column(12,uiOutput("binary_somatic_mutation_finish")),
                                        ####radiomics数据集
                                        column(12,uiOutput("radiomics_download_finish")),
                                        column(12,uiOutput("radiomics_finish"))
                                        
                                        # #有外部验证数据
                                        # conditionalPanel(condition="input.external_validation==1",
                                        #                  ##使用内置数据
                                        #                  conditionalPanel(condition="input.data_sources==1",
                                        # 
                                        #                                   ###创建用户
                                        #                                   column(12,uiOutput("user_create")),
                                        #                                   ###TCGA
                                        #                                   column(12,uiOutput("internal_process_finish")),
                                        #                                   column(12,uiOutput("datasets_integration_finish")),
                                        #                                   ###使用整合外部数据(Validation)
                                        #                                   conditionalPanel(condition="input.validation_use_integrated_datasets==1",
                                        # 
                                        #                                                    ####Validation文件夹
                                        #                                                    column(12,uiOutput("validation_create")),
                                        #                                                    ####Validation数据
                                        #                                                    column(12,uiOutput("validation_use_integrated_datasets_finish"))
                                        #                                   ),
                                        #                                   ###不使用整合外部数据(Validation)
                                        #                                   conditionalPanel(condition="input.validation_use_integrated_datasets==2",
                                        # 
                                        #                                                    ####Validation文件夹
                                        #                                                    column(12,uiOutput("validation_create")),
                                        #                                                    ####Validation数据
                                        #                                                    #####临床生存数据集
                                        #                                                    column(12,uiOutput("validation_clisur_download_finish")),
                                        #                                                    #####RNA数据集
                                        #                                                    column(12,uiOutput("validation_RNA_download_finish")),
                                        #                                                    #####DNA methylation数据集
                                        #                                                    column(12,uiOutput("validation_DNA_methylation_download_finish")),
                                        #                                                    #####copy number alterations数据集
                                        #                                                    column(12,uiOutput("validation_copy_number_alterations_download_finish")),
                                        #                                                    #####binary somatic mutation数据集
                                        #                                                    column(12,uiOutput("validation_binary_somatic_mutation_download_finish")),
                                        # 
                                        #                                                    column(12,uiOutput("validation_datasets_integration_finish"))
                                        #                                   )
                                        #                  ),
                                        #                  ##使用外部数据
                                        #                  conditionalPanel(condition="input.data_sources==2",
                                        # 
                                        #                                   ###创建用户
                                        #                                   column(12,uiOutput("user_create")),
                                        #                                   ###TCGA
                                        #                                   ####使用整合外部数据
                                        #                                   conditionalPanel(condition="input.use_integrated_datasets==1",
                                        # 
                                        #                                                    column(12,uiOutput("use_integrated_datasets_finish"))
                                        #                                   ),
                                        #                                   ####不使用整合外部数据
                                        #                                   conditionalPanel(condition="input.use_integrated_datasets==2",
                                        # 
                                        #                                                    #####临床生存数据集
                                        #                                                    column(12,uiOutput("clisur_download_finish")),
                                        #                                                    #####RNA数据集
                                        #                                                    column(12,uiOutput("RNA_download_finish")),
                                        #                                                    #####DNA methylation数据集
                                        #                                                    column(12,uiOutput("DNA_methylation_download_finish")),
                                        #                                                    #####copy number alterations数据集
                                        #                                                    column(12,uiOutput("copy_number_alterations_download_finish")),
                                        #                                                    #####binary somatic mutation数据集
                                        #                                                    column(12,uiOutput("binary_somatic_mutation_download_finish")),
                                        # 
                                        #                                                    column(12,uiOutput("datasets_integration_finish"))
                                        #                                   ),
                                        # 
                                        #                                   ###使用整合外部数据(Validation)
                                        #                                   conditionalPanel(condition="input.validation_use_integrated_datasets==1",
                                        # 
                                        #                                                    ####Validation文件夹
                                        #                                                    column(12,uiOutput("validation_create")),
                                        #                                                    ####Validation数据
                                        #                                                    column(12,uiOutput("validation_use_integrated_datasets_finish"))
                                        #                                   ),
                                        #                                   ###不使用整合外部数据(Validation)
                                        #                                   conditionalPanel(condition="input.validation_use_integrated_datasets==2",
                                        # 
                                        #                                                    ####Validation文件夹
                                        #                                                    column(12,uiOutput("validation_create")),
                                        #                                                    ####Validation数据
                                        #                                                    column(12,uiOutput("validation_use_integrated_datasets_finish")),
                                        #                                                    #####临床生存数据集
                                        #                                                    column(12,uiOutput("validation_clisur_download_finish")),
                                        #                                                    #####RNA数据集
                                        #                                                    column(12,uiOutput("validation_RNA_download_finish")),
                                        #                                                    #####DNA methylation数据集
                                        #                                                    column(12,uiOutput("validation_DNA_methylation_download_finish")),
                                        #                                                    #####copy number alterations数据集
                                        #                                                    column(12,uiOutput("validation_copy_number_alterations_download_finish")),
                                        #                                                    #####binary somatic mutation数据集
                                        #                                                    column(12,uiOutput("validation_binary_somatic_mutation_download_finish")),
                                        # 
                                        #                                                    column(12,uiOutput("validation_datasets_integration_finish"))
                                        #                                   )
                                        #                  )
                                        # ),
                                        # #无外部验证数据
                                        # conditionalPanel(condition="input.external_validation==2",
                                        # 
                                        #                  ##使用内置数据
                                        #                  conditionalPanel(condition="input.data_sources==1",
                                        # 
                                        #                                   ###创建用户
                                        #                                   column(12,uiOutput("user_create")),
                                        #                                   ###TCGA
                                        #                                   column(12,uiOutput("internal_process_finish")),
                                        #                                   column(12,uiOutput("datasets_integration_finish")),
                                        #                  ),
                                        #                  ##使用外部数据
                                        #                  conditionalPanel(condition="input.data_sources==2",
                                        # 
                                        #                                   ###创建用户
                                        #                                   column(12,uiOutput("user_create")),
                                        #                                   ###TCGA
                                        #                                   ####使用整合外部数据
                                        #                                   conditionalPanel(condition="input.use_integrated_datasets==1",
                                        # 
                                        #                                                    column(12,uiOutput("use_integrated_datasets_finish"))
                                        #                                   ),
                                        #                                   ####不使用整合外部数据
                                        #                                   conditionalPanel(condition="input.use_integrated_datasets==2",
                                        # 
                                        #                                                    #####临床生存数据集
                                        #                                                    column(12,uiOutput("clisur_download_finish")),
                                        #                                                    #####RNA数据集
                                        #                                                    column(12,uiOutput("RNA_download_finish")),
                                        #                                                    #####DNA methylation数据集
                                        #                                                    column(12,uiOutput("DNA_methylation_download_finish")),
                                        #                                                    #####copy number alterations数据集
                                        #                                                    column(12,uiOutput("copy_number_alterations_download_finish")),
                                        #                                                    #####binary somatic mutation数据集
                                        #                                                    column(12,uiOutput("binary_somatic_mutation_download_finish")),
                                        # 
                                        #                                                    column(12,uiOutput("datasets_integration_finish"))
                                        #                                   )
                                        #                  )
                                        # )
                                        # # #文件上传控件测试代码
                                        # # tableOutput("files"),
                                        # # column(12,uiOutput("files_show"))
                               ),
                               
                               #Data Integration (TCGA)
                               tabPanel(tags$h5(tags$strong("Data Integration (TCGA)"),align="center",style="color:#FFA500"),br(),
                                        
                                        ####外部数据整合
                                        column(12,uiOutput("datasets_integration_finish"))
                                        # column(12,uiOutput("datasets_integration_finish")),
                                        # br(),br(),
                                        # ####整合数据展示
                                        # #####mRNA(tpm)
                                        # uiOutput("table.integrated_mRNA_tpm_TCGA_Title"), #标题
                                        # column(12,br()),
                                        # dataTableOutput("table.integrated_mRNA_tpm_TCGA"), #表格
                                        # # uiOutput("table.integrated_mRNA_tpm_TCGA_Note"), #脚注
                                        # br(),br(),
                                        # #####lncRNA(tpm)
                                        # uiOutput("table.integrated_lncRNA_tpm_TCGA_Title"), #标题
                                        # column(12,br()),
                                        # dataTableOutput("table.integrated_lncRNA_tpm_TCGA"), #表格
                                        # # uiOutput("table.integrated_lncRNA_tpm_TCGA_Note"), #脚注
                                        # br(),br(),
                                        # #####DNA methylation
                                        # uiOutput("table.integrated_DNA_methylation_TCGA_Title"), #标题
                                        # column(12,br()),
                                        # dataTableOutput("table.integrated_DNA_methylation_TCGA"), #表格
                                        # # uiOutput("table.integrated_DNA_methylation_TCGA_Note"), #脚注
                                        # br(),br(),
                                        # #####copy number alterations
                                        # uiOutput("table.integrated_copy_number_alterations_TCGA_Title"), #标题
                                        # column(12,br()),
                                        # dataTableOutput("table.integrated_copy_number_alterations_TCGA"), #表格
                                        # # uiOutput("table.integrated_copy_number_alterations_TCGA_Note"), #脚注
                                        # br(),br(),
                                        # #####binary somatic mutation
                                        # uiOutput("table.integrated_binary_somatic_mutation_TCGA_Title"), #标题
                                        # column(12,br()),
                                        # dataTableOutput("table.integrated_binary_somatic_mutation_TCGA"), #表格
                                        # # uiOutput("table.integrated_binary_somatic_mutation_TCGA_Note"), #脚注
                                        # br(),br(),
                                        # #####radiomics
                                        # uiOutput("table.integrated_radiomics_TCGA_Title"), #标题
                                        # column(12,br()),
                                        # dataTableOutput("table.integrated_radiomics_TCGA"), #表格
                                        # # uiOutput("table.integrated_radiomics_TCGA_Note"), #脚注
                                        # br(),br(),
                                        # #####count
                                        # uiOutput("table.integrated_count_TCGA_Title"), #标题
                                        # column(12,br()),
                                        # dataTableOutput("table.integrated_count_TCGA"), #表格
                                        # # uiOutput("table.integrated_count_TCGA_Note"), #脚注
                                        # br(),br(),
                                        # #####mRNA(fpkm)
                                        # uiOutput("table.integrated_mRNA_fpkm_TCGA_Title"), #标题
                                        # column(12,br()),
                                        # dataTableOutput("table.integrated_mRNA_fpkm_TCGA"), #表格
                                        # # uiOutput("table.integrated_mRNA_fpkm_TCGA_Note"), #脚注
                                        # br(),br(),
                                        # #####lncRNA(fpkm)
                                        # uiOutput("table.integrated_lncRNA_fpkm_TCGA_Title"), #标题
                                        # column(12,br()),
                                        # dataTableOutput("table.integrated_lncRNA_fpkm_TCGA"), #表格
                                        # # uiOutput("table.integrated_lncRNA_fpkm_TCGA_Note"), #脚注
                                        # br(),br(),
                                        # #####segment
                                        # uiOutput("table.integrated_segment_TCGA_Title"), #标题
                                        # column(12,br()),
                                        # dataTableOutput("table.integrated_segment_TCGA"), #表格
                                        # # uiOutput("table.integrated_segment_TCGA_Note"), #脚注
                                        # br(),br(),
                                        # #####maf
                                        # uiOutput("table.integrated_maf_TCGA_Title"), #标题
                                        # column(12,br()),
                                        # dataTableOutput("table.integrated_maf_TCGA"), #表格
                                        # # uiOutput("table.integrated_maf_TCGA_Note"), #脚注
                                        # br(),br(),
                                        # #####clinSurv
                                        # uiOutput("table.integrated_clinSurv_TCGA_Title"), #标题
                                        # column(12,br()),
                                        # dataTableOutput("table.integrated_clinSurv_TCGA") #表格
                                        # # uiOutput("table.integrated_clinSurv_TCGA_Note") #脚注
                               ),
                               
                               #Data Preparation (Validation)
                               tabPanel(tags$h5(tags$strong("Data Preparation (Validation)"),align="center",style="color:#FFA500"),br(),
                                        
                                        ###Validation
                                        ####Validation文件夹
                                        column(12,uiOutput("validation_create")),
                                        ####使用整合外部数据
                                        column(12,uiOutput("validation_use_integrated_datasets_finish")),
                                        ####临床生存数据集
                                        column(12,uiOutput("validation_clisur_download_finish")),
                                        column(12,uiOutput("validation_clisur_finish")),
                                        ####RNA数据集
                                        column(12,uiOutput("validation_RNA_download_finish")),
                                        column(12,uiOutput("validation_RNA_finish")),
                                        ####DNA methylation数据集
                                        column(12,uiOutput("validation_DNA_methylation_download_finish")),
                                        column(12,uiOutput("validation_DNA_methylation_finish")),
                                        ####copy number alterations数据集
                                        column(12,uiOutput("validation_copy_number_alterations_download_finish")),
                                        column(12,uiOutput("validation_copy_number_alterations_finish")),
                                        ####binary somatic mutation数据集
                                        column(12,uiOutput("validation_binary_somatic_mutation_download_finish")),
                                        column(12,uiOutput("validation_binary_somatic_mutation_finish")),
                                        ####radiomics数据集
                                        column(12,uiOutput("validation_radiomics_download_finish")),
                                        column(12,uiOutput("validation_radiomics_finish"))
                                        
                                        # #有外部验证数据
                                        # conditionalPanel(condition="input.external_validation==1",
                                        #                  ##使用内置数据
                                        #                  conditionalPanel(condition="input.data_sources==1",
                                        # 
                                        #                                   ###创建用户
                                        #                                   column(12,uiOutput("user_create")),
                                        #                                   ###TCGA
                                        #                                   column(12,uiOutput("internal_process_finish")),
                                        #                                   column(12,uiOutput("datasets_integration_finish")),
                                        #                                   ###使用整合外部数据(Validation)
                                        #                                   conditionalPanel(condition="input.validation_use_integrated_datasets==1",
                                        # 
                                        #                                                    ####Validation文件夹
                                        #                                                    column(12,uiOutput("validation_create")),
                                        #                                                    ####Validation数据
                                        #                                                    column(12,uiOutput("validation_use_integrated_datasets_finish"))
                                        #                                   ),
                                        #                                   ###不使用整合外部数据(Validation)
                                        #                                   conditionalPanel(condition="input.validation_use_integrated_datasets==2",
                                        # 
                                        #                                                    ####Validation文件夹
                                        #                                                    column(12,uiOutput("validation_create")),
                                        #                                                    ####Validation数据
                                        #                                                    #####临床生存数据集
                                        #                                                    column(12,uiOutput("validation_clisur_download_finish")),
                                        #                                                    #####RNA数据集
                                        #                                                    column(12,uiOutput("validation_RNA_download_finish")),
                                        #                                                    #####DNA methylation数据集
                                        #                                                    column(12,uiOutput("validation_DNA_methylation_download_finish")),
                                        #                                                    #####copy number alterations数据集
                                        #                                                    column(12,uiOutput("validation_copy_number_alterations_download_finish")),
                                        #                                                    #####binary somatic mutation数据集
                                        #                                                    column(12,uiOutput("validation_binary_somatic_mutation_download_finish")),
                                        # 
                                        #                                                    column(12,uiOutput("validation_datasets_integration_finish"))
                                        #                                   )
                                        #                  ),
                                        #                  ##使用外部数据
                                        #                  conditionalPanel(condition="input.data_sources==2",
                                        # 
                                        #                                   ###创建用户
                                        #                                   column(12,uiOutput("user_create")),
                                        #                                   ###TCGA
                                        #                                   ####使用整合外部数据
                                        #                                   conditionalPanel(condition="input.use_integrated_datasets==1",
                                        # 
                                        #                                                    column(12,uiOutput("use_integrated_datasets_finish"))
                                        #                                   ),
                                        #                                   ####不使用整合外部数据
                                        #                                   conditionalPanel(condition="input.use_integrated_datasets==2",
                                        # 
                                        #                                                    #####临床生存数据集
                                        #                                                    column(12,uiOutput("clisur_download_finish")),
                                        #                                                    #####RNA数据集
                                        #                                                    column(12,uiOutput("RNA_download_finish")),
                                        #                                                    #####DNA methylation数据集
                                        #                                                    column(12,uiOutput("DNA_methylation_download_finish")),
                                        #                                                    #####copy number alterations数据集
                                        #                                                    column(12,uiOutput("copy_number_alterations_download_finish")),
                                        #                                                    #####binary somatic mutation数据集
                                        #                                                    column(12,uiOutput("binary_somatic_mutation_download_finish")),
                                        # 
                                        #                                                    column(12,uiOutput("datasets_integration_finish"))
                                        #                                   ),
                                        # 
                                        #                                   ###使用整合外部数据(Validation)
                                        #                                   conditionalPanel(condition="input.validation_use_integrated_datasets==1",
                                        # 
                                        #                                                    ####Validation文件夹
                                        #                                                    column(12,uiOutput("validation_create")),
                                        #                                                    ####Validation数据
                                        #                                                    column(12,uiOutput("validation_use_integrated_datasets_finish"))
                                        #                                   ),
                                        #                                   ###不使用整合外部数据(Validation)
                                        #                                   conditionalPanel(condition="input.validation_use_integrated_datasets==2",
                                        # 
                                        #                                                    ####Validation文件夹
                                        #                                                    column(12,uiOutput("validation_create")),
                                        #                                                    ####Validation数据
                                        #                                                    column(12,uiOutput("validation_use_integrated_datasets_finish")),
                                        #                                                    #####临床生存数据集
                                        #                                                    column(12,uiOutput("validation_clisur_download_finish")),
                                        #                                                    #####RNA数据集
                                        #                                                    column(12,uiOutput("validation_RNA_download_finish")),
                                        #                                                    #####DNA methylation数据集
                                        #                                                    column(12,uiOutput("validation_DNA_methylation_download_finish")),
                                        #                                                    #####copy number alterations数据集
                                        #                                                    column(12,uiOutput("validation_copy_number_alterations_download_finish")),
                                        #                                                    #####binary somatic mutation数据集
                                        #                                                    column(12,uiOutput("validation_binary_somatic_mutation_download_finish")),
                                        # 
                                        #                                                    column(12,uiOutput("validation_datasets_integration_finish"))
                                        #                                   )
                                        #                  )
                                        # ),
                                        # #无外部验证数据
                                        # conditionalPanel(condition="input.external_validation==2",
                                        # 
                                        #                  ##使用内置数据
                                        #                  conditionalPanel(condition="input.data_sources==1",
                                        # 
                                        #                                   ###创建用户
                                        #                                   column(12,uiOutput("user_create")),
                                        #                                   ###TCGA
                                        #                                   column(12,uiOutput("internal_process_finish")),
                                        #                                   column(12,uiOutput("datasets_integration_finish")),
                                        #                  ),
                                        #                  ##使用外部数据
                                        #                  conditionalPanel(condition="input.data_sources==2",
                                        # 
                                        #                                   ###创建用户
                                        #                                   column(12,uiOutput("user_create")),
                                        #                                   ###TCGA
                                        #                                   ####使用整合外部数据
                                        #                                   conditionalPanel(condition="input.use_integrated_datasets==1",
                                        # 
                                        #                                                    column(12,uiOutput("use_integrated_datasets_finish"))
                                        #                                   ),
                                        #                                   ####不使用整合外部数据
                                        #                                   conditionalPanel(condition="input.use_integrated_datasets==2",
                                        # 
                                        #                                                    #####临床生存数据集
                                        #                                                    column(12,uiOutput("clisur_download_finish")),
                                        #                                                    #####RNA数据集
                                        #                                                    column(12,uiOutput("RNA_download_finish")),
                                        #                                                    #####DNA methylation数据集
                                        #                                                    column(12,uiOutput("DNA_methylation_download_finish")),
                                        #                                                    #####copy number alterations数据集
                                        #                                                    column(12,uiOutput("copy_number_alterations_download_finish")),
                                        #                                                    #####binary somatic mutation数据集
                                        #                                                    column(12,uiOutput("binary_somatic_mutation_download_finish")),
                                        # 
                                        #                                                    column(12,uiOutput("datasets_integration_finish"))
                                        #                                   )
                                        #                  )
                                        # )
                                        # # #文件上传控件测试代码
                                        # # tableOutput("files"),
                                        # # column(12,uiOutput("files_show"))
                               ),
                               
                               #Data Integration (Validation)
                               tabPanel(tags$h5(tags$strong("Data Integration (Validation)"),align="center",style="color:#FFA500"),br(),
                                        
                                        ####外部数据整合
                                        column(12,uiOutput("validation_datasets_integration_finish"))
                                        # column(12,uiOutput("validation_datasets_integration_finish")),
                                        # br(),br(),
                                        # ####整合数据展示
                                        # #####mRNA(tpm)
                                        # uiOutput("table.integrated_mRNA_tpm_Validation_Title"), #标题
                                        # column(12,br()),
                                        # dataTableOutput("table.integrated_mRNA_tpm_Validation"), #表格
                                        # # uiOutput("table.integrated_mRNA_tpm_Validation_Note"), #脚注
                                        # br(),br(),
                                        # #####lncRNA(tpm)
                                        # uiOutput("table.integrated_lncRNA_tpm_Validation_Title"), #标题
                                        # column(12,br()),
                                        # dataTableOutput("table.integrated_lncRNA_tpm_Validation"), #表格
                                        # # uiOutput("table.integrated_lncRNA_tpm_Validation_Note"), #脚注
                                        # br(),br(),
                                        # #####DNA methylation
                                        # uiOutput("table.integrated_DNA_methylation_Validation_Title"), #标题
                                        # column(12,br()),
                                        # dataTableOutput("table.integrated_DNA_methylation_Validation"), #表格
                                        # # uiOutput("table.integrated_DNA_methylation_Validation_Note"), #脚注
                                        # br(),br(),
                                        # #####copy number alterations
                                        # uiOutput("table.integrated_copy_number_alterations_Validation_Title"), #标题
                                        # column(12,br()),
                                        # dataTableOutput("table.integrated_copy_number_alterations_Validation"), #表格
                                        # # uiOutput("table.integrated_copy_number_alterations_Validation_Note"), #脚注
                                        # br(),br(),
                                        # #####binary somatic mutation
                                        # uiOutput("table.integrated_binary_somatic_mutation_Validation_Title"), #标题
                                        # column(12,br()),
                                        # dataTableOutput("table.integrated_binary_somatic_mutation_Validation"), #表格
                                        # # uiOutput("table.integrated_binary_somatic_mutation_Validation_Note"), #脚注
                                        # br(),br(),
                                        # #####radiomics
                                        # uiOutput("table.integrated_radiomics_Validation_Title"), #标题
                                        # column(12,br()),
                                        # dataTableOutput("table.integrated_radiomics_Validation"), #表格
                                        # # uiOutput("table.integrated_radiomics_Validation_Note"), #脚注
                                        # br(),br(),
                                        # #####count
                                        # uiOutput("table.integrated_count_Validation_Title"), #标题
                                        # column(12,br()),
                                        # dataTableOutput("table.integrated_count_Validation"), #表格
                                        # # uiOutput("table.integrated_count_Validation_Note"), #脚注
                                        # br(),br(),
                                        # #####mRNA(fpkm)
                                        # uiOutput("table.integrated_mRNA_fpkm_Validation_Title"), #标题
                                        # column(12,br()),
                                        # dataTableOutput("table.integrated_mRNA_fpkm_Validation"), #表格
                                        # # uiOutput("table.integrated_mRNA_fpkm_Validation_Note"), #脚注
                                        # br(),br(),
                                        # #####lncRNA(fpkm)
                                        # uiOutput("table.integrated_lncRNA_fpkm_Validation_Title"), #标题
                                        # column(12,br()),
                                        # dataTableOutput("table.integrated_lncRNA_fpkm_Validation"), #表格
                                        # # uiOutput("table.integrated_lncRNA_fpkm_Validation_Note"), #脚注
                                        # br(),br(),
                                        # #####segment
                                        # uiOutput("table.integrated_segment_Validation_Title"), #标题
                                        # column(12,br()),
                                        # dataTableOutput("table.integrated_segment_Validation"), #表格
                                        # # uiOutput("table.integrated_segment_Validation_Note"), #脚注
                                        # br(),br(),
                                        # #####maf
                                        # uiOutput("table.integrated_maf_Validation_Title"), #标题
                                        # column(12,br()),
                                        # dataTableOutput("table.integrated_maf_Validation"), #表格
                                        # # uiOutput("table.integrated_maf_Validation_Note"), #脚注
                                        # br(),br(),
                                        # #####clinSurv
                                        # uiOutput("table.integrated_clinSurv_Validation_Title"), #标题
                                        # column(12,br()),
                                        # dataTableOutput("table.integrated_clinSurv_Validation") #表格
                                        # # uiOutput("table.integrated_clinSurv_Validation_Note") #脚注
                               )
                           )
                    )
                    
                ),
                
                #GET Module
                tabPanel(
                    title=tags$p(tags$h4(tags$strong("GET Module"),style="color:purple")),
                    value="GET Module",
                    
                    column(5,
                           
                           #1.总的步骤选项,对应不同的参数列表
                           wellPanel(style="background-color:#47AEE9;",
                                     
                                     fixedRow(
                                         #设置步骤选项,每个步骤对应不同的参数
                                         column(12,radioButtons("getModuleStep", h5(strong("Steps")), inline=FALSE,
                                                                choices = c("Get Elites"=1,"Get Clustering Number"=2,"Consensus Clustering"=3,"Silhouette"=4,"Multi-omics Heatmaps"=5),selected=1))
                                     )
                           ),
                           
                           #2.Get Elites
                           conditionalPanel(condition="input.getModuleStep==1",
                                wellPanel(style="background-color:#47AEE9",
                                    fixedRow(
                                        column(12,tags$h4(strong("Get Elites"),align="center",style="color:#834005")),
                                        column(12,offset=0,tags$p(h5(strong("In this step, we will filter out features that meet some stringent requirements as well as handle missing values. Now let's choose the datasets for 'Get Elites' first:",align="left")))),
                                        #Get elites on tcga datasets or validation datasets
                                        column(12,radioButtons("tcga_vali_getElites", h5(strong("Get elites on tcga datasets or validation datasets")), inline=TRUE, choices = c("TCGA"=1,"Validation"=2),selected=1)),
                                        
                                        conditionalPanel(condition="input.tcga_vali_getElites==1",
                                                         
                                                         conditionalPanel(condition="input.multiOmics.includes('1')",
                                                                          column(12,radioButtons("mRNA_getElites", h5(strong("Get Elites for mRNA dataset or not")), inline=TRUE, choices = c("Yes"=1,"No"=2),selected=1))
                                                         ),
                                                         conditionalPanel(condition="input.multiOmics.includes('2')",
                                                                          column(12,radioButtons("lncRNA_getElites", h5(strong("Get Elites for lncRNA dataset or not")), inline=TRUE, choices = c("Yes"=1,"No"=2),selected=1))
                                                         ),
                                                         conditionalPanel(condition="input.multiOmics.includes('3')",
                                                                          column(12,radioButtons("DNA_methylation_getElites", h5(strong("Get Elites for DNA methylation dataset or not")), inline=TRUE, choices = c("Yes"=1,"No"=2),selected=1))
                                                         ),
                                                         conditionalPanel(condition="input.multiOmics.includes('4')",
                                                                          column(12,radioButtons("copy_number_alterations_getElites", h5(strong("Get Elites for copy number alterations dataset or not")), inline=TRUE, choices = c("Yes"=1,"No"=2),selected=1))
                                                         ),
                                                         conditionalPanel(condition="input.multiOmics.includes('5')",
                                                                          column(12,radioButtons("binary_somatic_mutation_getElites", h5(strong("Get Elites for binary somatic mutation dataset or not")), inline=TRUE, choices = c("Yes"=1,"No"=2),selected=1))
                                                         ),
                                                         conditionalPanel(condition="input.multiOmics.includes('6')",
                                                                          column(12,radioButtons("radiomics_getElites", h5(strong("Get Elites for radiomics dataset or not")), inline=TRUE, choices = c("Yes"=1,"No"=2),selected=1))
                                                         )
                                        ),
                                        conditionalPanel(condition="input.tcga_vali_getElites==2",
                                                         
                                                         conditionalPanel(condition="input.validation_multiOmics.includes('1')",
                                                                          column(12,radioButtons("validation_mRNA_getElites", h5(strong("Get Elites for mRNA dataset or not")), inline=TRUE, choices = c("Yes"=1,"No"=2),selected=1))
                                                         ),
                                                         conditionalPanel(condition="input.validation_multiOmics.includes('2')",
                                                                          column(12,radioButtons("validation_lncRNA_getElites", h5(strong("Get Elites for lncRNA dataset or not")), inline=TRUE, choices = c("Yes"=1,"No"=2),selected=1))
                                                         ),
                                                         conditionalPanel(condition="input.validation_multiOmics.includes('3')",
                                                                          column(12,radioButtons("validation_DNA_methylation_getElites", h5(strong("Get Elites for DNA methylation dataset or not")), inline=TRUE, choices = c("Yes"=1,"No"=2),selected=1))
                                                         ),
                                                         conditionalPanel(condition="input.validation_multiOmics.includes('4')",
                                                                          column(12,radioButtons("validation_copy_number_alterations_getElites", h5(strong("Get Elites for copy number alterations dataset or not")), inline=TRUE, choices = c("Yes"=1,"No"=2),selected=1))
                                                         ),
                                                         conditionalPanel(condition="input.validation_multiOmics.includes('5')",
                                                                          column(12,radioButtons("validation_binary_somatic_mutation_getElites", h5(strong("Get Elites for binary somatic mutation dataset or not")), inline=TRUE, choices = c("Yes"=1,"No"=2),selected=1))
                                                         ),
                                                         conditionalPanel(condition="input.validation_multiOmics.includes('6')",
                                                                          column(12,radioButtons("validation_radiomics_getElites", h5(strong("Get Elites for radiomics dataset or not")), inline=TRUE, choices = c("Yes"=1,"No"=2),selected=1))
                                                         )
                                        )
                                    )
                                ),
                                
                                #TCGA
                                conditionalPanel(condition="input.tcga_vali_getElites==1",
                                                 
                                                 #Get Elites for mRNA dataset
                                                 conditionalPanel(condition="input.multiOmics.includes('1')",
                                                                  conditionalPanel(condition="input.mRNA_getElites==1",
                                                                                   wellPanel(style="background-color:#47AEE9",
                                                                                             fixedRow(
                                                                                                 column(12,tags$h4(strong("Get Elites settings for mRNA dataset"),align="center",style="color:#834005")),
                                                                                                 column(12,radioButtons("mRNA_getElites_survival", h5(strong("Use survival information or not")), inline=TRUE, choices = c("Yes"=1,"No"=2),selected=2)),
                                                                                                 column(12,radioButtons("mRNA_getElites_na", h5(strong("NA value action")), inline=TRUE, choices = c("Remove directly"=1,"KNN imputation"=2,"No action"=3),selected=1)),
                                                                                                 column(12,radioButtons("mRNA_getElites_log2", h5(strong("Perform log2 transformation for data before calculating statistics or not")), inline=TRUE, choices = c("Yes"=1,"No"=2),selected=2)),
                                                                                                 column(12,radioButtons("mRNA_getElites_lowpct", h5(strong("Set a numeric cutoff for removing low expression features or not")), inline=TRUE, choices = c("Yes"=1,"No"=2),selected=2)),
                                                                                                 conditionalPanel(condition = "input.mRNA_getElites_lowpct==1",
                                                                                                                  column(12,numericInput("mRNA_getElites_lowpct_value",label="The cutoff for removing low expression features",value=0.1,min=0,max=1,step=0.05))
                                                                                                 ),
                                                                                                 conditionalPanel(condition = "input.mRNA_getElites_survival==1",
                                                                                                                  column(12,numericInput("mRNA_getElites_p_cutoff",label="A numeric cutoff for nominal p value derived from univariate Cox proportional hazards regression",value=0.05,min=0,max=1,step=0.01))
                                                                                                 ),
                                                                                                 conditionalPanel(condition = "input.mRNA_getElites_survival==2",
                                                                                                                  column(12,radioButtons("mRNA_getElites_method", h5(strong("Choose a Get Elites method for mRNA dataset")), inline=TRUE, choices = c("mad"=1,"sd"=2,"pca"=3),selected=1)),
                                                                                                                  conditionalPanel(condition = "input.mRNA_getElites_method==1",
                                                                                                                                   column(12,radioButtons("mRNA_getElites_filtration1", h5(strong("Choose a filtration method for mRNA dataset")), inline=TRUE, choices = c("elite.num"=1,"elite.pct"=2),selected=1)),
                                                                                                                                   conditionalPanel(condition = "input.mRNA_getElites_filtration1==1",column(12,numericInput("mRNA_getElites_elite_num1",label="An integer cutoff of exact number for selecting top elites",value=1000,min=0,step=1))),
                                                                                                                                   conditionalPanel(condition = "input.mRNA_getElites_filtration1==2",column(12,numericInput("mRNA_getElites_elite_pct1",label="A numeric cutoff of percentage for selecting top elites",value=0.1,min=0,max=1,step=0.05)))
                                                                                                                  ),
                                                                                                                  conditionalPanel(condition = "input.mRNA_getElites_method==2",
                                                                                                                                   column(12,radioButtons("mRNA_getElites_filtration2", h5(strong("Choose a filtration method for mRNA dataset")), inline=TRUE, choices = c("elite.num"=1,"elite.pct"=2),selected=1)),
                                                                                                                                   conditionalPanel(condition = "input.mRNA_getElites_filtration2==1",column(12,numericInput("mRNA_getElites_elite_num2",label="An integer cutoff of exact number for selecting top elites",value=1000,min=0,step=1))),
                                                                                                                                   conditionalPanel(condition = "input.mRNA_getElites_filtration2==2",column(12,numericInput("mRNA_getElites_elite_pct2",label="A numeric cutoff of percentage for selecting top elites",value=0.1,min=0,max=1,step=0.05)))
                                                                                                                  ),
                                                                                                                  conditionalPanel(condition = "input.mRNA_getElites_method==3",
                                                                                                                                   column(12,numericInput("mRNA_getElites_pca_ratio",label="A numeric value which represents the ratio of principal components",value=0.9,min=0,max=1,step=0.05))
                                                                                                                  )
                                                                                                 ),
                                                                                                 column(12,radioButtons("mRNA_getElites_scaleFlag", h5(strong("Scaling the data after filtering or not")), inline=TRUE, choices = c("Yes"=1,"No"=2),selected=2)),
                                                                                                 column(12,radioButtons("mRNA_getElites_centerFlag", h5(strong("Centering the data after filtering or not")), inline=TRUE, choices = c("Yes"=1,"No"=2),selected=2))
                                                                                             )
                                                                                   )
                                                                  )           
                                                 ),
                                                 
                                                 #Get Elites for lncRNA dataset
                                                 conditionalPanel(condition="input.multiOmics.includes('2')",
                                                                  conditionalPanel(condition="input.lncRNA_getElites==1",
                                                                                   wellPanel(style="background-color:#47AEE9",
                                                                                             fixedRow(
                                                                                                 column(12,tags$h4(strong("Get Elites settings for lncRNA dataset"),align="center",style="color:#834005")),
                                                                                                 column(12,radioButtons("lncRNA_getElites_survival", h5(strong("Use survival information or not")), inline=TRUE, choices = c("Yes"=1,"No"=2),selected=2)),
                                                                                                 column(12,radioButtons("lncRNA_getElites_na", h5(strong("NA value action")), inline=TRUE, choices = c("Remove directly"=1,"KNN imputation"=2,"No action"=3),selected=1)),
                                                                                                 column(12,radioButtons("lncRNA_getElites_log2", h5(strong("Perform log2 transformation for data before calculating statistics or not")), inline=TRUE, choices = c("Yes"=1,"No"=2),selected=2)),
                                                                                                 column(12,radioButtons("lncRNA_getElites_lowpct", h5(strong("Set a numeric cutoff for removing low expression features or not")), inline=TRUE, choices = c("Yes"=1,"No"=2),selected=2)),
                                                                                                 conditionalPanel(condition = "input.lncRNA_getElites_lowpct==1",
                                                                                                                  column(12,numericInput("lncRNA_getElites_lowpct_value",label="The cutoff for removing low expression features",value=0.1,min=0,max=1,step=0.05))
                                                                                                 ),
                                                                                                 conditionalPanel(condition = "input.lncRNA_getElites_survival==1",
                                                                                                                  column(12,numericInput("lncRNA_getElites_p_cutoff",label="A numeric cutoff for nominal p value derived from univariate Cox proportional hazards regression",value=0.05,min=0,max=1,step=0.01))
                                                                                                 ),
                                                                                                 conditionalPanel(condition = "input.lncRNA_getElites_survival==2",
                                                                                                                  column(12,radioButtons("lncRNA_getElites_method", h5(strong("Choose a Get Elites method for lncRNA dataset")), inline=TRUE, choices = c("mad"=1,"sd"=2,"pca"=3),selected=1)),
                                                                                                                  conditionalPanel(condition = "input.lncRNA_getElites_method==1",
                                                                                                                                   column(12,radioButtons("lncRNA_getElites_filtration1", h5(strong("Choose a filtration method for lncRNA dataset")), inline=TRUE, choices = c("elite.num"=1,"elite.pct"=2),selected=1)),
                                                                                                                                   conditionalPanel(condition = "input.lncRNA_getElites_filtration1==1",column(12,numericInput("lncRNA_getElites_elite_num1",label="An integer cutoff of exact number for selecting top elites",value=1000,min=0,step=1))),
                                                                                                                                   conditionalPanel(condition = "input.lncRNA_getElites_filtration1==2",column(12,numericInput("lncRNA_getElites_elite_pct1",label="A numeric cutoff of percentage for selecting top elites",value=0.1,min=0,max=1,step=0.05)))
                                                                                                                  ),
                                                                                                                  conditionalPanel(condition = "input.lncRNA_getElites_method==2",
                                                                                                                                   column(12,radioButtons("lncRNA_getElites_filtration2", h5(strong("Choose a filtration method for lncRNA dataset")), inline=TRUE, choices = c("elite.num"=1,"elite.pct"=2),selected=1)),
                                                                                                                                   conditionalPanel(condition = "input.lncRNA_getElites_filtration2==1",column(12,numericInput("lncRNA_getElites_elite_num2",label="An integer cutoff of exact number for selecting top elites",value=1000,min=0,step=1))),
                                                                                                                                   conditionalPanel(condition = "input.lncRNA_getElites_filtration2==2",column(12,numericInput("lncRNA_getElites_elite_pct2",label="A numeric cutoff of percentage for selecting top elites",value=0.1,min=0,max=1,step=0.05)))
                                                                                                                  ),
                                                                                                                  conditionalPanel(condition = "input.lncRNA_getElites_method==3",
                                                                                                                                   column(12,numericInput("lncRNA_getElites_pca_ratio",label="A numeric value which represents the ratio of principal components",value=0.9,min=0,max=1,step=0.05))
                                                                                                                  )
                                                                                                 ),
                                                                                                 column(12,radioButtons("lncRNA_getElites_scaleFlag", h5(strong("Scaling the data after filtering or not")), inline=TRUE, choices = c("Yes"=1,"No"=2),selected=2)),
                                                                                                 column(12,radioButtons("lncRNA_getElites_centerFlag", h5(strong("Centering the data after filtering or not")), inline=TRUE, choices = c("Yes"=1,"No"=2),selected=2))                   
                                                                                             )
                                                                                   )
                                                                  )            
                                                 ),
                                                 
                                                 #Get Elites for DNA methylation dataset
                                                 conditionalPanel(condition="input.multiOmics.includes('3')",
                                                                  conditionalPanel(condition="input.DNA_methylation_getElites==1",
                                                                                   wellPanel(style="background-color:#47AEE9",
                                                                                             fixedRow(
                                                                                                 column(12,tags$h4(strong("Get Elites settings for DNA methylation dataset"),align="center",style="color:#834005")),
                                                                                                 column(12,radioButtons("DNA_methylation_getElites_survival", h5(strong("Use survival information or not")), inline=TRUE, choices = c("Yes"=1,"No"=2),selected=2)),
                                                                                                 column(12,radioButtons("DNA_methylation_getElites_na", h5(strong("NA value action")), inline=TRUE, choices = c("Remove directly"=1,"KNN imputation"=2,"No action"=3),selected=1)),
                                                                                                 column(12,radioButtons("DNA_methylation_getElites_log2", h5(strong("Perform log2 transformation for data before calculating statistics or not")), inline=TRUE, choices = c("Yes"=1,"No"=2),selected=2)),
                                                                                                 column(12,radioButtons("DNA_methylation_getElites_lowpct", h5(strong("Set a numeric cutoff for removing low expression features or not")), inline=TRUE, choices = c("Yes"=1,"No"=2),selected=2)),
                                                                                                 conditionalPanel(condition = "input.DNA_methylation_getElites_lowpct==1",
                                                                                                                  column(12,numericInput("DNA_methylation_getElites_lowpct_value",label="The cutoff for removing low expression features",value=0.1,min=0,max=1,step=0.05))
                                                                                                 ),
                                                                                                 conditionalPanel(condition = "input.DNA_methylation_getElites_survival==1",
                                                                                                                  column(12,numericInput("DNA_methylation_getElites_p_cutoff",label="A numeric cutoff for nominal p value derived from univariate Cox proportional hazards regression",value=0.05,min=0,max=1,step=0.01))
                                                                                                 ),
                                                                                                 conditionalPanel(condition = "input.DNA_methylation_getElites_survival==2",
                                                                                                                  column(12,radioButtons("DNA_methylation_getElites_method", h5(strong("Choose a Get Elites method for DNA methylation dataset")), inline=TRUE, choices = c("mad"=1,"sd"=2,"pca"=3),selected=1)),
                                                                                                                  conditionalPanel(condition = "input.DNA_methylation_getElites_method==1",
                                                                                                                                   column(12,radioButtons("DNA_methylation_getElites_filtration1", h5(strong("Choose a filtration method for DNA methylation dataset")), inline=TRUE, choices = c("elite.num"=1,"elite.pct"=2),selected=1)),
                                                                                                                                   conditionalPanel(condition = "input.DNA_methylation_getElites_filtration1==1",column(12,numericInput("DNA_methylation_getElites_elite_num1",label="An integer cutoff of exact number for selecting top elites",value=1000,min=0,step=1))),
                                                                                                                                   conditionalPanel(condition = "input.DNA_methylation_getElites_filtration1==2",column(12,numericInput("DNA_methylation_getElites_elite_pct1",label="A numeric cutoff of percentage for selecting top elites",value=0.1,min=0,max=1,step=0.05)))
                                                                                                                  ),
                                                                                                                  conditionalPanel(condition = "input.DNA_methylation_getElites_method==2",
                                                                                                                                   column(12,radioButtons("DNA_methylation_getElites_filtration2", h5(strong("Choose a filtration method for DNA methylation dataset")), inline=TRUE, choices = c("elite.num"=1,"elite.pct"=2),selected=1)),
                                                                                                                                   conditionalPanel(condition = "input.DNA_methylation_getElites_filtration2==1",column(12,numericInput("DNA_methylation_getElites_elite_num2",label="An integer cutoff of exact number for selecting top elites",value=1000,min=0,step=1))),
                                                                                                                                   conditionalPanel(condition = "input.DNA_methylation_getElites_filtration2==2",column(12,numericInput("DNA_methylation_getElites_elite_pct2",label="A numeric cutoff of percentage for selecting top elites",value=0.1,min=0,max=1,step=0.05)))
                                                                                                                  ),
                                                                                                                  conditionalPanel(condition = "input.DNA_methylation_getElites_method==3",
                                                                                                                                   column(12,numericInput("DNA_methylation_getElites_pca_ratio",label="A numeric value which represents the ratio of principal components",value=0.9,min=0,max=1,step=0.05))
                                                                                                                  )
                                                                                                 ),
                                                                                                 column(12,radioButtons("DNA_methylation_getElites_scaleFlag", h5(strong("Scaling the data after filtering or not")), inline=TRUE, choices = c("Yes"=1,"No"=2),selected=2)),
                                                                                                 column(12,radioButtons("DNA_methylation_getElites_centerFlag", h5(strong("Centering the data after filtering or not")), inline=TRUE, choices = c("Yes"=1,"No"=2),selected=2))                                     
                                                                                             )
                                                                                   )
                                                                  )
                                                 ),
                                                 
                                                 #Get Elites for copy number alterations dataset
                                                 conditionalPanel(condition="input.multiOmics.includes('4')",
                                                                  conditionalPanel(condition="input.copy_number_alterations_getElites==1",
                                                                                   wellPanel(style="background-color:#47AEE9",
                                                                                             fixedRow(
                                                                                                 column(12,tags$h4(strong("Get Elites settings for copy number alterations dataset"),align="center",style="color:#834005")),
                                                                                                 column(12,radioButtons("copy_number_alterations_getElites_survival", h5(strong("Use survival information or not")), inline=TRUE, choices = c("Yes"=1,"No"=2),selected=2)),
                                                                                                 column(12,radioButtons("copy_number_alterations_getElites_na", h5(strong("NA value action")), inline=TRUE, choices = c("Remove directly"=1,"KNN imputation"=2,"No action"=3),selected=1)),
                                                                                                 column(12,radioButtons("copy_number_alterations_getElites_log2", h5(strong("Perform log2 transformation for data before calculating statistics or not")), inline=TRUE, choices = c("Yes"=1,"No"=2),selected=2)),
                                                                                                 column(12,radioButtons("copy_number_alterations_getElites_lowpct", h5(strong("Set a numeric cutoff for removing low expression features or not")), inline=TRUE, choices = c("Yes"=1,"No"=2),selected=2)),
                                                                                                 conditionalPanel(condition = "input.copy_number_alterations_getElites_lowpct==1",
                                                                                                                  column(12,numericInput("copy_number_alterations_getElites_lowpct_value",label="The cutoff for removing low expression features",value=0.1,min=0,max=1,step=0.05))
                                                                                                 ),
                                                                                                 conditionalPanel(condition = "input.copy_number_alterations_getElites_survival==1",
                                                                                                                  column(12,numericInput("copy_number_alterations_getElites_p_cutoff",label="A numeric cutoff for nominal p value derived from univariate Cox proportional hazards regression",value=0.05,min=0,max=1,step=0.01))
                                                                                                 ),
                                                                                                 conditionalPanel(condition = "input.copy_number_alterations_getElites_survival==2",
                                                                                                                  column(12,radioButtons("copy_number_alterations_getElites_method", h5(strong("Choose a Get Elites method for copy number alterations dataset")), inline=TRUE, choices = c("mad"=1,"sd"=2,"pca"=3),selected=1)),
                                                                                                                  conditionalPanel(condition = "input.copy_number_alterations_getElites_method==1",
                                                                                                                                   column(12,radioButtons("copy_number_alterations_getElites_filtration1", h5(strong("Choose a filtration method for copy number alterations dataset")), inline=TRUE, choices = c("elite.num"=1,"elite.pct"=2),selected=1)),
                                                                                                                                   conditionalPanel(condition = "input.copy_number_alterations_getElites_filtration1==1",column(12,numericInput("copy_number_alterations_getElites_elite_num1",label="An integer cutoff of exact number for selecting top elites",value=1000,min=0,step=1))),
                                                                                                                                   conditionalPanel(condition = "input.copy_number_alterations_getElites_filtration1==2",column(12,numericInput("copy_number_alterations_getElites_elite_pct1",label="A numeric cutoff of percentage for selecting top elites",value=0.1,min=0,max=1,step=0.05)))
                                                                                                                  ),
                                                                                                                  conditionalPanel(condition = "input.copy_number_alterations_getElites_method==2",
                                                                                                                                   column(12,radioButtons("copy_number_alterations_getElites_filtration2", h5(strong("Choose a filtration method for copy number alterations dataset")), inline=TRUE, choices = c("elite.num"=1,"elite.pct"=2),selected=1)),
                                                                                                                                   conditionalPanel(condition = "input.copy_number_alterations_getElites_filtration2==1",column(12,numericInput("copy_number_alterations_getElites_elite_num2",label="An integer cutoff of exact number for selecting top elites",value=1000,min=0,step=1))),
                                                                                                                                   conditionalPanel(condition = "input.copy_number_alterations_getElites_filtration2==2",column(12,numericInput("copy_number_alterations_getElites_elite_pct2",label="A numeric cutoff of percentage for selecting top elites",value=0.1,min=0,max=1,step=0.05)))
                                                                                                                  ),
                                                                                                                  conditionalPanel(condition = "input.copy_number_alterations_getElites_method==3",
                                                                                                                                   column(12,numericInput("copy_number_alterations_getElites_pca_ratio",label="A numeric value which represents the ratio of principal components",value=0.9,min=0,max=1,step=0.05))
                                                                                                                  )
                                                                                                 ),
                                                                                                 column(12,radioButtons("copy_number_alterations_getElites_scaleFlag", h5(strong("Scaling the data after filtering or not")), inline=TRUE, choices = c("Yes"=1,"No"=2),selected=2)),
                                                                                                 column(12,radioButtons("copy_number_alterations_getElites_centerFlag", h5(strong("Centering the data after filtering or not")), inline=TRUE, choices = c("Yes"=1,"No"=2),selected=2))                                                       
                                                                                             )
                                                                                   )
                                                                  )            
                                                 ),
                                                 
                                                 #Get Elites for binary somatic mutation dataset
                                                 conditionalPanel(condition="input.multiOmics.includes('5')",
                                                                  conditionalPanel(condition="input.binary_somatic_mutation_getElites==1",
                                                                                   wellPanel(style="background-color:#47AEE9",
                                                                                             fixedRow(
                                                                                                 column(12,tags$h4(strong("Get Elites settings for binary somatic mutation dataset"),align="center",style="color:#834005")),
                                                                                                 column(12,radioButtons("binary_somatic_mutation_getElites_survival", h5(strong("Use survival information or not")), inline=TRUE, choices = c("Yes"=1,"No"=2),selected=2)),
                                                                                                 column(12,radioButtons("binary_somatic_mutation_getElites_na", h5(strong("NA value action")), inline=TRUE, choices = c("Remove directly"=1,"KNN imputation"=2,"No action"=3),selected=1)),
                                                                                                 column(12,radioButtons("binary_somatic_mutation_getElites_log2", h5(strong("Perform log2 transformation for data before calculating statistics or not")), inline=TRUE, choices = c("Yes"=1,"No"=2),selected=2)),
                                                                                                 # column(12,radioButtons("binary_somatic_mutation_getElites_lowpct", h5(strong("Set a numeric cutoff for removing low expression features or not")), inline=TRUE, choices = c("Yes"=1,"No"=2),selected=2)),
                                                                                                 # conditionalPanel(condition = "input.binary_somatic_mutation_getElites_lowpct==1",
                                                                                                 #                  column(12,numericInput("binary_somatic_mutation_getElites_lowpct_value",label="The cutoff for removing low expression features",value=0.1,min=0,max=1,step=0.05))
                                                                                                 # ),
                                                                                                 conditionalPanel(condition = "input.binary_somatic_mutation_getElites_survival==1",
                                                                                                                  column(12,numericInput("binary_somatic_mutation_getElites_p_cutoff",label="A numeric cutoff for nominal p value derived from univariate Cox proportional hazards regression",value=0.05,min=0,max=1,step=0.01))
                                                                                                 ),
                                                                                                 conditionalPanel(condition = "input.binary_somatic_mutation_getElites_survival==2",
                                                                                                                  column(12,radioButtons("binary_somatic_mutation_getElites_filtration", h5(strong("Choose a filtration method for binary somatic mutation dataset")), inline=TRUE, choices = c("elite.num"=1,"elite.pct"=2),selected=2)),
                                                                                                                  conditionalPanel(condition = "input.binary_somatic_mutation_getElites_filtration==1",column(12,numericInput("binary_somatic_mutation_getElites_elite_num",label="An integer cutoff of mutation frequency for selecting elites",value=100,min=0,step=1))),
                                                                                                                  conditionalPanel(condition = "input.binary_somatic_mutation_getElites_filtration==2",column(12,numericInput("binary_somatic_mutation_getElites_elite_pct",label="A numeric cutoff of 'mutation / sample' frequency for selecting elites",value=0.1,min=0,max=1,step=0.05)))                 
                                                                                                 ),
                                                                                                 column(12,radioButtons("binary_somatic_mutation_getElites_scaleFlag", h5(strong("Scaling the data after filtering or not")), inline=TRUE, choices = c("Yes"=1,"No"=2),selected=2)),
                                                                                                 column(12,radioButtons("binary_somatic_mutation_getElites_centerFlag", h5(strong("Centering the data after filtering or not")), inline=TRUE, choices = c("Yes"=1,"No"=2),selected=2))                                                                         
                                                                                             )
                                                                                   )
                                                                  )
                                                 ),
                                                 
                                                 #Get Elites for radiomics dataset
                                                 conditionalPanel(condition="input.multiOmics.includes('6')",
                                                                  conditionalPanel(condition="input.radiomics_getElites==1",
                                                                                   wellPanel(style="background-color:#47AEE9",
                                                                                             fixedRow(
                                                                                                 column(12,tags$h4(strong("Get Elites settings for radiomics dataset"),align="center",style="color:#834005")),
                                                                                                 column(12,radioButtons("radiomics_getElites_survival", h5(strong("Use survival information or not")), inline=TRUE, choices = c("Yes"=1,"No"=2),selected=2)),
                                                                                                 column(12,radioButtons("radiomics_getElites_na", h5(strong("NA value action")), inline=TRUE, choices = c("Remove directly"=1,"KNN imputation"=2,"No action"=3),selected=1)),
                                                                                                 column(12,radioButtons("radiomics_getElites_log2", h5(strong("Perform log2 transformation for data before calculating statistics or not")), inline=TRUE, choices = c("Yes"=1,"No"=2),selected=2)),
                                                                                                 column(12,radioButtons("radiomics_getElites_lowpct", h5(strong("Set a numeric cutoff for removing low expression features or not")), inline=TRUE, choices = c("Yes"=1,"No"=2),selected=2)),
                                                                                                 conditionalPanel(condition = "input.radiomics_getElites_lowpct==1",
                                                                                                                  column(12,numericInput("radiomics_getElites_lowpct_value",label="The cutoff for removing low expression features",value=0.1,min=0,max=1,step=0.05))
                                                                                                 ),
                                                                                                 conditionalPanel(condition = "input.radiomics_getElites_survival==1",
                                                                                                                  column(12,numericInput("radiomics_getElites_p_cutoff",label="A numeric cutoff for nominal p value derived from univariate Cox proportional hazards regression",value=0.05,min=0,max=1,step=0.01))
                                                                                                 ),
                                                                                                 conditionalPanel(condition = "input.radiomics_getElites_survival==2",
                                                                                                                  column(12,radioButtons("radiomics_getElites_method", h5(strong("Choose a Get Elites method for radiomics dataset")), inline=TRUE, choices = c("mad"=1,"sd"=2,"pca"=3),selected=1)),
                                                                                                                  conditionalPanel(condition = "input.radiomics_getElites_method==1",
                                                                                                                                   column(12,radioButtons("radiomics_getElites_filtration1", h5(strong("Choose a filtration method for radiomics dataset")), inline=TRUE, choices = c("elite.num"=1,"elite.pct"=2),selected=1)),
                                                                                                                                   conditionalPanel(condition = "input.radiomics_getElites_filtration1==1",column(12,numericInput("radiomics_getElites_elite_num1",label="An integer cutoff of exact number for selecting top elites",value=1000,min=0,step=1))),
                                                                                                                                   conditionalPanel(condition = "input.radiomics_getElites_filtration1==2",column(12,numericInput("radiomics_getElites_elite_pct1",label="A numeric cutoff of percentage for selecting top elites",value=0.1,min=0,max=1,step=0.05)))
                                                                                                                  ),
                                                                                                                  conditionalPanel(condition = "input.radiomics_getElites_method==2",
                                                                                                                                   column(12,radioButtons("radiomics_getElites_filtration2", h5(strong("Choose a filtration method for radiomics dataset")), inline=TRUE, choices = c("elite.num"=1,"elite.pct"=2),selected=1)),
                                                                                                                                   conditionalPanel(condition = "input.radiomics_getElites_filtration2==1",column(12,numericInput("radiomics_getElites_elite_num2",label="An integer cutoff of exact number for selecting top elites",value=1000,min=0,step=1))),
                                                                                                                                   conditionalPanel(condition = "input.radiomics_getElites_filtration2==2",column(12,numericInput("radiomics_getElites_elite_pct2",label="A numeric cutoff of percentage for selecting top elites",value=0.1,min=0,max=1,step=0.05)))
                                                                                                                  ),
                                                                                                                  conditionalPanel(condition = "input.radiomics_getElites_method==3",
                                                                                                                                   column(12,numericInput("radiomics_getElites_pca_ratio",label="A numeric value which represents the ratio of principal components",value=0.9,min=0,max=1,step=0.05))
                                                                                                                  )
                                                                                                 ),
                                                                                                 column(12,radioButtons("radiomics_getElites_scaleFlag", h5(strong("Scaling the data after filtering or not")), inline=TRUE, choices = c("Yes"=1,"No"=2),selected=2)),
                                                                                                 column(12,radioButtons("radiomics_getElites_centerFlag", h5(strong("Centering the data after filtering or not")), inline=TRUE, choices = c("Yes"=1,"No"=2),selected=2))
                                                                                             )
                                                                                   )
                                                                  )           
                                                 )
                                ),
                                
                                #Validation
                                conditionalPanel(condition="input.tcga_vali_getElites==2",
                                             
                                                 #Get Elites for mRNA dataset
                                                 conditionalPanel(condition="input.validation_multiOmics.includes('1')",
                                                                  conditionalPanel(condition="input.validation_mRNA_getElites==1",
                                                                                   wellPanel(style="background-color:#47AEE9",
                                                                                             fixedRow(
                                                                                                 column(12,tags$h4(strong("Get Elites settings for mRNA dataset"),align="center",style="color:#834005")),
                                                                                                 column(12,radioButtons("validation_mRNA_getElites_survival", h5(strong("Use survival information or not")), inline=TRUE, choices = c("Yes"=1,"No"=2),selected=2)),
                                                                                                 column(12,radioButtons("validation_mRNA_getElites_na", h5(strong("NA value action")), inline=TRUE, choices = c("Remove directly"=1,"KNN imputation"=2,"No action"=3),selected=1)),
                                                                                                 column(12,radioButtons("validation_mRNA_getElites_log2", h5(strong("Perform log2 transformation for data before calculating statistics or not")), inline=TRUE, choices = c("Yes"=1,"No"=2),selected=2)),
                                                                                                 column(12,radioButtons("validation_mRNA_getElites_lowpct", h5(strong("Set a numeric cutoff for removing low expression features or not")), inline=TRUE, choices = c("Yes"=1,"No"=2),selected=2)),
                                                                                                 conditionalPanel(condition = "input.validation_mRNA_getElites_lowpct==1",
                                                                                                                  column(12,numericInput("validation_mRNA_getElites_lowpct_value",label="The cutoff for removing low expression features",value=0.1,min=0,max=1,step=0.05))
                                                                                                 ),
                                                                                                 conditionalPanel(condition = "input.validation_mRNA_getElites_survival==1",
                                                                                                                  column(12,numericInput("validation_mRNA_getElites_p_cutoff",label="A numeric cutoff for nominal p value derived from univariate Cox proportional hazards regression",value=0.05,min=0,max=1,step=0.01))
                                                                                                 ),
                                                                                                 conditionalPanel(condition = "input.validation_mRNA_getElites_survival==2",
                                                                                                                  column(12,radioButtons("validation_mRNA_getElites_method", h5(strong("Choose a Get Elites method for mRNA dataset")), inline=TRUE, choices = c("mad"=1,"sd"=2,"pca"=3),selected=1)),
                                                                                                                  conditionalPanel(condition = "input.validation_mRNA_getElites_method==1",
                                                                                                                                   column(12,radioButtons("validation_mRNA_getElites_filtration1", h5(strong("Choose a filtration method for mRNA dataset")), inline=TRUE, choices = c("elite.num"=1,"elite.pct"=2),selected=1)),
                                                                                                                                   conditionalPanel(condition = "input.validation_mRNA_getElites_filtration1==1",column(12,numericInput("validation_mRNA_getElites_elite_num1",label="An integer cutoff of exact number for selecting top elites",value=1000,min=0,step=1))),
                                                                                                                                   conditionalPanel(condition = "input.validation_mRNA_getElites_filtration1==2",column(12,numericInput("validation_mRNA_getElites_elite_pct1",label="A numeric cutoff of percentage for selecting top elites",value=0.1,min=0,max=1,step=0.05)))
                                                                                                                  ),
                                                                                                                  conditionalPanel(condition = "input.validation_mRNA_getElites_method==2",
                                                                                                                                   column(12,radioButtons("validation_mRNA_getElites_filtration2", h5(strong("Choose a filtration method for mRNA dataset")), inline=TRUE, choices = c("elite.num"=1,"elite.pct"=2),selected=1)),
                                                                                                                                   conditionalPanel(condition = "input.validation_mRNA_getElites_filtration2==1",column(12,numericInput("validation_mRNA_getElites_elite_num2",label="An integer cutoff of exact number for selecting top elites",value=1000,min=0,step=1))),
                                                                                                                                   conditionalPanel(condition = "input.validation_mRNA_getElites_filtration2==2",column(12,numericInput("validation_mRNA_getElites_elite_pct2",label="A numeric cutoff of percentage for selecting top elites",value=0.1,min=0,max=1,step=0.05)))
                                                                                                                  ),
                                                                                                                  conditionalPanel(condition = "input.validation_mRNA_getElites_method==3",
                                                                                                                                   column(12,numericInput("validation_mRNA_getElites_pca_ratio",label="A numeric value which represents the ratio of principal components",value=0.9,min=0,max=1,step=0.05))
                                                                                                                  )
                                                                                                 ),
                                                                                                 column(12,radioButtons("validation_mRNA_getElites_scaleFlag", h5(strong("Scaling the data after filtering or not")), inline=TRUE, choices = c("Yes"=1,"No"=2),selected=2)),
                                                                                                 column(12,radioButtons("validation_mRNA_getElites_centerFlag", h5(strong("Centering the data after filtering or not")), inline=TRUE, choices = c("Yes"=1,"No"=2),selected=2))
                                                                                             )
                                                                                   )
                                                                  )           
                                                 ),
                                                 
                                                 #Get Elites for lncRNA dataset
                                                 conditionalPanel(condition="input.validation_multiOmics.includes('2')",
                                                                  conditionalPanel(condition="input.validation_lncRNA_getElites==1",
                                                                                   wellPanel(style="background-color:#47AEE9",
                                                                                             fixedRow(
                                                                                                 column(12,tags$h4(strong("Get Elites settings for lncRNA dataset"),align="center",style="color:#834005")),
                                                                                                 column(12,radioButtons("validation_lncRNA_getElites_survival", h5(strong("Use survival information or not")), inline=TRUE, choices = c("Yes"=1,"No"=2),selected=2)),
                                                                                                 column(12,radioButtons("validation_lncRNA_getElites_na", h5(strong("NA value action")), inline=TRUE, choices = c("Remove directly"=1,"KNN imputation"=2,"No action"=3),selected=1)),
                                                                                                 column(12,radioButtons("validation_lncRNA_getElites_log2", h5(strong("Perform log2 transformation for data before calculating statistics or not")), inline=TRUE, choices = c("Yes"=1,"No"=2),selected=2)),
                                                                                                 column(12,radioButtons("validation_lncRNA_getElites_lowpct", h5(strong("Set a numeric cutoff for removing low expression features or not")), inline=TRUE, choices = c("Yes"=1,"No"=2),selected=2)),
                                                                                                 conditionalPanel(condition = "input.validation_lncRNA_getElites_lowpct==1",
                                                                                                                  column(12,numericInput("validation_lncRNA_getElites_lowpct_value",label="The cutoff for removing low expression features",value=0.1,min=0,max=1,step=0.05))
                                                                                                 ),
                                                                                                 conditionalPanel(condition = "input.validation_lncRNA_getElites_survival==1",
                                                                                                                  column(12,numericInput("validation_lncRNA_getElites_p_cutoff",label="A numeric cutoff for nominal p value derived from univariate Cox proportional hazards regression",value=0.05,min=0,max=1,step=0.01))
                                                                                                 ),
                                                                                                 conditionalPanel(condition = "input.validation_lncRNA_getElites_survival==2",
                                                                                                                  column(12,radioButtons("validation_lncRNA_getElites_method", h5(strong("Choose a Get Elites method for lncRNA dataset")), inline=TRUE, choices = c("mad"=1,"sd"=2,"pca"=3),selected=1)),
                                                                                                                  conditionalPanel(condition = "input.validation_lncRNA_getElites_method==1",
                                                                                                                                   column(12,radioButtons("validation_lncRNA_getElites_filtration1", h5(strong("Choose a filtration method for lncRNA dataset")), inline=TRUE, choices = c("elite.num"=1,"elite.pct"=2),selected=1)),
                                                                                                                                   conditionalPanel(condition = "input.validation_lncRNA_getElites_filtration1==1",column(12,numericInput("validation_lncRNA_getElites_elite_num1",label="An integer cutoff of exact number for selecting top elites",value=1000,min=0,step=1))),
                                                                                                                                   conditionalPanel(condition = "input.validation_lncRNA_getElites_filtration1==2",column(12,numericInput("validation_lncRNA_getElites_elite_pct1",label="A numeric cutoff of percentage for selecting top elites",value=0.1,min=0,max=1,step=0.05)))
                                                                                                                  ),
                                                                                                                  conditionalPanel(condition = "input.validation_lncRNA_getElites_method==2",
                                                                                                                                   column(12,radioButtons("validation_lncRNA_getElites_filtration2", h5(strong("Choose a filtration method for lncRNA dataset")), inline=TRUE, choices = c("elite.num"=1,"elite.pct"=2),selected=1)),
                                                                                                                                   conditionalPanel(condition = "input.validation_lncRNA_getElites_filtration2==1",column(12,numericInput("validation_lncRNA_getElites_elite_num2",label="An integer cutoff of exact number for selecting top elites",value=1000,min=0,step=1))),
                                                                                                                                   conditionalPanel(condition = "input.validation_lncRNA_getElites_filtration2==2",column(12,numericInput("validation_lncRNA_getElites_elite_pct2",label="A numeric cutoff of percentage for selecting top elites",value=0.1,min=0,max=1,step=0.05)))
                                                                                                                  ),
                                                                                                                  conditionalPanel(condition = "input.validation_lncRNA_getElites_method==3",
                                                                                                                                   column(12,numericInput("validation_lncRNA_getElites_pca_ratio",label="A numeric value which represents the ratio of principal components",value=0.9,min=0,max=1,step=0.05))
                                                                                                                  )
                                                                                                 ),
                                                                                                 column(12,radioButtons("validation_lncRNA_getElites_scaleFlag", h5(strong("Scaling the data after filtering or not")), inline=TRUE, choices = c("Yes"=1,"No"=2),selected=2)),
                                                                                                 column(12,radioButtons("validation_lncRNA_getElites_centerFlag", h5(strong("Centering the data after filtering or not")), inline=TRUE, choices = c("Yes"=1,"No"=2),selected=2))                   
                                                                                             )
                                                                                   )
                                                                  )            
                                                 ),
                                                 
                                                 #Get Elites for DNA methylation dataset
                                                 conditionalPanel(condition="input.validation_multiOmics.includes('3')",
                                                                  conditionalPanel(condition="input.validation_DNA_methylation_getElites==1",
                                                                                   wellPanel(style="background-color:#47AEE9",
                                                                                             fixedRow(
                                                                                                 column(12,tags$h4(strong("Get Elites settings for DNA methylation dataset"),align="center",style="color:#834005")),
                                                                                                 column(12,radioButtons("validation_DNA_methylation_getElites_survival", h5(strong("Use survival information or not")), inline=TRUE, choices = c("Yes"=1,"No"=2),selected=2)),
                                                                                                 column(12,radioButtons("validation_DNA_methylation_getElites_na", h5(strong("NA value action")), inline=TRUE, choices = c("Remove directly"=1,"KNN imputation"=2,"No action"=3),selected=1)),
                                                                                                 column(12,radioButtons("validation_DNA_methylation_getElites_log2", h5(strong("Perform log2 transformation for data before calculating statistics or not")), inline=TRUE, choices = c("Yes"=1,"No"=2),selected=2)),
                                                                                                 column(12,radioButtons("validation_DNA_methylation_getElites_lowpct", h5(strong("Set a numeric cutoff for removing low expression features or not")), inline=TRUE, choices = c("Yes"=1,"No"=2),selected=2)),
                                                                                                 conditionalPanel(condition = "input.validation_DNA_methylation_getElites_lowpct==1",
                                                                                                                  column(12,numericInput("validation_DNA_methylation_getElites_lowpct_value",label="The cutoff for removing low expression features",value=0.1,min=0,max=1,step=0.05))
                                                                                                 ),
                                                                                                 conditionalPanel(condition = "input.validation_DNA_methylation_getElites_survival==1",
                                                                                                                  column(12,numericInput("validation_DNA_methylation_getElites_p_cutoff",label="A numeric cutoff for nominal p value derived from univariate Cox proportional hazards regression",value=0.05,min=0,max=1,step=0.01))
                                                                                                 ),
                                                                                                 conditionalPanel(condition = "input.validation_DNA_methylation_getElites_survival==2",
                                                                                                                  column(12,radioButtons("validation_DNA_methylation_getElites_method", h5(strong("Choose a Get Elites method for DNA methylation dataset")), inline=TRUE, choices = c("mad"=1,"sd"=2,"pca"=3),selected=1)),
                                                                                                                  conditionalPanel(condition = "input.validation_DNA_methylation_getElites_method==1",
                                                                                                                                   column(12,radioButtons("validation_DNA_methylation_getElites_filtration1", h5(strong("Choose a filtration method for DNA methylation dataset")), inline=TRUE, choices = c("elite.num"=1,"elite.pct"=2),selected=1)),
                                                                                                                                   conditionalPanel(condition = "input.validation_DNA_methylation_getElites_filtration1==1",column(12,numericInput("validation_DNA_methylation_getElites_elite_num1",label="An integer cutoff of exact number for selecting top elites",value=1000,min=0,step=1))),
                                                                                                                                   conditionalPanel(condition = "input.validation_DNA_methylation_getElites_filtration1==2",column(12,numericInput("validation_DNA_methylation_getElites_elite_pct1",label="A numeric cutoff of percentage for selecting top elites",value=0.1,min=0,max=1,step=0.05)))
                                                                                                                  ),
                                                                                                                  conditionalPanel(condition = "input.validation_DNA_methylation_getElites_method==2",
                                                                                                                                   column(12,radioButtons("validation_DNA_methylation_getElites_filtration2", h5(strong("Choose a filtration method for DNA methylation dataset")), inline=TRUE, choices = c("elite.num"=1,"elite.pct"=2),selected=1)),
                                                                                                                                   conditionalPanel(condition = "input.validation_DNA_methylation_getElites_filtration2==1",column(12,numericInput("validation_DNA_methylation_getElites_elite_num2",label="An integer cutoff of exact number for selecting top elites",value=1000,min=0,step=1))),
                                                                                                                                   conditionalPanel(condition = "input.validation_DNA_methylation_getElites_filtration2==2",column(12,numericInput("validation_DNA_methylation_getElites_elite_pct2",label="A numeric cutoff of percentage for selecting top elites",value=0.1,min=0,max=1,step=0.05)))
                                                                                                                  ),
                                                                                                                  conditionalPanel(condition = "input.validation_DNA_methylation_getElites_method==3",
                                                                                                                                   column(12,numericInput("validation_DNA_methylation_getElites_pca_ratio",label="A numeric value which represents the ratio of principal components",value=0.9,min=0,max=1,step=0.05))
                                                                                                                  )
                                                                                                 ),
                                                                                                 column(12,radioButtons("validation_DNA_methylation_getElites_scaleFlag", h5(strong("Scaling the data after filtering or not")), inline=TRUE, choices = c("Yes"=1,"No"=2),selected=2)),
                                                                                                 column(12,radioButtons("validation_DNA_methylation_getElites_centerFlag", h5(strong("Centering the data after filtering or not")), inline=TRUE, choices = c("Yes"=1,"No"=2),selected=2))                                     
                                                                                             )
                                                                                   )
                                                                  )
                                                 ),
                                                 
                                                 #Get Elites for copy number alterations dataset
                                                 conditionalPanel(condition="input.validation_multiOmics.includes('4')",
                                                                  conditionalPanel(condition="input.validation_copy_number_alterations_getElites==1",
                                                                                   wellPanel(style="background-color:#47AEE9",
                                                                                             fixedRow(
                                                                                                 column(12,tags$h4(strong("Get Elites settings for copy number alterations dataset"),align="center",style="color:#834005")),
                                                                                                 column(12,radioButtons("validation_copy_number_alterations_getElites_survival", h5(strong("Use survival information or not")), inline=TRUE, choices = c("Yes"=1,"No"=2),selected=2)),
                                                                                                 column(12,radioButtons("validation_copy_number_alterations_getElites_na", h5(strong("NA value action")), inline=TRUE, choices = c("Remove directly"=1,"KNN imputation"=2,"No action"=3),selected=1)),
                                                                                                 column(12,radioButtons("validation_copy_number_alterations_getElites_log2", h5(strong("Perform log2 transformation for data before calculating statistics or not")), inline=TRUE, choices = c("Yes"=1,"No"=2),selected=2)),
                                                                                                 column(12,radioButtons("validation_copy_number_alterations_getElites_lowpct", h5(strong("Set a numeric cutoff for removing low expression features or not")), inline=TRUE, choices = c("Yes"=1,"No"=2),selected=2)),
                                                                                                 conditionalPanel(condition = "input.validation_copy_number_alterations_getElites_lowpct==1",
                                                                                                                  column(12,numericInput("validation_copy_number_alterations_getElites_lowpct_value",label="The cutoff for removing low expression features",value=0.1,min=0,max=1,step=0.05))
                                                                                                 ),
                                                                                                 conditionalPanel(condition = "input.validation_copy_number_alterations_getElites_survival==1",
                                                                                                                  column(12,numericInput("validation_copy_number_alterations_getElites_p_cutoff",label="A numeric cutoff for nominal p value derived from univariate Cox proportional hazards regression",value=0.05,min=0,max=1,step=0.01))
                                                                                                 ),
                                                                                                 conditionalPanel(condition = "input.validation_copy_number_alterations_getElites_survival==2",
                                                                                                                  column(12,radioButtons("validation_copy_number_alterations_getElites_method", h5(strong("Choose a Get Elites method for copy number alterations dataset")), inline=TRUE, choices = c("mad"=1,"sd"=2,"pca"=3),selected=1)),
                                                                                                                  conditionalPanel(condition = "input.validation_copy_number_alterations_getElites_method==1",
                                                                                                                                   column(12,radioButtons("validation_copy_number_alterations_getElites_filtration1", h5(strong("Choose a filtration method for copy number alterations dataset")), inline=TRUE, choices = c("elite.num"=1,"elite.pct"=2),selected=1)),
                                                                                                                                   conditionalPanel(condition = "input.validation_copy_number_alterations_getElites_filtration1==1",column(12,numericInput("validation_copy_number_alterations_getElites_elite_num1",label="An integer cutoff of exact number for selecting top elites",value=1000,min=0,step=1))),
                                                                                                                                   conditionalPanel(condition = "input.validation_copy_number_alterations_getElites_filtration1==2",column(12,numericInput("validation_copy_number_alterations_getElites_elite_pct1",label="A numeric cutoff of percentage for selecting top elites",value=0.1,min=0,max=1,step=0.05)))
                                                                                                                  ),
                                                                                                                  conditionalPanel(condition = "input.validation_copy_number_alterations_getElites_method==2",
                                                                                                                                   column(12,radioButtons("validation_copy_number_alterations_getElites_filtration2", h5(strong("Choose a filtration method for copy number alterations dataset")), inline=TRUE, choices = c("elite.num"=1,"elite.pct"=2),selected=1)),
                                                                                                                                   conditionalPanel(condition = "input.validation_copy_number_alterations_getElites_filtration2==1",column(12,numericInput("validation_copy_number_alterations_getElites_elite_num2",label="An integer cutoff of exact number for selecting top elites",value=1000,min=0,step=1))),
                                                                                                                                   conditionalPanel(condition = "input.validation_copy_number_alterations_getElites_filtration2==2",column(12,numericInput("validation_copy_number_alterations_getElites_elite_pct2",label="A numeric cutoff of percentage for selecting top elites",value=0.1,min=0,max=1,step=0.05)))
                                                                                                                  ),
                                                                                                                  conditionalPanel(condition = "input.validation_copy_number_alterations_getElites_method==3",
                                                                                                                                   column(12,numericInput("validation_copy_number_alterations_getElites_pca_ratio",label="A numeric value which represents the ratio of principal components",value=0.9,min=0,max=1,step=0.05))
                                                                                                                  )
                                                                                                 ),
                                                                                                 column(12,radioButtons("validation_copy_number_alterations_getElites_scaleFlag", h5(strong("Scaling the data after filtering or not")), inline=TRUE, choices = c("Yes"=1,"No"=2),selected=2)),
                                                                                                 column(12,radioButtons("validation_copy_number_alterations_getElites_centerFlag", h5(strong("Centering the data after filtering or not")), inline=TRUE, choices = c("Yes"=1,"No"=2),selected=2))                                                       
                                                                                             )
                                                                                   )
                                                                  )            
                                                 ),
                                                 
                                                 #Get Elites for binary somatic mutation dataset
                                                 conditionalPanel(condition="input.validation_multiOmics.includes('5')",
                                                                  conditionalPanel(condition="input.validation_binary_somatic_mutation_getElites==1",
                                                                                   wellPanel(style="background-color:#47AEE9",
                                                                                             fixedRow(
                                                                                                 column(12,tags$h4(strong("Get Elites settings for binary somatic mutation dataset"),align="center",style="color:#834005")),
                                                                                                 column(12,radioButtons("validation_binary_somatic_mutation_getElites_survival", h5(strong("Use survival information or not")), inline=TRUE, choices = c("Yes"=1,"No"=2),selected=2)),
                                                                                                 column(12,radioButtons("validation_binary_somatic_mutation_getElites_na", h5(strong("NA value action")), inline=TRUE, choices = c("Remove directly"=1,"KNN imputation"=2,"No action"=3),selected=1)),
                                                                                                 column(12,radioButtons("validation_binary_somatic_mutation_getElites_log2", h5(strong("Perform log2 transformation for data before calculating statistics or not")), inline=TRUE, choices = c("Yes"=1,"No"=2),selected=2)),
                                                                                                 # column(12,radioButtons("validation_binary_somatic_mutation_getElites_lowpct", h5(strong("Set a numeric cutoff for removing low expression features or not")), inline=TRUE, choices = c("Yes"=1,"No"=2),selected=2)),
                                                                                                 # conditionalPanel(condition = "input.validation_binary_somatic_mutation_getElites_lowpct==1",
                                                                                                 #                  column(12,numericInput("validation_binary_somatic_mutation_getElites_lowpct_value",label="The cutoff for removing low expression features",value=0.1,min=0,max=1,step=0.05))
                                                                                                 # ),
                                                                                                 conditionalPanel(condition = "input.validation_binary_somatic_mutation_getElites_survival==1",
                                                                                                                  column(12,numericInput("validation_binary_somatic_mutation_getElites_p_cutoff",label="A numeric cutoff for nominal p value derived from univariate Cox proportional hazards regression",value=0.05,min=0,max=1,step=0.01))
                                                                                                 ),
                                                                                                 conditionalPanel(condition = "input.validation_binary_somatic_mutation_getElites_survival==2",
                                                                                                                  column(12,radioButtons("validation_binary_somatic_mutation_getElites_filtration", h5(strong("Choose a filtration method for binary somatic mutation dataset")), inline=TRUE, choices = c("elite.num"=1,"elite.pct"=2),selected=2)),
                                                                                                                  conditionalPanel(condition = "input.validation_binary_somatic_mutation_getElites_filtration==1",column(12,numericInput("validation_binary_somatic_mutation_getElites_elite_num",label="An integer cutoff of mutation frequency for selecting elites",value=100,min=0,step=1))),
                                                                                                                  conditionalPanel(condition = "input.validation_binary_somatic_mutation_getElites_filtration==2",column(12,numericInput("validation_binary_somatic_mutation_getElites_elite_pct",label="A numeric cutoff of 'mutation / sample' frequency for selecting elites",value=0.1,min=0,max=1,step=0.05)))                 
                                                                                                 ),
                                                                                                 column(12,radioButtons("validation_binary_somatic_mutation_getElites_scaleFlag", h5(strong("Scaling the data after filtering or not")), inline=TRUE, choices = c("Yes"=1,"No"=2),selected=2)),
                                                                                                 column(12,radioButtons("validation_binary_somatic_mutation_getElites_centerFlag", h5(strong("Centering the data after filtering or not")), inline=TRUE, choices = c("Yes"=1,"No"=2),selected=2))                                                                         
                                                                                             )
                                                                                   )
                                                                  )
                                                 ),
                                                 
                                                 #Get Elites for radiomics dataset
                                                 conditionalPanel(condition="input.validation_multiOmics.includes('6')",
                                                                  conditionalPanel(condition="input.validation_radiomics_getElites==1",
                                                                                   wellPanel(style="background-color:#47AEE9",
                                                                                             fixedRow(
                                                                                                 column(12,tags$h4(strong("Get Elites settings for radiomics dataset"),align="center",style="color:#834005")),
                                                                                                 column(12,radioButtons("validation_radiomics_getElites_survival", h5(strong("Use survival information or not")), inline=TRUE, choices = c("Yes"=1,"No"=2),selected=2)),
                                                                                                 column(12,radioButtons("validation_radiomics_getElites_na", h5(strong("NA value action")), inline=TRUE, choices = c("Remove directly"=1,"KNN imputation"=2,"No action"=3),selected=1)),
                                                                                                 column(12,radioButtons("validation_radiomics_getElites_log2", h5(strong("Perform log2 transformation for data before calculating statistics or not")), inline=TRUE, choices = c("Yes"=1,"No"=2),selected=2)),
                                                                                                 column(12,radioButtons("validation_radiomics_getElites_lowpct", h5(strong("Set a numeric cutoff for removing low expression features or not")), inline=TRUE, choices = c("Yes"=1,"No"=2),selected=2)),
                                                                                                 conditionalPanel(condition = "input.validation_radiomics_getElites_lowpct==1",
                                                                                                                  column(12,numericInput("validation_radiomics_getElites_lowpct_value",label="The cutoff for removing low expression features",value=0.1,min=0,max=1,step=0.05))
                                                                                                 ),
                                                                                                 conditionalPanel(condition = "input.validation_radiomics_getElites_survival==1",
                                                                                                                  column(12,numericInput("validation_radiomics_getElites_p_cutoff",label="A numeric cutoff for nominal p value derived from univariate Cox proportional hazards regression",value=0.05,min=0,max=1,step=0.01))
                                                                                                 ),
                                                                                                 conditionalPanel(condition = "input.validation_radiomics_getElites_survival==2",
                                                                                                                  column(12,radioButtons("validation_radiomics_getElites_method", h5(strong("Choose a Get Elites method for radiomics dataset")), inline=TRUE, choices = c("mad"=1,"sd"=2,"pca"=3),selected=1)),
                                                                                                                  conditionalPanel(condition = "input.validation_radiomics_getElites_method==1",
                                                                                                                                   column(12,radioButtons("validation_radiomics_getElites_filtration1", h5(strong("Choose a filtration method for radiomics dataset")), inline=TRUE, choices = c("elite.num"=1,"elite.pct"=2),selected=1)),
                                                                                                                                   conditionalPanel(condition = "input.validation_radiomics_getElites_filtration1==1",column(12,numericInput("validation_radiomics_getElites_elite_num1",label="An integer cutoff of exact number for selecting top elites",value=1000,min=0,step=1))),
                                                                                                                                   conditionalPanel(condition = "input.validation_radiomics_getElites_filtration1==2",column(12,numericInput("validation_radiomics_getElites_elite_pct1",label="A numeric cutoff of percentage for selecting top elites",value=0.1,min=0,max=1,step=0.05)))
                                                                                                                  ),
                                                                                                                  conditionalPanel(condition = "input.validation_radiomics_getElites_method==2",
                                                                                                                                   column(12,radioButtons("validation_radiomics_getElites_filtration2", h5(strong("Choose a filtration method for radiomics dataset")), inline=TRUE, choices = c("elite.num"=1,"elite.pct"=2),selected=1)),
                                                                                                                                   conditionalPanel(condition = "input.validation_radiomics_getElites_filtration2==1",column(12,numericInput("validation_radiomics_getElites_elite_num2",label="An integer cutoff of exact number for selecting top elites",value=1000,min=0,step=1))),
                                                                                                                                   conditionalPanel(condition = "input.validation_radiomics_getElites_filtration2==2",column(12,numericInput("validation_radiomics_getElites_elite_pct2",label="A numeric cutoff of percentage for selecting top elites",value=0.1,min=0,max=1,step=0.05)))
                                                                                                                  ),
                                                                                                                  conditionalPanel(condition = "input.validation_radiomics_getElites_method==3",
                                                                                                                                   column(12,numericInput("validation_radiomics_getElites_pca_ratio",label="A numeric value which represents the ratio of principal components",value=0.9,min=0,max=1,step=0.05))
                                                                                                                  )
                                                                                                 ),
                                                                                                 column(12,radioButtons("validation_radiomics_getElites_scaleFlag", h5(strong("Scaling the data after filtering or not")), inline=TRUE, choices = c("Yes"=1,"No"=2),selected=2)),
                                                                                                 column(12,radioButtons("validation_radiomics_getElites_centerFlag", h5(strong("Centering the data after filtering or not")), inline=TRUE, choices = c("Yes"=1,"No"=2),selected=2))
                                                                                             )
                                                                                   )
                                                                  )           
                                                 )
                                ),
                          
                                #Process Get Elites and integrate datasets
                                wellPanel(style="background-color:#47AEE9;",
                                          
                                          fixedRow(
                                              column(12,tags$h4(strong("Process Get Elites"),align="center",style="color:#834005")),
                                              column(12,offset=0,tags$p(h5("Now click the 'Process' button below to process 'Get Elites' based on the settings above, and then integrate datasets for following steps.",align="left"))),
                                              column(12,actionButton("Get_Elites_process","Process",width="100%",class="btn btn-primary"))
                                          )
                                )
                           ),
                           
                           #2.Get Clustering Number
                           conditionalPanel(condition="input.getModuleStep==2",
                                wellPanel(style="background-color:#47AEE9",
                                    fixedRow(
                                        column(12,tags$h4(strong("Get Clustering Number"),align="center",style="color:#834005")),
                                        column(12,offset=0,tags$p(h5(strong("This step aims to search the optimal number for multi-omics integrative clustering determined by 'clustering prediction index' (CPI) and 'Gap-statistics' (Gapk).",align="left")))),
                                        column(12,offset=0,tags$p(h5(strong("Please indicate the range of clustering number first:",align="left")))),
                                        column(12,numericInput("getClusteringNumber_min",label="Minimum",value=2,min=1,step=1)),
                                        column(12,numericInput("getClusteringNumber_max",label="Maximum",value=8,step=1)),
                                        column(12,radioButtons("getClusteringNumber_center", h5(strong("Centering the data or not")), inline=TRUE, choices = c("Yes"=1,"No"=2),selected=1)),
                                        column(12,radioButtons("getClusteringNumber_scale", h5(strong("Scaling the data or not")), inline=TRUE, choices = c("Yes"=1,"No"=2),selected=1)),
                                        column(12,textAreaInput(inputId = "getClusteringNumber_FigureName",
                                                                label = "Figure Name",
                                                                value = "optimal_number_cluster",
                                                                width = '100%', height = '35px')
                                        ),
                                        column(12,actionButton("getClusteringNumber_process","Process",width="100%",class="btn btn-primary"))
                                    )
                                )
                          ),
                          
                          #3.Consensus Clustering
                          conditionalPanel(condition="input.getModuleStep==3",
                               wellPanel(style="background-color:#47AEE9",
                                   fixedRow(
                                       column(12,tags$h4(strong("Consensus Clustering"),align="center",style="color:#834005")),
                                       column(12,offset=0,tags$p(h5(strong("This step aims to perform multi-omics integrative clustering by specifying one or more algorithms at once.",align="left")))),
                                       column(12,offset=0,tags$p(h5(strong("Now let's choose clustering algorithms at first. If you want to get consensus results from different algorithms, please choose at least two algorithms!",align="left")))),
                                       column(12,checkboxGroupInput(inputId = "multiClusteringAlgorithms",
                                                                    label = "Clustering algorithms (Choose at least two types if you want to get consensus results)",
                                                                    choices = list("iClusterBayes","SNF","PINSPlus","NEMO","COCA","LRAcluster","ConsensusClustering","IntNMF","CIMLR","MoCluster"),
                                                                    selected = list("iClusterBayes","SNF","PINSPlus","NEMO","COCA","LRAcluster","ConsensusClustering","IntNMF","CIMLR","MoCluster")
                                       ))
                                   )
                               ),
                               
                               #只选择了一种聚类方法
                               conditionalPanel(condition="input.multiClusteringAlgorithms.length==1",
                                    
                                    #1.iClusterBayes
                                    conditionalPanel(condition="input.multiClusteringAlgorithms.includes('iClusterBayes')",
                                        wellPanel(style="background-color:#47AEE9",
                                            fixedRow(
                                                column(12,tags$h4(strong("iClusterBayes"),align="center",style="color:#834005")),
                                                
                                                #The number of clusters
                                                column(12,radioButtons("iClusterBayes_N_clust", h5(strong("The number of clusters")), inline=TRUE, choices = c("System optimal"=1,"User defined"=2),selected=1)),
                                                conditionalPanel(condition="input.iClusterBayes_N_clust==2",
                                                   column(12,numericInput("iClusterBayes_N_clusts",label="User defined cluster number",value=2,min=2,step=1))
                                                ),
                                                
                                                #The prior probability for the indicator variable gamma of each subdataset
                                                ##mRNA
                                                conditionalPanel(condition="input.multiOmics.includes('1')",
                                                    column(12,numericInput("iClusterBayes_mRNA_prior_gamma",label="The prior probability for the indicator variable gamma of mRNA subdataset",value=0.5,min=0,max=1,step=0.05))
                                                ),
                                                
                                                ##lncRNA
                                                conditionalPanel(condition="input.multiOmics.includes('2')",
                                                    column(12,numericInput("iClusterBayes_lncRNA_prior_gamma",label="The prior probability for the indicator variable gamma of lncRNA subdataset",value=0.5,min=0,max=1,step=0.05))             
                                                ),
                                                
                                                ##DNA methylation
                                                conditionalPanel(condition="input.multiOmics.includes('3')",
                                                    column(12,numericInput("iClusterBayes_DNA_methylation_prior_gamma",label="The prior probability for the indicator variable gamma of DNA methylation subdataset",value=0.5,min=0,max=1,step=0.05))                         
                                                ),
                                                
                                                ##copy number alterations
                                                conditionalPanel(condition="input.multiOmics.includes('4')",
                                                    column(12,numericInput("iClusterBayes_copy_number_alterations_prior_gamma",label="The prior probability for the indicator variable gamma of copy number alterations subdataset",value=0.5,min=0,max=1,step=0.05))            
                                                ),
                                                
                                                ##binary somatic mutation
                                                conditionalPanel(condition="input.multiOmics.includes('5')",
                                                    column(12,numericInput("iClusterBayes_binary_somatic_mutation_prior_gamma",label="The prior probability for the indicator variable gamma of binary somatic mutation subdataset",value=0.5,min=0,max=1,step=0.05))                         
                                                ),
                                                
                                                ##radiomics
                                                conditionalPanel(condition="input.multiOmics.includes('6')",
                                                    column(12,numericInput("iClusterBayes_radiomics_prior_gamma",label="The prior probability for the indicator variable gamma of radiomics subdataset",value=0.5,min=0,max=1,step=0.05))                         
                                                ),
                                                
                                                #The number of MCMC burnin
                                                column(12,numericInput("iClusterBayes_burnin",label="The number of MCMC burnin",value=18000,min=0,step=1000)),
                                                
                                                #The number of MCMC draw
                                                column(12,numericInput("iClusterBayes_draw",label="The number of MCMC draw",value=12000,min=0,step=1000)),
                                                
                                                #The standard deviation of random walk proposal for the latent variable
                                                column(12,numericInput("iClusterBayes_sdev",label="The standard deviation of random walk proposal for the latent variable",value=0.05,min=0,step=0.05)),
                                                
                                                #A numerical value to thin the MCMC chain in order to reduce autocorrelation
                                                column(12,numericInput("iClusterBayes_thin",label="A numerical value to thin the MCMC chain in order to reduce autocorrelation",value=3,step=1)),
                                                
                                                column(12,actionButton("iClusterBayes_process","Process",width="100%",class="btn btn-primary"))
                                            )
                                        )
                                    ),
                                    
                                    #2.SNF
                                    conditionalPanel(condition="input.multiClusteringAlgorithms.includes('SNF')",
                                        wellPanel(style="background-color:#47AEE9",
                                            fixedRow(
                                                column(12,tags$h4(strong("SNF"),align="center",style="color:#834005")),
                                                
                                                #The number of clusters
                                                column(12,radioButtons("SNF_N_clust", h5(strong("The number of clusters")), inline=TRUE, choices = c("System optimal"=1,"User defined"=2),selected=1)),
                                                conditionalPanel(condition="input.SNF_N_clust==2",
                                                    column(12,numericInput("SNF_N_clusts",label="User defined cluster number",value=2,min=2,step=1))
                                                ),
                                                
                                                #The number of neighbors in K-nearest neighbors part of the algorithm
                                                column(12,numericInput("SNF_K",label="The number of neighbors in K-nearest neighbors part of the algorithm",value=30,min=0,step=1)),
                                                
                                                #The number of interations for the diffusion process
                                                column(12,numericInput("SNF_t",label="The number of interations for the diffusion process",value=20,min=0,step=1)),
                                                
                                                #The variance for local model
                                                column(12,numericInput("SNF_sigma",label="The variance for local model",value=0.5,min=0,step=0.05)),
                                                
                                                column(12,actionButton("SNF_process","Process",width="100%",class="btn btn-primary"))
                                            )
                                        )
                                    ),
                                    
                                    #3.PINSPlus
                                    conditionalPanel(condition="input.multiClusteringAlgorithms.includes('PINSPlus')",
                                        wellPanel(style="background-color:#47AEE9",
                                            fixedRow(
                                                column(12,tags$h4(strong("PINSPlus"),align="center",style="color:#834005")),
                                                
                                                #The number of clusters
                                                column(12,radioButtons("PINSPlus_N_clust", h5(strong("The number of clusters")), inline=TRUE, choices = c("System optimal"=1,"User defined"=2),selected=1)),
                                                conditionalPanel(condition="input.PINSPlus_N_clust==2",
                                                    column(12,numericInput("PINSPlus_N_clusts",label="User defined cluster number",value=2,min=2,step=1))
                                                ),
                                                
                                                #The normalization method for consensus clustering
                                                column(12,radioButtons("PINSPlus_norMethod", h5(strong("The normalization method for consensus clustering")), inline=FALSE, choices = c("median-centered"=1,"mean-centered"=2,"z-score"=3,"none"=4),selected=4)),
                                                
                                                #Built-in clustering algorithm that PerturbationClustering will use
                                                column(12,radioButtons("PINSPlus_clusteringMethod", h5(strong("Built-in clustering algorithm that PerturbationClustering will use")), inline=FALSE, choices = c("kmeans"=1,"pam"=2,"hclust"=3),selected=1)),
                                                
                                                #The minimum number of iterations
                                                column(12,numericInput("PINSPlus_iterMin",label="The minimum number of iterations",value=50,min=0,step=10)),
                                                
                                                #The maximum number of iterations
                                                column(12,numericInput("PINSPlus_iterMax",label="The maximum number of iterations",value=500,min=0,step=10)),
                                                
                                                column(12,actionButton("PINSPlus_process","Process",width="100%",class="btn btn-primary"))
                                            )
                                        )
                                    ),
                                    
                                    #4.NEMO
                                    conditionalPanel(condition="input.multiClusteringAlgorithms.includes('NEMO')",
                                        wellPanel(style="background-color:#47AEE9",
                                            fixedRow(
                                                column(12,tags$h4(strong("NEMO"),align="center",style="color:#834005")),
                                                
                                                #The number of clusters
                                                column(12,radioButtons("NEMO_N_clust", h5(strong("The number of clusters")), inline=TRUE, choices = c("System optimal"=1,"User defined"=2),selected=1)),
                                                conditionalPanel(condition="input.NEMO_N_clust==2",
                                                    column(12,numericInput("NEMO_N_clusts",label="User defined cluster number",value=2,min=2,step=1))
                                                ),
                                                
                                                #The number of neighbors to use for each omic
                                                column(12,radioButtons("NEMO_numNeighbors", h5(strong("The number of neighbors to use for each omic")), inline=FALSE, choices = c("The number of neighbors used for all omics"=1,"A list of numbers for each omic"=2,"The number of samples divided by NUM.NEIGHBORS.RATIO"=3),selected=3)),
                                                conditionalPanel(condition="input.NEMO_numNeighbors==1",
                                                    column(12,numericInput("NEMO_numNeighbors1",label="The number of neighbors used for all omics",value=30,min=0,step=1))            
                                                ),
                                                conditionalPanel(condition="input.NEMO_numNeighbors==2",
                                                    column(12,offset=0,tags$p(h5(strong("Input the number of neighbors for each omic",align="left")))),
                                                    column(10,rHandsontableOutput("NEMO_numNeighbors2")),
                                                    column(12,offset=0,tags$p(h5(strong("",align="left"))))  #空一行
                                                ),
                                                
                                                column(12,actionButton("NEMO_process","Process",width="100%",class="btn btn-primary"))
                                            )
                                        )
                                    ),
                                    
                                    #5.COCA
                                    conditionalPanel(condition="input.multiClusteringAlgorithms.includes('COCA')",
                                        wellPanel(style="background-color:#47AEE9",
                                            fixedRow(
                                                column(12,tags$h4(strong("COCA"),align="center",style="color:#834005")),
                                                
                                                #The number of clusters
                                                column(12,radioButtons("COCA_N_clust", h5(strong("The number of clusters")), inline=TRUE, choices = c("System optimal"=1,"User defined"=2),selected=1)),
                                                conditionalPanel(condition="input.COCA_N_clust==2",
                                                    column(12,numericInput("COCA_N_clusts",label="User defined cluster number",value=2,min=2,step=1))
                                                ),
                                                
                                                #Clustering methods to be used to cluster the observations in each subdataset
                                                column(12,textInput(inputId="COCA_methods",label="Clustering methods to be used to cluster the observations in each subdataset",value="hclust")),
                                                
                                                #Distances to be used in the clustering step for each subdataset
                                                column(12,textInput(inputId="COCA_distances",label="Distances to be used in the clustering step for each subdataset",value="euclidean")),
                                                
                                                column(12,actionButton("COCA_process","Process",width="100%",class="btn btn-primary"))
                                            )
                                        )
                                    ),
                                    
                                    #6.LRAcluster
                                    conditionalPanel(condition="input.multiClusteringAlgorithms.includes('LRAcluster')",
                                        wellPanel(style="background-color:#47AEE9",
                                            fixedRow(
                                                column(12,tags$h4(strong("LRAcluster"),align="center",style="color:#834005")),
                                                
                                                #The number of clusters
                                                column(12,radioButtons("LRAcluster_N_clust", h5(strong("The number of clusters")), inline=TRUE, choices = c("System optimal"=1,"User defined"=2),selected=1)),
                                                conditionalPanel(condition="input.LRAcluster_N_clust==2",
                                                    column(12,numericInput("LRAcluster_N_clusts",label="User defined cluster number",value=2,min=2,step=1))
                                                ),
                                                
                                                #The cluster algorithm for similarity matrix
                                                column(12,textInput(inputId="LRAcluster_clusterAlg",label="The cluster algorithm for similarity matrix",value="ward.D")),
                                                
                                                column(12,actionButton("LRAcluster_process","Process",width="100%",class="btn btn-primary"))
                                            )
                                        )
                                    ),
                                    
                                    #7.ConsensusClustering
                                    conditionalPanel(condition="input.multiClusteringAlgorithms.includes('ConsensusClustering')",
                                        wellPanel(style="background-color:#47AEE9",
                                            fixedRow(
                                                column(12,tags$h4(strong("ConsensusClustering"),align="center",style="color:#834005")),
                                                
                                                #The number of clusters
                                                column(12,radioButtons("ConsensusClustering_N_clust", h5(strong("The number of clusters")), inline=TRUE, choices = c("System optimal"=1,"User defined"=2),selected=1)),
                                                conditionalPanel(condition="input.ConsensusClustering_N_clust==2",
                                                    column(12,numericInput("ConsensusClustering_N_clusts",label="User defined cluster number",value=2,min=2,step=1))
                                                ),
                                                
                                                #The normalization method for consensus clustering
                                                column(12,radioButtons("ConsensusClustering_norMethod", h5(strong("The normalization method for consensus clustering")), inline=FALSE, choices = c("median-centered"=1,"mean-centered"=2,"z-score"=3,"none"=4),selected=4)),
                                                
                                                #The number of subsamples
                                                column(12,numericInput("ConsensusClustering_reps",label="The number of subsamples",value=500,min=0,step=50)),
                                                
                                                #The proportion of items to sample
                                                column(12,numericInput("ConsensusClustering_pItem",label="The proportion of items to sample",value=0.8,min=0,max=1,step=0.05)),
                                                
                                                #The proportion of features to sample
                                                column(12,numericInput("ConsensusClustering_pFeature",label="The proportion of features to sample",value=0.8,min=0,max=1,step=0.05)),
                                                
                                                #The cluster algorithm
                                                column(12,textInput(inputId="ConsensusClustering_clusterAlg",label="The cluster algorithm",value="hc")),
                                                
                                                #The heirachical linakge method for subsampling
                                                column(12,textInput(inputId="ConsensusClustering_innerLinkage",label="The heirachical linakge method for subsampling",value="ward.D")),
                                                
                                                #The heirarchical method for consensus matrix
                                                column(12,textInput(inputId="ConsensusClustering_finalLinkage",label="The heirarchical method for consensus matrix",value="ward.D")),
                                                
                                                #The distance function
                                                column(12,textInput(inputId="ConsensusClustering_distance",label="The distance function",value="pearson")),
                                                
                                                #The output format for heatmap
                                                column(12,radioButtons("ConsensusClustering_plot", h5(strong("The output format for heatmap")), inline=FALSE, choices = c("pdf"=1,"png"=2,"pngBMP"=3,"none"=4),selected=4)),
                                                
                                                #Write output and log to csv or not
                                                column(12,radioButtons("ConsensusClustering_writeTable", h5(strong("Write output and log to csv or not")), inline=TRUE, choices = c("Yes"=1,"No"=2),selected=2)),
                                                
                                                #The name of output directory
                                                conditionalPanel(condition="(input.ConsensusClustering_plot != 4) | (input.ConsensusClustering_writeTable == 1)",
                                                    column(12,textInput(inputId="ConsensusClustering_title",label="The name of output directory",value="consensuscluster"))
                                                ),
                                                
                                                #A random seed for reproducible results
                                                column(12,numericInput("ConsensusClustering_seed",label="A random seed for reproducible results",value=123456,step=1)),
                                                
                                                column(12,actionButton("ConsensusClustering_process","Process",width="100%",class="btn btn-primary"))
                                            )
                                        )
                                    ),
                                    
                                    #8.IntNMF
                                    conditionalPanel(condition="input.multiClusteringAlgorithms.includes('IntNMF')",
                                        wellPanel(style="background-color:#47AEE9",
                                            fixedRow(
                                                column(12,tags$h4(strong("IntNMF"),align="center",style="color:#834005")),
                                                # column(12,offset=0,tags$p(h5(strong("You don't need to indicate any parameters, just click the 'process' button below.",align="left")))),
                                                
                                                #The number of clusters
                                                column(12,radioButtons("IntNMF_N_clust", h5(strong("The number of clusters")), inline=TRUE, choices = c("System optimal"=1,"User defined"=2),selected=1)),
                                                conditionalPanel(condition="input.IntNMF_N_clust==2",
                                                    column(12,numericInput("IntNMF_N_clusts",label="User defined cluster number",value=2,min=2,step=1))
                                                ),
                                                
                                                column(12,actionButton("IntNMF_process","Process",width="100%",class="btn btn-primary"))
                                            )
                                        )
                                    ),
                                    
                                    #9.CIMLR
                                    conditionalPanel(condition="input.multiClusteringAlgorithms.includes('CIMLR')",
                                        wellPanel(style="background-color:#47AEE9",
                                            fixedRow(
                                                column(12,tags$h4(strong("CIMLR"),align="center",style="color:#834005")),
                                                
                                                #The number of clusters
                                                column(12,radioButtons("CIMLR_N_clust", h5(strong("The number of clusters")), inline=TRUE, choices = c("System optimal"=1,"User defined"=2),selected=1)),
                                                conditionalPanel(condition="input.CIMLR_N_clust==2",
                                                    column(12,numericInput("CIMLR_N_clusts",label="User defined cluster number",value=2,min=2,step=1))
                                                ),
                                                
                                                #Ratio of the number of cores to be used when computing the multi-kernel
                                                column(12,numericInput("CIMLR_cores_ratio",label="Ratio of the number of cores to be used when computing the multi-kernel",value=0,min=0,max=1,step=0.05)),
                                                
                                                column(12,actionButton("CIMLR_process","Process",width="100%",class="btn btn-primary"))
                                            )
                                        )
                                    ),
                                    
                                    #10.MoCluster
                                    conditionalPanel(condition="input.multiClusteringAlgorithms.includes('MoCluster')",
                                        wellPanel(style="background-color:#47AEE9",
                                            fixedRow(
                                                column(12,tags$h4(strong("MoCluster"),align="center",style="color:#834005")),
                                                
                                                #The number of clusters
                                                column(12,radioButtons("MoCluster_N_clust", h5(strong("The number of clusters")), inline=TRUE, choices = c("System optimal"=1,"User defined"=2),selected=1)),
                                                conditionalPanel(condition="input.MoCluster_N_clust==2",
                                                    column(12,numericInput("MoCluster_N_clusts",label="User defined cluster number",value=2,min=2,step=1))
                                                ),
                                                
                                                #The number of components to calculate
                                                column(12,numericInput("MoCluster_ncomp",label="The number of components to calculate",value=10,min=0,step=1)),
                                                
                                                #Method
                                                column(12,radioButtons("MoCluster_method", h5(strong("Method")), inline=FALSE, choices = c("CPCA"=1,"GCCA"=2,"MCIA"=3),selected=1)),
                                                
                                                #Normalization of different matrices
                                                column(12,radioButtons("MoCluster_option", h5(strong("Normalization of different matrices")), inline=FALSE, choices = c("lambda1"=1,"inertia"=2,"uniform"=3),selected=1)),
                                                
                                                #A numeric value to indicate the absolute number (if k >= 1) or the proportion (if 0 < k < 1) of non-zero coefficients for the variable loading vectors
                                                column(12,offset=0,tags$p(h5(strong("A numeric value to indicate the absolute number (if k >= 1) or the proportion (if 0 < k < 1) of non-zero coefficients for the variable loading vectors.",align="left")))),
                                                ##The format
                                                column(12,radioButtons("MoCluster_k_format", h5(strong("Format")), inline=FALSE, choices = c("Absolute number"=1,"Proportion"=2,"all"=3),selected=1)),
                                                ##Setting
                                                conditionalPanel(condition="input.MoCluster_k_format==1",
                                                    column(12,radioButtons("MoCluster_k_setting1", h5(strong("Setting")), inline=TRUE, choices = c("For all omics"=1,"For each omic"=2),selected=1)),
                                                    conditionalPanel(condition="input.MoCluster_k_setting1==1",
                                                        column(12,numericInput("MoCluster_k1",label="The absolute number",value=10,min=1,step=1))
                                                    ),
                                                    conditionalPanel(condition="input.MoCluster_k_setting1==2",
                                                        column(12,offset=0,tags$p(h5(strong("Input the absolute number for each omic",align="left")))),
                                                        column(10,rHandsontableOutput("table.MoCluster_k1")),
                                                        column(12,offset=0,tags$p(h5(strong("",align="left"))))  #空一行
                                                    )
                                                ),
                                                conditionalPanel(condition="input.MoCluster_k_format==2",
                                                    column(12,radioButtons("MoCluster_k_setting2", h5(strong("Setting")), inline=TRUE, choices = c("For all omics"=1,"For each omic"=2),selected=1)),
                                                    conditionalPanel(condition="input.MoCluster_k_setting2==1",
                                                        column(12,numericInput("MoCluster_k2",label="The proportion",value=0.1,min=0,max=1,step=0.05))
                                                    ),
                                                    conditionalPanel(condition="input.MoCluster_k_setting2==2",
                                                        column(12,offset=0,tags$p(h5(strong("Input the proportion for each omic",align="left")))),
                                                        column(10,rHandsontableOutput("table.MoCluster_k2")),
                                                        column(12,offset=0,tags$p(h5(strong("",align="left"))))  #空一行
                                                    )
                                                ),
                                                
                                                #The variables should be centered or not
                                                column(12,radioButtons("MoCluster_center", h5(strong("The variables should be centered or not")), inline=TRUE, choices = c("Yes"=1,"No"=2),selected=1)),
                                                
                                                #The variables should be scaled or not
                                                column(12,radioButtons("MoCluster_scale", h5(strong("The variables should be scaled or not")), inline=TRUE, choices = c("Yes"=1,"No"=2),selected=1)),
                                                
                                                #The cluster algorithm for distance
                                                column(12,textInput(inputId="MoCluster_clusterAlg",label="The cluster algorithm for distance",value="ward.D")),
                                                
                                                column(12,actionButton("MoCluster_process","Process",width="100%",class="btn btn-primary"))
                                            )
                                        )
                                    )
                               ),
                               
                               #选择了至少两种聚类方法
                               conditionalPanel(condition="input.multiClusteringAlgorithms.length>1",
                                                wellPanel(style="background-color:#47AEE9",
                                                          fixedRow(
                                                              column(12,tags$h4(strong("Consensus Clustering"),align="center",style="color:#834005")),
                                                              column(12,offset=0,tags$p(h5(strong("You choose more than 1 algorithm and all of them shall be run with parameters by default.",align="left")))),
                                                              
                                                              #The number of clusters
                                                              column(12,radioButtons("Consensus_N_clust", h5(strong("The number of clusters")), inline=TRUE, choices = c("System optimal"=1,"User defined"=2),selected=1)),
                                                              conditionalPanel(condition="input.Consensus_N_clust==2",
                                                                    column(12,numericInput("Consensus_N_clusts",label="User defined cluster number",value=2,min=2,step=1))
                                                              ),
                                                              
                                                              column(12,actionButton("consensus_clustering_process1","Process",width="100%",class="btn btn-primary"))
                                                          )
                                                ),
                                                wellPanel(style="background-color:#47AEE9",
                                                          fixedRow(
                                                              column(12,tags$h4(strong("Consensus Heatmap"),align="center",style="color:#834005")),
                                                              #Distance measurement for hierarchical clustering
                                                              column(12,textInput(inputId="consensus_clustering_distance",label="Distance measurement for hierarchical clustering",value="euclidean")),
                                                              #Clustering method for hierarchical clustering
                                                              column(12,textInput(inputId="consensus_clustering_linkage",label="Clustering method for hierarchical clustering",value="ward.D")),
                                                              #Heatmap mapping color
                                                              column(12,radioButtons("consensus_clustering_mapcolor", h5(strong("Heatmap mapping color settings")), inline=TRUE, choices = c("System default"=1,"User defined"=2),selected=1)),
                                                              conditionalPanel(condition="input.consensus_clustering_mapcolor==2",
                                                                  column(12,offset=0,tags$p(h5(strong("Input heatmap mapping colors (use hex color format, e.g. #000004FF)",align="left")))),
                                                                  column(9,rHandsontableOutput("table.consensus_clustering_mapcolor")),
                                                                  column(12,offset=0,tags$p(h5(strong("",align="left"))))  #空一行
                                                              ),
                                                              #Colors for annotating each cluster at the top of heatmap
                                                              column(12,radioButtons("consensus_clustering_clustcolor", h5(strong("Clustering subtypes color settings")), inline=TRUE, choices = c("System default"=1,"User defined"=2),selected=1)),
                                                              conditionalPanel(condition="input.consensus_clustering_clustcolor==2",
                                                                  column(12,offset=0,tags$p(h5(strong("Input colors for clustering subtypes (use hex color format, e.g. #2EC4B6)",align="left")))),
                                                                  column(9,rHandsontableOutput("table.consensus_clustering_clustcolor")),
                                                                  column(12,offset=0,tags$p(h5(strong("",align="left"))))  #空一行             
                                                              ),
                                                              #Show sample ID or not
                                                              column(12,radioButtons("consensus_clustering_showID", h5(strong("Show sample ID or not")), inline=TRUE, choices = c("Yes"=1,"No"=2),selected=2)),
                                                              #The name of the consensus heatmap
                                                              column(12,textAreaInput(inputId = "consensus_clustering_FigureName",
                                                                                      label = "Figure Name",
                                                                                      value = "consensusheatmap",
                                                                                      width = '100%', height = '35px')
                                                              ),
                                                              #The width of output figure
                                                              column(12,numericInput("consensus_clustering_FigureWidth",label="The width of output figure",value=5.5,min=0.5,step=0.5)),
                                                              #The height of output figure
                                                              column(12,numericInput("consensus_clustering_FigureHeight",label="The height of output figure",value=5,min=0.5,step=0.5)),
                                                              
                                                              column(12,actionButton("consensus_clustering_process2","Process",width="100%",class="btn btn-primary"))
                                                          )
                                                )
                               )
                          ),
                          
                          #4.Silhouette
                          conditionalPanel(condition="input.getModuleStep==4",
                              
                              wellPanel(style="background-color:#47AEE9",
                                        fixedRow(
                                            column(12,tags$h4(strong("Silhouette"),align="center",style="color:#834005")),
                                            column(12,offset=0,tags$p(h5(strong("This step aims to visualize silhouette information from consensus clustering.",align="left")))),
                                            #Colors for annotating each cluster
                                            column(12,radioButtons("silhouette_clustcolor", h5(strong("Colors for annotating each cluster")), inline=TRUE, choices = c("System default"=1,"User defined"=2),selected=1)),
                                            conditionalPanel(condition="input.silhouette_clustcolor==2",
                                                             column(12,offset=0,tags$p(h5(strong("Input colors for annotating each cluster (use hex color format, e.g. #2EC4B6)",align="left")))),
                                                             column(9,rHandsontableOutput("table.silhouette_clustcolor")),
                                                             column(12,offset=0,tags$p(h5(strong("",align="left"))))  #空一行             
                                            ),
                                            #The name of the silhouette plot
                                            column(12,textAreaInput(inputId = "silhouette_FigureName",
                                                                    label = "Figure Name",
                                                                    value = "silhouette",
                                                                    width = '100%', height = '35px')
                                            ),
                                            #The width of output figure
                                            column(12,numericInput("silhouette_FigureWidth",label="The width of output figure",value=5.5,min=0.5,step=0.5)),
                                            #The height of output figure
                                            column(12,numericInput("silhouette_FigureHeight",label="The height of output figure",value=5,min=0.5,step=0.5)),
                                            
                                            column(12,actionButton("silhouette_process","Process",width="100%",class="btn btn-primary"))
                                        )
                              )
                          ),
                          
                          #5.Multi-omics Heatmaps
                          conditionalPanel(condition="input.getModuleStep==5",
                            
                                           ##基本设置和标准化
                                           wellPanel(style="background-color:#47AEE9",
                                                     
                                                     fixedRow(
                                                         
                                                         ##聚类方法数为1
                                                         conditionalPanel(condition="input.multiClusteringAlgorithms.length==1",
                                                                          #iClusterBayes
                                                                          conditionalPanel(condition="input.multiClusteringAlgorithms.includes('iClusterBayes')",
                                                                                           
                                                                                           column(12,tags$h4(strong("Multi-omics Heatmaps (iClusterBayes)"),align="center",style="color:#834005"))        
                                                                          ),
                                                                          
                                                                          #SNF
                                                                          conditionalPanel(condition="input.multiClusteringAlgorithms.includes('SNF')",
                                                                                           
                                                                                           column(12,tags$h4(strong("Multi-omics Heatmaps (SNF)"),align="center",style="color:#834005"))        
                                                                          ),
                                                                          
                                                                          #PINSPlus
                                                                          conditionalPanel(condition="input.multiClusteringAlgorithms.includes('PINSPlus')",
                                                                                           
                                                                                           column(12,tags$h4(strong("Multi-omics Heatmaps (PINSPlus)"),align="center",style="color:#834005"))        
                                                                          ),
                                                                          
                                                                          #NEMO
                                                                          conditionalPanel(condition="input.multiClusteringAlgorithms.includes('NEMO')",
                                                                                           
                                                                                           column(12,tags$h4(strong("Multi-omics Heatmaps (NEMO)"),align="center",style="color:#834005"))        
                                                                          ),
                                                                          
                                                                          #COCA
                                                                          conditionalPanel(condition="input.multiClusteringAlgorithms.includes('COCA')",
                                                                                           
                                                                                           column(12,tags$h4(strong("Multi-omics Heatmaps (COCA)"),align="center",style="color:#834005"))        
                                                                          ),
                                                                          
                                                                          #LRAcluster
                                                                          conditionalPanel(condition="input.multiClusteringAlgorithms.includes('LRAcluster')",
                                                                                           
                                                                                           column(12,tags$h4(strong("Multi-omics Heatmaps (LRAcluster)"),align="center",style="color:#834005"))        
                                                                          ),
                                                                          
                                                                          #ConsensusClustering
                                                                          conditionalPanel(condition="input.multiClusteringAlgorithms.includes('ConsensusClustering')",
                                                                                           
                                                                                           column(12,tags$h4(strong("Multi-omics Heatmaps (ConsensusClustering)"),align="center",style="color:#834005"))        
                                                                          ),
                                                                          
                                                                          #IntNMF
                                                                          conditionalPanel(condition="input.multiClusteringAlgorithms.includes('IntNMF')",
                                                                                           
                                                                                           column(12,tags$h4(strong("Multi-omics Heatmaps (IntNMF)"),align="center",style="color:#834005"))        
                                                                          ),
                                                                          
                                                                          #CIMLR
                                                                          conditionalPanel(condition="input.multiClusteringAlgorithms.includes('CIMLR')",
                                                                                           
                                                                                           column(12,tags$h4(strong("Multi-omics Heatmaps (CIMLR)"),align="center",style="color:#834005"))        
                                                                          ),
                                                                          
                                                                          #MoCluster
                                                                          conditionalPanel(condition="input.multiClusteringAlgorithms.includes('MoCluster')",
                                                                                           
                                                                                           column(12,tags$h4(strong("Multi-omics Heatmaps (MoCluster)"),align="center",style="color:#834005"))        
                                                                          )
                                                         ),
                                                         
                                                         ##聚类方法数大于1
                                                         conditionalPanel(condition="input.multiClusteringAlgorithms.length>1",
                                                                          
                                                                          column(12,tags$h4(strong("Multi-omics Heatmaps (multiple algorithms)"),align="center",style="color:#834005"))
                                                         ),
                                                         
                                                         column(12,offset=0,tags$p(h5(strong("This step aims to vertically concatenate multiple heatmap derived from each omics data combined with clustering results and other annotation information.",align="left")))),
                                                         
                                                         ##聚类方法数大于1
                                                         conditionalPanel(condition="input.multiClusteringAlgorithms.length>1",
                                                                          
                                                                          column(12,radioButtons("multiple_algorithms_heatmap", h5(strong("Heatmaps for consensus clustering results of multiple clustering algorithms or not")), inline=TRUE, choices = c("Yes"=1,"No"=2),selected=1)),

                                                                          conditionalPanel(condition="input.multiple_algorithms_heatmap==2",
                                                                                           
                                                                                           #Choose the clustering results derived from specified algorithm for multi-omics heatmaps
                                                                                           column(12,textInput(inputId="specified_algorithm_heatmap",label="Indicate the clustering results derived from specified algorithm for multi-omics heatmaps (input one of the clustering algorithm listed in 'Consensus Clustering' step)",value=""))
                                                                          )
                                                         ),
                                                         
                                                         ##数据标准化设置
                                                         column(12,offset=0,tags$p(h5(strong("First let's get standardized multi-omics data:",align="left")))),
                                                         
                                                         ###Indicate if each omic data should be centered
                                                         column(12,offset=0,tags$p(h5(strong("Indicate if each omic data should be centered",align="left")))),
                                                         column(10,rHandsontableOutput("table.centerFlag_heatmap")),
                                                         column(12,offset=0,tags$p(h5(strong("",align="left")))),  #空一行
                                                         
                                                         ###Indicate if each omic data should be scaled
                                                         column(12,offset=0,tags$p(h5(strong("Indicate if each omic data should be scaled",align="left")))),
                                                         column(10,rHandsontableOutput("table.scaleFlag_heatmap")),
                                                         column(12,offset=0,tags$p(h5(strong("",align="left")))),  #空一行
                                                         
                                                         ###Assign truncating values for extreme values in continuous normalized multi-omics data
                                                         column(12,offset=0,tags$p(h5(strong("Assign truncating values for extreme values in continuous normalized multi-omics data (normalized values that exceed the truncating values will be replaced by truncating values, which is useful to map colors in heatmap)",align="left")))),
                                                         column(10,rHandsontableOutput("table.halfwidth_heatmap")),
                                                         column(12,offset=0,tags$p(h5(strong("",align="left"))))  #空一行
                                                     )
                                           ),
                                           
                                           ##热图参数
                                           wellPanel(style="background-color:#47AEE9",
                                                     
                                                    fixedRow(
                                                        
                                                        ##聚类方法数为1
                                                        conditionalPanel(condition="input.multiClusteringAlgorithms.length==1",
                                                                         #iClusterBayes
                                                                         conditionalPanel(condition="input.multiClusteringAlgorithms.includes('iClusterBayes')",
                                                                                          
                                                                                          column(12,tags$h4(strong("Multi-omics Heatmaps (iClusterBayes)"),align="center",style="color:#834005"))        
                                                                         ),
                                                                         
                                                                         #SNF
                                                                         conditionalPanel(condition="input.multiClusteringAlgorithms.includes('SNF')",
                                                                                          
                                                                                          column(12,tags$h4(strong("Multi-omics Heatmaps (SNF)"),align="center",style="color:#834005"))        
                                                                         ),
                                                                         
                                                                         #PINSPlus
                                                                         conditionalPanel(condition="input.multiClusteringAlgorithms.includes('PINSPlus')",
                                                                                          
                                                                                          column(12,tags$h4(strong("Multi-omics Heatmaps (PINSPlus)"),align="center",style="color:#834005"))        
                                                                         ),
                                                                         
                                                                         #NEMO
                                                                         conditionalPanel(condition="input.multiClusteringAlgorithms.includes('NEMO')",
                                                                                          
                                                                                          column(12,tags$h4(strong("Multi-omics Heatmaps (NEMO)"),align="center",style="color:#834005"))        
                                                                         ),
                                                                         
                                                                         #COCA
                                                                         conditionalPanel(condition="input.multiClusteringAlgorithms.includes('COCA')",
                                                                                          
                                                                                          column(12,tags$h4(strong("Multi-omics Heatmaps (COCA)"),align="center",style="color:#834005"))        
                                                                         ),
                                                                         
                                                                         #LRAcluster
                                                                         conditionalPanel(condition="input.multiClusteringAlgorithms.includes('LRAcluster')",
                                                                                          
                                                                                          column(12,tags$h4(strong("Multi-omics Heatmaps (LRAcluster)"),align="center",style="color:#834005"))        
                                                                         ),
                                                                         
                                                                         #ConsensusClustering
                                                                         conditionalPanel(condition="input.multiClusteringAlgorithms.includes('ConsensusClustering')",
                                                                                          
                                                                                          column(12,tags$h4(strong("Multi-omics Heatmaps (ConsensusClustering)"),align="center",style="color:#834005"))        
                                                                         ),
                                                                         
                                                                         #IntNMF
                                                                         conditionalPanel(condition="input.multiClusteringAlgorithms.includes('IntNMF')",
                                                                                          
                                                                                          column(12,tags$h4(strong("Multi-omics Heatmaps (IntNMF)"),align="center",style="color:#834005"))        
                                                                         ),
                                                                         
                                                                         #CIMLR
                                                                         conditionalPanel(condition="input.multiClusteringAlgorithms.includes('CIMLR')",
                                                                                          
                                                                                          column(12,tags$h4(strong("Multi-omics Heatmaps (CIMLR)"),align="center",style="color:#834005"))        
                                                                         ),
                                                                         
                                                                         #MoCluster
                                                                         conditionalPanel(condition="input.multiClusteringAlgorithms.includes('MoCluster')",
                                                                                          
                                                                                          column(12,tags$h4(strong("Multi-omics Heatmaps (MoCluster)"),align="center",style="color:#834005"))        
                                                                         )
                                                        ),
                                                        
                                                        ##聚类方法数大于1
                                                        conditionalPanel(condition="input.multiClusteringAlgorithms.length>1",
                                                                         
                                                                         column(12,tags$h4(strong("Multi-omics Heatmaps (multiple algorithms)"),align="center",style="color:#834005"))
                                                        ),
                                                        
                                                        ##热图参数设置
                                                        column(12,offset=0,tags$p(h5(strong("Then let's set the parameters of the heatmap:",align="left")))),
                                                        
                                                        ##Row title for each omic data
                                                        column(12,radioButtons("row_title_heatmap", h5(strong("Row title settings for each omic data")), inline=TRUE, choices = c("System default"=1,"User defined"=2),selected=1)),
                                                        ###User defined
                                                        conditionalPanel(condition="input.row_title_heatmap==2",
                                                                         column(12,offset=0,tags$p(h5(strong("Input row title for each omic data",align="left")))),
                                                                         column(10,rHandsontableOutput("table.row_title_heatmap")),
                                                                         column(12,offset=0,tags$p(h5(strong("",align="left"))))  #空一行             
                                                        ),
                                                        
                                                        ##Legend title for each omic data
                                                        column(12,radioButtons("legend_title_heatmap", h5(strong("Legend title settings for each omic data")), inline=TRUE, choices = c("System default"=1,"User defined"=2),selected=1)),
                                                        ###User defined
                                                        conditionalPanel(condition="input.legend_title_heatmap==2",
                                                                         column(12,offset=0,tags$p(h5(strong("Input legend title for each omic data",align="left")))),
                                                                         column(10,rHandsontableOutput("table.legend_title_heatmap")),
                                                                         column(12,offset=0,tags$p(h5(strong("",align="left"))))  #空一行             
                                                        ),
                                                        
                                                        ##Show the dendrogram for columns at the top of heatmap or not
                                                        conditionalPanel(condition="(input.multiClusteringAlgorithms.length==1 && (input.multiClusteringAlgorithms.includes('COCA') | input.multiClusteringAlgorithms.includes('LRAcluster') | input.multiClusteringAlgorithms.includes('ConsensusClustering') | input.multiClusteringAlgorithms.includes('MoCluster'))) | (input.multiClusteringAlgorithms.length>1 && (input.multiple_algorithms_heatmap==1 | (input.multiple_algorithms_heatmap==2 && (input.specified_algorithm_heatmap.trim()=='COCA' | input.specified_algorithm_heatmap.trim()=='LRAcluster' | input.specified_algorithm_heatmap.trim()=='ConsensusClustering' | input.specified_algorithm_heatmap.trim()=='MoCluster'))))",
                                                                         
                                                                         column(12,radioButtons("show_col_dend_heatmap", h5(strong("Show the dendrogram for columns at the top of heatmap or not")), inline=TRUE, choices = c("Yes"=1,"No"=2),selected=1))         
                                                        ),
                                                        
                                                        ##Show the sample names for columns at the bottom of heatmap or not
                                                        column(12,radioButtons("show_colnames_heatmap", h5(strong("Show the sample names for columns at the bottom of heatmap or not")), inline=TRUE, choices = c("Yes"=1,"No"=2),selected=2)),
                                                        
                                                        ##Show the feature names for rows of each omic data
                                                        column(12,offset=0,tags$p(h5(strong("Click the boxes below to indicate whether show the feature names for rows of each omic data",align="left")))),
                                                        column(10,rHandsontableOutput("table.show_rownames_heatmap")),
                                                        column(12,offset=0,tags$p(h5(strong("",align="left")))),  #空一行 
                                                        
                                                        ##Distance method for clustering each omic data at feature dimension
                                                        column(12,offset=0,tags$p(h5(strong("Input distance method for clustering each omic data at feature dimension",align="left")))),
                                                        column(10,rHandsontableOutput("table.clust_dist_row_heatmap")),
                                                        column(12,offset=0,tags$p(h5(strong("",align="left")))),  #空一行 
                                                        
                                                        ##Clustering method for each omic data at feature dimension
                                                        column(12,offset=0,tags$p(h5(strong("Input clustering method for each omic data at feature dimension",align="left")))),
                                                        column(10,rHandsontableOutput("table.clust_method_row_heatmap")),
                                                        column(12,offset=0,tags$p(h5(strong("",align="left")))),  #空一行
                                                        
                                                        ##Show dendrogram for rows of each omic data
                                                        column(12,offset=0,tags$p(h5(strong("Click the boxes below to indicate whether show dendrogram for rows of each omic data",align="left")))),
                                                        column(10,rHandsontableOutput("table.show_row_dend_heatmap")),
                                                        column(12,offset=0,tags$p(h5(strong("",align="left")))),  #空一行 
                                                        
                                                        ##Colors for annotating each cluster at the top of heatmap
                                                        column(12,radioButtons("heatmap_clustcolor", h5(strong("Colors for annotating each cluster at the top of heatmap")), inline=TRUE, choices = c("System default"=1,"User defined"=2),selected=1)),
                                                        conditionalPanel(condition="input.heatmap_clustcolor==2",
                                                                         column(12,offset=0,tags$p(h5(strong("Input colors for annotating each cluster (use hex color format, e.g. #2EC4B6)",align="left")))),
                                                                         column(9,rHandsontableOutput("table.heatmap_clustcolor")),
                                                                         column(12,offset=0,tags$p(h5(strong("",align="left"))))  #空一行             
                                                        ),
                                                        
                                                        ##Colors for heatmap of each omic data
                                                        column(12,radioButtons("heatmap_color", h5(strong("Colors for heatmap of each omic data")), inline=TRUE, choices = c("System default"=1,"User defined"=2),selected=1)),
                                                        conditionalPanel(condition="input.heatmap_color==2",
                                                                         
                                                                         ###mRNA
                                                                         conditionalPanel(condition="input.multiOmics.includes('1')",
                                                                                          
                                                                                          #The number of colors for mRNA dataset heatmap
                                                                                          column(12,numericInput("mRNA_heatmap_color_number",label="The number of colors for mRNA dataset heatmap",value=3,min=3,step=1)),
                                                                                          
                                                                                          #Input colors for mRNA dataset heatmap
                                                                                          column(12,offset=0,tags$p(h5(strong("Input colors for mRNA dataset heatmap (use hex color format, e.g. #00FF00)",align="left")))),
                                                                                          column(9,rHandsontableOutput("table.mRNA_heatmap_color")),
                                                                                          column(12,offset=0,tags$p(h5(strong("",align="left"))))  #空一行             
                                                                         ),
                                                                         
                                                                         ###lncRNA
                                                                         conditionalPanel(condition="input.multiOmics.includes('2')",
                                                                                          
                                                                                          #The number of colors for lncRNA dataset heatmap
                                                                                          column(12,numericInput("lncRNA_heatmap_color_number",label="The number of colors for lncRNA dataset heatmap",value=3,min=3,step=1)),
                                                                                          
                                                                                          #Input colors for lncRNA dataset heatmap
                                                                                          column(12,offset=0,tags$p(h5(strong("Input colors for lncRNA dataset heatmap (use hex color format, e.g. #6699CC)",align="left")))),
                                                                                          column(9,rHandsontableOutput("table.lncRNA_heatmap_color")),
                                                                                          column(12,offset=0,tags$p(h5(strong("",align="left"))))  #空一行             
                                                                         ),
                                                                         
                                                                         ###DNA methylation
                                                                         conditionalPanel(condition="input.multiOmics.includes('3')",
                                                                                          
                                                                                          #The number of colors for DNA methylation dataset heatmap
                                                                                          column(12,numericInput("DNA_methylation_heatmap_color_number",label="The number of colors for DNA methylation dataset heatmap",value=3,min=3,step=1)),
                                                                                          
                                                                                          #Input colors for DNA methylation dataset heatmap
                                                                                          column(12,offset=0,tags$p(h5(strong("Input colors for DNA methylation dataset heatmap (use hex color format, e.g. #0074FE)",align="left")))),
                                                                                          column(9,rHandsontableOutput("table.DNA_methylation_heatmap_color")),
                                                                                          column(12,offset=0,tags$p(h5(strong("",align="left"))))  #空一行             
                                                                         ),
                                                                         
                                                                         ###copy number alterations
                                                                         conditionalPanel(condition="input.multiOmics.includes('4')",
                                                                                          
                                                                                          #The number of colors for copy number alterations dataset heatmap
                                                                                          column(12,numericInput("copy_number_alterations_heatmap_color_number",label="The number of colors for copy number alterations dataset heatmap",value=3,min=3,step=1)),
                                                                                          
                                                                                          #Input colors for copy number alterations dataset heatmap
                                                                                          column(12,offset=0,tags$p(h5(strong("Input colors for copy number alterations dataset heatmap (use hex color format, e.g. #00FFFF)",align="left")))),
                                                                                          column(9,rHandsontableOutput("table.copy_number_alterations_heatmap_color")),
                                                                                          column(12,offset=0,tags$p(h5(strong("",align="left"))))  #空一行             
                                                                         ),
                                                                         
                                                                         ###binary somatic mutation
                                                                         conditionalPanel(condition="input.multiOmics.includes('5')",
                                                                                          
                                                                                          #Input colors for binary somatic mutation dataset heatmap
                                                                                          column(12,offset=0,tags$p(h5(strong("Input colors for binary somatic mutation dataset heatmap (use hex color format, e.g. #E5E5E5)",align="left")))),
                                                                                          column(9,rHandsontableOutput("table.binary_somatic_mutation_heatmap_color")),
                                                                                          column(12,offset=0,tags$p(h5(strong("",align="left"))))  #空一行             
                                                                         ),
                                                                         
                                                                         ###radiomics
                                                                         conditionalPanel(condition="input.multiOmics.includes('6')",
                                                                                          
                                                                                          #The number of colors for radiomics dataset heatmap
                                                                                          column(12,numericInput("radiomics_heatmap_color_number",label="The number of colors for radiomics dataset heatmap",value=3,min=3,step=1)),
                                                                                          
                                                                                          #Input colors for radiomics dataset heatmap
                                                                                          column(12,offset=0,tags$p(h5(strong("Input colors for radiomics dataset heatmap (use hex color format, e.g. #A52A2A)",align="left")))),
                                                                                          column(9,rHandsontableOutput("table.radiomics_heatmap_color")),
                                                                                          column(12,offset=0,tags$p(h5(strong("",align="left"))))  #空一行             
                                                                         )
                                                        ),
                                                        
                                                        ##Features should be annotated specifically in each omic data
                                                        conditionalPanel(condition="(input.multiClusteringAlgorithms.length==1 && (input.multiClusteringAlgorithms.includes('iClusterBayes') | input.multiClusteringAlgorithms.includes('CIMLR') | input.multiClusteringAlgorithms.includes('MoCluster'))) | (input.multiClusteringAlgorithms.length>1 && input.multiple_algorithms_heatmap==2 && (input.specified_algorithm_heatmap.trim()=='iClusterBayes' | input.specified_algorithm_heatmap.trim()=='CIMLR' | input.specified_algorithm_heatmap.trim()=='MoCluster'))",
                                                                         
                                                                         ###mRNA
                                                                         conditionalPanel(condition="input.multiOmics.includes('1')",
                                                                                          
                                                                                          #Features should be annotated specifically in mRNA dataset heatmap or not
                                                                                          column(12,radioButtons("mRNA_heatmap_annRow", h5(strong("Features should be annotated specifically in mRNA dataset heatmap or not")), inline=TRUE, choices = c("Yes"=1,"No"=2),selected=2)),
                                                                                          
                                                                                          conditionalPanel(condition="input.mRNA_heatmap_annRow==1",
                                                                                                           
                                                                                                           #Indicate the number of features to be annotated
                                                                                                           column(12,numericInput("mRNA_heatmap_annRow_number",label="The number of features to be annotated for mRNA dataset heatmap",value=10,min=1,step=1)),
                                                                                                           
                                                                                                           #Settings for features annotation in mRNA dataset heatmap
                                                                                                           column(12,radioButtons("mRNA_heatmap_annRow_setting", h5(strong("Settings for features annotation in mRNA dataset heatmap")), inline=TRUE, choices = c("Annotate top number of features"=1,"User defined"=2),selected=1)),
                                                                                                           
                                                                                                           #User defined features to be annotated
                                                                                                           conditionalPanel(condition="input.mRNA_heatmap_annRow_setting==2",
                                                                                                                            column(12,offset=0,tags$p(h5(strong("Input features to be annotated for mRNA dataset heatmap",align="left")))),
                                                                                                                            column(9,rHandsontableOutput("table.mRNA_heatmap_annRow")),
                                                                                                                            column(12,offset=0,tags$p(h5(strong("",align="left"))))  #空一行 
                                                                                                           )
                                                                                          )
                                                                                          
                                                                         ),
                                                                         
                                                                         ###lncRNA
                                                                         conditionalPanel(condition="input.multiOmics.includes('2')",
                                                                                          
                                                                                          #Features should be annotated specifically in lncRNA dataset heatmap or not
                                                                                          column(12,radioButtons("lncRNA_heatmap_annRow", h5(strong("Features should be annotated specifically in lncRNA dataset heatmap or not")), inline=TRUE, choices = c("Yes"=1,"No"=2),selected=2)),
                                                                                          
                                                                                          conditionalPanel(condition="input.lncRNA_heatmap_annRow==1",
                                                                                                           
                                                                                                           #Indicate the number of features to be annotated
                                                                                                           column(12,numericInput("lncRNA_heatmap_annRow_number",label="The number of features to be annotated for lncRNA dataset heatmap",value=10,min=1,step=1)),
                                                                                                           
                                                                                                           #Settings for features annotation in lncRNA dataset heatmap
                                                                                                           column(12,radioButtons("lncRNA_heatmap_annRow_setting", h5(strong("Settings for features annotation in lncRNA dataset heatmap")), inline=TRUE, choices = c("Annotate top number of features"=1,"User defined"=2),selected=1)),
                                                                                                           
                                                                                                           #User defined features to be annotated
                                                                                                           conditionalPanel(condition="input.lncRNA_heatmap_annRow_setting==2",
                                                                                                                            column(12,offset=0,tags$p(h5(strong("Input features to be annotated for lncRNA dataset heatmap",align="left")))),
                                                                                                                            column(9,rHandsontableOutput("table.lncRNA_heatmap_annRow")),
                                                                                                                            column(12,offset=0,tags$p(h5(strong("",align="left"))))  #空一行 
                                                                                                           )
                                                                                          )
                                                                         ),
                                                                         
                                                                         ###DNA methylation
                                                                         conditionalPanel(condition="input.multiOmics.includes('3')",
                                                                                          
                                                                                          #Features should be annotated specifically in DNA methylation dataset heatmap or not
                                                                                          column(12,radioButtons("DNA_methylation_heatmap_annRow", h5(strong("Features should be annotated specifically in DNA methylation dataset heatmap or not")), inline=TRUE, choices = c("Yes"=1,"No"=2),selected=2)),
                                                                                          
                                                                                          conditionalPanel(condition="input.DNA_methylation_heatmap_annRow==1",
                                                                                                           
                                                                                                           #Indicate the number of features to be annotated
                                                                                                           column(12,numericInput("DNA_methylation_heatmap_annRow_number",label="The number of features to be annotated for DNA methylation dataset heatmap",value=10,min=1,step=1)),
                                                                                                           
                                                                                                           #Settings for features annotation in DNA methylation dataset heatmap
                                                                                                           column(12,radioButtons("DNA_methylation_heatmap_annRow_setting", h5(strong("Settings for features annotation in DNA methylation dataset heatmap")), inline=TRUE, choices = c("Annotate top number of features"=1,"User defined"=2),selected=1)),
                                                                                                           
                                                                                                           #User defined features to be annotated
                                                                                                           conditionalPanel(condition="input.DNA_methylation_heatmap_annRow_setting==2",
                                                                                                                            column(12,offset=0,tags$p(h5(strong("Input features to be annotated for DNA methylation dataset heatmap",align="left")))),
                                                                                                                            column(9,rHandsontableOutput("table.DNA_methylation_heatmap_annRow")),
                                                                                                                            column(12,offset=0,tags$p(h5(strong("",align="left"))))  #空一行 
                                                                                                           )
                                                                                          )
                                                                         ),
                                                                         
                                                                         ###copy number alterations
                                                                         conditionalPanel(condition="input.multiOmics.includes('4')",
                                                                                          
                                                                                          #Features should be annotated specifically in copy number alterations dataset heatmap or not
                                                                                          column(12,radioButtons("copy_number_alterations_heatmap_annRow", h5(strong("Features should be annotated specifically in copy number alterations dataset heatmap or not")), inline=TRUE, choices = c("Yes"=1,"No"=2),selected=2)),
                                                                                          
                                                                                          conditionalPanel(condition="input.copy_number_alterations_heatmap_annRow==1",
                                                                                                           
                                                                                                           #Indicate the number of features to be annotated
                                                                                                           column(12,numericInput("copy_number_alterations_heatmap_annRow_number",label="The number of features to be annotated for copy number alterations dataset heatmap",value=10,min=1,step=1)),
                                                                                                           
                                                                                                           #Settings for features annotation in copy number alterations dataset heatmap
                                                                                                           column(12,radioButtons("copy_number_alterations_heatmap_annRow_setting", h5(strong("Settings for features annotation in copy number alterations dataset heatmap")), inline=TRUE, choices = c("Annotate top number of features"=1,"User defined"=2),selected=1)),
                                                                                                           
                                                                                                           #User defined features to be annotated
                                                                                                           conditionalPanel(condition="input.copy_number_alterations_heatmap_annRow_setting==2",
                                                                                                                            column(12,offset=0,tags$p(h5(strong("Input features to be annotated for copy number alterations dataset heatmap",align="left")))),
                                                                                                                            column(9,rHandsontableOutput("table.copy_number_alterations_heatmap_annRow")),
                                                                                                                            column(12,offset=0,tags$p(h5(strong("",align="left"))))  #空一行 
                                                                                                           )
                                                                                          )
                                                                                          
                                                                         ),
                                                                         
                                                                         ###binary somatic mutation
                                                                         conditionalPanel(condition="input.multiOmics.includes('5')",
                                                                                          
                                                                                          #Features should be annotated specifically in binary somatic mutation dataset heatmap or not
                                                                                          column(12,radioButtons("binary_somatic_mutation_heatmap_annRow", h5(strong("Features should be annotated specifically in binary somatic mutation dataset heatmap or not")), inline=TRUE, choices = c("Yes"=1,"No"=2),selected=2)),
                                                                                          
                                                                                          conditionalPanel(condition="input.binary_somatic_mutation_heatmap_annRow==1",
                                                                                                           
                                                                                                           #Indicate the number of features to be annotated
                                                                                                           column(12,numericInput("binary_somatic_mutation_heatmap_annRow_number",label="The number of features to be annotated for binary somatic mutation dataset heatmap",value=10,min=1,step=1)),
                                                                                                           
                                                                                                           #Settings for features annotation in binary somatic mutation dataset heatmap
                                                                                                           column(12,radioButtons("binary_somatic_mutation_heatmap_annRow_setting", h5(strong("Settings for features annotation in binary somatic mutation dataset heatmap")), inline=TRUE, choices = c("Annotate top number of features"=1,"User defined"=2),selected=1)),
                                                                                                           
                                                                                                           #User defined features to be annotated
                                                                                                           conditionalPanel(condition="input.binary_somatic_mutation_heatmap_annRow_setting==2",
                                                                                                                            column(12,offset=0,tags$p(h5(strong("Input features to be annotated for binary somatic mutation dataset heatmap",align="left")))),
                                                                                                                            column(9,rHandsontableOutput("table.binary_somatic_mutation_heatmap_annRow")),
                                                                                                                            column(12,offset=0,tags$p(h5(strong("",align="left"))))  #空一行 
                                                                                                           )
                                                                                          )
                                                                         ),
                                                                         
                                                                         ###radiomics
                                                                         conditionalPanel(condition="input.multiOmics.includes('6')",
                                                                                          
                                                                                          #Features should be annotated specifically in radiomics dataset heatmap or not
                                                                                          column(12,radioButtons("radiomics_heatmap_annRow", h5(strong("Features should be annotated specifically in radiomics dataset heatmap or not")), inline=TRUE, choices = c("Yes"=1,"No"=2),selected=2)),
                                                                                          
                                                                                          conditionalPanel(condition="input.radiomics_heatmap_annRow==1",
                                                                                                           
                                                                                                           #Indicate the number of features to be annotated
                                                                                                           column(12,numericInput("radiomics_heatmap_annRow_number",label="The number of features to be annotated for radiomics dataset heatmap",value=10,min=1,step=1)),
                                                                                                           
                                                                                                           #Settings for features annotation in mRNA dataset heatmap
                                                                                                           column(12,radioButtons("radiomics_heatmap_annRow_setting", h5(strong("Settings for features annotation in radiomics dataset heatmap")), inline=TRUE, choices = c("Annotate top number of features"=1,"User defined"=2),selected=1)),
                                                                                                           
                                                                                                           #User defined features to be annotated
                                                                                                           conditionalPanel(condition="input.radiomics_heatmap_annRow_setting==2",
                                                                                                                            column(12,offset=0,tags$p(h5(strong("Input features to be annotated for radiomics dataset heatmap",align="left")))),
                                                                                                                            column(9,rHandsontableOutput("table.radiomics_heatmap_annRow")),
                                                                                                                            column(12,offset=0,tags$p(h5(strong("",align="left"))))  #空一行 
                                                                                                           )
                                                                                          )
                                                                                          
                                                                         )
                                                        ),
                                                        
                                                        ##Sample annotations from survival information for heatmap
                                                        ###Sample annotations from survival information for heatmap or not
                                                        column(12,radioButtons("heatmap_annCol", h5(strong("Sample annotations from survival information for heatmap or not")), inline=TRUE, choices = c("Yes"=1,"No"=2),selected=2)),
                                                        
                                                        conditionalPanel(condition="input.heatmap_annCol==1",
                                                                         
                                                                         #The number of sample annotations from survival information for heatmap
                                                                         column(12,numericInput("heatmap_annCol_number",label="The number of sample annotations from survival information for heatmap",value=3,min=1,step=1)),
                                                                         
                                                                         #Input sample annotations from survival information for heatmap
                                                                         column(12,offset=0,tags$p(h5(strong("Input sample annotations from survival information for heatmap",align="left")))),
                                                                         column(12,offset=0,tags$p(h5("First line: Please input the sample annotation variables from survival information",align="left"))),
                                                                         column(12,offset=0,tags$p(h5("Second line: Please input 'Continuous' or 'Categorical' to indicate the type of each sample annotation variable",align="left"))),
                                                                         column(12,offset=0,tags$p(h5("Last line: Please input the colors for each sample annotation variable (use hex color format, e.g. #000004FF and English semicolons should be used to separate the input colors)",align="left"))),
                                                                         column(12,offset=0,tags$p(h5("Note1: If the sample annotation variable is continuous, the number of indicated colors should be equal to 3, which represents the minimum, median and maximum value of this variable)",align="left"))),
                                                                         column(12,offset=0,tags$p(h5("Note2: If the sample annotation variable is categorical, the number of indicated colors should be equal to the number of categories for this variable)",align="left"))),
                                                                         column(10,rHandsontableOutput("table.heatmap_annCol")),
                                                                         column(12,offset=0,tags$p(h5(strong("",align="left"))))  #空一行
                                                        ),
                                                        
                                                        ##The width for each heatmap
                                                        column(12,numericInput("width_heatmap",label="The width for each heatmap",value=6,min=0.5,step=0.5)),
                                                        
                                                        ##The height for each heatmap
                                                        column(12,numericInput("height_heatmap",label="The height for each heatmap",value=2,min=0.5,step=0.5)),
                                                        
                                                        ##The name of the heatmap
                                                        column(12,textAreaInput(inputId = "heatmap_FigureName",
                                                                                label = "Figure Name",
                                                                                value = "moheatmap",
                                                                                width = '100%', height = '35px')
                                                        ),
                                                        
                                                        column(12,actionButton("heatmap_process","Process",width="100%",class="btn btn-primary"))
                                                    )
                                           )
                          )
                    ),
                    br(),
                    column(7,
                           tabsetPanel(
                               tabPanel(tags$h5(tags$strong("Get Elites"),align="center",style="color:#FFA500"),br(),
                                        column(12,uiOutput("getElites_finish"))
                                        # column(12,uiOutput("getElites_finish")),
                                        # br(),br(),
                                        # ####getElites后数据展示
                                        # #####mRNA(tpm)
                                        # uiOutput("table.getElites_mRNA_tpm_Title"), #标题
                                        # column(12,br()),
                                        # dataTableOutput("table.getElites_mRNA_tpm"), #表格
                                        # # uiOutput("table.getElites_mRNA_tpm_Note"), #脚注
                                        # br(),br(),
                                        # #####lncRNA(tpm)
                                        # uiOutput("table.getElites_lncRNA_tpm_Title"), #标题
                                        # column(12,br()),
                                        # dataTableOutput("table.getElites_lncRNA_tpm"), #表格
                                        # # uiOutput("table.getElites_lncRNA_tpm_Note"), #脚注
                                        # br(),br(),
                                        # #####DNA methylation
                                        # uiOutput("table.getElites_DNA_methylation_Title"), #标题
                                        # column(12,br()),
                                        # dataTableOutput("table.getElites_DNA_methylation"), #表格
                                        # # uiOutput("table.getElites_DNA_methylation_Note"), #脚注
                                        # br(),br(),
                                        # #####copy number alterations
                                        # uiOutput("table.getElites_copy_number_alterations_Title"), #标题
                                        # column(12,br()),
                                        # dataTableOutput("table.getElites_copy_number_alterations"), #表格
                                        # # uiOutput("table.getElites_copy_number_alterations_Note"), #脚注
                                        # br(),br(),
                                        # #####binary somatic mutation
                                        # uiOutput("table.getElites_binary_somatic_mutation_Title"), #标题
                                        # column(12,br()),
                                        # dataTableOutput("table.getElites_binary_somatic_mutation"), #表格
                                        # # uiOutput("table.getElites_binary_somatic_mutation_Note"), #脚注
                                        # br(),br(),
                                        # #####radiomics
                                        # uiOutput("table.getElites_radiomics_Title"), #标题
                                        # column(12,br()),
                                        # dataTableOutput("table.getElites_radiomics") #表格
                                        # # uiOutput("table.getElites_radiomics_Note") #脚注
                               ),
                               tabPanel(tags$h5(tags$strong("Get Clustering Number"),align="center",style="color:#FFA500"),br(),
                                        # column(12,h3(plotOutput("getClusteringNumber_figure"))),
                                        uiOutput("getClusteringNumber_figure_unit"),
                                        br(),br(),
                                        column(12,uiOutput("getClusteringNumber_finish"))
                               ),
                               tabPanel(tags$h5(tags$strong("Consensus Clustering"),align="center",style="color:#FFA500"),br(),
                                        column(12,uiOutput("iClusterBayes_finish")),
                                        column(12,uiOutput("SNF_finish")),
                                        column(12,uiOutput("PINSPlus_finish")),
                                        column(12,uiOutput("NEMO_finish")),
                                        column(12,uiOutput("COCA_finish")),
                                        column(12,uiOutput("LRAcluster_finish")),
                                        column(12,uiOutput("ConsensusClustering_finish")),
                                        column(12,uiOutput("IntNMF_finish")),
                                        column(12,uiOutput("CIMLR_finish")),
                                        column(12,uiOutput("MoCluster_finish")),
                                        column(12,uiOutput("MOIC_finish")),
                                        br(),br(),
                                        # column(12,h3(plotOutput("consensus_heatmap_figure"))),
                                        column(12,uiOutput("consensus_heatmap_figure_unit")),
                                        br(),br(),
                                        column(12,uiOutput("cmoic_finish"))
                               ),
                               tabPanel(tags$h5(tags$strong("Silhouette"),align="center",style="color:#FFA500"),br(),
                                        # column(12,h3(plotOutput("silhouette_figure"))),
                                        uiOutput("silhouette_figure_unit"),
                                        br(),br(),
                                        column(12,uiOutput("silhouette_finish"))
                               ),
                               tabPanel(tags$h5(tags$strong("Multi-omics Heatmaps"),align="center",style="color:#FFA500"),br(),
                                        # column(12,h3(plotOutput("heatmap_figure"))),
                                        uiOutput("heatmap_figure_unit"),
                                        br(),br(),
                                        column(12,uiOutput("heatmap_finish"))
                               )
                           )
                    )
                ),
                
                #COMP Module
                tabPanel(
                    title=tags$p(tags$h4(tags$strong("COMP Module"),style="color:purple")),
                    value="COMP Module",
                    
                    column(5,
                           
                           #1.总的步骤选项,对应不同的参数列表
                           wellPanel(style="background-color:#47AEE9;",
                                     
                                     fixedRow(
                                         #设置步骤选项,每个步骤对应不同的参数
                                         column(12,radioButtons("compModuleStep", h5(strong("Steps")), inline=FALSE,
                                                                choices = c("Compare survival outcome"=1,"Compare clinical features"=2,"Compare mutational frequency"=3,"Compare total mutation burden"=4,"Compare fraction genome altered"=5,"Compare drug sensitivity"=6,"Compare agreement with other subtypes"=7),selected=1))
                                     )
                           ),
                           
                           #2.Compare survival outcome
                           conditionalPanel(condition="input.compModuleStep==1",
                                            wellPanel(style="background-color:#47AEE9",
                                                      fixedRow(
                                                          column(12,tags$h4(strong("Compare survival outcome"),align="center",style="color:#834005")),
                                                          column(12,offset=0,tags$p(h5(strong("In this step, we will compare the prognosis of different subtypes based on the clustering results from 'GET Module' by Kaplan-Meier survival curve.",align="left")))),
                                                          column(12,offset=0,tags$p(h5(strong("Pay attention: the format of survival time should be days and the values of survival status should be 0 or 1 (0: censoring; 1: event). Please make sure you provide the correct survival information first.",align="left")))),
                                                          #Compare survival outcome on tcga datasets or validation datasets
                                                          column(12,radioButtons("tcga_vali_compSurv", h5(strong("Compare survival outcome on TCGA datasets or Validation datasets")), inline=TRUE, choices = c("TCGA"=1,"Validation"=2),selected=1)),
                                                          #Model-free approaches for subtype prediction in validation cohort
                                                          conditionalPanel(condition="input.tcga_vali_compSurv==2",
                                                                           column(12,radioButtons("validation_method_compSurv", h5(strong("Model-free approaches for subtype prediction in validation cohort")), inline=TRUE, choices = c("NTP"=1,"PAM"=2),selected=1))
                                                          ),
                                                          
                                                          #Format convertion of the survival time
                                                          column(12,radioButtons("convt_time_compSurv", h5(strong("Format conversion of the survival time")), inline=TRUE, choices = c("Years"=1,"Months"=2, "No conversion"=3),selected=3)),
                                                          
                                                          #Setting for the x-axis cutoff for showing the maximal survival time
                                                          column(12,radioButtons("surv_cut_compSurv", h5(strong("Setting for the x-axis cutoff for showing the maximal survival time")), inline=TRUE, choices = c("System default"=1,"User defined"=2),selected=1)),
                                                          
                                                          #The x-axis cutoff for showing the maximal survival time
                                                          ##Year
                                                          conditionalPanel(condition="input.convt_time_compSurv==1 & input.surv_cut_compSurv==2",
                                                                           column(12,numericInput("surv_cut_value_compSurv",label="The x-axis cutoff for showing the maximal survival time (Unit: year)",value=10,step=1))             
                                                          ),
                                                          ##Month
                                                          conditionalPanel(condition="input.convt_time_compSurv==2 & input.surv_cut_compSurv==2",
                                                                           column(12,numericInput("surv_cut_value_compSurv",label="The x-axis cutoff for showing the maximal survival time (Unit: month)",value=120,step=12))             
                                                          ),
                                                          ##Day
                                                          conditionalPanel(condition="input.convt_time_compSurv==3 & input.surv_cut_compSurv==2",
                                                                           column(12,numericInput("surv_cut_value_compSurv",label="The x-axis cutoff for showing the maximal survival time (Unit: day)",value=3650,step=365))             
                                                          ),
                                                          
                                                          #Estimate probability of surviving beyond a certain number (x) of years (Estimating x-year survival)
                                                          column(12,radioButtons("xyrs_est_compSurv", h5(strong("Estimate probability of surviving beyond a certain number of years (Estimating x-year survival) or not")), inline=TRUE, choices = c("Yes"=1,"No"=2),selected=2)),
                                                          
                                                          conditionalPanel(condition="input.xyrs_est_compSurv==1",
                                                                           
                                                                           #The number of years for surviving probability estimation
                                                                           column(12,numericInput("xyrs_est_number_compSurv",label="The number of years for surviving probability estimation",value=2,min=1,max=10,step=1)),
                                                                           
                                                                           #Input the year
                                                                           column(12,offset=0,tags$p(h5(strong("Input the year",align="left")))),
                                                                           column(9,rHandsontableOutput("table.xyrs_est_compSurv")),
                                                                           column(12,offset=0,tags$p(h5(strong("",align="left"))))  #空一行
                                                          ),
                                                          
                                                          #Colors for each subtype
                                                          column(12,radioButtons("clustcolor_compSurv", h5(strong("Setting for colors of clustering subtypes")), inline=TRUE, choices = c("System default"=1,"User defined"=2),selected=1)),
                                                          
                                                          conditionalPanel(condition="input.clustcolor_compSurv==2",
                                                                           
                                                                           #Input color for each subtype
                                                                           column(12,offset=0,tags$p(h5(strong("Input color for each subtype (use hex color format, e.g. #2EC4B6)",align="left")))),
                                                                           column(9,rHandsontableOutput("table.compSurv_clustcolor")),
                                                                           column(12,offset=0,tags$p(h5(strong("",align="left"))))  #空一行             
                                                          ),
                                                          
                                                          #Method for adjusting p values
                                                          column(12,radioButtons("p_adjust_method_compSurv", h5(strong("Method for adjusting p values")), inline=FALSE, choices = c('holm'=1, 'hochberg'=2, 'hommel'=3, 'bonferroni'=4, 'BH'=5, 'BY'=6, 'fdr'=7, 'none'=8),selected=5)),
                                                          
                                                          #The way for drawing a horizontal/vertical line at median survival
                                                          column(12,radioButtons("surv_median_line_compSurv", h5(strong("The way for drawing a horizontal/vertical line at median survival")), inline=FALSE, choices = c('Horizontal line'=1, 'Vertical line'=2, 'Horizontal & Vertical'=3, 'No lines'=4),selected=4)),
                                                          
                                                          #The name of the kaplan-meier curve
                                                          conditionalPanel(condition="input.tcga_vali_compSurv==1",
                                                                           column(12,textAreaInput(inputId = "compSurv_FigureName_tcga",
                                                                                                   label = "Figure Name",
                                                                                                   value = "KAPLAN-MEIER_CURVE(TCGA)",
                                                                                                   width = '100%', height = '35px')
                                                                           )             
                                                          ),
                                                          conditionalPanel(condition="input.tcga_vali_compSurv==2",
                                                                           column(12,textAreaInput(inputId = "compSurv_FigureName_validation",
                                                                                                   label = "Figure Name",
                                                                                                   value = "KAPLAN-MEIER_CURVE(Validation)",
                                                                                                   width = '100%', height = '35px')
                                                                           )             
                                                          ),
                                                          
                                                          column(12,actionButton("compSurv_process","Process",width="100%",class="btn btn-primary"))
                                                      )
                                            )
                           ),
                           
                           #3.Compare clinical features
                           conditionalPanel(condition="input.compModuleStep==2",
                                            wellPanel(style="background-color:#47AEE9",
                                                      fixedRow(
                                                          column(12,tags$h4(strong("Compare clinical features"),align="center",style="color:#834005")),
                                                          column(12,offset=0,tags$p(h5(strong("In this step, a table that is easy to use in medical research papers will be created to summarize the specified baseline variables (continuous & categorical) stratified by specified categorical variable and then perform statistical tests.",align="left")))),
                                                          #Compare clinical features on tcga datasets or validation datasets
                                                          column(12,radioButtons("tcga_vali_compClinvar", h5(strong("Compare clinical features on tcga datasets or validation datasets")), inline=TRUE, choices = c("TCGA"=1,"Validation"=2),selected=1)),
                                                          #Model-free approaches for subtype prediction in validation cohort
                                                          conditionalPanel(condition="input.tcga_vali_compClinvar==2",
                                                                           column(12,radioButtons("validation_method_compClinvar", h5(strong("Model-free approaches for subtype prediction in validation cohort")), inline=TRUE, choices = c("NTP"=1,"PAM"=2),selected=1))
                                                          ),
                                                          
                                                          #Specify variables in survival and clinical information for summary and statistical tests
                                                          ##The number of variables in survival and clinical information chosen for summary and statistical tests
                                                          column(12,numericInput("variables_number_compClinvar",label="The number of variables in survival and clinical information chosen for summary and statistical tests",min=1,value=5,step=1)),
                                                          ##Input the variables
                                                          column(12,offset=0,tags$p(h5(strong("Input the variables",align="left")))),
                                                          column(9,rHandsontableOutput("table.variables_compClinvar")),
                                                          column(12,offset=0,tags$p(h5(strong("",align="left")))),  #空一行
                                                          
                                                          #User defined stratifying variable or not
                                                          column(12,radioButtons("strata_compClinvar", h5(strong("User defined stratifying variable or not")), inline=TRUE, choices = c('Yes'=1, 'No'=2),selected=2)),
                                                          conditionalPanel(condition="input.strata_compClinvar==1",
                                                                           column(12,textInput(inputId="strata_variable_compClinvar",label="Input the stratifying variable",value=""))             
                                                          ),
                                                          
                                                          #Indicate the categorical variables
                                                          ##The number of categorical variables in chosen variables
                                                          column(12,numericInput("categorical_number_compClinvar",label="The number of the categorical variables in chosen variables",min=0,value=3,step=1)),
                                                          ##Input the categorical variables
                                                          conditionalPanel(condition="input.categorical_number_compClinvar>0",
                                                                           column(12,offset=0,tags$p(h5(strong("Input the categorical variables",align="left")))),
                                                                           column(9,rHandsontableOutput("table.categorical_compClinvar")),
                                                                           column(12,offset=0,tags$p(h5(strong("",align="left"))))  #空一行                                                                           
                                                          ),
                                                          
                                                          #Specify the variables for which the p-values should be those of nonparametric tests
                                                          ##The number of variables for which the p-values should be those of nonparametric tests
                                                          column(12,numericInput("nonparametric_number_compClinvar",label="The number of variables for which the p-values should be those of nonparametric tests",min=0,value=0,step=1)),
                                                          ##Input the variables for which the p-values should be those of nonparametric tests
                                                          conditionalPanel(condition="input.nonparametric_number_compClinvar>0",
                                                                           column(12,offset=0,tags$p(h5(strong("Input variables for which the p-values should be those of nonparametric tests",align="left")))),
                                                                           column(9,rHandsontableOutput("table.nonparametric_compClinvar")),
                                                                           column(12,offset=0,tags$p(h5(strong("",align="left"))))  #空一行                                                                           
                                                          ),
                                                          
                                                          #Specify the variables for which the p-values should be those of exact tests
                                                          ##The number of variables for which the p-values should be those of exact tests
                                                          column(12,numericInput("exact_number_compClinvar",label="The number of variables for which the p-values should be those of exact tests",min=0,value=0,step=1)),
                                                          ##Input the variables for which the p-values should be those of exact tests
                                                          conditionalPanel(condition="input.exact_number_compClinvar>0",
                                                                           column(12,offset=0,tags$p(h5(strong("Input variables for which the p-values should be those of exact tests",align="left")))),
                                                                           column(9,rHandsontableOutput("table.exact_compClinvar")),
                                                                           column(12,offset=0,tags$p(h5(strong("",align="left"))))  #空一行                                                                           
                                                          ),
                                                          
                                                          #Whether NA should be handled as a regular factor level rather than missing value
                                                          column(12,radioButtons("includeNA_compClinvar", h5(strong("Whether NA should be handled as a regular factor level rather than missing value")), inline=TRUE, choices = c('Yes'=1, 'No'=2),selected=2)),
                                                          
                                                          #Transform the '.txt' output file to a '.docx' WORD file or not ('.txt' file will be also kept)
                                                          column(12,radioButtons("doWord_compClinvar", h5(strong("Transform the '.txt' output file to a '.docx' WORD file or not ('.txt' file will be also kept)")), inline=TRUE, choices = c('Yes'=1, 'No'=2),selected=1)),
                                                          
                                                          #The name of the output table
                                                          conditionalPanel(condition="input.tcga_vali_compClinvar==1",
                                                                           conditionalPanel(condition="input.strata_compClinvar==1",
                                                                                            column(12,textAreaInput(inputId = "compClinvar_TableName_tcga",
                                                                                                                    label = "The name of the output table",
                                                                                                                    value = "Summarization_of_clinical_variables_stratified_by_user_defined_stratifying_variable(TCGA)",
                                                                                                                    width = '100%', height = '35px')
                                                                                            )           
                                                                           ),
                                                                           conditionalPanel(condition="input.strata_compClinvar==2",
                                                                                            column(12,textAreaInput(inputId = "compClinvar_TableName_tcga",
                                                                                                                    label = "The name of the output table",
                                                                                                                    value = "Summarization_of_clinical_variables_stratified_by_current_subtypes(TCGA)",
                                                                                                                    width = '100%', height = '35px')
                                                                                            )           
                                                                           )
                                                                                  
                                                          ),
                                                          conditionalPanel(condition="input.tcga_vali_compClinvar==2",
                                                                           conditionalPanel(condition="input.strata_compClinvar==1",
                                                                                            column(12,textAreaInput(inputId = "compClinvar_TableName_validation",
                                                                                                                    label = "The name of the output table",
                                                                                                                    value = "Summarization_of_clinical_variables_stratified_by_user_defined_stratifying_variable(Validation)",
                                                                                                                    width = '100%', height = '35px')
                                                                                            )           
                                                                           ),
                                                                           conditionalPanel(condition="input.strata_compClinvar==2",
                                                                                            column(12,textAreaInput(inputId = "compClinvar_TableName_validation",
                                                                                                                    label = "The name of the output table",
                                                                                                                    value = "Summarization_of_clinical_variables_stratified_by_current_subtypes(Validation)",
                                                                                                                    width = '100%', height = '35px')
                                                                                            )           
                                                                           )     
                                                          ),
                                                         
                                                          column(12,actionButton("compClinvar_process","Process",width="100%",class="btn btn-primary"))
                                                      )
                                            )
                           ),
                           
                           #4.Compare mutational frequency
                           conditionalPanel(condition="input.compModuleStep==3",
                                            wellPanel(style="background-color:#47AEE9",
                                                      fixedRow(
                                                          column(12,tags$h4(strong("Compare mutational frequency"),align="center",style="color:#834005")),
                                                          column(12,offset=0,tags$p(h5(strong("In this step, a table and an oncoprint will be generated to compare mutational frequency among different multi-omics integerative clusters to test the independency between subtypes and mutational status.",align="left")))),
                                                          #Compare mutational frequency on tcga datasets or validation datasets
                                                          column(12,radioButtons("tcga_vali_compMut", h5(strong("Compare mutational frequency on tcga datasets or validation datasets")), inline=TRUE, choices = c("TCGA"=1,"Validation"=2),selected=1)),
                                                          #Model-free approaches for subtype prediction in validation cohort
                                                          conditionalPanel(condition="input.tcga_vali_compMut==2",
                                                                           column(12,radioButtons("validation_method_compMut", h5(strong("Model-free approaches for subtype prediction in validation cohort")), inline=TRUE, choices = c("NTP"=1,"PAM"=2),selected=1))
                                                          ),
                                                          
                                                          #The frequency cutoff for mutation data (Only features that mutated in over than such proportion would be included in testing)
                                                          column(12,numericInput("freq_cutoff_compMut",label="The frequency cutoff for mutation data (Only features that mutated in over than such proportion would be included in testing)",min=0,max=1,value=0.05,step=0.05)),
                                                          
                                                          #Statistical method for independence testing
                                                          column(12,radioButtons("test_method_compMut", h5(strong("Statistical method for independence testing")), inline=TRUE, choices = c('fisher'=1, 'chisq'=2),selected=1)),
                                                          
                                                          #The correction method for multiple comparison
                                                          column(12,radioButtons("p_adj_method_compMut", h5(strong("The correction method for multiple comparison")), inline=FALSE, choices = c('holm'=1, 'hochberg'=2, 'hommel'=3, 'bonferroni'=4, 'BH'=5, 'BY'=6, 'fdr'=7),selected=5)),
                                                          
                                                          #Transform the '.txt' output file to a '.docx' WORD file or not ('.txt' file will be also kept)
                                                          column(12,radioButtons("doWord_compMut", h5(strong("Transform the '.txt' output file to a '.docx' WORD file or not ('.txt' file will be also kept)")), inline=TRUE, choices = c('Yes'=1, 'No'=2),selected=1)),
                                                          
                                                          #The name of the output table
                                                          conditionalPanel(condition="input.tcga_vali_compMut==1",
                                                                           column(12,textAreaInput(inputId = "compMut_TableName_tcga",
                                                                                                   label = "The name of the output table",
                                                                                                   value = "Independent_test_between_subtype_and_mutation(TCGA)",
                                                                                                   width = '100%', height = '35px')
                                                                           )      
                                                          ),
                                                          
                                                          conditionalPanel(condition="input.tcga_vali_compMut==2",
                                                                           column(12,textAreaInput(inputId = "compMut_TableName_validation",
                                                                                                   label = "The name of the output table",
                                                                                                   value = "Independent_test_between_subtype_and_mutation(Validation)",
                                                                                                   width = '100%', height = '35px')
                                                                           )      
                                                          ),
                                                          
                                                          #Perform clustering within each subtype for oncoprint or not
                                                          column(12,radioButtons("innerclust_compMut", h5(strong("Perform clustering within each subtype for oncoprint or not")), inline=TRUE, choices = c('Yes'=1, 'No'=2),selected=1)),
                                                          
                                                          #The nominal p value cutoff for significant mutations shown in oncoprint
                                                          column(12,numericInput("p_cutoff_compMut",label="The nominal p value cutoff for significant mutations shown in oncoprint",min=0,max=1,value=0.05,step=0.05)),
                                                          
                                                          #The adjusted p value cutoff for significant mutations shown in oncoprint
                                                          column(12,numericInput("p_adj_cutoff_compMut",label="The adjusted p value cutoff for significant mutations shown in oncoprint",min=0,max=1,value=0.05,step=0.05)),
                                                          
                                                          #User defined mutation color for oncoprint or not
                                                          column(12,radioButtons("mut_col_compMut", h5(strong("User defined mutation color for oncoprint or not")), inline=TRUE, choices = c('Yes'=1, 'No'=2),selected=2)),
                                                          conditionalPanel(condition="input.mut_col_compMut==1",
                                                                           column(12,textInput(inputId="mut_col_value_compMut",label="Input the mutation color for oncoprint (use hex color format, e.g. #21498D)",value=""))       
                                                          ),
                                                          
                                                          #User defined background color for oncoprint or not
                                                          column(12,radioButtons("bg_col_compMut", h5(strong("User defined background color for oncoprint or not")), inline=TRUE, choices = c('Yes'=1, 'No'=2),selected=2)),
                                                          conditionalPanel(condition="input.bg_col_compMut==1",
                                                                           column(12,textInput(inputId="bg_col_value_compMut",label="Input the background color for oncoprint (use hex color format, e.g. #dcddde)",value=""))       
                                                          ),
                                                          
                                                          #Colors for each subtype
                                                          column(12,radioButtons("clustcolor_compMut", h5(strong("Setting for colors to annotate each subtype in oncoprint")), inline=TRUE, choices = c("System default"=1,"User defined"=2),selected=1)),
                                                          
                                                          conditionalPanel(condition="input.clustcolor_compMut==2",
                                                                           
                                                                           #Input color for each subtype
                                                                           column(12,offset=0,tags$p(h5(strong("Input color for each subtype (use hex color format, e.g. #2EC4B6)",align="left")))),
                                                                           column(9,rHandsontableOutput("table.compMut_clustcolor")),
                                                                           column(12,offset=0,tags$p(h5(strong("",align="left"))))  #空一行             
                                                          ),
                                                          
                                                          #Sample annotations from survival information for oncoprint
                                                          column(12,radioButtons("compMut_annCol", h5(strong("Sample annotations from survival information for oncoprint or not")), inline=TRUE, choices = c("Yes"=1,"No"=2),selected=2)),
                                                          
                                                          conditionalPanel(condition="input.compMut_annCol==1",
                                                                           
                                                                           #The number of sample annotations from survival information for oncoprint
                                                                           column(12,numericInput("compMut_annCol_number",label="The number of sample annotations from survival information for oncoprint",value=3,min=1,step=1)),
                                                                           
                                                                           #Input sample annotations from survival information for oncoprint
                                                                           column(12,offset=0,tags$p(h5(strong("Input sample annotations from survival information for oncoprint",align="left")))),
                                                                           column(12,offset=0,tags$p(h5("First line: Please input the sample annotation variables from survival information",align="left"))),
                                                                           column(12,offset=0,tags$p(h5("Second line: Please input 'Continuous' or 'Categorical' to indicate the type of each sample annotation variable",align="left"))),
                                                                           column(12,offset=0,tags$p(h5("Last line: Please input the colors for each sample annotation variable (use hex color format, e.g. #000004FF and English semicolons should be used to separate the input colors)",align="left"))),
                                                                           column(12,offset=0,tags$p(h5("Note1: If the sample annotation variable is continuous, the number of indicated colors should be equal to 3, which represents the minimum, median and maximum value of this variable)",align="left"))),
                                                                           column(12,offset=0,tags$p(h5("Note2: If the sample annotation variable is categorical, the number of indicated colors should be equal to the number of categories for this variable)",align="left"))),
                                                                           column(10,rHandsontableOutput("table.compMut_annCol")),
                                                                           column(12,offset=0,tags$p(h5(strong("",align="left"))))  #空一行
                                                          ),
                                                          
                                                          #The width of output figure
                                                          column(12,numericInput("width_compMut",label="The width of output figure",value=8,min=0.5,step=0.5)),
                                                          
                                                          #The height of output figure
                                                          column(12,numericInput("height_compMut",label="The height of output figure",value=4,min=0.5,step=0.5)),
                                                          
                                                          #The name of output figure
                                                          conditionalPanel(condition="input.tcga_vali_compMut==1",
                                                                           column(12,textAreaInput(inputId = "compMut_FigureName_tcga",
                                                                                                   label = "Figure Name",
                                                                                                   value = "oncoprint_for_mutations_with_frequency_over_than_the_setting_cutoff(TCGA)",
                                                                                                   width = '100%', height = '35px')
                                                                           )   
                                                          ),
                                                          conditionalPanel(condition="input.tcga_vali_compMut==2",
                                                                           column(12,textAreaInput(inputId = "compMut_FigureName_validation",
                                                                                                   label = "Figure Name",
                                                                                                   value = "oncoprint_for_mutations_with_frequency_over_than_the_setting_cutoff(Validation)",
                                                                                                   width = '100%', height = '35px')
                                                                           ) 
                                                          ),
                                                          
                                                          column(12,actionButton("compMut_process","Process",width="100%",class="btn btn-primary"))
                                                      )
                                            )
                           ),
                           
                           #5.Compare total mutation burden
                           conditionalPanel(condition="input.compModuleStep==4",
                                            wellPanel(style="background-color:#47AEE9",
                                                      fixedRow(
                                                          column(12,tags$h4(strong("Compare total mutation burden"),align="center",style="color:#834005")),
                                                          column(12,offset=0,tags$p(h5(strong("In this step, we will calculate Total Mutation Burden (TMB) and compare them among current subtypes.",align="left")))),
                                                          #Compare total mutation burden on tcga datasets or validation datasets
                                                          column(12,radioButtons("tcga_vali_compTMB", h5(strong("Compare total mutation burden on tcga datasets or validation datasets")), inline=TRUE, choices = c("TCGA"=1,"Validation"=2),selected=1)),
                                                          #Model-free approaches for subtype prediction in validation cohort
                                                          conditionalPanel(condition="input.tcga_vali_compTMB==2",
                                                                           column(12,radioButtons("validation_method_compTMB", h5(strong("Model-free approaches for subtype prediction in validation cohort")), inline=TRUE, choices = c("NTP"=1,"PAM"=2),selected=1))
                                                          ),
                                                          
                                                          #Whether remove repeated variants in a particuar sample, mapped to multiple transcripts of same gene
                                                          column(12,radioButtons("rmDup_compTMB", h5(strong("Whether remove repeated variants in a particuar sample, mapped to multiple transcripts of same gene")), inline=TRUE, choices = c('Yes'=1, 'No'=2),selected=1)),
                                                          
                                                          #Remove possible FLAGS
                                                          column(12,radioButtons("rmFLAGS_compTMB", h5(strong("Remove possible FLAGS (These FLAGS genes are often non-pathogenic and passengers, but are frequently mutated in most of the public exome studies, some of which are fishy. Examples of such genes include TTN, MUC16, etc)")), inline=TRUE, choices = c('Yes'=1, 'No'=2),selected=2)),
                                                          
                                                          #Whether user defined a list of variant classifications that should be considered as non-synonymous and the rest will be considered synonymous
                                                          column(12,radioButtons("nonSyn_compTMB", h5(strong("Whether user defined a list of variant classifications that should be considered as non-synonymous and the rest will be considered synonymous")), inline=TRUE, choices = c('Yes'=1, 'No'=2),selected=2)),
                                                          conditionalPanel(condition="input.nonSyn_compTMB==1",
                                                              
                                                                           #The number of variant classifications that should be considered as non-synonymous
                                                                           column(12,numericInput("nonSyn_compTMB_number",label="The number of variant classifications that should be considered as non-synonymous",value=3,min=1,step=1)),
                                                                           
                                                                           #Input variant classification
                                                                           column(12,offset=0,tags$p(h5(strong("Input variant classification",align="left")))),
                                                                           column(9,rHandsontableOutput("table.nonSyn_compTMB")),
                                                                           column(12,offset=0,tags$p(h5(strong("",align="left"))))  #空一行             
                                                          ),
                                                          
                                                          #The estimation of exome size
                                                          column(12,numericInput("exome_size_compTMB",label="The estimation of exome size",min=0,value=38,step=1)),
                                                          
                                                          #The method for statistical testing
                                                          column(12,radioButtons("test_method_compTMB", h5(strong("The method for statistical testing")), inline=TRUE, choices = c('parametric'=1, 'nonparametric'=2),selected=2)),
                                                          
                                                          #Show the sample size within each subtype at the top of the figure or not
                                                          column(12,radioButtons("show_size_compTMB", h5(strong("Show the sample size within each subtype at the top of the figure or not")), inline=TRUE, choices = c('Yes'=1, 'No'=2),selected=1)),
                                                          
                                                          #Colors for each subtype
                                                          column(12,radioButtons("clustcolor_compTMB", h5(strong("Setting for colors to annotate each subtype")), inline=TRUE, choices = c("System default"=1,"User defined"=2),selected=1)),
                                                          
                                                          conditionalPanel(condition="input.clustcolor_compTMB==2",
                                                                           
                                                                           #Input color for each subtype
                                                                           column(12,offset=0,tags$p(h5(strong("Input color for each subtype (use hex color format, e.g. #2EC4B6)",align="left")))),
                                                                           column(9,rHandsontableOutput("table.compTMB_clustcolor")),
                                                                           column(12,offset=0,tags$p(h5(strong("",align="left"))))  #空一行             
                                                          ),
                                                          
                                                          #The width of boxviolin plot
                                                          column(12,numericInput("width_compTMB",label="The width of boxviolin plot",value=6,min=0.5,step=0.5)),
                                                          
                                                          #The height of boxviolin plot
                                                          column(12,numericInput("height_compTMB",label="The height of boxviolin plot",value=6,min=0.5,step=0.5)),
                                                          
                                                          #The name of output figure
                                                          conditionalPanel(condition="input.tcga_vali_compTMB==1",
                                                                           column(12,textAreaInput(inputId = "compTMB_FigureName_tcga",
                                                                                                   label = "Figure Name",
                                                                                                   value = "distribution_of_TMB_and_titv(TCGA)",
                                                                                                   width = '100%', height = '35px')
                                                                           )
                                                          ),
                                                          conditionalPanel(condition="input.tcga_vali_compTMB==2",
                                                                           column(12,textAreaInput(inputId = "compTMB_FigureName_validation",
                                                                                                   label = "Figure Name",
                                                                                                   value = "distribution_of_TMB_and_titv(Validation)",
                                                                                                   width = '100%', height = '35px')
                                                                           )
                                                          ),
                                                         
                                                          column(12,actionButton("compTMB_process","Process",width="100%",class="btn btn-primary"))
                                                      )
                                            )
                           ),
                           
                           #6.Compare fraction genome altered
                           conditionalPanel(condition="input.compModuleStep==5",
                                            wellPanel(style="background-color:#47AEE9",
                                                      fixedRow(
                                                          column(12,tags$h4(strong("Compare fraction genome altered"),align="center",style="color:#834005")),
                                                          column(12,offset=0,tags$p(h5(strong("This step calculates Fraction Genome Altered (FGA), Fraction Genome Gained (FGG) as well as Fraction Genome Lost (FGL) seperately, and compares them among current subtypes.",align="left")))),
                                                          #Compare fraction genome altered on tcga datasets or validation datasets
                                                          column(12,radioButtons("tcga_vali_compFGA", h5(strong("Compare fraction genome altered on tcga datasets or validation datasets")), inline=TRUE, choices = c("TCGA"=1,"Validation"=2),selected=1)),
                                                          #Model-free approaches for subtype prediction in validation cohort
                                                          conditionalPanel(condition="input.tcga_vali_compFGA==2",
                                                                           column(12,radioButtons("validation_method_compFGA", h5(strong("Model-free approaches for subtype prediction in validation cohort")), inline=TRUE, choices = c("NTP"=1,"PAM"=2),selected=1))
                                                          ),
                                                          
                                                          #The type of the 'value' column in segmented copy number dataset
                                                          column(12,radioButtons("iscopynumber_compFGA", h5(strong("The type of the 'value' column in segmented copy number dataset")), inline=TRUE, choices = c('copy-number value (original)'=1, 'segments value (segment_mean)'=2),selected=2)),
                                                          
                                                          #The cutoff for identifying copy-number gain or loss
                                                          column(12,numericInput("cnathreshold_compFGA",label="The cutoff for identifying copy-number gain or loss",value=0.2,min=0,step=0.05)),
                                                          
                                                          #The method for statistical testing
                                                          column(12,radioButtons("test_method_compFGA", h5(strong("The method for statistical testing")), inline=TRUE, choices = c('parametric'=1, 'nonparametric'=2),selected=2)),
                                                          
                                                          #The mapping colors for bars of FGA, FGG and FGL
                                                          column(12,radioButtons("barcolor_compFGA", h5(strong("Setting for mapping colors for bars of FGA, FGG and FGL")), inline=TRUE, choices = c("System default"=1,"User defined"=2),selected=1)),
                                                          
                                                          conditionalPanel(condition="input.barcolor_compFGA==2",
                                                              
                                                                           #Input colors for bars of FGA, FGG and FGL
                                                                           column(12,offset=0,tags$p(h5(strong("The mapping colors for bars of FGA, FGG and FGL",align="left")))),
                                                                           column(9,rHandsontableOutput("table.barcolor_compFGA")),
                                                                           column(12,offset=0,tags$p(h5(strong("",align="left"))))  #空一行
                                                          ),
                                                         
                                                          #Colors for each subtype
                                                          column(12,radioButtons("clustcolor_compFGA", h5(strong("Setting for colors to annotate each subtype")), inline=TRUE, choices = c("System default"=1,"User defined"=2),selected=1)),
                                                          
                                                          conditionalPanel(condition="input.clustcolor_compFGA==2",
                                                                           
                                                                           #Input color for each subtype
                                                                           column(12,offset=0,tags$p(h5(strong("Input color for each subtype (use hex color format, e.g. #2EC4B6)",align="left")))),
                                                                           column(9,rHandsontableOutput("table.compFGA_clustcolor")),
                                                                           column(12,offset=0,tags$p(h5(strong("",align="left"))))  #空一行             
                                                          ),
                                                          
                                                          #The width of barplot
                                                          column(12,numericInput("width_compFGA",label="The width of barplot",value=8,min=0.5,step=0.5)),
                                                          
                                                          #The height of barplot
                                                          column(12,numericInput("height_compFGA",label="The height of barplot",value=4,min=0.5,step=0.5)),
                                                          
                                                          #The name of the barplot
                                                          conditionalPanel(condition="input.tcga_vali_compFGA==1",
                                                                           column(12,textAreaInput(inputId = "compFGA_FigureName_tcga",
                                                                                                   label = "Figure Name",
                                                                                                   value = "barplot_of_fraction_genome_altered(TCGA)",
                                                                                                   width = '100%', height = '35px')
                                                                           )
                                                          ),
                                                          conditionalPanel(condition="input.tcga_vali_compFGA==2",
                                                                           column(12,textAreaInput(inputId = "compFGA_FigureName_validation",
                                                                                                   label = "Figure Name",
                                                                                                   value = "barplot_of_fraction_genome_altered(Validation)",
                                                                                                   width = '100%', height = '35px')
                                                                           )
                                                          ),
                                                          
                                                          column(12,actionButton("compFGA_process","Process",width="100%",class="btn btn-primary"))
                                                      )
                                            )
                           ),
                           
                           #7.Compare drug sensitivity
                           conditionalPanel(condition="input.compModuleStep==6",
                                            wellPanel(style="background-color:#47AEE9",
                                                      fixedRow(
                                                          column(12,tags$h4(strong("Compare drug sensitivity"),align="center",style="color:#834005")),
                                                          column(12,offset=0,tags$p(h5(strong("This step estimates the IC50 of specific drug for each subtype by developing a ridge regression predictive model based on all/specific cell lines derived from Genomics of Drug Sensitivity in Cancer (GDSC) and compares the IC50 among current subtypes.",align="left")))),
                                                          #Compare drug sensitivity on tcga datasets or validation datasets
                                                          column(12,radioButtons("tcga_vali_compDrugsen", h5(strong("Compare drug sensitivity on tcga datasets or validation datasets")), inline=TRUE, choices = c("TCGA"=1,"Validation"=2),selected=1)),
                                                          #Model-free approaches for subtype prediction in validation cohort
                                                          conditionalPanel(condition="input.tcga_vali_compDrugsen==2",
                                                                           column(12,radioButtons("validation_method_compDrugsen", h5(strong("Model-free approaches for subtype prediction in validation cohort")), inline=TRUE, choices = c("NTP"=1,"PAM"=2),selected=1))
                                                          ),
                                                          
                                                          #The name of the drug from GDSC for which you would like to predict sensitivity
                                                          column(12,textInput(inputId="drugs_compDrugsen",label="The name of the drug from GDSC for which you would like to predict sensitivity",value="Cisplatin")),
                                                          
                                                          #Train the models on only a subset of the CGP cell lines or not
                                                          column(12,radioButtons("tissueType_compDrugsen", h5(strong("Train the models on only a subset of the CGP cell lines or not")), inline=TRUE, choices = c("Yes"=1,"No"=2),selected=2)),
                                                          
                                                          #The tissue type from which the cell lines originated
                                                          conditionalPanel(condition="input.tissueType_compDrugsen==1",
                                                                           
                                                                column(12,radioButtons("tissueTypeValue_compDrugsen", h5(strong("The tissue type from which the cell lines originated")), inline=FALSE, choices = c("aero_digestive_tract"=1,"blood"=2,"bone"=3,"breast"=4,"digestive_system"=5,"lung"=6,"nervous_system"=7,"skin"=8,"urogenital_system"=9),selected=6))
                                                          ),
                                                          
                                                          #The method for statistical testing
                                                          column(12,radioButtons("test_method_compDrugsen", h5(strong("The method for statistical testing")), inline=TRUE, choices = c('parametric'=1, 'nonparametric'=2),selected=2)),
                                                          
                                                          #Colors for each subtype
                                                          column(12,radioButtons("clustcolor_compDrugsen", h5(strong("Setting for colors to annotate each subtype")), inline=TRUE, choices = c("System default"=1,"User defined"=2),selected=1)),
                                                          
                                                          conditionalPanel(condition="input.clustcolor_compDrugsen==2",
                                                                           
                                                                           #Input color for each subtype
                                                                           column(12,offset=0,tags$p(h5(strong("Input color for each subtype (use hex color format, e.g. #2EC4B6)",align="left")))),
                                                                           column(9,rHandsontableOutput("table.compDrugsen_clustcolor")),
                                                                           column(12,offset=0,tags$p(h5(strong("",align="left"))))  #空一行             
                                                          ),
                                                          
                                                          #The seed for reproducing the result of comparing drug sensitivity
                                                          column(12,numericInput("seed_compDrugsen",label="The seed for reproducing the result of comparing drug sensitivity",value=123456,step=1)),
                                                          
                                                          #The width of boxviolin plot
                                                          column(12,numericInput("width_compDrugsen",label="The width of boxviolin plot",value=5,min=0.5,step=0.5)),
                                                          
                                                          #The height of boxviolin plot
                                                          column(12,numericInput("height_compDrugsen",label="The height of boxviolin plot",value=5,min=0.5,step=0.5)),
                                                          
                                                          #The prefix for the name of output boxviolin plot
                                                          conditionalPanel(condition="input.tcga_vali_compDrugsen==1",
                                                                           column(12,textAreaInput(inputId = "compDrugsen_prefix_tcga",
                                                                                                   label = "The prefix for the name of output boxviolin plot (the name of output boxviolin plot will be defined as 'prefix+the name of the indicated drug')",
                                                                                                   value = "boxviolin_of_estimated_IC50(TCGA)",
                                                                                                   width = '100%', height = '35px')
                                                                           )
                                                          ),
                                                          conditionalPanel(condition="input.tcga_vali_compDrugsen==2",
                                                                           column(12,textAreaInput(inputId = "compDrugsen_prefix_validation",
                                                                                                   label = "The prefix for the name of output boxviolin plot (the name of output boxviolin plot will be defined as 'prefix+the name of the indicated drug')",
                                                                                                   value = "boxviolin_of_estimated_IC50(Validation)",
                                                                                                   width = '100%', height = '35px')
                                                                           )
                                                          ),
                                                          
                                                          column(12,actionButton("compDrugsen_process","Process",width="100%",class="btn btn-primary"))
                                                      )
                                            )
                           ),
                           
                           #8.Compare agreement with other subtypes
                           conditionalPanel(condition="input.compModuleStep==7",
                                            wellPanel(style="background-color:#47AEE9",
                                                      fixedRow(
                                                          column(12,tags$h4(strong("Compare agreement with other subtypes"),align="center",style="color:#834005")),
                                                          column(12,offset=0,tags$p(h5(strong("This step aims to compute four evaluation indicators, including Rand Index, Jaccard Index, Fowlkes-Mallows, and Normalized Mutual Information for agreement of two partitions, then generate a barplot and an alluvial diagram for visualization.",align="left")))),
                                                          #Compare agreement with other subtypes on tcga datasets or validation datasets
                                                          column(12,radioButtons("tcga_vali_compAgree", h5(strong("Compare agreement with other subtypes on tcga datasets or validation datasets")), inline=TRUE, choices = c("TCGA"=1,"Validation"=2),selected=1)),
                                                          #Model-free approaches for subtype prediction in validation cohort
                                                          conditionalPanel(condition="input.tcga_vali_compAgree==2",
                                                                           column(12,radioButtons("validation_method_compAgree", h5(strong("Model-free approaches for subtype prediction in validation cohort")), inline=TRUE, choices = c("NTP"=1,"PAM"=2),selected=1))
                                                          ),
                                                          
                                                          #The number of the traditional subtypes for comparison (1-6)
                                                          column(12,numericInput("subt2comp_number_compAgree",label="The number of the traditional subtypes for comparison (1-6)",max=6,min=1,value=1,step=1)),
                                                          
                                                          #Input the variable name of traditional subtypes in survival and clinical information for comparison
                                                          column(12,offset=0,tags$p(h5(strong("Input the variable name of traditional subtypes in survival and clinical information for comparison",align="left")))),
                                                          column(9,rHandsontableOutput("table.subt2comp_compAgree")),
                                                          column(12,offset=0,tags$p(h5(strong("",align="left")))),  #空一行
                                                          
                                                          #Colors for each subtype
                                                          column(12,radioButtons("clustcolor_compAgree", h5(strong("Setting for colors to annotate each subtype")), inline=TRUE, choices = c("System default"=1,"User defined"=2),selected=1)),
                                                          
                                                          conditionalPanel(condition="input.clustcolor_compAgree==2",
                                                                           
                                                                           #Input color for each subtype
                                                                           column(12,offset=0,tags$p(h5(strong("Input color for each subtype (use hex color format, e.g. #2EC4B6)",align="left")))),
                                                                           column(9,rHandsontableOutput("table.compAgree_clustcolor")),
                                                                           column(12,offset=0,tags$p(h5(strong("",align="left"))))  #空一行             
                                                          ),
                                                          
                                                          #The width for box in alluvial diagram
                                                          column(12,numericInput("box_width_compAgree",label="The width for box in alluvial diagram",min=0.05,value=0.1,step=0.05)),
                                                          
                                                          #The width of the figure
                                                          column(12,numericInput("width_compAgree",label="The width of the figure",value=6,min=0.5,step=0.5)),
                                                          
                                                          #The height of the figure
                                                          column(12,numericInput("height_compAgree",label="The height of the figure",value=5,min=0.5,step=0.5)),
                                                          
                                                          #The name of the figure
                                                          conditionalPanel(condition="input.tcga_vali_compAgree==1",
                                                                           column(12,textAreaInput(inputId = "compAgree_FigureName_tcga",
                                                                                                   label = "Figure Name",
                                                                                                   value = "agreement_between_current_subtype_and_other_classifications(TCGA)",
                                                                                                   width = '100%', height = '35px')
                                                                           )
                                                          ),
                                                          conditionalPanel(condition="input.tcga_vali_compAgree==2",
                                                                           column(12,textAreaInput(inputId = "compAgree_FigureName_validation",
                                                                                                   label = "Figure Name",
                                                                                                   value = "agreement_between_current_subtype_and_other_classifications(Validation)",
                                                                                                   width = '100%', height = '35px')
                                                                           )
                                                          ),
                                                          
                                                          column(12,actionButton("compAgree_process","Process",width="100%",class="btn btn-primary"))
                                                      )
                                            )
                           )
                           
                    ),
                    br(),
                    column(7,
                           tabsetPanel(
                               tabPanel(tags$h5(tags$strong("Survival"),align="center",style="color:#FFA500"),br(),
                                        
                                    # column(12,h3(plotOutput("compSurv_figure"))),
                                    uiOutput("compSurv_figure_unit"), #Kaplan-Meier curve
                                    br(),br(),
                                    column(12,uiOutput("compSurv_finish"))    
                               ),
                               tabPanel(tags$h5(tags$strong("Clinical Features"),align="center",style="color:#FFA500"),br(),
                                    
                                    uiOutput("table.compClinvarTitle"), #标题
                                    column(12,br()),
                                    dataTableOutput("table.compClinvar"), #表格
                                    # uiOutput("table.compClinvarNote"), #脚注
                                    br(),br(),
                                    column(12,uiOutput("compClinvar_finish"))    
                               ),
                               tabPanel(tags$h5(tags$strong("Mutational Frequency"),align="center",style="color:#FFA500"),br(),
                                    
                                    # column(12,h3(plotOutput("compMut_figure"))),
                                    uiOutput("compMut_figure_unit"), #oncoprint
                                    br(),br(),
                                    uiOutput("table.compMutTitle"), #标题
                                    column(12,br()),
                                    dataTableOutput("table.compMut"), #表格
                                    # uiOutput("table.compMutNote"), #脚注
                                    br(),br(),
                                    column(12,uiOutput("compMut_finish"))      
                               ),
                               tabPanel(tags$h5(tags$strong("TMB"),align="center",style="color:#FFA500"),br(),
                                    
                                    # column(12,h3(plotOutput("compTMB_figure"))),
                                    uiOutput("compTMB_figure_unit"), #boxviolin plot    
                                    br(),br(),
                                    uiOutput("table.compTMBdatTitle"), #标题
                                    column(12,br()),
                                    dataTableOutput("table.compTMBdat"), #表格
                                    # uiOutput("table.compTMBdatNote"), #脚注
                                    br(),br(),
                                    column(12,uiOutput("compTMB_finish"))
                               ),
                               tabPanel(tags$h5(tags$strong("FGA"),align="center",style="color:#FFA500"),br(),
                                     
                                        # column(12,h3(plotOutput("compFGA_figure"))),
                                        uiOutput("compFGA_figure_unit"), #barplot    
                                        br(),br(),
                                        uiOutput("table.compFGAsummaryTitle"), #标题
                                        column(12,br()),
                                        dataTableOutput("table.compFGAsummary"), #表格
                                        # uiOutput("table.compFGAsummaryNote"), #脚注
                                        br(),br(),
                                        column(12,uiOutput("compFGA_finish"))
                               ),
                               tabPanel(tags$h5(tags$strong("Drug Sensitivity"),align="center",style="color:#FFA500"),br(),
                                        
                                        # column(12,h3(plotOutput("compDrugsen_figure"))),
                                        uiOutput("compDrugsen_figure_unit"), #barplot    
                                        br(),br(),
                                        uiOutput("table.compDrugsenpredictedBoxdatTitle"), #标题
                                        column(12,br()),
                                        dataTableOutput("table.compDrugsenpredictedBoxdat"), #表格
                                        # uiOutput("table.compDrugsenpredictedBoxdatNote"), #脚注
                                        br(),br(),
                                        column(12,uiOutput("compDrugsen_finish")) 
                               ),
                               tabPanel(tags$h5(tags$strong("Agreement"),align="center",style="color:#FFA500"),br(),
                                        
                                        # column(12,h3(plotOutput("compAgree_figure"))),
                                        uiOutput("compAgree_figure_unit"), #barplot+alluvial diagram    
                                        br(),br(),
                                        uiOutput("table.compAgreeTitle"), #标题
                                        column(12,br()),
                                        dataTableOutput("table.compAgree"), #表格
                                        # uiOutput("table.compAgreeNote"), #脚注
                                        br(),br(),
                                        column(12,uiOutput("compAgree_finish"))
                               )
                           )
                    )
                ),
                
                #RUN Module
                tabPanel(
                    title=tags$p(tags$h4(tags$strong("RUN Module"),style="color:purple")),
                    value="RUN Module",
                    
                    column(5,
                           
                           #1.总的步骤选项,对应不同的参数列表
                           wellPanel(style="background-color:#47AEE9",
                            
                                     fixedRow(
                                         #设置步骤选项,每个步骤对应不同的参数
                                         column(12,radioButtons("runModuleStep", h5(strong("Steps")), inline=FALSE,
                                                                choices = c("Run differential expression analysis"=1,"Run biomarker identification procedure"=2,"Run gene set enrichment analysis"=3,"Run gene set variation analysis"=4,"Run nearest template prediction"=5,"Run partition around medoids classifier"=6,"Run consistency evaluation using Kappa statistics"=7),selected=1))
                                     )
                           ),
                           
                           #2.Run differential expression analysis
                           conditionalPanel(condition="input.runModuleStep==1",
                                            wellPanel(style="background-color:#47AEE9",
                                                      fixedRow(
                                                          column(12,tags$h4(strong("Run differential expression analysis"),align="center",style="color:#834005")),
                                                          column(12,offset=0,tags$p(h5(strong("In this step, we will perform differential expression analysis using chosen algorithm (deseq2 or edger or limma) between two classes identified by multi-omics clustering process.",align="left")))),
                                                          
                                                          #Choose the algorithm for differential expression analysis
                                                          column(12,radioButtons("dea_method_runDEA", h5(strong("Choose the algorithm for differential expression analysis")), inline=TRUE, choices = c("deseq2"=1,"edger"=2,"limma"=3),selected=1)),
                                                          
                                                          #Indicate the prefix of output file (e.g. if the name of output file is 'consensusMOIC_nsclc_limma_test_result.CS1_vs_Others.txt', the prefix would be 'nsclc'. If you do not want a prefix, ignore this parameter and type nothing)
                                                          column(12,textAreaInput(inputId = "prefix_runDEA",
                                                                                  label = "Indicate the prefix of output file (e.g. if the name of output file is 'consensusMOIC_nsclc_limma_test_result.CS1_vs_Others.txt', the prefix would be 'nsclc'. If you do not want a prefix, ignore this parameter and type nothing)",
                                                                                  value = "",
                                                                                  width = '100%', height = '35px')
                                                          ),
                                                          
                                                          #The result already exists, overwrite it or skip this step directly
                                                          column(12,radioButtons("overwt_runDEA", h5(strong("The result already exists, overwrite it or skip this step directly")), inline=TRUE, choices = c("Overwrite"=1,"Skip"=2),selected=1)),
                                                          
                                                          #Whether sort adjusted p value in ascending order for output table
                                                          column(12,radioButtons("sort_p_runDEA", h5(strong("Whether sort adjusted p value in ascending order for output table")), inline=TRUE, choices = c("Yes"=1,"No"=2),selected=1)),
                                                          
                                                          #Only select 'id', 'log2fc', 'pvalue' and 'padj' columns for output table or not
                                                          column(12,radioButtons("verbose_runDEA", h5(strong("Only select 'id', 'log2fc', 'pvalue' and 'padj' columns for output table or not")), inline=TRUE, choices = c("Yes"=1,"No"=2),selected=1)),
                                                          
                                                          column(12,actionButton("runDEA_process","Process",width="100%",class="btn btn-primary"))
                                                      )
                                            )
                           ),
                           
                           #3.Run biomarker identification procedure
                           conditionalPanel(condition="input.runModuleStep==2",
                                            wellPanel(style="background-color:#47AEE9",
                                                      fixedRow(
                                                          column(12,tags$h4(strong("Run biomarker identification procedure"),align="center",style="color:#834005")),
                                                          column(12,offset=0,tags$p(h5(strong("This step aims to identify uniquely and significantly expressed (overexpressed or downexpressed) biomarkers for each subtype identified by multi-omics clustering process. A template including top markers will be genearated for subtype external verification and a heatmap will also be generated.",align="left")))),
                                                          
                                                          #Indicate the algorithm for completed differential expression analysis
                                                          column(12,radioButtons("dea_method_runMarker", h5(strong("Indicate the algorithm for completed differential expression analysis")), inline=TRUE, choices = c("deseq2"=1,"edger"=2,"limma"=3),selected=1)),
                                                          
                                                          #Indicate the nominal p value for identifying significant markers
                                                          column(12,numericInput("p_cutoff_runMarker",label="Indicate the nominal p value for identifying significant markers",max=1,min=0,value=0.05,step=0.05)),
                                                          
                                                          #Indicate the adjusted p value for identifying significant markers
                                                          column(12,numericInput("p_adj_cutoff_runMarker",label="Indicate the adjusted p value for identifying significant markers",max=1,min=0,value=0.05,step=0.05)),
                                                          
                                                          #Indicate the direction of identifying significant marker
                                                          column(12,radioButtons("dirct_runMarker", h5(strong("Indicate the direction of identifying significant marker")), inline=TRUE, choices = c("up-regulated"=1,"down-regulated"=2),selected=1)),
                                                          
                                                          #Indicate the number of top markers sorted by log2fc should be identified for each subtype
                                                          column(12,numericInput("n_marker_runMarker",label="Indicate the number of top markers sorted by log2fc should be identified for each subtype",min=0,value=200,step=1)),
                                                          
                                                          ##Sample annotations from survival information for heatmap
                                                          ###Sample annotations from survival information for heatmap or not
                                                          column(12,radioButtons("runMarker_annCol", h5(strong("Sample annotations from survival information for heatmap or not")), inline=TRUE, choices = c("Yes"=1,"No"=2),selected=2)),
                                                          
                                                          conditionalPanel(condition="input.runMarker_annCol==1",
                                                                           
                                                                           #The number of sample annotations from survival information for heatmap
                                                                           column(12,numericInput("runMarker_annCol_number",label="The number of sample annotations from survival information for heatmap",value=3,min=1,step=1)),
                                                                           
                                                                           #Input sample annotations from survival information for heatmap
                                                                           column(12,offset=0,tags$p(h5(strong("Input sample annotations from survival information for heatmap",align="left")))),
                                                                           column(12,offset=0,tags$p(h5("First line: Please input the sample annotation variables from survival information",align="left"))),
                                                                           column(12,offset=0,tags$p(h5("Second line: Please input 'Continuous' or 'Categorical' to indicate the type of each sample annotation variable",align="left"))),
                                                                           column(12,offset=0,tags$p(h5("Last line: Please input the colors for each sample annotation variable (use hex color format, e.g. #000004FF and English semicolons should be used to separate the input colors)",align="left"))),
                                                                           column(12,offset=0,tags$p(h5("Note1: If the sample annotation variable is continuous, the number of indicated colors should be equal to 3, which represents the minimum, median and maximum value of this variable)",align="left"))),
                                                                           column(12,offset=0,tags$p(h5("Note2: If the sample annotation variable is categorical, the number of indicated colors should be equal to the number of categories for this variable)",align="left"))),
                                                                           column(10,rHandsontableOutput("table.runMarker_annCol")),
                                                                           column(12,offset=0,tags$p(h5(strong("",align="left"))))  #空一行
                                                          ),
                                                          
                                                          ##Colors for annotating each cluster at the top of heatmap
                                                          column(12,radioButtons("runMarker_clustcolor", h5(strong("Colors for annotating each cluster at the top of heatmap")), inline=TRUE, choices = c("System default"=1,"User defined"=2),selected=1)),
                                                          conditionalPanel(condition="input.runMarker_clustcolor==2",
                                                                           column(12,offset=0,tags$p(h5(strong("Input colors for annotating each cluster (use hex color format, e.g. #2EC4B6)",align="left")))),
                                                                           column(9,rHandsontableOutput("table.runMarker_clustcolor")),
                                                                           column(12,offset=0,tags$p(h5(strong("",align="left"))))  #空一行             
                                                          ),
                                                          
                                                          ##Indicate if expression data should be centered
                                                          column(12,radioButtons("runMarker_centerFlag", h5(strong("Indicate if expression data should be centered")), inline=TRUE, choices = c("Yes"=1,"No"=2),selected=1)),
                                                         
                                                          ##Indicate if expression data should be scaled
                                                          column(12,radioButtons("runMarker_scaleFlag", h5(strong("Indicate if expression data should be scaled")), inline=TRUE, choices = c("Yes"=1,"No"=2),selected=1)),
                                                          
                                                          #Assign marginal cutoff for truncating values in data or not
                                                          column(12,radioButtons("halfwidth_runMarker", h5(strong("Assign marginal cutoff for truncating values in data or not")), inline=TRUE, choices = c("Yes"=1,"No"=2),selected=1)),
                                                          
                                                          #Marginal cutoff for truncating values in data
                                                          column(12,numericInput("halfwidth_value_runMarker",label="Marginal cutoff for truncating values in data",min=0.5,value=3,step=0.5)),
                                                          
                                                          ##Show rownames (feature names) in heatmap or not
                                                          column(12,radioButtons("runMarker_show_rownames", h5(strong("Show rownames (feature names) in heatmap or not")), inline=TRUE, choices = c("Yes"=1,"No"=2),selected=2)),
                                                          
                                                          ##Show colnames (sample ID) in heatmap or not
                                                          column(12,radioButtons("runMarker_show_colnames", h5(strong("Show colnames (sample ID) in heatmap or not")), inline=TRUE, choices = c("Yes"=1,"No"=2),selected=2)),
                                                          
                                                          ##Colors for heatmap
                                                          column(12,radioButtons("runMarker_heatmap_color", h5(strong("Colors for heatmap")), inline=TRUE, choices = c("System default"=1,"User defined"=2),selected=1)),
                                                          conditionalPanel(condition="input.runMarker_heatmap_color==2",
                                                                           
                                                                           #The number of colors for heatmap
                                                                           column(12,numericInput("runMarker_heatmap_color_number",label="The number of colors for heatmap",value=3,min=3,step=1)),
                                                                           
                                                                           #Input colors for heatmap
                                                                           column(12,offset=0,tags$p(h5(strong("Input colors for heatmap (use hex color format, e.g. #00FF00)",align="left")))),
                                                                           column(9,rHandsontableOutput("table.runMarker_heatmap_color")),
                                                                           column(12,offset=0,tags$p(h5(strong("",align="left"))))  #空一行             
                                                          ),
                                                          
                                                          #The width of output figure
                                                          column(12,numericInput("width_runMarker",label="The width of output figure",value=8,min=0.5,step=0.5)),
                                                          
                                                          #The height of output figure
                                                          column(12,numericInput("height_runMarker",label="The height of output figure",value=8,min=0.5,step=0.5)),
                                                          
                                                          #Indicate the prefix of output figure (e.g. if the name of output figure is 'markerheatmap_using_upregulated_genes.pdf', the prefix would be 'markerheatmap'. If you do not want a prefix but using the system default prefix 'markerheatmap', ignore this parameter and type nothing)
                                                          column(12,textAreaInput(inputId = "runMarker_FigureName",
                                                                                  label = "Indicate the prefix of output figure (e.g. if the name of output figure is 'markerheatmap_using_upregulated_genes.pdf', the prefix would be 'markerheatmap'. If you do not want a prefix but using the system default prefix 'markerheatmap', ignore this parameter and type nothing)",
                                                                                  value = "",
                                                                                  width = '100%', height = '35px')
                                                          ),
                                                          
                                                          column(12,actionButton("runMarker_process","Process",width="100%",class="btn btn-primary"))
                                                      )
                                            )
                           ),
                           
                           #4.Run gene set enrichment analysis
                           conditionalPanel(condition="input.runModuleStep==3",
                                            
                                            ##1.Prepare a gene set background file
                                            wellPanel(style="background-color:#47AEE9",
                                                      fixedRow(
                                                          column(12,tags$h4(strong("Prepare a gene set background file"),align="center",style="color:#834005")),
                                                          column(12,offset=0,tags$p(h5(strong("First we need to prepare a gene set background file for gene set enrichment analysis.",align="left")))),
                                                          
                                                          #The manner for gene set background file preparation
                                                          column(12,radioButtons("gsea_geneset_background_preparation", h5(strong("The manner for gene set background file preparation")), inline=TRUE,choices = c("System default"=1,"Download from specified url"=2,"Upload manually"=3),selected=1)),
                                                          
                                                          conditionalPanel(condition="input.gsea_geneset_background_preparation==2",
                                                                          
                                                                                            column(12,textAreaInput(inputId = "gsea_geneset_background_preparation_downloadUrl",
                                                                                                                    label = "Download url for gene set background file",
                                                                                                                    value = "",
                                                                                                                    width = '100%', height = '35px')),
                                                                           
                                                                           column(12,actionButton("gsea_geneset_background_preparation_process2","Process",width="100%",class="btn btn-primary"))
                                                          ),
                                                          conditionalPanel(condition="input.gsea_geneset_background_preparation==3",
                                                                           
                                                                           column(12,fileInput('gsea_geneset_background_preparation_fileUpload',tags$p(tags$b("Upload the gene set background file")),
                                                                                               accept=c('text/csv', 'text/comma-separated-values,text/plain', '.csv', '.xls', '.xlsx', '.gmt', '.gz', '.tar.gz', '.rar', '.zip', '.json') )),
                                                                          
                                                                           column(12,actionButton("gsea_geneset_background_preparation_process3","Process",width="100%",class="btn btn-primary"))
                                                          )
                                                          
                                                      )
                                            ),
                                            
                                            ##2.Perform gene set enrichment analysis
                                            wellPanel(style="background-color:#47AEE9",
                                                      fixedRow(
                                                          column(12,tags$h4(strong("Run gene set enrichment analysis"),align="center",style="color:#834005")),
                                                          column(12,offset=0,tags$p(h5(strong("This step aims to perform gene set enrichment analysis using a background file to identify subtype-specific (overexpressed or downexpressed) functional pathways for each subtype.",align="left")))),
                                                          
                                                          #Indicate the algorithm for completed differential expression analysis
                                                          column(12,radioButtons("dea_method_runGSEA", h5(strong("Indicate the algorithm for completed differential expression analysis")), inline=TRUE, choices = c("deseq2"=1,"edger"=2,"limma"=3),selected=1)),
                                                          
                                                          #Indicate the direction of identifying significant pathway
                                                          column(12,radioButtons("dirct_runGSEA", h5(strong("Indicate the direction of identifying significant pathway")), inline=TRUE, choices = c("up-regulated"=1,"down-regulated"=2),selected=1)),
                                                          
                                                          #Indicate the number of top pathways sorted by NES should be identified for each subtype
                                                          column(12,numericInput("n_path_runGSEA",label="Indicate the number of top pathways sorted by NES should be identified for each subtype",min=0,value=10,step=1)),
                                                          
                                                          #Indicate the number of permutations for gene set enrichment analysis (1000 by default and 10000 will be better for reproducibility)
                                                          column(12,numericInput("nPerm_runGSEA",label="Indicate the number of permutations for gene set enrichment analysis (1000 by default and 10000 will be better for reproducibility)",min=0,value=1000,step=1)),
                                                          
                                                          #Indicate minimal size of each gene set for analysis
                                                          column(12,numericInput("minGSSize_runGSEA",label="Indicate minimal size of each gene set for analysis",min=0,value=10,step=1)),
                                                          
                                                          #Indicate maximal size of each gene set for analysis
                                                          column(12,numericInput("maxGSSize_runGSEA",label="Indicate maximal size of each gene set for analysis",min=0,value=500,step=1)),
                                                          
                                                          #Indicate the nominal p value for identifying significant pathways
                                                          column(12,numericInput("p_cutoff_runGSEA",label="Indicate the nominal p value for identifying significant pathways",max=1,min=0,value=0.05,step=0.05)),
                                                          
                                                          #Indicate the adjusted p value for identifying significant pathways
                                                          column(12,numericInput("p_adj_cutoff_runGSEA",label="Indicate the adjusted p value for identifying significant pathways",max=1,min=0,value=0.25,step=0.05)),
                                                          
                                                          #Indicate the method to employ in the estimation of gene set enrichment scores per sample
                                                          column(12,radioButtons("gsva_method_runGSEA", h5(strong("Indicate the method to employ in the estimation of gene set enrichment scores per sample")), inline=TRUE, choices = c("gsva"=1,"ssgsea"=2,"zscore"=3,"plage"=4),selected=1)),
                                                          
                                                          #Indicate the method to calculate subtype-specific pathway enrichment scores
                                                          column(12,radioButtons("norm_method_runGSEA", h5(strong("Indicate the method to calculate subtype-specific pathway enrichment scores")), inline=TRUE, choices = c("mean"=1,"median"=2),selected=1)),
                                                          
                                                          ##Colors for annotating each cluster at the top of heatmap
                                                          column(12,radioButtons("runGSEA_clustcolor", h5(strong("Colors for annotating each cluster at the top of heatmap")), inline=TRUE, choices = c("System default"=1,"User defined"=2),selected=1)),
                                                          conditionalPanel(condition="input.runGSEA_clustcolor==2",
                                                                           column(12,offset=0,tags$p(h5(strong("Input colors for annotating each cluster (use hex color format, e.g. #2EC4B6)",align="left")))),
                                                                           column(9,rHandsontableOutput("table.runGSEA_clustcolor")),
                                                                           column(12,offset=0,tags$p(h5(strong("",align="left"))))  #空一行             
                                                          ),
                                                          
                                                          ##Colors for heatmap
                                                          column(12,radioButtons("runGSEA_heatmap_color", h5(strong("Colors for heatmap")), inline=TRUE, choices = c("System default"=1,"User defined"=2),selected=1)),
                                                          conditionalPanel(condition="input.runGSEA_heatmap_color==2",
                                                                           
                                                                           #Input colors for heatmap
                                                                           column(12,offset=0,tags$p(h5(strong("Input colors for heatmap (use hex color format, e.g. #0000FF)",align="left")))),
                                                                           column(9,rHandsontableOutput("table.runGSEA_heatmap_color")),
                                                                           column(12,offset=0,tags$p(h5(strong("",align="left"))))  #空一行             
                                                          ),
                                                          
                                                          ##The width of output figure
                                                          column(12,numericInput("width_runGSEA",label="The width of output figure",value=15,min=0.5,step=0.5)),
                                                          
                                                          ##The height of output figure
                                                          column(12,numericInput("height_runGSEA",label="The height of output figure",value=10,min=0.5,step=0.5)),
                                                          
                                                          ##Indicate the prefix of output figure (e.g. if the name of output figure is 'gseaheatmap_using_upregulated_pathways.pdf', the prefix would be 'gseaheatmap'. If you do not want a prefix but using the system default prefix 'gseaheatmap', ignore this parameter and type nothing)
                                                          column(12,textAreaInput(inputId = "runGSEA_FigureName",
                                                                                  label = "Indicate the prefix of output figure (e.g. if the name of output figure is 'gseaheatmap_using_upregulated_pathways.pdf', the prefix would be 'gseaheatmap'. If you do not want a prefix but using the system default prefix 'gseaheatmap', ignore this parameter and type nothing)",
                                                                                  value = "",
                                                                                  width = '100%', height = '35px')
                                                          ),
                                                          
                                                          column(12,actionButton("runGSEA_process","Process",width="100%",class="btn btn-primary"))
                                                      )
                                            )
                           ),
                           
                           #5.Run gene set variation analysis
                           conditionalPanel(condition="input.runModuleStep==4",
                                            
                                            ##1.Prepare a gene set list of interest
                                            wellPanel(style="background-color:#47AEE9",
                                                      fixedRow(
                                                          column(12,tags$h4(strong("Prepare a gene set list of interest"),align="center",style="color:#834005")),
                                                          column(12,offset=0,tags$p(h5(strong("First we need to prepare a gene set list of interest for gene set variation analysis.",align="left")))),
                                                          
                                                          #The manner for the gene set list of interest preparation
                                                          column(12,radioButtons("gsva_geneset_interest_preparation", h5(strong("The manner for the gene set list of interest preparation")), inline=TRUE,choices = c("System default"=1,"Download from specified url"=2,"Upload manually"=3),selected=1)),
                                                          
                                                          conditionalPanel(condition="input.gsva_geneset_interest_preparation==2",
                                                                           
                                                                           column(12,textAreaInput(inputId = "gsva_geneset_interest_preparation_downloadUrl",
                                                                                                   label = "Download url for the gene set list of interest",
                                                                                                   value = "",
                                                                                                   width = '100%', height = '35px')),
                                                                           
                                                                           column(12,actionButton("gsva_geneset_interest_preparation_process2","Process",width="100%",class="btn btn-primary"))
                                                          ),
                                                          conditionalPanel(condition="input.gsva_geneset_interest_preparation==3",
                                                                           
                                                                           column(12,fileInput('gsva_geneset_interest_preparation_fileUpload',tags$p(tags$b("Upload the gene set list of interest")),
                                                                                               accept=c('text/csv', 'text/comma-separated-values,text/plain', '.csv', '.xls', '.xlsx', '.gmt', '.gz', '.tar.gz', '.rar', '.zip', '.json') )),
                                                                           
                                                                           column(12,actionButton("gsva_geneset_interest_preparation_process3","Process",width="100%",class="btn btn-primary"))
                                                          )
                                                      )
                                            ),
                                            
                                            ##2.Perform gene set variation analysis
                                            wellPanel(style="background-color:#47AEE9",
                                                      fixedRow(
                                                          column(12,tags$h4(strong("Run gene set variation analysis"),align="center",style="color:#834005")),
                                                          column(12,offset=0,tags$p(h5(strong("This step aims to use gene set variation analysis to calculate enrichment score of each sample in each subtype based on a given gene set list of interest.",align="left")))),
                                                          
                                                          #Indicate the method to employ in the estimation of gene set enrichment scores per sample
                                                          column(12,radioButtons("gsva_method_runGSVA", h5(strong("Indicate the method to employ in the estimation of gene set enrichment scores per sample")), inline=TRUE, choices = c("gsva"=1,"ssgsea"=2,"zscore"=3,"plage"=4),selected=1)),
                                                          
                                                          #Indicate if enrichment scores should be centered or not
                                                          column(12,radioButtons("centerFlag_runGSVA", h5(strong("Indicate if enrichment scores should be centered or not")), inline=TRUE, choices = c("Yes"=1,"No"=2),selected=1)),
                                                          
                                                          #Indicate if enrichment scores should be scaled or not
                                                          column(12,radioButtons("scaleFlag_runGSVA", h5(strong("Indicate if enrichment scores should be scaled or not")), inline=TRUE, choices = c("Yes"=1,"No"=2),selected=1)),
                                                          
                                                          #Assign marginal cutoff for truncating enrichment scores or not
                                                          column(12,radioButtons("halfwidth_runGSVA", h5(strong("Assign marginal cutoff for truncating enrichment scores or not")), inline=TRUE, choices = c("Yes"=1,"No"=2),selected=1)),
                                                          
                                                          #Marginal cutoff for truncating enrichment scores
                                                          column(12,numericInput("halfwidth_value_runGSVA",label="Marginal cutoff for truncating enrichment scores",min=0.5,value=1,step=0.5)),
                                                          
                                                          ##Sample annotations from survival information for heatmap
                                                          ###Sample annotations from survival information for heatmap or not
                                                          column(12,radioButtons("runGSVA_annCol", h5(strong("Sample annotations from survival information for heatmap or not")), inline=TRUE, choices = c("Yes"=1,"No"=2),selected=2)),
                                                          
                                                          conditionalPanel(condition="input.runGSVA_annCol==1",
                                                                           
                                                                           #The number of sample annotations from survival information for heatmap
                                                                           column(12,numericInput("runGSVA_annCol_number",label="The number of sample annotations from survival information for heatmap",value=3,min=1,step=1)),
                                                                           
                                                                           #Input sample annotations from survival information for heatmap
                                                                           column(12,offset=0,tags$p(h5(strong("Input sample annotations from survival information for heatmap",align="left")))),
                                                                           column(12,offset=0,tags$p(h5("First line: Please input the sample annotation variables from survival information",align="left"))),
                                                                           column(12,offset=0,tags$p(h5("Second line: Please input 'Continuous' or 'Categorical' to indicate the type of each sample annotation variable",align="left"))),
                                                                           column(12,offset=0,tags$p(h5("Last line: Please input the colors for each sample annotation variable (use hex color format, e.g. #000004FF and English semicolons should be used to separate the input colors)",align="left"))),
                                                                           column(12,offset=0,tags$p(h5("Note1: If the sample annotation variable is continuous, the number of indicated colors should be equal to 3, which represents the minimum, median and maximum value of this variable)",align="left"))),
                                                                           column(12,offset=0,tags$p(h5("Note2: If the sample annotation variable is categorical, the number of indicated colors should be equal to the number of categories for this variable)",align="left"))),
                                                                           column(10,rHandsontableOutput("table.runGSVA_annCol")),
                                                                           column(12,offset=0,tags$p(h5(strong("",align="left"))))  #空一行
                                                          ),
                                                          
                                                          ##Colors for annotating each cluster at the top of heatmap
                                                          column(12,radioButtons("runGSVA_clustcolor", h5(strong("Colors for annotating each cluster at the top of heatmap")), inline=TRUE, choices = c("System default"=1,"User defined"=2),selected=1)),
                                                          conditionalPanel(condition="input.runGSVA_clustcolor==2",
                                                                           column(12,offset=0,tags$p(h5(strong("Input colors for annotating each cluster (use hex color format, e.g. #2EC4B6)",align="left")))),
                                                                           column(9,rHandsontableOutput("table.runGSVA_clustcolor")),
                                                                           column(12,offset=0,tags$p(h5(strong("",align="left"))))  #空一行             
                                                          ),
                                                          
                                                          ##Colors for heatmap
                                                          column(12,radioButtons("runGSVA_heatmap_color", h5(strong("Colors for heatmap")), inline=TRUE, choices = c("System default"=1,"User defined"=2),selected=1)),
                                                          conditionalPanel(condition="input.runGSVA_heatmap_color==2",
                                                                           
                                                                           #The number of colors for heatmap
                                                                           column(12,numericInput("runGSVA_heatmap_color_number",label="The number of colors for heatmap",value=5,min=3,step=1)),
                                                                           
                                                                           #Input colors for heatmap
                                                                           column(12,offset=0,tags$p(h5(strong("Input colors for heatmap (use hex color format, e.g. #00FF00)",align="left")))),
                                                                           column(9,rHandsontableOutput("table.runGSVA_heatmap_color")),
                                                                           column(12,offset=0,tags$p(h5(strong("",align="left"))))  #空一行             
                                                          ),
                                                          
                                                          #Distance measurement for hierarchical clustering
                                                          column(12,textInput(inputId="runGSVA_distance",label="Distance measurement for hierarchical clustering",value="euclidean")),
                                                          
                                                          #Clustering method for hierarchical clustering
                                                          column(12,textInput(inputId="runGSVA_linkage",label="Clustering method for hierarchical clustering",value="ward.D")),
                                                          
                                                          #Show rownames (feature names) in heatmap or not
                                                          column(12,radioButtons("show_rownames_runGSVA", h5(strong("Show rownames (feature names) in heatmap or not")), inline=TRUE, choices = c("Yes"=1,"No"=2),selected=1)),
                                                          
                                                          #Show colnames (sample ID) in heatmap or not
                                                          column(12,radioButtons("show_colnames_runGSVA", h5(strong("Show colnames (sample ID) in heatmap or not")), inline=TRUE, choices = c("Yes"=1,"No"=2),selected=2)),
                                                          
                                                          #The width of output figure
                                                          column(12,numericInput("width_runGSVA",label="The width of output figure",value=8,min=0.5,step=0.5)),
                                                          
                                                          #The height of output figure
                                                          column(12,numericInput("height_runGSVA",label="The height of output figure",value=8,min=0.5,step=0.5)),
                                                          
                                                          #Indicate the prefix of output figure (e.g. if the name of output figure is 'enrichment_heatmap_using_gsva.pdf', the prefix would be 'enrichment_heatmap_using'. If you do not want a prefix but using the system default prefix 'enrichment_heatmap_using', ignore this parameter and type nothing)
                                                          column(12,textAreaInput(inputId = "runGSVA_FigureName",
                                                                                  label = "Indicate the prefix of output figure (e.g. if the name of output figure is 'enrichment_heatmap_using_gsva.pdf', the prefix would be 'enrichment_heatmap_using'. If you do not want a prefix but using the system default prefix 'enrichment_heatmap_using', ignore this parameter and type nothing)",
                                                                                  value = "",
                                                                                  width = '100%', height = '35px')
                                                          ),
                                                          
                                                          column(12,actionButton("runGSVA_process","Process",width="100%",class="btn btn-primary"))
                                                      )
                                            )
                           ),
                           
                           #6.Run nearest template prediction
                           conditionalPanel(condition="input.runModuleStep==5",
                                            
                                            wellPanel(style="background-color:#47AEE9",
                                                      fixedRow(
                                                          column(12,tags$h4(strong("Run nearest template prediction"),align="center",style="color:#834005")),
                                                          column(12,offset=0,tags$p(h5(strong("This step aims to assign potential subtype labels on tcga or validation cohort using Nearest Template Prediction (NTP) based on predefined templates derived from current identified subtypes.",align="left")))),
                                                          
                                                          #Run nearest template prediction on tcga or validation cohort
                                                          column(12,radioButtons("tcga_vali_runNTP", h5(strong("Run nearest template prediction on tcga or validation cohort")), inline=TRUE, choices = c("TCGA"=1,"Validation"=2),selected=2)),
                                                          
                                                          #Choose template from 'Run biomarker identification procedure'
                                                          column(12,offset=0,tags$p(h5(strong("Choose template from 'Run biomarker identification procedure':",align="left")))),
                                                          
                                                          #Indicate the algorithm
                                                          column(12,radioButtons("dea_method_runNTP", h5(strong("Indicate the algorithm")), inline=TRUE, choices = c("deseq2"=1,"edger"=2,"limma"=3),selected=2)),
                                                          
                                                          #Indicate the direction
                                                          column(12,radioButtons("dirct_runNTP", h5(strong("Indicate the direction")), inline=TRUE, choices = c("up-regulated"=1,"down-regulated"=2),selected=1)),
                                                          
                                                          #Indicate if the expression data should be further scaled
                                                          column(12,radioButtons("scaleFlag_runNTP", h5(strong("Indicate if the expression data should be further scaled")), inline=TRUE, choices = c("Yes"=1,"No"=2),selected=1)),
                                                          
                                                          #Indicate if the expression data should be further centered
                                                          column(12,radioButtons("centerFlag_runNTP", h5(strong("Indicate if the expression data should be further centered")), inline=TRUE, choices = c("Yes"=1,"No"=2),selected=1)),
                                                          
                                                          #Indicate the permutations for p-value estimation
                                                          column(12,numericInput("nPerm_runNTP",label="Indicate the permutations for p-value estimation",min=0,value=1000,step=1)),
                                                          
                                                          #Indicate the distance measurement
                                                          column(12,radioButtons("distance_runNTP", h5(strong("Indicate the distance measurement")), inline=TRUE, choices = c("cosine"=1,"pearson"=2,"spearman"=3,"kendall"=4),selected=1)),
                                                          
                                                          #Input an integer value for p-value reproducibility
                                                          column(12,numericInput("seed_runNTP",label="Input an integer value for p-value reproducibility",value=123456,step=1)),
                                                          
                                                          #The width of output figure
                                                          column(12,numericInput("width_runNTP",label="The width of output figure",value=5,min=0.5,step=0.5)),
                                                          
                                                          #The height of output figure
                                                          column(12,numericInput("height_runNTP",label="The height of output figure",value=5,min=0.5,step=0.5)),
                                                          
                                                          #The name of the nearest template prediction heatmap
                                                          conditionalPanel(condition="input.tcga_vali_runNTP==1",
                                                                           
                                                                           column(12,textAreaInput(inputId = "runNTP_FigureName_tcga",
                                                                                                   label = "The name of the nearest template prediction heatmap",
                                                                                                   value = "ntpheatmap(TCGA)",
                                                                                                   width = '100%', height = '35px')
                                                                           )
                                                          ),
                                                          conditionalPanel(condition="input.tcga_vali_runNTP==2",
                                                                           
                                                                           column(12,textAreaInput(inputId = "runNTP_FigureName_validation",
                                                                                                   label = "The name of the nearest template prediction heatmap",
                                                                                                   value = "ntpheatmap(Validation)",
                                                                                                   width = '100%', height = '35px')
                                                                           )
                                                          ),
                                                          
                                                          column(12,actionButton("runNTP_process","Process",width="100%",class="btn btn-primary"))
                                                      )
                                            )
                           ),
                           
                           #7.Run partition around medoids classifier
                           conditionalPanel(condition="input.runModuleStep==6",
                                            
                                            wellPanel(style="background-color:#47AEE9",
                                                      fixedRow(
                                                          column(12,tags$h4(strong("Run partition around medoids classifier"),align="center",style="color:#834005")),
                                                          column(12,offset=0,tags$p(h5(strong("This step aims to use partition around medoids (PAM) classifier to predict potential subtype labels on tcga or validation cohort and calculate in-group proportions (IGP) statistics.",align="left")))),
                                                          
                                                          #Run partition around medoids classifier on tcga or validation cohort
                                                          column(12,radioButtons("tcga_vali_runPAM", h5(strong("Run partition around medoids classifier on tcga or validation cohort")), inline=TRUE, choices = c("TCGA"=1,"Validation"=2),selected=2)),
                                                          
                                                          #Indicate action for NA values in normalized expression training data
                                                          column(12,radioButtons("runPAM_training_na", h5(strong("Indicate action for NA values in normalized expression training data")), inline=TRUE, choices = c("Remove directly"=1,"KNN imputation"=2),selected=1)),
                                                          
                                                          #Indicate action for NA values in normalized expression testing data
                                                          column(12,radioButtons("runPAM_testing_na", h5(strong("Indicate action for NA values in normalized expression testing data")), inline=TRUE, choices = c("Remove directly"=1,"KNN imputation"=2),selected=1)),
                                                          
                                                          #Whether indicate a subset of genes to be used
                                                          column(12,radioButtons("gene_subset_runPAM", h5(strong("Whether indicate a subset of genes to be used")), inline=TRUE, choices = c("Yes"=1,"No"=2),selected=2)),
                                                          
                                                          #Indicate a subset of genes to be used (Use english commas to separate gene names)
                                                          conditionalPanel(condition="input.gene_subset_runPAM==1",
                                                                
                                                                           column(12,textAreaInput(inputId = "gene_subset_list_runPAM",
                                                                                                   label = "Indicate a subset of genes to be used (use english commas to separate gene names)",
                                                                                                   value = "",
                                                                                                   width = '100%', height = '35px')
                                                                           )            
                                                          ),
                                                          
                                                          column(12,actionButton("runPAM_process","Process",width="100%",class="btn btn-primary"))
                                                      )
                                            )
                           ),
                           
                           #8.Run consistency evaluation using Kappa statistics
                           conditionalPanel(condition="input.runModuleStep==7",
                                            
                                            #1.Choose the results you have obtained on tcga or validation cohort using NTP or PAM subtype prediction method
                                            wellPanel(style="background-color:#47AEE9",
                                                      fixedRow(
                                                          column(12,tags$h4(strong("Run consistency evaluation using Kappa statistics"),align="center",style="color:#834005")),
                                                          column(12,offset=0,tags$p(h5(strong("This step aims to calculate Kappa statistic to measure the consistency between two appraisements.",align="left")))),
                                                          # column(12,offset=0,tags$p(h5(strong("First let's choose the results you have obtained on tcga or validation cohort using NTP or PAM subtype prediction method:",align="left")))),
                                                          
                                                          #Choose the results you have obtained in 'Run nearest template prediction' and 'Run partition around medoids classifier' procedures
                                                          column(12,checkboxGroupInput(inputId = "ntp_pam_runKappa",
                                                                                       label = "Choose the results you have obtained in 'Run nearest template prediction' and 'Run partition around medoids classifier' procedures",
                                                                                       choices = list("Run nearest template prediction on tcga cohort"=1,"Run partition around medoids classifier on tcga cohort"=2,"Run nearest template prediction on validation cohort"=3,"Run partition around medoids classifier on validation cohort"=4)
                                                          ))
                                                      )
                                            ),
                                            
                                            #2.Run consistency evaluation between current subtypes derived from multi-omics clustering and NTP-predicted subtypes on tcga cohort
                                            conditionalPanel(condition="input.ntp_pam_runKappa.includes('1')",
                                                
                                                wellPanel(style="background-color:#47AEE9",
                                                        fixedRow(
                                                            column(12,tags$h4(strong("Run Kappa on tcga cohort (CMOIC vs NTP)"),align="center",style="color:#834005")),
                                                            
                                                            #Indicate the label of the first subtype
                                                            column(12,textInput(inputId="subt1_lab_runKappa1",label="Indicate the label of the first subtype",value="CMOIC_TCGA")),
                                                            
                                                            #Indicate the label of the second subtype
                                                            column(12,textInput(inputId="subt2_lab_runKappa1",label="Indicate the label of the second subtype",value="NTP_TCGA")),
                                                            
                                                            #The width of output figure
                                                            column(12,numericInput("width_runKappa1",label="The width of output figure",value=5,min=0.5,step=0.5)),
                                                            
                                                            #The height of output figure
                                                            column(12,numericInput("height_runKappa1",label="The height of output figure",value=5,min=0.5,step=0.5)),
                                                            
                                                            #The name of the consistency heatmap
                                                            column(12,textAreaInput(inputId = "runKappa1_FigureName",
                                                                                    label = "The name of the consistency heatmap",
                                                                                    value = "constheatmap(CMOIC_VS_NTP_TCGA)",
                                                                                    width = '100%', height = '35px')
                                                            ),
                                                            
                                                            column(12,actionButton("runKappa1_process","Process",width="100%",class="btn btn-primary"))
                                                        )
                                                )
                                            ),
                                            
                                            #3.Run consistency evaluation between current subtypes derived from multi-omics clustering and PAM-predicted subtypes on tcga cohort
                                            conditionalPanel(condition="input.ntp_pam_runKappa.includes('2')",
                                                 
                                                             wellPanel(style="background-color:#47AEE9",
                                                                       fixedRow(
                                                                           column(12,tags$h4(strong("Run Kappa on tcga cohort (CMOIC vs PAM)"),align="center",style="color:#834005")),
                                                                           
                                                                           #Indicate the label of the first subtype
                                                                           column(12,textInput(inputId="subt1_lab_runKappa2",label="Indicate the label of the first subtype",value="CMOIC_TCGA")),
                                                                           
                                                                           #Indicate the label of the second subtype
                                                                           column(12,textInput(inputId="subt2_lab_runKappa2",label="Indicate the label of the second subtype",value="PAM_TCGA")),
                                                                           
                                                                           #The width of output figure
                                                                           column(12,numericInput("width_runKappa2",label="The width of output figure",value=5,min=0.5,step=0.5)),
                                                                           
                                                                           #The height of output figure
                                                                           column(12,numericInput("height_runKappa2",label="The height of output figure",value=5,min=0.5,step=0.5)),
                                                                           
                                                                           #The name of the consistency heatmap
                                                                           column(12,textAreaInput(inputId = "runKappa2_FigureName",
                                                                                                   label = "The name of the consistency heatmap",
                                                                                                   value = "constheatmap(CMOIC_VS_PAM_TCGA)",
                                                                                                   width = '100%', height = '35px')
                                                                           ),
                                                                           
                                                                           column(12,actionButton("runKappa2_process","Process",width="100%",class="btn btn-primary"))            
                                                                       )
                                                             )
                                            ),
                                            
                                            #4.Run consistency evaluation between NTP-predicted subtypes and PAM-predicted subtypes on tcga cohort
                                            conditionalPanel(condition="input.ntp_pam_runKappa.includes('1') & input.ntp_pam_runKappa.includes('2')",
                                                             
                                                             wellPanel(style="background-color:#47AEE9",
                                                                       fixedRow(
                                                                           column(12,tags$h4(strong("Run Kappa on tcga cohort (NTP vs PAM)"),align="center",style="color:#834005")),
                                                                           
                                                                           #Indicate the label of the first subtype
                                                                           column(12,textInput(inputId="subt1_lab_runKappa3",label="Indicate the label of the first subtype",value="NTP_TCGA")),
                                                                           
                                                                           #Indicate the label of the second subtype
                                                                           column(12,textInput(inputId="subt2_lab_runKappa3",label="Indicate the label of the second subtype",value="PAM_TCGA")),
                                                                           
                                                                           #The width of output figure
                                                                           column(12,numericInput("width_runKappa3",label="The width of output figure",value=5,min=0.5,step=0.5)),
                                                                           
                                                                           #The height of output figure
                                                                           column(12,numericInput("height_runKappa3",label="The height of output figure",value=5,min=0.5,step=0.5)),
                                                                           
                                                                           #The name of the consistency heatmap
                                                                           column(12,textAreaInput(inputId = "runKappa3_FigureName",
                                                                                                   label = "The name of the consistency heatmap",
                                                                                                   value = "constheatmap(NTP_VS_PAM_TCGA)",
                                                                                                   width = '100%', height = '35px')
                                                                           ),
                                                                           
                                                                           column(12,actionButton("runKappa3_process","Process",width="100%",class="btn btn-primary"))            
                                                                       )
                                                             )
                                            ),
                                            
                                            #5.Run consistency evaluation between NTP-predicted subtypes and PAM-predicted subtypes on validation cohort
                                            conditionalPanel(condition="input.ntp_pam_runKappa.includes('3') & input.ntp_pam_runKappa.includes('4')",
                                                             
                                                             wellPanel(style="background-color:#47AEE9",
                                                                       fixedRow(
                                                                           column(12,tags$h4(strong("Run Kappa on validation cohort (NTP vs PAM)"),align="center",style="color:#834005")),
                                                                           
                                                                           #Indicate the label of the first subtype
                                                                           column(12,textInput(inputId="subt1_lab_runKappa4",label="Indicate the label of the first subtype",value="NTP_Validation")),
                                                                           
                                                                           #Indicate the label of the second subtype
                                                                           column(12,textInput(inputId="subt2_lab_runKappa4",label="Indicate the label of the second subtype",value="PAM_Validation")),
                                                                           
                                                                           #The width of output figure
                                                                           column(12,numericInput("width_runKappa4",label="The width of output figure",value=5,min=0.5,step=0.5)),
                                                                           
                                                                           #The height of output figure
                                                                           column(12,numericInput("height_runKappa4",label="The height of output figure",value=5,min=0.5,step=0.5)),
                                                                           
                                                                           #The name of the consistency heatmap
                                                                           column(12,textAreaInput(inputId = "runKappa4_FigureName",
                                                                                                   label = "The name of the consistency heatmap",
                                                                                                   value = "constheatmap(NTP_VS_PAM_Validation)",
                                                                                                   width = '100%', height = '35px')
                                                                           ),
                                                                           
                                                                           column(12,actionButton("runKappa4_process","Process",width="100%",class="btn btn-primary"))            
                                                                       )
                                                             )
                                            ),
                                            
                                            #All steps in MOVICS have been finished
                                            wellPanel(style="background-color:#47AEE9",
                                                      fixedRow(
                                                          column(12,tags$h4(strong("Finish"),align="center",style="color:#834005")),
                                                          column(12,offset=0,tags$p(h5(strong("All steps in MOVICS have been finished!",align="left")))),
                                                          column(12,actionButton("MOVICS_finish_process","Finish",width="100%",class="btn btn-primary"))            
                                                      )
                                            )
                           )
                    ),
                    br(),
                    column(7,
                           tabsetPanel(
                               tabPanel(tags$h5(tags$strong("DEA"),align="center",style="color:#FFA500"),br(),
                                        
                                        uiOutput("table.runDEATitle"), #标题
                                        column(12,br()),
                                        dataTableOutput("table.runDEA"), #表格
                                        # uiOutput("table.runDEANote"), #脚注
                                        br(),br(),
                                        column(12,uiOutput("runDEA_finish"))        
                               ),
                               tabPanel(tags$h5(tags$strong("Biomarkers Identification"),align="center",style="color:#FFA500"),br(),
                                        
                                        # column(12,h3(plotOutput("runMarker_figure"))),
                                        uiOutput("runMarker_figure_unit"), #heatmap  
                                        br(),br(),
                                        uiOutput("table.runMarkerTitle"), #标题
                                        column(12,br()),
                                        dataTableOutput("table.runMarker"), #表格
                                        # uiOutput("table.runMarkerNote"), #脚注
                                        br(),br(),
                                        column(12,uiOutput("runMarker_finish"))
                               ),
                               tabPanel(tags$h5(tags$strong("GSEA"),align="center",style="color:#FFA500"),br(),
                                        
                                        column(12,uiOutput("gsea_geneset_background_preparation_finish2")),
                                        column(12,uiOutput("gsea_geneset_background_preparation_finish3")),
                                        br(),br(),
                                        # column(12,h3(plotOutput("runGSEA_figure"))),
                                        uiOutput("runGSEA_figure_unit"), #heatmap  
                                        br(),br(),
                                        uiOutput("table.gseaResultTitle"), #标题1
                                        column(12,br()),
                                        dataTableOutput("table.gseaResult"), #表格1
                                        # uiOutput("table.gseaResultNote"), #脚注1
                                        br(),br(),
                                        uiOutput("table.enrichmentScoresTitle"), #标题2
                                        column(12,br()),
                                        dataTableOutput("table.enrichmentScores"), #表格2
                                        # uiOutput("table.enrichmentScoresNote"), #脚注2
                                        br(),br(),
                                        column(12,uiOutput("runGSEA_finish"))
                               ),
                               tabPanel(tags$h5(tags$strong("GSVA"),align="center",style="color:#FFA500"),br(),
                                        
                                        column(12,uiOutput("gsva_geneset_interest_preparation_finish2")),
                                        column(12,uiOutput("gsva_geneset_interest_preparation_finish3")),
                                        br(),br(),
                                        # column(12,h3(plotOutput("runGSVA_figure"))),
                                        uiOutput("runGSVA_figure_unit"), #heatmap  
                                        br(),br(),
                                        uiOutput("table.gsvaRawTitle"), #标题1
                                        column(12,br()),
                                        dataTableOutput("table.gsvaRaw"), #表格1
                                        # uiOutput("table.gsvaRawNote"), #脚注1
                                        br(),br(),
                                        uiOutput("table.gsvaZscoredTitle"), #标题2
                                        column(12,br()),
                                        dataTableOutput("table.gsvaZscored"), #表格2
                                        # uiOutput("table.gsvaZscoredNote"), #脚注2
                                        br(),br(),
                                        column(12,uiOutput("runGSVA_finish"))
                               ),
                               tabPanel(tags$h5(tags$strong("NTP"),align="center",style="color:#FFA500"),br(),
                                        
                                        # column(12,h3(plotOutput("runNTP_figure"))),
                                        uiOutput("runNTP_figure_unit"), #heatmap  
                                        br(),br(),
                                        uiOutput("table.ntpTitle"), #标题
                                        column(12,br()),
                                        dataTableOutput("table.ntp"), #表格
                                        # uiOutput("table.ntpNote"), #脚注
                                        br(),br(),
                                        column(12,uiOutput("runNTP_finish"))
                               ),
                               tabPanel(tags$h5(tags$strong("PAM"),align="center",style="color:#FFA500"),br(),
                                        
                                        uiOutput("table.pamIGPTitle"), #标题1
                                        column(12,br()),
                                        dataTableOutput("table.pamIGP"), #表格1
                                        # uiOutput("table.pamIGPNote"), #脚注1
                                        br(),br(),
                                        uiOutput("table.pamPredictionTitle"), #标题2
                                        column(12,br()),
                                        dataTableOutput("table.pamPrediction"), #表格2
                                        # uiOutput("table.pamPredictionNote"), #脚注2
                                        br(),br(),
                                        column(12,uiOutput("runPAM_finish"))
                               ),
                               tabPanel(tags$h5(tags$strong("Kappa Statistics"),align="center",style="color:#FFA500"),br(),
                                        
                                        # column(12,h3(plotOutput("runKappa_figure1"))),
                                        uiOutput("runKappa_figure_unit1"), #consistency heatmap(CMOIC_VS_NTP_TCGA)  
                                        br(),br(),
                                        column(12,uiOutput("runKappa_finish1")),
                                        br(),br(),
                                        # column(12,h3(plotOutput("runKappa_figure2"))),
                                        uiOutput("runKappa_figure_unit2"), #consistency heatmap(CMOIC_VS_PAM_TCGA)  
                                        br(),br(),
                                        column(12,uiOutput("runKappa_finish2")),
                                        br(),br(),
                                        # column(12,h3(plotOutput("runKappa_figure3"))),
                                        uiOutput("runKappa_figure_unit3"), #consistency heatmap(NTP_VS_PAM_TCGA)  
                                        br(),br(),
                                        column(12,uiOutput("runKappa_finish3")),
                                        br(),br(),
                                        # column(12,h3(plotOutput("runKappa_figure4"))),
                                        uiOutput("runKappa_figure_unit4"), #consistency heatmap(NTP_VS_PAM_Validation)  
                                        br(),br(),
                                        column(12,uiOutput("runKappa_finish4")),
                                        br(),br(),
                                        column(12,uiOutput("MOVICS_finish"))
                               )
                           )
                    )
                ),
                
                #Users Guide Section
                tabPanel(
                    title=tags$p(tags$h4(tags$strong("Users Guide"),style="color:purple")),
                    value="Users Guide",
                    
                    column(12,
                           uiOutput("usersGuide")
                    )
                    
                ),
                
                #About
                tabPanel(
                    title=tags$p(tags$h4(tags$strong("About"),style="color:purple")),
                    value="About",
                    
                    column(12,
                           uiOutput("about")
                    )
                    
                )
    )
    
    
)

