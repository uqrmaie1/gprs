shinyUI(fluidPage(
  
	titlePanel(fluidRow(column(11, span(paste0("GPRS, an online tool for Genomic Profile Risk Score analyses version 0.9.0"), title='')),
             column(1, actionButton("showTutorial", label = "?"))), windowTitle='GPRS'),
	
	# loading message & conditionalPanel for the loading message
		tags$head(tags$style(type="text/css", "
	             #loadmessage {
	               position: fixed;
	               top: 0px;
	               left: 0px;
	               width: 100%;
	               padding: 5px 0px 5px 0px;
	               text-align: center;
	               font-weight: bold;
	               font-size: 100%;
	               color: black;
	               background-color: grey;
	               z-index: 105;
	             }
	          "),
		          tags$script(HTML('
      Shiny.addCustomMessageHandler("jsCode",
        function(message) {
          console.log(message)
          eval(message.code);
        }
      );
    '))
              ),
		  conditionalPanel(condition="$('html').hasClass('shiny-busy')",
	                   tags$div("Loading...",id="loadmessage")),

	# tutorialUI call & cleanup button
		fluidRow(
			#actionButton("showTutorial", label = "Show/hide Tutorial"),
			#actionButton("cleanupButton", label = "Cleanup"),
			uiOutput("tutorialUI")
		),
	hr(),

	fluidRow(
		column(4,
			wellPanel(# 1. Upload Input
				h3("1. Upload Input"),
					fileInput('file1', p('Score Files -', a(paste("Example"), href="example_profile_scores.zip", target="_blank")), multiple = "T"),
					fileInput('file2', 'Optional: Covariate Files', multiple="T"),
					fileInput('file3', p('Label File -', a("Example", href="examples/disease1_label.csv", target="_blank"))),
        h4('   or'),
				actionButton("exampleButton", label = "Load examples")
				),

			wellPanel(# 2. Choose Covariates
				h3("2. Choose Covariates"),
					p("Unselected covariates are excluded from the analysis"),
					uiOutput("covariatesChecklistUI")
			),

			wellPanel(# 3. Run GPRS Regression Script
				h3("3. Run GPRS Regression Script"),
					#p("The button below runs the GPRS regression script, which takes information from the score and label files, analyzes them, and stores the output in two files called 'GPRS.csv' and 'GPRS_deciles.csv'. The tables and graphs below depend on these files."),
					actionButton("action", label = "Run Main Script"),
					actionButton("GxEaction", label = "Run GxE (experimental)"),
				  #checkboxInput('gxeMethod', 'alternative gxe method', value = FALSE),
					br()
					#br(),
					#p("May take 30-60 seconds. Any object that depends on the output of this script will turn translucent while the script is running.")
			)#,

# 			wellPanel(# 4. Run AEGP Script (Optional)
# 				h3("4. Run AEGP Script (Optional)"),
# 					#p("The buttons below run one of two versions of the AEGP script, which takes mainly summary statistics from the label file and 'GPRS.csv', and stores the output in either 'Vg2FromP.csv' or 'CorrFromP.csv' depending on which version is run. You can view these files in a tabular format in the tab 'Output - AEGP'."),
# 					#actionButton("actionVg2", label = "Predict Vg2 from P"),
# 					#actionButton("actionCorr", label = "Predict Corr from P"),
# 					#br(),
# 					br(),
# 					p('Dudbridge, F. (2013). "Power and predictive accuracy of polygenic risk scores." ',
# 						tags$u('Plos Genetics'),
# 						strong('9'),
# 						'(3): e1003348. ', 
# 						a('PMID: 23555274', href='http://www.ncbi.nlm.nih.gov/pubmed/23555274', target='_blank')
# 					)
# 			)
		),
    column(8,
           fluidRow(uiOutput("sf")),
           fluidRow(
             uiOutput("paneltest")
       )
	)
  
	  
		
)))
