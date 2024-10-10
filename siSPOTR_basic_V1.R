library(shiny)
library(DT)
library(shinythemes)
library(Biostrings)
library(data.table)
library(stringr)
library(shinyjs)

# UI
ui <- fluidPage(
  theme = shinytheme("united"),
  titlePanel("siSPOTR App"),
  navbarPage(
    "I want to...",
    id = "nav",
    tabPanel(
      "check my gene(s) for best guide RNAs",
      sidebarLayout(
        sidebarPanel(
          selectInput("species", "Select Organism",
                      choices = c("human", "mouse")),
          numericInput("gsize", "guide RNA size (nt)",
                      min = 19,
                      max = 22,
                      value = 22),
          fileInput("fastaFile1", "Upload FASTA File"),
          actionButton("runButton1", "Go!")
        ),
        mainPanel(
          #actionButton("runButton1", "Run"),
          dataTableOutput("resultTable1")
        )
      )
    ),
    tabPanel(
      "assess the off-target potential of my guide RNA(s)",
      sidebarLayout(
        sidebarPanel(
          selectInput("species", "Select Organism",
                      choices = c("human", "mouse")),
          checkboxInput("revcomp", "These are target sites (anti-sense)", value = TRUE),
          fileInput("fastaFile2", "Upload FASTA File"),
          actionButton("runButton2", "Go!")
        ),
        mainPanel(
          #actionButton("runButton2", "Run"),
          dataTableOutput("resultTable2")
        )
      )
    )
  )
)

# Server
server <- function(input, output, session) {
#Function 1 takes gene sequence and spits out guide RNAs

  processFASTA1 <- function(fastaFile1) {
    target.mRNA <- readDNAStringSet(fastaFile1$datapath)
    names(target.mRNA) <- make.names(names(target.mRNA), unique = TRUE)
    POTS.file <- "POTS.txt"
    POTS.dt <- fread(POTS.file)
    guide.size <- input$gsize #22 # Probably limit to 19-22. Or warn if different.
    target.species <- input$species # For seed off-target POTS score. Only Human and Mouse right now.
    # allow.seed.gu.wobble <- TRUE 
    gc.min <- 0.2
    gc.max <- 0.7
    n.gc.pass <- 2 # Hard-code this for interface. 
    is.pol3 <- TRUE # Will need to do the pol3 check once guides are put into scaffolds (to-do)
    
    # Functions ----
    # Splits target mRNA into guide-sized windows
    get.all.target.sites.f <- function(seqs, size=22){
      seq.list <- lapply(seq_along(seqs),
                         function(i){
                           t.seq <- seqs[[i]]
                           seq.len <- nchar(t.seq)
                           nm <- names(seqs)[[i]]
                           v <- trim(Views(t.seq, start = seq(seq.len), width = size))
                           v <- v[width(v)==size]
                           names(v) <- paste0(nm, "_p_", start(v), "_", end(v))
                           return(RNAStringSet(v))
                         })
      do.call(c, seq.list)
    }
    # Takes the target sites and generates appropriate guide and passenger,injecting G:U wobbles for biasing, where appropriate.
    # Also adds, POTS, gc% and flag for passing filters.
    characterize.and.bias.guides.f <- function(target.seq, POTS=POTS.dt, n.pass.gc=n.gc.pass, species=target.species){
      # Help bias guides first
      guides.rna <- reverseComplement(target.seq)
      guide.names <- names(target.seq)
      guides.rna.right <- as.character(subseq(guides.rna, start=2, end=width(guides.rna)))
      guides.rna.1 <- subseq(guides.rna, start=1, width=1) %>% ifelse(.=="C", "u", .)
      guide.rna.out <- paste0(guides.rna.1, guides.rna.right)
      # Now passenger across from the 5' end of guide
      passenger.rna <- reverseComplement(subseq(guides.rna, 1, width=guide.size-2))
      pass.width <- width(passenger.rna)
      passenger.rna.guideP1P2 <- as.character(subseq(passenger.rna, 
                                                     start=pass.width-1, 
                                                     end=pass.width)) %>% gsub("C", "u", .)
      passenger.rna.out <- paste0(subseq(passenger.rna, 1, pass.width-2),
                                  passenger.rna.guideP1P2,
                                  "uu")
      dt <- data.table(seq.id=guide.names,
                       target.seq=as.character(target.seq),
                       guide.seq=guide.rna.out,
                       pass.seq=passenger.rna.out)
      # Calc GC%, and other features from user-specified 
      pass.5p.search <- ifelse(n.pass.gc<2, 2, n.pass.gc) # Will search the first 2 bases of the passenger strand for GC, unless user specifies more than 2
      dt[, `:=`(Seed=subseq(guide.seq, start=2, end=8),
                GC=stringr::str_count(guide.seq, "[GC]")/nchar(guide.seq),
                passenger.bias=stringr::str_count(subseq(pass.seq, 1, width=pass.5p.search), "[GC]"))]
      dt.merged <- merge.data.table(dt, POTS, by="Seed")
      dt.merged[, c("target.name", "target.pos.tmp"):=tstrsplit(seq.id, split="_p_")]
      dt.merged[, c("start.target", "end.target"):=tstrsplit(target.pos.tmp, split="_", type.convert = TRUE)]
      
      output.order <- c("seq.id", "guide.seq", "pass.seq", "target.name", "start.target", "end.target", "GC", "passenger.bias", "POTS.HSA", "POTS.MMU", "SPS")
      # Clean and output
      if(tolower(species)=="mouse" | tolower(species)=="mmu"){
        setorder(dt.merged, POTS.MMU)
      } else{
        setorder(dt.merged, POTS.HSA)
      }
      return(dt.merged[, .SD, .SDcols=output.order])
    }
    
    data1 <- get.all.target.sites.f(target.mRNA) %>% 
      characterize.and.bias.guides.f(target.seq=.)
    
    # Pass target.sites all to shiny. Will set filters here, but may want to be reactive ----
    
    #target.sites.all[between(GC, gc.min, gc.max) & passenger.bias>=n.gc.pass]
    return(data1)
    
  }
# Function 2 takes guide RNAs and assesses off-target potential
  processFASTA2 <- function(fastaFile2) {
    POTS.file <- "POTS.txt"
    POTS.dt <- fread(POTS.file)
    # User-supplied:
    is.guide <- TRUE # Have user select if the inputs are Guides (Antisense) or target sites
    # Guide input:
    # THis is a user-interface thing. I'd recommend allowing file upload, or
    # direct pasting of 22-mers. Pasting could be one per line, or fasta
    # For now. Hard-coded and example set of guides
    guide.fl <- "data/test_guides.fa" 
    guide.in <- readBStringSet(fastaFile2$datapath)
    if(is.guide==FALSE){
      guides <- RNAStringSet(reverseComplement(guide.in))
    } else{
      guides <- RNAStringSet(guide.in)
    }
    if(length(names(guides))==0){
      names(guides) <- paste0("Guide_", seq_along(guides))
    }
    seeds <- subseq(guides, start = 2, end = 8)
    guides.dt <- data.table(seq.id=names(guides),
                            seq.in=as.character(guide.in),
                            guide.seq=as.character(guides),
                            Seed=as.character(seeds))
    data2 <- merge.data.table(guides.dt, POTS.dt, on="Seed", all.x=TRUE)
    #guide.dt.annot
    return(data2)
  }
  
  # Store the processed data as reactive values
  data1 <- reactiveValues()
  data2 <- reactiveValues()
  
  # Event handler for Analysis 1 Run button
  observeEvent(input$runButton1, {
    if (!is.null(input$fastaFile1)) {
      data1$table <- processFASTA1(input$fastaFile1)
    }
  })
  
  # Event handler for Analysis 2 Run button
  observeEvent(input$runButton2, {
    if (!is.null(input$fastaFile2)) {
      data2$table <- processFASTA2(input$fastaFile2)
    }
  })
  
  # Render the Analysis 1 data table in the UI
  output$resultTable1 <- renderDataTable({
    if (!is.null(data1$table)) {
      data1$table
    }
  })
  
  # Render the Analysis 2 data table in the UI
  output$resultTable2 <- renderDataTable({
    if (!is.null(data2$table)) {
      data2$table
    }
  })
}

# Run the Shiny app
shinyApp(ui, server)
