###################################################################################
#
# FALCONER SHINY APP
#
# Description : The aim of this App is to show how the combination of gene action 
# and allele frequencies at causal loci translate to genetic variance and genetic
# variance components for a complex trait. Although the theory underlying the App
# is more than a century old, it is highly relevant in the current era of 
# genome-wide association studies (GWAS). It highlights the specific definition of
# the effect size estimates by GWAS and the variation it generates in the population,
# i.e. how locus-specific effects leads to individual differences. The App can be used
# to demonstrate the relationship between a SNP effect size estimated from GWAS and the 
# variation the SNP generates in the population, i.e. how locus-specific effects lead 
# to individual differences. In addition, it can also be used to demonstrate how within
# and between locus interactions (dominance and epistasis, respectively) usually do not
# lead to a large amount of non-additive variance relative to additive variance, and
# therefore that these interactions usually do not explain individual differences in 
# a population.
#  
# The three models described mainly illustrate the Chapters 7 and 8 of Falconer 
# and Mackay (1996) and the Chapter 5 of Lynch and Walsh (1998).
# 
# Authors: Valentin Hivert, Naomi Wray and Peter Visscher
# Date: 03 Feb 2021

# Version 1.0
###################################################################################

##Libraries
require(shiny)
require(plot3D)
require(htmlTable)
require(pBrackets)
require(shinyalert)
require(rhandsontable)
require(htmltools)

### Function to derive components of genetic variance #############
Compute_GeneticVariances <- function(GV, p, q){
  
  #GV matrix of Genotypic values in the format row for locus B (B2B2, B1B2, B1B1) and column for locus A (A2A2, A1A2, A1A1)
  
  freqA=c( (1-p)**2, 2*p*(1-p), p**2 )
  freqB=c( (1-q)**2, 2*q*(1-q), q**2 )
  Geno_freq=freqB%*%t(freqA)
  
  #Population mean
  M=sum(GV*Geno_freq)
  
  ############################################################
  #Additive effects
  CondQ.freq=c(1-q,q)%*%t(freqA)#Conditionnal frequencies to Q
  CondP.freq=freqB%*%t(c(1-p,p))#Conditionnal frequencies to P
  Cond.Mean=Alpha=matrix(0,nrow = 2,ncol = 2)
  #To optimize with matrix calc.
  for(i in 1:2){
    Cond.Mean[1,i]=sum(GV[,c(i,i+1)]*CondP.freq)#Conditionnal Mean for Locus 1 (Locus A)
    Cond.Mean[2,i]=sum(GV[c(i,i+1),]*CondQ.freq)#Conditionnal Mean for Locus 2 (Locus B)
  }
  Alpha=Cond.Mean-M #Average effects
  Expected.Geno=cbind(2*Alpha[,1],rowSums(Alpha),2*Alpha[,2])
  
  ############################################################
  #Dominance effects
  
  D_effects=matrix(NA,nrow = 2,ncol = 3)#3*3 matrix, row 1 = Locus A, row 2 = locus B. 1 col per genotype
  D_effects[1,]=freqB%*%GV-M-Expected.Geno[1,]#Locus A dominance effects
  D_effects[2,]=freqA%*%t(GV)-M-Expected.Geno[2,]#Locus B dominance effects
  
  ############################################################
  #Additive X Additive effects
  
  #Conditional mean of A1.B1. etc...
  i=j=c(1,2)
  tmp=c(1-q,q)%*%t(c(1-p,p))
  aa=do.call(what = rbind,lapply(i,FUN = function(x,y,idx,tmp){res=rep(0,2);for(j in idx){res[j]=sum(y[c(x,x+1),c(j,j+1)]*tmp)};return(res)},y=GV,idx=j,tmp=tmp))#For each allele at locus A
  AA_effects=aa-M-rbind(Alpha[2,1]+Alpha[1,],Alpha[2,2]+Alpha[1,])
  aa=sum(AA_effects[,1]-AA_effects[,2])
  rm(tmp);rm(i);rm(j)
  
  #Additive X Dominance effects
  AD_effects=matrix(NA,nrow = 4,ncol = 3)
  for(i in 1:2){
    AD_effects[i,]=colSums(c(1-p,p)*t(GV[,c(i,i+1)]))-M-D_effects[2,]-Expected.Geno[2,]-Alpha[1,i]-c(2*AA_effects[1,i],sum(AA_effects[,i]),2*AA_effects[2,i])#Conditionnal on locus A ad_A1.B1B1 etc...
    AD_effects[i+2,]=colSums(c(1-q,q)*GV[c(i,i+1),])-M-D_effects[1,]-Expected.Geno[1,]-Alpha[2,i]-c(2*AA_effects[i,1],sum(AA_effects[i,]),2*AA_effects[i,2])
  }
  
  ############################################################
  #Dominance X Dominance effects
  BV=t(outer(Expected.Geno[1,],Expected.Geno[2,],FUN = "+"))#Breeding values
  D.Geno=t(outer(D_effects[1,],D_effects[2,],FUN = "+"))
  AA.Geno=t(apply(AA_effects,MARGIN = 1,FUN = function(x){return(c(2*x[1],sum(x),2*x[2]))}))
  AA.Geno=rbind(2*AA.Geno[1,],AA.Geno[1,]+AA.Geno[2,],2*AA.Geno[2,])
  AD.Geno=t(rbind(2*AD_effects[1,],colSums(AD_effects[1:2,]),2*AD_effects[2,])) + rbind(2*AD_effects[3,],colSums(AD_effects[3:4,]),2*AD_effects[4,])#AD effect for each genotypes
  
  DD_effects=GV-M-BV-D.Geno-AA.Geno-AD.Geno
  
  ############################################################
  #VARIANCE COMPONENTS
  
  #Total genotypic variance
  VG=sum(Geno_freq*(GV**2))-M**2
  
  #Additive variance
  VA_locus=2*rowSums(matrix(c(1-p,p,1-q,q),ncol = 2,byrow = T)*Alpha**2)
  names(VA_locus)=c("Locus_A","Locus_B")
  VA=sum(VA_locus)
  # Check from the breeding values : weighted.mean(expected.Geno**2,w = Geno_freq)
  
  #Dominance variance
  VD_locus=rowSums(rbind(freqA,freqB)*D_effects**2)
  names(VD_locus)=c("Locus_A","Locus_B")
  VD=sum(VD_locus)
  
  #Additive-by-Additive variance
  VAA=4*weighted.mean(AA_effects**2,w=matrix(c((1-p)*(1-q),(1-q)*p,q*(1-p),q*p),2,2,byrow = T))
  VAA/sqrt(2*p*(1-p)*2*q*(1-q))
  
  
  #Additive-by-Dominance variance
  VAD=4*weighted.mean(AD_effects**2,w = matrix(c((1-p)*freqB,p*freqB,(1-q)*freqA,q*freqA),ncol = 3,byrow = T))
  
  #Dominance X Dominance variance
  VDD=weighted.mean(DD_effects**2,w = Geno_freq)
  
  
  #Summary
  # VG
  # VA+VD+VAA+VAD+VDD
  
  res.var=list(VA=list(VA_locus,VA), VD=list(VD_locus,VD), VAA=VAA, VAD=VAD, VDD=VDD, VG=VG)
  return(list(mu=M,alpha=Alpha[,2]-Alpha[,1],var=res.var))
}
###################################################

### For user defined GV model, initial matrix with example values from Lynch and Walsh 1998  #############
GV=matrix(c(18,40.9,61.1,54.6,47.6,66.5,47.8,83.6,101.7),ncol = 3,byrow = T,dimnames = list(c("B\u2082B\u2082","B\u2081B\u2082","B\u2081B\u2081"),c("A\u2082A\u2082","A\u2081A\u2082","A\u2081A\u2081")))

###################################################

ui <- fluidPage(
  useShinyalert(),
  
  titlePanel(title = "The Falconer ShinyApp"),
  HTML("<b>Choose a model to display:</b>"),
  navbarPage(tags$b("Model"), id="model",
             tabPanel("Single-locus Additive and Dominance", value="AD"),
             tabPanel("Two-locus Additive and Additive-by-Additive", value="AA"),
             tabPanel("General two-locus model", value="Perso")
  ),
  conditionalPanel(
    condition = "input.model == 'AD'",
    div(HTML("This model uses the notation of Falconer and Mackay (1996) Chapters 7 & 8, and uses their equations to generate the population specific variance components given the allele frequency p and the genotypic values <i>a</i> and <i>d</i>.<br><br>")
        ,style="text-align: left;", style = "background-color:#f2f2f2;", style = "font-size:110%;")
  ),
  conditionalPanel(
    condition = "input.model == 'AA'",
    div(HTML("This model uses the notation of Falconer and Mackay (1996) Chapters 7 & 8 as well as Mäki-Tanila and Hill (2014). It uses their equations to generate the population specific variance components given the allele frequencies <i>p</i> (locus A) and <i>q</i> (locus B), the genotypic values of each locus <i>a</i><sub>A</sub> and <i>a</i><sub>B</sub>, as well as the aditive-by-additive effect <i>a</i><sub>AB</sub>.<br><br>")
        ,style="text-align: left;", style = "background-color:#f2f2f2;", style = "font-size:110%;")
    
  ),
  conditionalPanel(
    condition = "input.model == 'Perso'",
    div(HTML("This model uses the general least square model described in Lynch and Walsh (1998) Chapter 5 to derive the different variance components given the genotypic values input by the user as well as the chosen allele frequencies <i>p</i> (locus A) and <i>q</i> (locus B).<br><br>")
        ,style="text-align: left;", style = "background-color:#f2f2f2;", style = "font-size:110%;")
    
  ),
  div(HTML("To obtain help and details about the application, please click on the README button: "),actionButton("preview", "README",icon = icon("question-circle")), style="text-align: left;", style = "background-color:#f2f2f2;", style = "font-size:110%;"),
  HTML("<br>"),
  sidebarLayout(
    sidebarPanel(
      sliderInput("sliderP", HTML("A<sub>1</sub> allele frequency (<i>p</i>)"),
                  min = 0.001, max = 0.999, value = 0.3,step = 0.01),
      
      conditionalPanel(
        condition = "input.model == 'AA' || input.model == 'Perso'",
        sliderInput("sliderQ", HTML("B<sub>1</sub> allele frequency (<i>q</i>)"),
                    min = 0.001, max = 0.999, value = 0.3,step = 0.01)
      ),
      
      conditionalPanel(
        condition = "input.model == 'AD' || input.model == 'AA'",
        uiOutput("SliderA1")
      ),
     
      conditionalPanel(
        condition = "input.model == 'AD'",
        sliderInput("sliderD", HTML("Genotypic value <i>d</i>"),
                    min = -10.0, max = 10.0, value = 4.0,step = 1)
      ),  
      
      conditionalPanel(
        condition = "input.model == 'AA'",
        sliderInput("sliderA2", HTML("Locus B genotypic value <i>a</i><sub>B</sub>"),
                    min = -10.0, max = 10.0, value = 4.0),
        sliderInput("sliderAA", HTML("Additive-by-Additive interaction effect (<i>a</i><sub>AB</sub>)"),
                    min = -10.0, max = 10.0, value = 2.0)
      ),  
      
      conditionalPanel(
        condition = "input.model == 'Perso'",
        HTML("<b>Input your genotypic values in the table:</b><br><br>"),
        rHandsontableOutput("hot"),
        HTML("<br>")
      )
    ),
    mainPanel(
      conditionalPanel(
        condition = "input.model == 'AD'",
        div(img(src='AD_model.png', width="50%"), br(), tags$b("Arbitrarily assigned genotypic values (i.e. trait means per genotype class)"), br(), br(), style="text-align: center;", style = "background-color:#f2f2f2;")
      ),
      conditionalPanel(
        condition = "input.model == 'AA'",
        div(img(src='AAA_model.png', width="40%"), br(), tags$b("Arbitrarily assigned genotypic values (i.e. trait means per genotype class)"), br(), br(), style="text-align: center;", style = "background-color:#f2f2f2;")
      ),
      
      splitLayout(
        align="left",
        cellWidths = c("40%","60%"),
        verticalLayout(htmlOutput("table"),
                       HTML("<b>Variance components</b>"),
                       htmlOutput("varianceA"),
                       conditionalPanel(
                         condition = "input.model == 'AD' || input.model == 'Perso'",
                         htmlOutput("varianceD")
                       ),
                       conditionalPanel(
                         condition = "input.model == 'AA' || input.model == 'Perso'",
                         htmlOutput("varianceAA")
                       ),
                       conditionalPanel(
                         condition = "input.model == 'Perso'",
                         htmlOutput("varianceAD")
                       ),
                       conditionalPanel(
                         condition = "input.model == 'Perso'",
                         htmlOutput("varianceDD")
                       ),
                       htmlOutput("varianceG"),
                       htmlOutput("varianceRatio")
                       ),
        plotOutput("plot"),
        cellArgs = list(style = "vertical-align: top;")
      ),
      
      fluidRow(

        column(4,HTML("")
        ),

        column(8,
               conditionalPanel(
                 condition = "input.model == 'AA' || input.model == 'Perso'",
                 HTML("Graphical representation of genotypic values (closed circles) at two biallelic loci A and B. The horizontal scale show  the number of A<sub>1</sub> alleles in the genotype. The different genotypes at locus B are depicted in different colors. For each of the locus B genotypes, the linear regression line between the number of A<sub>1</sub> allele and the genotypic value is drawn.")
               ),
               conditionalPanel(
                 condition = "input.model == 'AD'",
                 HTML("Reproduction of the Figure 7.2 of Falconer and Mackay (1996). Graphical representation of genotypic (closed blue circles) and breeding (open blue circles) values,
             of the genotypes for a locus with two alleles A<sub>1</sub> and A<sub>2</sub> at freqencies <i>p</i> and 1-<i>p</i>.
             Horizontal scale: number of A<sub>1</sub> allele in the genotype. Vertical scales of values are : on the left-arbitrarily assigned values (see figure at the top of the page);
             on right-deviation from the population mean (black cross). Each point size is weighted by its genotype frequency. A linear regression line between the number of A<sub>1</sub> allele and the genotypic values is fitted by weighted least square.")
               )
        )
      ),
      
      conditionalPanel(
        condition = "input.model == 'AA'",
        plotOutput(outputId="contourplot")
      ),
      conditionalPanel(
        condition = "input.model == 'AD'",
        plotOutput(outputId="plotDom")
      ),
      conditionalPanel(
        condition = "input.model == 'AA'",
        HTML("Additive (<i>V<sub>A</sub></i>), additive-by-additive (<i>V<sub>AA</sub></i>) and proportion of genotypic variance explained by additive variance (<i>V<sub>A</sub>/V<sub>G</sub></i>) as a function of the allele frequencies <i>p</i> and <i>q</i>. Current setting of <i>p</i> and <i>q</i> is depicted with a white cross.")
      ),
      conditionalPanel(
        condition = "input.model == 'AD'",
        HTML("Distributions of additive (<i>V<sub>A</sub></i>), dominance (<i>V<sub>D</sub></i>) and total (<i>V<sub>G</sub></i>) genetic variance on the left, and proportion of genotypic variance explained by additive variance (<i>V<sub>A</sub>/V<sub>G</sub></i>) as a function of the allele frequency <i>p</i> on the right. The current user input allele frequency <i>p</i> is depicted by a vertical red solid line in the left panel and by a red cross in the right panel.")
      )
    )
  ),
  hr(),
  uiOutput("References"),
  hr(),
  HTML("<br><b>Authors:</b> The Falconer ShinyApp was written by Valentin Hivert, based on the previous versions that had input from Luke Lloyd-Jones, Alex Holloway and Matt Robinson.<br><br><b>Contact:</b> v.hivert@imb.uq.edu.au")
  )

server <- function(input, output, session) {
  
  v <- reactiveValues(mu = 0, beta1 = 0, beta2 = 0, Va = 0, Va1 = 0, Va2 = 0, Vd = 0, Vd1 = 0, Vd2 = 0, Vaa = 0, Vad = 0, Vdd = 0, Vg = 0)#Pop mean, betas and Variance components.They are updated in the main plot output
  values <- reactiveValues()#Reactive values for the user defined genotypic values matrix
  init <- reactiveValues(idx=0)
  ### Helper dialogue #############
  observeEvent(input$preview, {
    showModal(modalDialog(
      title = "README",
      HTML("  <h1>Welcome in the Falconer ShinyApp</h1><br><br>
The aim of this App is to show how the combination of gene action and allele frequencies at causal loci translate to genetic variance and genetic variance components for a complex trait. Although the theory underlying the App is more than a century old, it is highly relevant in the current era of genome-wide association studies (GWAS). The App can be used to demonstrate the relationship between a SNP effect size estimated from GWAS and the variation the SNP generates in the population, <i>i.e.</i> how locus-specific effects lead to individual differences. In addition, it can also be used to demonstrate how within and between locus interactions (dominance and epistasis, respectively) usually do not lead to a large amount of non-additive variance relative to additive variance, and therefore that these interactions usually do not explain individual differences in a population.<br><br>

The three models described below mainly illustrate the Chapters 7 and 8 of Falconer and Mackay (1996) and the Chapter 5 of Lynch and Walsh (1998).<br><br>

<h3><b>Single-locus Model with additive and dominance effect:</b></h3><br>

In this single-locus model, we consider a biallelic locus with allele A<sub>1</sub> and A<sub>2</sub> in frequencies <i>p</i> and 1-<i>p</i>. Under panmixia (<i>i.e.</i> random mating) and Hardy-Weinberg equilibrium, the expected genotype frequencies are (1-<i>p</i>)<sup>2</sup>,2<i>p</i> (1-<i>p</i>) and <i>p</i><sup>2</sup>, for A<sub>2</sub>A<sub>2</sub>, A<sub>1</sub>A<sub>2</sub> and A<sub>1</sub>A<sub>1</sub> respectively. 
We arbitrarily assign genotypic values (<i>i.e.</i> trait means per genotype class) -<i>a</i>, <i>d</i> and <i>a</i> to the three genotypes, <i>d</i> representing the dominance effect (within locus interaction, no interaction when <i>d</i> = 0) and 2<i>a</i> the difference between the two homozygotes. Under this model, the population mean is:<br><br>

M = (2<i>p</i>-1)<i>a</i> + 2<i>p</i>(1-<i>p</i>)<i>d</i><br><br>

<b>Average effect of gene (allele) substitution (also called additive effect in the literature)</b><br>

The transmission of value from parents to offspring occurs through their genes (alleles) and not their genotypes. The average effect of gene substitution (&#120572;) is defined as the average effect on the trait when substituting alleles at this locus in the population. It can also be defined as the mean value of genotypes produced by different gametes:<br><br>

&#120572; = <i>a</i> + (1-2<i>p</i>)<i>d</i> <br><br>

Importantly, &#120572; is also the slope of the linear regression of the genotype means, weighted by their frequency, on the A<sub>1</sub> allele dosage (0, 1 or 2). <br><br>

When performing a standard GWAS, individual phenotypes <i>y</i> are regressed on the number <i>x</i> (<i>x</i> = 0, 1, 2) of reference alleles at a given locus, , <i>i.e.</i>, the allelic “dosage”, where the reference allele for this dosage count is arbitrarily the major or minor allele (but this arbitrary choice is reflected in the sign of the regression coefficient &beta;:<br><br>

y = &mu; + &beta;x + e<br><br>

Where the residuals e include both the non-additive genetic effects at the locus, the genetic effects (additive and non-additive) at other loci and an environmental and/or chance (non-genetic) effect. The quantity of interest is the slope &beta; of the model (the effect size of the locus), which is the average effect of allele substitution, hence &beta; = &#120572;.<br><br>

<b>Additive (breeding) values and dominance deviations</b><br>

The breeding values are the expected genotypic values under additivity (the predictions from the linear model). Expressed as deviations from the population mean M, the breeding values of the 3 genotypes A<sub>2</sub>A<sub>2</sub>, A<sub>1</sub>A<sub>2</sub> and A<sub>1</sub>A<sub>1</sub> are
      -2<i>p</i>&#120572;, (1-2<i>p</i>)&#120572; and (2-2<i>p</i>)&#120572;. The residuals of the linear regression are the deviations due to the within locus interaction (dominance)<br><br>

<b>Genetic variance</b><br>

The total genotypic variance (<i>V<sub>G</sub></i>) is partitioned into Additive (<i>V<sub>A</sub></i>) and Dominance (<i>V<sub>D</sub></i>) variance.<br><br>

<i>V<sub>G</sub></i> = <i>V<sub>A</sub></i> + <i>V<sub>D</sub></i><br><br>

<b>Additive variance</b><br>

The Additive variance (<i>V<sub>A</sub></i>) is the variance of additive (breeding) values. When values are expressed in terms of deviation from the population mean, the variance simply become the mean of the squared values. Hence, <i>V<sub>A</sub></i> is obtained by squaring the additive (breeding) values described above, multiplying by the corresponding frequencies and summing over the 3 genotypes, leading to:<br><br>

<i>V<sub>A</sub></i> = 2<i>p</i>(1-<i>p</i>)&#120572;<sup>2</sup> = 2<i>p</i>(1-<i>p</i>)[<i>a</i>+<i>d</i>(1-2<i>p</i>)]<sup>2</sup> = <i>H</i>&#120572;<sup>2</sup>, <br><br>

with <i>H</i> the heterozygosity at the locus. Note that <i>V<sub>A</sub></i> is the variance explained by a SNP in GWAS (2<i>p</i>(1-<i>p</i>)&#120572;<sup>2</sup> = 2<i>p</i>(1-<i>p</i>)&beta;<sup>2</sup>) and contain both a term due to additivity (<i>a</i>) and dominance (<i>d</i>) through the average effect &#120572;.<br><br>

<b>Dominance variance</b><br>

Similarly, the dominance variance (<i>V<sub>D</sub></i>) is the variance of dominance deviations:<br><br>

<i>V<sub>D</sub></i> = (2<i>p</i> (1-<i>p</i>)<i>d</i>)<sup>2</sup> = <i>H</i><sup>2</sup><i>d</i><sup>2</sup><br><br>

Therefore, the dominance variance disproportionally depends on the locus heterozygosity compared to the additive variance (<i>H</i><sup>2</sup> versus <i>H</i>).<br><br>

<h3><b>Two-locus Model with additive and additive-by-additive effect:</b></h3><br>

We extend the one-locus to a two-locus model with additive and additive-by-additive epistatic interaction only, assuming no within loci dominance effects (<i>d</i> = 0 at both loci). We introduce a second (unlinked) locus with alleles B<sub>1</sub> and B<sub>2</sub> in frequencies <i>q</i> and 1-<i>q</i> respectively. The genotypic values and allele frequencies of the 9 genotypes are:<br>",
           
           htmlTable(matrix(data = c("<b>-<i>a</i><sub>A</sub>-<i>a</i><sub>B</sub>+<i>a</i><sub>AB</sub></b><br>&nbsp;&nbsp;&nbsp;(1-<i>p</i>)<sup>2</sup>(1-<i>q</i>)<sup>2</sup>&nbsp;&nbsp;&nbsp;","<b>-<i>a</i><sub>B</sub></b><br>2<i>p</i>(1-<i>p</i>)(1-<i>q</i>)<sup>2</sup>","<b><i>a</i><sub>A</sub>-<i>a</i><sub>B</sub>-<i>a</i><sub>AB</sub></b><br><i>p</i><sup>2</sup>(1-<i>q</i>)<sup>2</sup>","<b>-<i>a</i><sub>A</sub></b><br>(1-<i>p</i>)<sup>2</sup>2<i>q</i>(1-<i>q</i>)","<b>0</b><br>4<i>p</i>(1-<i>p</i>)<i>q</i>(1-<i>q</i>)","<b><i>a</i><sub>A</sub></b><br>2<i>p</i><sup>2</sup><i>q</i>(1-<i>q</i>)","<b>-<i>a</i><sub>A</sub>+<i>a</i><sub>B</sub>-<i>a</i><sub>AB</sub></b><br>(1-<i>p</i>)<i>q</i><sup>2</sup>","<b><i>a</i><sub>B</sub></b><br>2<i>p</i>(1-<i>p</i>)<i>q</i><sup>2</sup>","<b><i>a</i><sub>A</sub>+<i>a</i><sub>B</sub>+<i>a</i><sub>AB</sub></b><br><i>p</i><sup>2</sup><i>q</i><sup>2</sup>"),nrow = 3,ncol = 3,byrow = T),
                     header =  c("A<sub>2</sub>A<sub>2</sub>","A<sub>1</sub>A<sub>2</sub>","A<sub>1</sub>A<sub>1</sub>"),
                     rnames = c("<b>B<sub>2</sub>B<sub>2</sub></b>","<b>B<sub>1</sub>B<sub>2</sub></b>","<b>B<sub>1</sub>B<sub>1</sub></b>"),
                     align.header = c('l',rep('c',3)),
                     align='l|ccc'),
           
           "where <i>a</i><sub>A</sub> (<i>a</i><sub>B</sub>) is the genotypic value for the upper homozygote A<sub>1</sub>A<sub>1</sub> (B<sub>1</sub>B<sub>1</sub>) and <i>a</i><sub>AB</sub> is the additive-by-additive interaction effect. This is a re-parametrization of the model described by Mäki-Tanila and Hill (2014).<br><br>

<b>Population mean</b><br>

In our model, the mean of the genotypic values is:<br><br>

M = <i>a</i><sub>A</sub>(2<i>p</i>-1) + <i>a</i><sub>B</sub>(2<i>q</i>-1) + <i>a</i><sub>AB</sub>(1-2(<i>p</i>+<i>q</i>)+4<i>pq</i>)<br><br>

Note that the expression of M depends on the arbitrarily assigned genotypic values.<br><br>

<b>Average effect of gene (allele) substitution</b><br>

In this model, the locus specific average effects are:<br><br>

&#120572;<sub>A</sub> = <i>a</i><sub>A</sub> + 2<i>qa</i><sub>AB</sub><br>

&#120572;<sub>B</sub> = <i>a</i><sub>B</sub> + 2<i>pa</i><sub>AB</sub><br><br>

<b>Genetic variance</b><br>

The total genotypic variance (<i>V<sub>G</sub></i>) of the model is partitioned into Additive (<i>V<sub>A</sub></i>) and Additive-by-Additive (<i>V<sub>AA</sub></i>) variance.<br><br>

<i>V<sub>G</sub></i> = <i>V<sub>A</sub></i> + <i>V<sub>AA</sub></i><br><br>

<b>Additive variance</b><br>

The additive variance of the model is:<br><br>

<i>V<sub>A</sub></i> = &sum;<sub><i>i</i></sub><i>H<sub>i</sub></i>&#120572;<sub><i>i</i></sub><sup>2</sup>, with <i>H<sub>i</sub></i> the heterozygosity at locus <i>i</i> (<i>i</i> = A, B) and &#120572;<sub><i>i</i></sub> the average effect of locus <i>i</i>. Hence:<br><br>

<i>V<sub>A</sub></i> = 2<i>p</i>(1-<i>p</i>)[<i>a</i><sub>A</sub>+2<i>qa</i><sub>AB</sub>] + 2<i>q</i>(1-<i>q</i>)[<i>a</i><sub>B</sub>+2<i>pa</i><sub>AB</sub>] <br><br>

Note that <i>V<sub>A</sub></i> contains a term due pairwise additive-by-additive interaction between locus A and B (<i>a</i><sub>AB</sub>).<br><br>

<b>Additive-by-Additive variance</b><br>

The additive-by-additive variance of the model is:<br><br>

<i>V<sub>AA</sub></i> = &sum;<sub><i>i</i></sub>&sum;<sub><i>j>i</i></sub><i>H<sub>i</sub>H<sub>j</sub>a<sub>ij</sub></i><sup>2</sup>, with <i>H<sub>i</sub></i> the heterozygosity at locus <i>i</i> (<i>i</i> = A, B) and <i>a<sub>ij</sub></i> the additive-by-additive interaction effect between locus <i>i</i> and <i>j</i>. Hence:<br><br>

<i>V<sub>AA</sub></i> = 4<i>p</i>(1-<i>p</i>)<i>q</i>(1-<i>q</i>)<i>a</i><sub>AB</sub><br><br>

Therefore, the additive-by-additive variance disproportionally depends on the locus heterozygosity as compared to the additive variance.<br><br>

<h3><b>General two locus model:</b></h3><br>

Lastly, we use a generalized two-locus model where the user can provide all the genotypic values in an interactive table and choose the allele frequencies at the two loci (<i>p</i> and <i>q</i>). The genotypic values as well as the linear regressions are plotted as a function of the A<sub>1</sub> allelic dosage for the different genotypes at locus B, as well as the linear regression of the genotypic values weighted by their frequency on the A<sub>1</sub> allele dosage. The total genotypic variance (<i>V<sub>G</sub></i>) of this model is then partitioned in five components:<br><br>

<i>V<sub>G</sub></i> = <i>V<sub>A</sub></i> + <i>V<sub>D</sub></i> + <i>V<sub>AA</sub></i> + <i>V<sub>AD</sub></i> + <i>V<sub>DD</sub></i><br><br>

Where <i>V<sub>AD</sub></i> is the additive-by-dominance variance and <i>V<sub>DD</sub></i> the dominance-by-dominance variance. We use the least square approach described in Lynch and Walsh (1998) Chapter 5 to derive the different variance components and display their values.
  "),
      easyClose = TRUE,
      size="l",
      footer = tagList(
        modalButton("OK")
      )
    ))
  },ignoreNULL = FALSE)
  ###################################################
  
  ######################################################################################################
  #
  # Synchronized updates of sliders
  #
  
  output$SliderA1<- renderUI({
    if(isolate(init$idx)==0){
      init$idx=1
      a=4
    }else{a=input$sliderA1}
    label=ifelse(input$model=="AD","Genotypic value <i>a</i> ","Locus A genotypic value <i>a</i><sub>A</sub>")
    sliderInput("sliderA1",label =  HTML(label), min = -10, max = 10, step = 1,value=a)
  })
  
  ### Synchronize  Additive value's name to the model display ############# DEPRECATED
   # observeEvent(input$model, {
   #   label=ifelse(input$model=="AD","Genotypic value &alpha;","Locus A genotypic value &alpha;<sub>A</sub>")
   #   updateSliderInput(session, "sliderA1", label=HTML(label))
   # })
  
  ## Handsontable ####
  observe({
    if (!is.null(input$hot)) {
      values[["previous"]] <- isolate(values[["GV"]])
      GV = hot_to_r(input$hot)
    } else {
      if (is.null(values[["GV"]]))
        GV <- GV
      else
        GV <- values[["GV"]]
    }
    values[["GV"]] <- GV
  })
  
  output$hot <- renderRHandsontable({
    GV <- values[["GV"]]
    row.names(GV)=c("B\u2082B\u2082","B\u2081B\u2082","B\u2081B\u2081"); colnames(GV) = c("A\u2082A\u2082","A\u2081A\u2082","A\u2081A\u2081")
    if (!is.null(GV))
      rhandsontable(GV) #rowHeaders =  c("B2B2","B1B1","B1B1"), colHeaders = c("A2A2","A1A1","A1A1")
  })
  
  
  ### Output Var A  #############
  output$varianceA <- renderUI({
    if(input$model=="AD"){
      HTML(paste("Additive variance (<b><i>V<sub>A</sub></i></b>) :", format(v$Va, digits = 2, nsmall = 2)))
    }else{
      HTML(paste("Additive variance (<b><i>V<sub>A</sub></i></b>) : ", format(v$Va, digits = 2, nsmall = 2), " (locus A = ",format(v$Va1, digits = 2, nsmall = 2),", locus B = ", format(v$Va2, digits = 2, nsmall = 2),")",sep=""))
    }
  })
  ### Output Var D  #############
  output$varianceD <- renderUI({
    if(input$model=="AD"){
      HTML(paste("Dominance variance (<b><i>V<sub>D</sub></i></b>) :", format(v$Vd, digits = 2, nsmall = 2)))
    }else{
      HTML(paste("Dominance variance (<b><i>V<sub>D</sub></i></b>) : ", format(v$Vd, digits = 2, nsmall = 2), " (locus A = ",format(v$Vd1, digits = 2, nsmall = 2),", locus B = ", format(v$Vd2, digits = 2, nsmall = 2),")",sep=""))
    }
  })
  ### Output Var AA  #############
  output$varianceAA <- renderUI({
    HTML(paste("Additive-by-Additive variance (<b><i>V<sub>AA</sub></i></b>) :", format(v$Vaa, digits = 2, nsmall = 2)))
  })
  ### Output Var AD  #############
  output$varianceAD <- renderUI({
    HTML(paste("Additive-by-Dominance variance (<b><i>V<sub>AD</sub></i></b>) :", format(v$Vad, digits = 2, nsmall = 2)))
  })
  ### Output Var DD  #############
  output$varianceDD <- renderUI({
    HTML(paste("Dominance-by-Dominance variance (<b><i>V<sub>DD</sub></i></b>) :", format(v$Vdd, digits = 2, nsmall = 2)))
  })
  ### Output Var G  #############
  output$varianceG <- renderUI({
    if(input$model=="AA"){
      HTML(paste("Genotypic variance (<b><i>V<sub>G</sub> = V<sub>A</sub> + V<sub>AA</sub></i></b>) :", format(v$Vg, digits = 2, nsmall = 2)))
    }else if(input$model=="AD"){
      HTML(paste("Genotypic variance (<b><i>V<sub>G</sub> = V<sub>A</sub> + V<sub>D</sub></i></b>) :", format(v$Vg, digits = 2, nsmall = 2)))
    }else{
      HTML(paste("Genotypic variance (<b><i>V<sub>G</sub> = V<sub>A</sub> + V<sub>D</sub> + V<sub>I</sub></i></b>) :", format(v$Vg, digits = 2, nsmall = 2)))
    }
  })
  ### Output Var A  / Var G #############
  output$varianceRatio <- renderUI({
    HTML(paste("<b><i>V<sub>A</sub> / V<sub>G</sub></i></b> :", format(v$Va / v$Vg, digits = 2, nsmall = 2)))
  })
  
  ### Output Ref  #############
  output$References <- renderUI({
    if(input$model=="AD"){
      HTML("<b>References</b><br>Falconer, D.S., and Mackay T.F.C. (1996). Introduction to quantitative Genetics, Ed. 4th. Longmans Green, Harlow, Essex.")
    }else if(input$model=="AA"){
      HTML("<b>References</b><br>Falconer, D.S., and Mackay T.F.C. (1996). Introduction to quantitative Genetics, Ed. 4th. Longmans Green, Harlow, Essex.<br><br>
           Mäki-Tanila A., Hill W.G. (2014). Influence of gene interaction on complex trait variation with multi-locus models. <i>Genetics</i>, 198(1):355-367.")
    }else{
      HTML("<b>References</b><br>Falconer, D.S., and Mackay T.F.C. (1996). Introduction to quantitative Genetics, Ed. 4th. Longmans Green, Harlow, Essex.<br><br>
           Mäki-Tanila A., Hill W.G. (2014). Influence of gene interaction on complex trait variation with multi-locus models. <i>Genetics</i>, 198(1):355-367.<br><br>
           Lynch, M. and Walsh, B. (1998). Genetics and Analysis of Quantitative Traits. Sinauer Associates.")
    }
  })
  
  ### Output Main Plot  #############
  output$plot <- renderPlot({
    req(input$sliderA1)
    p  <- input$sliderP
    a1  <- input$sliderA1
    
    if(input$model=="AA" | input$model=="Perso"){
      
      par(mar=c(5.1, 4.1, 0.5, 2.1))
      
      q  <- input$sliderQ
      x <- c(-1,0,1)#c(0,1,2) #xa Genotype coding
      
      if(input$model=="AA"){
        a2  <- input$sliderA2
        aa  <- input$sliderAA
        #Falconer model (-a 0 a)
        v$mu <- round(a1*(p-(1-p))+a2*(q-(1-q))+aa*(1-2*p-2*q+4*p*q),digits = 3)
        v$beta1 <- beta1 <- a1+aa*(2*q-1)
        v$beta2 <- beta2 <- a2+aa*(2*p-1)
        H1 <- 2*p*(1-p)
        H2 <- 2*q*(1-q)
        # Variance components
        v$Va1 <- H1*(beta1^2)
        v$Va2 <- H2*(beta2^2)
        v$Va <- v$Va1 + v$Va2
        v$Vaa <- H1*H2*(aa^2)
        v$Vg <- v$Va + v$Vaa
        # For plotting
        X_AA <- x%*%t(x) * aa
        Geno_Mat <- X_AA + matrix(x*a1,nrow = 3,ncol = 3,byrow = T) + matrix(x*a2,nrow = 3,ncol = 3,byrow = F)
      }else if(input$model=="Perso"){
        Geno_Mat <- values[["GV"]]
        freqA <- c( (1-p)**2, 2*p*(1-p), p**2 )
        freqB <- c( (1-q)**2, 2*q*(1-q), q**2 )
        Geno_freq=freqB%*%t(freqA)
        glsq = Compute_GeneticVariances(GV = Geno_Mat,p = p,q = q)
        # Mu, Alpha and Variance components
        v$mu <- glsq$mu
        v$beta1 <- glsq$alpha[1]
        v$beta2 <- glsq$alpha[2]
        v$Va1 <- glsq$var$VA[[1]][1]
        v$Va2 <- glsq$var$VA[[1]][2]
        v$Va <- glsq$var$VA[[2]]
        v$Vd1 <- glsq$var$VD[[1]][1]
        v$Vd2 <- glsq$var$VD[[1]][2]
        v$Vd <- glsq$var$VD[[2]]
        v$Vaa <- glsq$var$VAA
        v$Vad <- glsq$var$VAD
        v$Vdd <- glsq$var$VDD
        v$Vg <- glsq$var$VG
      }
      
      w1 <- c((1-q)^2*(1-p)^2, (1-q)^2*2*p*(1-p),(1-q)^2*p^2) #With B2B2
      w2 <- c(2*q*(1-q)*(1-p)^2, 2*q*(1-q)*2*p*(1-p),2*q*(1-q)*p^2) #With B1B2
      w3 <- c(q^2*(1-p)^2, q^2*2*p*(1-p),q^2*p^2) #With B1B1
      
      cex.val1 <- (1 + w1)^2
      cex.val2 <- (1 + w2)^2
      cex.val3 <- (1 + w3)^2
      
      lm1 <- lm(Geno_Mat[1,] ~ x,weights = w1)#With B2B2
      lm2 <- lm(Geno_Mat[2,] ~ x,weights = w2)#With B1B2
      lm3 <- lm(Geno_Mat[3,] ~ x,weights = w3)#With B1B1
      
      ylim=c(min(c(Geno_Mat,as.vector(predict(lm1)),as.vector(predict(lm2)),as.vector(predict(lm3)))),max(c(Geno_Mat,as.vector(predict(lm1)),as.vector(predict(lm2)),as.vector(predict(lm3)))))
      plot(x, Geno_Mat[1,], cex = cex.val1, ylim = ylim, xaxt = "n", xlab = "Genotype (effect allele counts)", ylab = "Genotypic value",pch=16,cex.lab=1.3,col="black")
      axis(1, at = x, labels = c(expression(A[2]*A[2]~(0)), expression(A[1]*A[2]~(1)), expression(A[1]*A[1]~(2))),cex.lab=1.5)
      legend(x="topleft",legend = c(expression(B[1]*B[1]),expression(B[1]*B[2]),expression(B[2]*B[2])),col=c("blue","red","black"),lty=1,lwd=2,title=expression(bold("Locus B genotype")),bty="n")
      abline(lm1,col="black")
      points(x=x,y=Geno_Mat[2,],pch=16,col="red",cex = cex.val2,lwd=2) #with B1B2
      abline(lm2,col="red")
      points(x=x,y=Geno_Mat[3,],pch=16,col="blue",cex = cex.val3,lwd=2) #with B1B1
      abline(lm3,col="blue")
      
    }else{ #Dominance model AD
      
      par(mar=c(5.1, 4.1, 0.5, 5.1))
      
      d  <- input$sliderD
      v$mu <- mu <- a1*(2*p - 1) + 2*d*p*(1-p)#Population Mean
      v$beta1 <- beta1 <- a1 + d*(1-2*p)
      # Genotype values
      gAA  <- -a1
      gAB  <- d
      gBB  <- a1
      #varG <- (1-p)^2 * (gAA-mu)^2 + 2*p*(1-p)*(gAB-mu)^2 + p^2*(gBB-mu)^2
      # Breeding values for the 3 genotypes
      bvAA <- -2*p*beta1
      bvAB <- (1-2*p)*beta1
      bvBB <- 2*(1-p)*beta1
      # varBV <- (1-p)^2 * (bvAA)^2 + 2*p*(1-p)*(bvAB)^2 + p^2*(bvBB)^2
      # Variance components
      v$Va <- 2*p*(1-p)*beta1^2
      v$Vd <- (2*p*(1-p)*d)^2
      v$Vg <- v$Va + v$Vd
      # For plotting
      x <- c(0,1,2)
      y <- c(-a1,d,a1)
      w <- c((1-p)^2, 2*p*(1-p),p^2)
      lm <- lm(y ~ x,weights = w)
      pred <- as.vector(predict(lm))
      cex.val <- (1 + w)^2
      ylim=c(min(c(-10,pred)),max(c(10,pred)))
      plot(x, y, cex = cex.val, ylim = ylim, xaxt = "n", xlab = "Genotype (effect allele counts)", ylab = "Value",pch=16,col="blue",cex.lab=1.3)
      axis(1, at = c(0, 1, 2), labels = c(expression(A[2]*A[2]~(0)), expression(A[1]*A[2]~(1)), expression(A[1]*A[1]~(2))))
      points(y=pred,x=seq(0,2),col="blue",cex=2)
      abline(lm,col="blue")
      axis(4, at = c(pred,mu), labels = round(c(bvAA,bvAB,bvBB,0),digits = 2))
      mtext("Deviation from population mean", side = 4, line = 3, cex=1.3)
      points(y=mu,x=(mu-lm$coefficient[1])/lm$coefficient[2], pch=3, col="black",cex=2,lwd=4)
      legend(x=ifelse(a1>1,"topleft","bottomleft"),legend = c("Genotypic value", "Additive (breeding) value"),col=c("blue","blue"),pch=c(16,1),bty="n")
      
      if(d!=0){
        lines(x = c(0,0),y = c(pred[1],y[1]),lty=3)
        lines(x = c(1,1),y = c(pred[2],y[2]),lty=3)
        lines(x = c(2,2),y = c(pred[3],y[3]),lty=3)
      }
      
      if(beta1!=0){
        lines(x = c(0,1),y = c(pred[1],pred[1]),lty=2)
        lines(x = c(1,1),y = c(pred[1],pred[2]),lty=2)
        brackets(x1 = 1,x2 = 1,y1 = pred[1],y2=pred[2],type = 1,h=ifelse(beta1>0,-0.1,0.1))
        
        text(x = 1.1, y = pred[1]-((pred[1]-pred[2])/2), bquote(alpha == .(beta1)), srt = 0, cex = 1,pos = 4)
      }
    }
    
  })
  
  ### Output Plot Dominance  #############
  # output$plotDom <- renderPlot({ #Plot for dominance model
  #   req(input$sliderA1)
  #   p  <- seq(0,1,0.01)
  #   a1  <- input$sliderA1
  #   d  <- input$sliderD
  #   beta1 <- a1 + d*(1-2*p)
  #   Va <- 2*p*(1-p)*beta1^2
  #   Vd <- (2*p*(1-p)*d)^2
  #   # For plotting
  #   par(mar=c(5.1, 5.1, 1, 2.1))
  #   par(mfrow=c(1,3))
  #   plot(x = p,y = Va,type="l",col="black",lwd=2,xlab=bquote("Frequency of "~A[1]~(italic(p))),ylab=bquote(italic(V[A])),cex.lab=1.8,cex.axis=1.8)
  #   points(x=input$sliderP,y=v$Va,pch=3,col="red",cex=4,lwd=4)
  #   
  #   plot(x = p,y = Vd,type="l",col="black",lwd=2,xlab=bquote("Frequency of "~A[1]~(italic(p))),ylab=bquote(italic(V[D])),cex.lab=1.8,cex.axis=1.8)
  #   points(x=input$sliderP,y=v$Vd,pch=3,col="red",cex=4,lwd=4)
  #   
  #   plot(x = p,y = Va/(Va+Vd),type="l",col="black",lwd=2,xlab=bquote("Frequency of "~A[1]~(italic(p))),ylab=bquote(italic(V[A]/V[G])),cex.lab=1.8,cex.axis=1.8)
  #   points(x=input$sliderP,y=v$Va/(v$Va+v$Vd),pch=3,col="red",cex=4,lwd=4)
  #   
  # })
  
  output$plotDom <- renderPlot({ #Plot for dominance model
    req(input$sliderA1)
    p  <- seq(0,1,0.01)
    a1  <- input$sliderA1
    d  <- input$sliderD
    beta1 <- a1 + d*(1-2*p)
    Va <- 2*p*(1-p)*beta1^2
    Vd <- (2*p*(1-p)*d)^2
    # For plotting
    par(mar=c(5.1, 5.1, 2.1, 2.1))
    par(mfrow=c(1,2))
    plot(x = p,y = Va+Vd,type="l",col="black",lwd=2, lty=4,xlab=bquote("Frequency of "~A[1]~(italic(p))),ylab="Genetic variance",cex.lab=1.8,cex.axis=1.8)
    points(x=p,y=Va,type="l",col="black",lwd=2)
    points(x=p,y=Vd,type="l",col="black",lwd=2, lty=3)
    abline(v = input$sliderP,col="red",lwd=2)
    legend(x="topright",legend = c(expression(V[A]),expression(V[D]),expression(V[G])),lty = c(1,3,4),lwd = 2, bty = "n",cex=1.5)
    #title(main = "A.", font=3, adj = 0,mgp = c(2.25,1,0),cex.main = 2.5)
    
    plot(x = p,y = Va/(Va+Vd),type="l",col="black",lwd=2,xlab=bquote("Frequency of "~A[1]~(italic(p))),ylab=bquote(italic(V[A]/V[G])),cex.lab=1.8,cex.axis=1.8)
    points(x=input$sliderP,y=v$Va/(v$Va+v$Vd),pch=3,col="red",cex=4,lwd=4)
    #title(main = "B.", font=3, adj = 0,mgp = c(2.25,1,0),cex.main = 2.5)
  })
  
  ### Output ContourPlot AA  #############
  output$contourplot <- renderPlot({
    req(input$sliderA1)
    a1  <- input$sliderA1
    a2  <- input$sliderA2
    aa  <- input$sliderAA
    # For plotting
    x=y=seq(0,1,0.01)
    z=z2=matrix(0,nrow = length(x),ncol=length(y))
    for(i in 1:length(x)){
      for(j in 1:length(y)){
        p <- x[i]
        q <- y[j]
        H1=2*p*(1-p)
        H2=2*q*(1-q)
        
        beta1=a1+aa*(2*q-1)
        beta2=a2+aa*(2*p-1)
        
        z[i,j] <- H1*(beta1^2) + H2*(beta2^2)
        z2[i,j] <- H1*H2*(aa^2)
      }
    }
    par(mar=c(5.1, 5.1, 4.1, 2.1))
    par(mfrow=c(1,3))
    image2D(z,x,y,contour=T,rasterImage=T,clab=expression(bold(V[A])),xlab=bquote("Frequency of "~A[1]~(italic(p))),ylab=bquote("Frequency of "~B[1]~(italic(q))),cex.lab=1.8,cex.axis=1.5,colkey=list(cex.axis=1.5,cex.clab=1.8))
    points(x=input$sliderP,y=input$sliderQ,pch=3,col="white",cex=4,lwd=4)
    
    image2D(z2,x,y,contour=T,rasterImage=T,clab=expression(bold(V[AA])),xlab=bquote("Frequency of "~A[1]~(italic(p))),ylab=bquote("Frequency of "~B[1]~(italic(q))),cex.lab=1.8,cex.axis=1.5,colkey=list(cex.axis=1.5,cex.clab=1.8))
    points(x=input$sliderP,y=input$sliderQ,pch=3,col="white",cex=4,lwd=4)
    
    image2D(z/(z+z2),x,y,contour=T,rasterImage=T,clab=bquote(bold(frac(V[A],V[G]))),xlab=bquote("Frequency of "~A[1]~(italic(p))),ylab=bquote("Frequency of "~B[1]~(italic(q))),cex.lab=1.8,cex.axis=1.5,colkey=list(cex.axis=1.5,cex.clab=1.8))
    points(x=input$sliderP,y=input$sliderQ,pch=3,col="white",cex=4,lwd=4)
  })
  
  ############################################################################################
  #
  # TABLES
  # 
  ### Output Main Table  #############
  output$table <- renderUI({
    req(input$sliderA1)
    p <- input$sliderP
    
    if(input$model=="AD"){ #Additive Dominance model AD
      q <- 1-input$sliderP
      a  <- input$sliderA1
      d  <- input$sliderD
      beta <- v$beta1
      mu <- round(v$mu,digits = 3)
      rowNames <- c("Frequencies","Assigned values","Genotypic value","Additive (breeding) value","Dominance deviation")
      colNames <- c("A<sub>2</sub>A<sub>2</sub>","A<sub>1</sub>A<sub>2</sub>","A<sub>1</sub>A<sub>1</sub>")
      data <- round(c( q**2, 2*p*q, p**2,
                       -a, d, a,
                       -2*p*(a+q*d), a*(q-p)+d*(1-2*p*q), 2*q*(a-p*d),
                       -2*p*beta, (q-p)*beta, 2*q*beta,
                       -2*(p**2)*d, 2*p*q*d, -2*(q**2)*d
      ), digits = 2)
      
      HTML(
        htmlTable(matrix(data, 
                         ncol=3, byrow = TRUE),
                  header =  colNames,
                  rnames = rowNames,
                  rgroup = c("","Deviations from population mean:"),
                  n.rgroup = c(2,3),
                  cgroup = c("Genotypes"),
                  n.cgroup = c(3), 
                  caption=paste0("<br>Population mean: M = ",mu,"<br>Average effect of gene-substitution: &#120572; = ",beta)) 
      )
    }else if(input$model == "AA" | input$model == "Perso"){#Additive Additive-by-Additive model A-AA
      q <- input$sliderQ
      freq.1 <- c((1-p)**2,2*p*(1-p),p**2)
      freq.2 <- c((1-q)**2,2*q*(1-q),q**2)
      geno.freq <- round(freq.1%*%t(freq.2),digits = 2)
      beta1 <- round(v$beta1,digits = 2)
      beta2 <- round(v$beta2, digits = 2)
      mu <- round(v$mu,digits = 2)
      if(input$model == "AA"){
        a1  <- input$sliderA1
        a2  <- input$sliderA2
        aa  <- input$sliderAA
        #Maki-Tanila and Hill model (0 a 2a)
        # beta1 <- a1 + 2*q*aa
        # beta2 <- a2 + 2*p*aa
        
        #Falconer model (-a 0 a)
        geno.val <- sapply(X=c(-1,0,1),FUN=function(x,y,a){return(a[1]*x+y*a[2]+y*x*a[3])},y=c(-1,0,1),a=c(a1,a2,aa))
      }else{
        geno.val <- values[["GV"]]
      }
      
      
      geno.1=c("A<sub>2</sub>A<sub>2</sub>","A<sub>1</sub>A<sub>2</sub>","A<sub>1</sub>A<sub>1</sub>")
      geno.2=c("<b>B<sub>2</sub>B<sub>2</sub></b>","<b>B<sub>1</sub>B<sub>2</sub></b>","<b>B<sub>1</sub>B<sub>1</sub></b>")
      
      data <- matrix(paste(as.vector(t(geno.freq))," (", as.vector(t(geno.val)),")"),nrow = 3,ncol = 3,byrow = T)
      HTML(
        htmlTable(data,
                  header =  geno.1,
                  rnames = geno.2,
                  rowlabel = paste0("q=",q,"\\p=",p),
                  align.header = c('l',rep('c',3)),
                  caption=paste0("<br>Population mean : M = ",mu,"<br>Locus A average effect of gene-substitution: &#120572;<sub>A</sub> = ",beta1,"<br>Locus B average effect of gene-substitution: &#120572;<sub>B</sub> = ",beta2,"<br><br>Genotype frequencies (genotypic values) are :"))
      )
      ##########################################
    }
  })#End table UI
}#End server

#############################################
# Run the app
#
shinyApp(ui = ui, server = server)





