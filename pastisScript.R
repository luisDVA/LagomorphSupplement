# USING PASTIS TO INCORPORATE MISSING TAXA AT THE PHYLOGENETIC ASSEMBLY STAGE
library (ape)
library (pastis)

# "read.GenBank" function used to:
# download sequences using accesion list from Ge et al. (2013)
# download additional sequences
##  Lepus alleni                HQ596458  
##  Lepus arcticus	            HQ596461	
##  Lepus callotis            	HQ596469	
##  Lepus castroviejoi        	JN037350	
##  Lepus coreanus            	AB687532	
##  Lepus corsicanus          	JN037359	
##	Lepus flavigularis        	HQ596475	
##	Lepus granatensis	          HQ596476	
##	Lepus insularis           	HM222713	
##	Lepus othus	                HQ596479	
##	Lepus tolai	                AY649610	
##	Sylvilagus robustus       	HQ143449	
##	Sylvilagus transitionalis 	AF034256	

# aligned sequences using Jalview (Muscle alignment, nucleotide preset)

# read inputs
ge13phylo <- read.tree ("ge13full.txt") # from Ge et al. (2013)
taxonClade <- read.table ("taxonClade.txt", header=TRUE)
# references for clade assignment provided in "cladePlacement.pdf")
# sequence alignment file and missing clade file must be in workspace folder

# prepare PASTIS input
# PASTIS uses the phylogeny to create flexible constraints including placement of missing clades
LagPast <- read_input (ge13phylo, taxonClade, missing_clades="missingClades.txt", sequences="allSeqs.fas")

# modify output template to assign genes and model following Ge et al. (2013)
LagPast$output_template <- "#NEXUS
                            \n\n Begin DATA; 
                            \n\t Dimensions ntax=<ntax> nchar=<nchar>;
                            \n\t Format datatype=DNA gap=- missing=?; 
                            \n\t Matrix \n<sequences>   \n\n; 
                            \n\n\n begin MRBAYES; \n\n 
                            \n\n\t charset cytb = 1-1619;
                            \n\t charset ND4 = 1620-19485;
                            \n\t charset X12s = 19486-37397;
                            \n\n\t partition genes = 3: cytb,ND4,X12s;
                            \n\t set partition = genes;
                            \n lset applyto=(all) nst=6 rates=invgamma;                              
                            \n\nunlink shape=(all) tratio=(all) statefreq=(all) revmat=(all) pinvar=(all); 
                            \n<constraints>\n  \n  
                            \n  set usebeagle=no Beaglesse=no; 
                            \n\n prset brlenspr=clock:birthdeath; 
                            \n prset Extinctionpr = Fixed(0); 
                            \n prset Speciationpr=exponential(1); 
                            \n prset clockvarpr=ibr; 
                            \n prset ibrvarpr=exponential(10); 
                            \nmcmcp nruns=1 nchains=1 ngen=1400000 samplefreq=1000; \nmcmc; 
                            \n\nsumt filename=<outputfile> burnin=400000 contype=halfcompat;\n\nend; "

# create MrBayes input file
pastis_main (LagPast, output_file="lagomPastB.nex")
# execute input file in MrBayes 3.2 

# post-processing:
# dated the consensus tree using PATHd8 and the same dating information in Ge et al. (2013)
# read consensus tree into R
# final output provided as "lphyloPASTIS.nex"