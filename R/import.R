# import of mass spectrometry proteomics data analysis software reports.
import <-
function(ms_filenames, ms_filetype, ...) UseMethod("import")

import.default <- function(ms_filenames, ms_filetype, concentration_filename=NA, averageruns=FALSE, sumruns=FALSE, mprophet_cutoff=0.01, mprophet_protein_cutoff=NA, openswath_superimpose_identifications=FALSE, openswath_superimpose_protein_identifications=FALSE, openswath_replace_run_id=FALSE, openswath_filtertop=FALSE, openswath_removedecoys=TRUE, peptideprophet_cutoff=0.95, abacus_column="ADJNSAF", pepxml2csv_runsplit="~", diau_peptide_column="Top6tra", diau_protein_column="iBAQ", fasta=NA, ...) {
	transition_intensity <- peptide_intensity <- protein_intensity <- NULL

	if (!ms_filetype %in% c("openswath","mprophet","openmslfq","skyline","abacus","pepxml2csv","diau_peptide","diau_protein","maxquant_peptide")) {
		stop("Please select a valid filetype. Options:  \"openswath\", \"mprophet\", \"openmslfq\", \"skyline\", \"abacus\", \"pepxml2csv\", \"diau_peptide\", \"diau_protein\", \"maxquant_peptide\"")
	}
	
	# ms_filenames must be provided as vector and are converted to a list
	ms_filenames <- as.list(ms_filenames)

	# Skyline import is facilitated by the Transitions Results report
	if (ms_filetype=="skyline") {
		# Skyline specfic adapter
		data.ms <- skyline_converter.import(ms_filenames)
		if (averageruns==TRUE){
			data.ms <- averageruns.import(data.ms,target="transition_intensity")
			data.ms$run_id <- "averaged"
		}
		if (sumruns==TRUE){
			data.ms <- sumruns.import(data.ms,target="transition_intensity")
			data.ms$run_id <- "summed"
		}
		# Skyline exports all selected transitions. But we only want to use those with an associated MS2 intensity
		data <- subset(data.ms,is.finite(transition_intensity))
	}
	# OpenSWATH import is facilitated by using the output from either mProphet or the TOPPtool OpenSwathFeatureXMLToTSV
	else if (ms_filetype=="openswath") {
		# OpenSWATH specific adapter
		data.ms <- openswath_converter.import(ms_filenames,mprophet_cutoff=mprophet_cutoff,mprophet_protein_cutoff=mprophet_protein_cutoff,openswath_superimpose_identifications=openswath_superimpose_identifications, openswath_superimpose_protein_identifications=openswath_superimpose_protein_identifications, openswath_replace_run_id=openswath_replace_run_id, openswath_filtertop = openswath_filtertop, openswath_removedecoys=openswath_removedecoys)

		if (averageruns==TRUE){
			data.ms <- averageruns.import(data.ms,target="transition_intensity")
			data.ms$run_id <- "averaged"
		}
		if (sumruns==TRUE){
			data.ms <- sumruns.import(data.ms,target="transition_intensity")
			data.ms$run_id <- "summed"
		}
		data <- subset(data.ms,is.finite(transition_intensity))
	}
	# mProphet import is facilitated by using the all_peakgroups.xls output
	else if (ms_filetype=="mprophet") {
		# mProphet specific adapter
		data.ms <- mprophet_converter.import(ms_filenames,mprophet_cutoff=mprophet_cutoff, openswath_replace_run_id=openswath_replace_run_id, openswath_filtertop = openswath_filtertop, openswath_removedecoys=openswath_removedecoys)

		if (averageruns==TRUE){
			data.ms <- averageruns.import(data.ms,target="transition_intensity")
			data.ms$run_id <- "averaged"
		}
		if (sumruns==TRUE){
			data.ms <- sumruns.import(data.ms,target="transition_intensity")
			data.ms$run_id <- "summed"
		}
		data <- subset(data.ms,is.finite(transition_intensity))
	}
	# OpenMS LFQ import is facilitated from the TOPPtool ProteinQuantifier and the resulting peptides.csv file
	else if (ms_filetype=="openmslfq") {
		# OpenMS LFQ specific adaptor
		data.ms <- openmslfq_converter.import(ms_filenames)
		if (averageruns==TRUE){
			data.ms <- averageruns.import(data.ms,target="peptide_intensity")
			data.ms$run_id <- "averaged"
		}
		if (sumruns==TRUE){
			data.ms <- sumruns.import(data.ms,target="peptide_intensity")
			data.ms$run_id <- "summed"
		}
		data <- subset(data.ms,is.finite(peptide_intensity))
	}
	# Abacus protein import is facilitated from Abacus with default output settings
	else if (ms_filetype=="abacus") {
		data.ms <- abacus_converter.import(ms_filenames, peptideprophet_cutoff=peptideprophet_cutoff, abacus_column=abacus_column)
		if (averageruns==TRUE){
			data.ms <- averageruns.import(data.ms,target="protein_intensity")
			data.ms$run_id <- "averaged"
		}
		if (sumruns==TRUE){
			data.ms <- sumruns.import(data.ms,target="protein_intensity")
			data.ms$run_id <- "summed"
		}
		data <- subset(data.ms,is.finite(protein_intensity))
	}
	else if (ms_filetype=="pepxml2csv") {
		data.ms <- pepxml2csv_converter.import(ms_filenames, peptideprophet_cutoff=peptideprophet_cutoff, pepxml2csv_runsplit=pepxml2csv_runsplit)
		if (averageruns==TRUE){
			data.ms <- averageruns.import(data.ms,target="peptide_intensity")
			data.ms$run_id <- "averaged"
		}
		if (sumruns==TRUE){
			data.ms <- sumruns.import(data.ms,target="peptide_intensity")
			data.ms$run_id <- "summed"
		}
		data <- subset(data.ms,is.finite(peptide_intensity))
	}
	else if (ms_filetype=="diau_peptide") {
		data.ms <- diau_peptide_converter.import(ms_filenames, peptideprophet_cutoff=peptideprophet_cutoff, diau_peptide_column=diau_peptide_column)
		if (averageruns==TRUE){
			data.ms <- averageruns.import(data.ms,target="peptide_intensity")
			data.ms$run_id <- "averaged"
		}
		if (sumruns==TRUE){
			data.ms <- sumruns.import(data.ms,target="peptide_intensity")
			data.ms$run_id <- "summed"
		}
		data <- subset(data.ms,is.finite(peptide_intensity))
	}
	else if (ms_filetype=="diau_protein") {
		data.ms <- diau_protein_converter.import(ms_filenames, peptideprophet_cutoff=peptideprophet_cutoff, diau_protein_column=diau_protein_column)
		if (averageruns==TRUE){
			data.ms <- averageruns.import(data.ms,target="peptide_intensity")
			data.ms$run_id <- "averaged"
		}
		if (sumruns==TRUE){
			data.ms <- sumruns.import(data.ms,target="peptide_intensity")
			data.ms$run_id <- "summed"
		}
		data <- subset(data.ms,is.finite(protein_intensity))
	}
	# MaxQuant peptide import is facilitated from peptides.txt
	else if (ms_filetype=="maxquant_peptide") {
		data.ms <- maxquant_peptide_converter.import(ms_filenames)
		if (averageruns==TRUE){
			data.ms <- averageruns.import(data.ms,target="peptide_intensity")
			data.ms$run_id <- "averaged"
		}
		if (sumruns==TRUE){
			data.ms <- sumruns.import(data.ms,target="peptide_intensity")
			data.ms$run_id <- "summed"
		}
		data <- subset(data.ms,is.finite(peptide_intensity))
	}

	if (!is.na(fasta)) {
		proteins <- read.fasta(file = fasta, seqtype = "AA", as.string = TRUE, seqonly = FALSE, strip.desc = TRUE)
		sequences<-sapply(proteins,function(X){return(X[1][[1]])})
		index<-data.frame("peptide_sequence"=unique(data$peptide_sequence),stringsAsFactors=FALSE)
		index$protein_id<-sapply(index$peptide_sequence,function(X){mappeptides.import(X,sequences)})
		data<-merge(data[,!(names(data) %in% c("protein_id"))],index,by="peptide_sequence")
		data[which(data$protein_id=="0/"),]$protein_id<-NA
	}

	if ("peptide_sequence" %in% names(data)) {
		# peptide_sequence is stripped and can only contain the 20 natural AA without any modifications
		data$peptide_sequence<-stripsequence.import(data$peptide_sequence)
	}

	# if no concentration file is supplied, "?" is used to indicate absence
	if (!is.na(concentration_filename)) {
		data.conc <- read.csv(concentration_filename, stringsAsFactors=FALSE)
		data.conc$concentration<-as.numeric(data.conc$concentration)
		if ("protein_id" %in% names(data.conc)) {
			# if proteins were manually (e.g. with Skyline quantified)
			if ("run_id" %in% names(data.conc)) {
				data <- merge(data,data.conc, by.x=c("run_id", "protein_id"), by.y=c("run_id", "protein_id"), all.x = TRUE)
			}
			else {
				data <- merge(data,data.conc, by.x="protein_id", by.y="protein_id", all.x = TRUE)
			}
		}
		else if ("peptide_id" %in% names(data.conc)) {
			# if only the spiked-in peptides are known
			if ("run_id" %in% names(data.conc)) {
				data <- merge(data, data.conc, by.x=c("run_id", "peptide_id"), by.y=c("run_id", "peptide_id"), all.x = TRUE)
			}
			else {
				data <- merge(data, data.conc, by.x="peptide_id", by.y="peptide_id", all.x = TRUE)
			}
		}

		data$concentration[ is.na(data$concentration) ] <- "?"
	}
	else {
		data$concentration <- "?"
	}

	return(data)
}

skyline_converter.import <- function(ms_filenames, ...)  {
	data.files = lapply(ms_filenames, read.csv, stringsAsFactors=FALSE)
	data.import <- do.call("rbind", data.files)

	# Skyline headers
	data.ms <- data.import[,c("ReplicateName","ProteinName","PeptideSequence","FragmentIon","PeptideSequence","PrecursorCharge","Area")]
	
	# aLFQ headers
	names(data.ms) <- c("run_id","protein_id","peptide_id","transition_id","peptide_sequence","precursor_charge","transition_intensity")
	
	# replace Skyline NA with R NA values.
	data.ms <- replace(data.ms, data.ms=="#N/A", NA)
	
	data.ms$precursor_charge <- as.numeric(data.ms$precursor_charge)
	data.ms$transition_intensity <- as.numeric(data.ms$transition_intensity)

	return(data.ms)
}

openswath_converter.import <- function(ms_filenames, mprophet_cutoff, mprophet_protein_cutoff,openswath_superimpose_identifications, openswath_superimpose_protein_identifications, openswath_replace_run_id, openswath_filtertop, openswath_removedecoys,...)  {
	transition_code <- transition_code_intensity <- transition_id <- peptide_id <- transition_fragment <- m_score <- peak_group_rank <- decoy <- NULL

	data.files = lapply(ms_filenames, read_table.import, openswath_replace_run_id=openswath_replace_run_id)
	data.import <- data.table(do.call("rbind", data.files))


	if ("align_origfilename" %in% names(data.import)) {
		# OpenSWATH Feature Alignment headers
		data.ms <- data.import[,c("align_origfilename","ProteinName","FullPeptideName","m_score","peak_group_rank","Sequence","Charge","aggr_Fragment_Annotation","aggr_Peak_Area","decoy"), with = FALSE]

		# aLFQ headers
		setnames(data.ms,c("align_origfilename","ProteinName","FullPeptideName","m_score","peak_group_rank","Sequence","Charge","aggr_Fragment_Annotation","aggr_Peak_Area","decoy"),c("run_id","protein_id","peptide_id","m_score","peak_group_rank","peptide_sequence","precursor_charge","transition_code","transition_code_intensity","decoy"))
	}
	else {
		# OpenSWATH headers
		data.ms <- data.import[,c("filename","ProteinName","FullPeptideName","m_score","peak_group_rank","Sequence","Charge","aggr_Fragment_Annotation","aggr_Peak_Area","decoy"), with = FALSE]

		# aLFQ headers
		setnames(data.ms,c("filename","ProteinName","FullPeptideName","m_score","peak_group_rank","Sequence","Charge","aggr_Fragment_Annotation","aggr_Peak_Area","decoy"),c("run_id","protein_id","peptide_id","m_score","peak_group_rank","peptide_sequence","precursor_charge","transition_code","transition_code_intensity","decoy"))
	}
	
	setkeyv(data.ms,c("run_id","transition_code","protein_id","peptide_id","m_score","peak_group_rank","peptide_sequence","precursor_charge","decoy"))

	data.ms<-data.ms[, list("transition_fragment"=strsplit(transition_code,";")[[1]],"transition_intensity"=as.numeric(strsplit(transition_code_intensity,";")[[1]])), by=key(data.ms)]

	data.ms[, transition_id:=paste(peptide_id,transition_fragment)]

	# filter mprophet_cutoff on assay-level
	if (openswath_superimpose_identifications) {
		tcset <- unique(subset(data.ms, m_score <= mprophet_cutoff)$transition_code)
		data.ms <- subset(data.ms, transition_code %in% tcset)
	}
	else {
		data.ms <- subset(data.ms, m_score <= mprophet_cutoff)
	}

	# filter mprophet_cutoff on protein-level
	if (openswath_superimpose_protein_identifications) {
		protset <- unique(subset(data.ms, m_score <= mprophet_protein_cutoff)$protein_id)
		data.ms <- subset(data.ms, protein_id %in% protset)
	}

	# filter top peakgroup
	if (openswath_filtertop) {
		data.ms <- subset(data.ms, peak_group_rank == 1)
	}

	# filter decoys
	if (openswath_removedecoys) {
		data.ms <- subset(data.ms, decoy == FALSE)
	}

	return(data.frame(data.ms[,c("run_id","protein_id","peptide_id","transition_id","peptide_sequence","precursor_charge","transition_intensity"), with = FALSE]))
}

mprophet_converter.import <- function(ms_filenames, mprophet_cutoff, openswath_replace_run_id, openswath_filtertop, openswath_removedecoys, ...)  {
	transition_code <- transition_code_intensity <- transition_id <- peptide_id <- transition_fragment <- m_score <- peak_group_rank <- decoy <- NULL

	data.files = lapply(ms_filenames, read_table.import)
	data.import <- data.table(do.call("rbind", data.files))

	# mProphet headers
	data.ms <- data.import[,c("run_id","protein","transition_group_record","m_score","peak_group_rank","transition_group_pepseq","transition_group_charge","transition_code_target","abs_area_code_target","decoy"), with = FALSE]

	# aLFQ headers
	setnames(data.ms,c("run_id","protein","transition_group_record","m_score","peak_group_rank","transition_group_pepseq","transition_group_charge","transition_code_target","abs_area_code_target","decoy"),c("run_id","protein_id","peptide_id","m_score","peak_group_rank","peptide_sequence","precursor_charge","transition_code","transition_code_intensity","decoy"))
	
	setkeyv(data.ms,c("run_id","protein_id","peptide_id","m_score","peak_group_rank","peptide_sequence","precursor_charge","decoy"))

	data.ms<-data.ms[, list("transition_fragment"=strsplit(transition_code,",")[[1]],"transition_intensity"=as.numeric(strsplit(transition_code_intensity,",")[[1]])), by=key(data.ms)]

	data.ms[, transition_id:=paste(peptide_id,transition_fragment)]

	# filter mprophet_cutoff
	data.ms <- subset(data.ms, m_score <= mprophet_cutoff)

	# filter top peakgroup
	if (openswath_filtertop) {
		data.ms <- subset(data.ms, peak_group_rank == 1)
	}

	# filter decoys
	if (openswath_removedecoys) {
		data.ms <- subset(data.ms, decoy == FALSE)
	}

	return(data.frame(data.ms[,c("run_id","protein_id","peptide_id","transition_id","peptide_sequence","precursor_charge","transition_intensity"), with = FALSE]))
}

read_table.import <- function(filename, openswath_replace_run_id = FALSE, ...) {
	# helper function for mProphet to read runs. filenames are converted to the run_id.
	data <- as.data.frame(fread(filename))

	if (openswath_replace_run_id) {
		data$run_id <- filename
	}

	return(data)
}

openmslfq_converter.import <- function(ms_filenames, ...)  {
	# we need to skip the comments in the header
	data.files = lapply(ms_filenames, read.csv, sep = ",",  comment.char = "#", blank.lines.skip = TRUE, stringsAsFactors=FALSE)

	data.trans<-lapply(data.files,openmslfq_transform.import)

	data.ms <- do.call("rbind", data.trans)

	return(data.ms)
}
openmslfq_transform.import <- function(data, ...) {
	n_proteins <- NULL

	# converts short to long format
	# filter proteotypic peptides
	data <- subset(data, n_proteins == 1)
	
	data.trans <- data.frame(run_id = NA, protein_id = NA, peptide_id = NA, peptide_sequence = NA, precursor_charge = NA, peptide_intensity = NA, stringsAsFactors=FALSE)[0,]

	i <- length(names(data)) #-1
	while( i > 4 ) {
		data.trans <- rbind(data.trans,data.frame(run_id = names(data)[i], protein_id = data$protein, peptide_id = data$peptide, peptide_sequence = data$peptide, precursor_charge = data$charge, peptide_intensity = data[,i], stringsAsFactors=FALSE))
	
		i <- i-1
	}

	return(data.trans)
}

maxquant_peptide_converter.import <- function(ms_filenames, ...)  {
	# we need to skip the comments in the header
	data.files = lapply(ms_filenames, read.csv, sep = "\t",  comment.char = "#", blank.lines.skip = TRUE, stringsAsFactors=FALSE)

	data.trans<-lapply(data.files,maxquant_peptide_transform.import)

	data.ms <- do.call("rbind", data.trans)

	return(data.ms)
}
maxquant_peptide_transform.import <- function(data, ...) {
	# converts short to long format
	# filter proteotypic peptides	
	data.trans <- data.frame(run_id = NA, protein_id = NA, peptide_id = NA, peptide_sequence = NA, precursor_charge = NA, peptide_intensity = NA, stringsAsFactors=FALSE)[0,]

	run_intensities<-names(data)[as.vector(sapply(names(data),function(X){substr(X,1,10)=="Intensity."}))]
	for (i in 1:length(run_intensities)) {
		data.trans <- rbind(data.trans,data.frame(run_id = gsub("Intensity.", "", run_intensities[i]), protein_id = data$Proteins, peptide_id = data$Sequence, peptide_sequence = data$Sequence, precursor_charge = data$Charges[1], peptide_intensity = data[,run_intensities[i]], stringsAsFactors=FALSE))	
	}

	return(data.trans)
}

abacus_converter.import <- function(ms_filenames, peptideprophet_cutoff, abacus_column, ...)  {
	data.files = lapply(ms_filenames, read.csv, sep="\t", stringsAsFactors=FALSE)

	data.trans<-lapply(data.files,function(X){abacus_transform.import(X,peptideprophet_cutoff,abacus_column)})

	data.ms <- do.call("rbind", data.trans)

	return(data.ms)
}
abacus_transform.import <- function(data, peptideprophet_cutoff, abacus_column, ...)  {
	PW <- ISFWD <- NULL
	data.re<-ldply(lapply(unlist(strsplit(names(data)[which(substr(names(data),1,3)!="ALL")][grepl("_ID", names(data)[which(substr(names(data),1,3)!="ALL")])],"_ID")),function(X){data.frame("RUN"=X,"PROTID"=data$PROTID,"PROTLEN"=data$PROTLEN,"ISFWD"=data$ISFWD,"PW"=data[,paste(X,"PW",sep="_")],"NUMSPECSTOT"=data[,paste(X,"NUMSPECSTOT",sep="_")],"TOTNSAF"=data[,paste(X,"TOTNSAF",sep="_")],"NUMSPECSUNIQ"=data[,paste(X,"NUMSPECSUNIQ",sep="_")],"UNIQNSAF"=data[,paste(X,"UNIQNSAF",sep="_")],"NUMSPECSADJ"=data[,paste(X,"NUMSPECSADJ",sep="_")],"ADJNSAF"=data[,paste(X,"ADJNSAF",sep="_")])}),data.frame)

	data.filt <- subset(data.re, PW >= peptideprophet_cutoff & ISFWD==1)

	data.ms <- data.filt[,c("RUN","PROTID",abacus_column)]

	names(data.ms) <- c("run_id","protein_id","protein_intensity")
	
	return(data.ms)
}

pepxml2csv_converter.import <- function(ms_filenames, peptideprophet_cutoff, pepxml2csv_runsplit = "~", ...)  {
	probability <- peptide_intensity <- NULL
	data.files = lapply(ms_filenames, read_table.import)
	data.import <- do.call("rbind", data.files)
	data.import[which(is.na(data.import$modified_peptide)),]$modified_peptide<-data.import[which(is.na(data.import$modified_peptide)),]$peptide
	data.import$run_id<-data.frame(do.call(rbind, strsplit(as.vector(data.import$spectrum), split = pepxml2csv_runsplit)),stringsAsFactors=FALSE)[,1]

	# pepxml2csv headers
	if ("probability_ip" %in% names(data.import)) {
		data.ms <- data.table(data.import)[,c("run_id","assumed_charge","peptide","modified_peptide","protein","probability_ip"), with = FALSE]
		setnames(data.ms,c("run_id","assumed_charge","peptide","modified_peptide","protein","probability_ip"),c("run_id","precursor_charge","peptide_sequence","peptide_id","protein_id","probability"))
	}
	else {
		data.ms <- data.table(data.import)[,c("run_id","assumed_charge","peptide","modified_peptide","protein","probability_pp"), with = FALSE]
		setnames(data.ms,c("run_id","assumed_charge","peptide","modified_peptide","protein","probability_pp"),c("run_id","precursor_charge","peptide_sequence","peptide_id","protein_id","probability"))
	}

	data.ms<-subset(data.ms,probability >= peptideprophet_cutoff)
	
	setkeyv(data.ms,c("run_id","protein_id","peptide_id","peptide_sequence","precursor_charge"))
	data.ms[, peptide_intensity:=1]

	data.ms<-data.ms[, list("peptide_intensity"=sum(peptide_intensity)), by=key(data.ms)]

	return(as.data.frame(data.ms))
}

diau_peptide_converter.import <- function(ms_filenames, peptideprophet_cutoff, diau_peptide_column, ...)  {
	data.files = lapply(ms_filenames, read.csv, sep="\t", stringsAsFactors=FALSE)

	data.trans<-lapply(data.files,function(X){diau_peptide_transform.import(X,peptideprophet_cutoff,diau_peptide_column)})

	data.ms <- do.call("rbind", data.trans)

	return(data.ms)
}
diau_peptide_transform.import <- function(data, peptideprophet_cutoff, diau_peptide_column, ...)  {
	data.re<-ldply(lapply(unlist(strsplit(names(data)[grepl("_Spec_Centric_Prob", names(data))],"_Spec_Centric_Prob")),function(X){data.frame("run"=X,"peptide_id"=data$Peptide.Key,"precursor_charge"=data$Charge,"Spec_Centric_Prob"=data[,paste(X,"Spec_Centric_Prob",sep="_")],"Pep_Centric_Prob"=data[,paste(X,"Pep_Centric_Prob",sep="_")],"RT"=data[,paste(X,"RT",sep="_")],"MS1"=data[,paste(X,"MS1",sep="_")],"Top6pepTop6tra.Freq.0.5."=data[,paste(X,"Top6tra.Freq.0.5.",sep="_")],"Top6tra"=data[,paste(X,"Top6tra",sep="_")])}),data.frame)

	data.filt <- subset(data.re, Pep_Centric_Prob >= peptideprophet_cutoff)

	data.ms <- data.filt[,c("run","peptide_id","precursor_charge",diau_peptide_column)]

	names(data.ms) <- c("run_id","peptide_id","precursor_charge","peptide_intensity")

	data.ms$peptide_sequence<-unlist(sapply(data.ms$peptide_id,function(X){gsub("[^[:alpha:]]","",gsub("\\[[^\\]]*\\]", "", X, perl=TRUE))}))
	
	return(data.ms)
}
diau_protein_converter.import <- function(ms_filenames, peptideprophet_cutoff, diau_protein_column, ...)  {
	data.files = lapply(ms_filenames, read.csv, sep="\t", stringsAsFactors=FALSE)

	data.trans<-lapply(data.files,function(X){diau_protein_transform.import(X,peptideprophet_cutoff,diau_protein_column)})

	data.ms <- do.call("rbind", data.trans)

	return(data.ms)
}
diau_protein_transform.import <- function(data, peptideprophet_cutoff, diau_protein_column, ...)  {
	data.re<-ldply(lapply(unlist(strsplit(names(data)[grepl("_Prob", names(data))],"_Prob")),function(X){data.frame("run"=X,"protein_id"=data$Protein.Key,"Prob"=data[,paste(X,"Prob",sep="_")],"peptides"=data[,paste(X,"Peptides",sep="_")],"Top6pepTop6tra.Freq.0.5."=data[,paste(X,"Top6pepTop6tra.Freq.0.5.",sep="_")],"Top6pepTop6tra"=data[,paste(X,"Top6pepTop6tra",sep="_")],"iBAQ"=data[,paste(X,"iBAQ",sep="_")])}),data.frame)

	data.filt <- subset(data.re, Prob >= peptideprophet_cutoff)

	data.ms <- data.filt[,c("run","protein_id",diau_protein_column)]

	names(data.ms) <- c("run_id","protein_id","protein_intensity")
	
	return(data.ms)
}

stripsequence.import <- function(X, ...) {
	return(gsub('[a-z]','',gsub('\\[[^\\[\\]]*\\]','',gsub('\\[[^\\[\\]]*\\]','',gsub('\\([^()]*\\)','',gsub('\\([^()]*\\)','',gsub('\\{[^{}]*\\}','',gsub('\\{[^{}]*\\}','',X)))),perl=TRUE),perl=TRUE),perl=TRUE))
}

mappeptides.import <- function(X, sequences, ...) {
	protein_ids<-names(sequences)[grep(X,sequences)]
	return(paste(length(protein_ids),paste(protein_ids,collapse="/"),sep="/"))
}

averageruns.import <- function(data,target="transition_intensity", ...)  {
	protein_id <- peptide_id <- transition_id <- peptide_sequence <- precursor_charge <- NULL

	if (target=="transition_intensity") {
		data.averaged <- ddply(data, .(protein_id,peptide_id,transition_id,peptide_sequence,precursor_charge), function(X) {data.frame("transition_intensity" = mean(X$transition_intensity, na.rm = TRUE), stringsAsFactors=FALSE)})
	}
	else if (target=="peptide_intensity") {
		data.averaged <- ddply(data, .(protein_id,peptide_id,peptide_sequence,precursor_charge), function(X) {data.frame("peptide_intensity" = mean(X$peptide_intensity, na.rm = TRUE, stringsAsFactors=FALSE))})
	}
	else if (target=="protein_intensity") {
		data.averaged <- ddply(data, .(protein_id), function(X) {data.frame("protein_intensity" = mean(X$protein_intensity, na.rm = TRUE, stringsAsFactors=FALSE))})
	}
	
	return(data.averaged)
}

sumruns.import <- function(data,target="transition_intensity", ...)  {
	protein_id <- peptide_id <- transition_id <- peptide_sequence <- precursor_charge <- NULL

	if (target=="transition_intensity") {
		data.summed <- ddply(data, .(protein_id,peptide_id,transition_id,peptide_sequence,precursor_charge), function(X) {data.frame("transition_intensity" = sum(X$transition_intensity, na.rm = TRUE), stringsAsFactors=FALSE)})
	}
	else if (target=="peptide_intensity") {
		data.summed <- ddply(data, .(protein_id,peptide_id,peptide_sequence,precursor_charge), function(X) {data.frame("peptide_intensity" = sum(X$peptide_intensity, na.rm = TRUE, stringsAsFactors=FALSE))})
	}
	else if (target=="protein_intensity") {
		data.summed <- ddply(data, .(protein_id), function(X) {data.frame("protein_intensity" = sum(X$protein_intensity, na.rm = TRUE, stringsAsFactors=FALSE))})
	}
	
	return(data.summed)
}
