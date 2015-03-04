normalizeBetweenRuns <- function(x, ...)  UseMethod("normalizeBetweenRuns")

normalizeBetweenRuns.default <- function(data, spikedin_peptide_sequences=c("LGGNEQVTR","GAGSSEPVTGLDAK","VEATFGVDESNAK","YILAGVENSK","TPVISGGPYEYR","TPVITGAPYEYR","DGLDAASYYAPVR","ADVTPADFSEWSK","GTFIIDPGGVIR","GTFIIDPAAVIR","LFLQFGAQGSPFLK"), span=0.7, iterations = 3, method = "fast", ...) {
    if (!inherits(data, "data.frame")) stop("Is not a data.frame")

	original_names<-names(data)

	if ("transition_intensity" %in% names(data)) {
		if (dim(unique(data[,c("run_id","protein_id","peptide_id","transition_id")]))[1] != dim(data[,c("run_id","protein_id","peptide_id","transition_id")])[1]) {
			stop("normalizeBetweenRuns: run_id / protein_id / peptide_id / transition_id are not unique among datasets!")
		}

		data$unit_id<-paste(data$protein_id,data$peptide_id,data$transition_id,sep="_")

		temp<-data[,c("run_id","unit_id","transition_intensity")]
		temp$transition_intensity<-log(temp$transition_intensity)
		ttemp<-acast(temp, unit_id ~ run_id, value.var="transition_intensity", fill=NA)

		peptide_weights<-data.frame("unit_id"=row.names(ttemp))
		peptide_weights<-merge(peptide_weights,unique(data[,c('unit_id','peptide_sequence')]))
		peptide_weights$weight<-0
		peptide_weights[which(peptide_weights$peptide_sequence %in% spikedin_peptide_sequences),"weight"]<-1

		nttemp<-normalizeCyclicLoess(ttemp, weights = peptide_weights$weights, span = span, iterations = iterations, method = method)
		ntemp<-melt(nttemp)[,c("Var1","Var2","value")]
		names(ntemp)<-c("unit_id","run_id","transition_intensity")
		ndata<-merge(data[,! names(data) %in% c("transition_intensity")],subset(ntemp,!is.na(transition_intensity)),by=c("run_id","unit_id"))[,original_names]
		ndata$transition_intensity<-exp(ndata$transition_intensity)
	}
	else if ("peptide_intensity" %in% names(data)) {
		if (dim(unique(data[,c("run_id","protein_id","peptide_id")]))[1] != dim(data[,c("run_id","protein_id","peptide_id")])[1]) {
			stop("normalizeBetweenRuns: run_id / protein_id / peptide_id are not unique among datasets!")
		}
		data$unit_id<-paste(data$protein_id,data$peptide_id,sep="_")

		temp<-data[,c("run_id","unit_id","peptide_intensity")]
		temp$peptide_intensity<-log(temp$peptide_intensity)
		ttemp<-acast(temp, unit_id ~ run_id, value.var="peptide_intensity", fill=NA)

		peptide_weights<-data.frame("unit_id"=row.names(ttemp))
		peptide_weights<-merge(peptide_weights,unique(data[,c('unit_id','peptide_sequence')]))
		peptide_weights$weight<-0
		peptide_weights[which(peptide_weights$peptide_sequence %in% spikedin_peptide_sequences),"weight"]<-1

		nttemp<-normalizeCyclicLoess(ttemp, weights = peptide_weights$weights, span = span, iterations = iterations, method = method)
		ntemp<-melt(nttemp)[,c("Var1","Var2","value")]
		names(ntemp)<-c("unit_id","run_id","peptide_intensity")
		ndata<-merge(data[,! names(data) %in% c("peptide_intensity")],subset(ntemp,!is.na(peptide_intensity)),by=c("run_id","unit_id"))[,original_names]
		ndata$peptide_intensity<-exp(ndata$peptide_intensity)
	}
	else {
		stop("Is not a generic aLFQ input data.frame")
	}

	return(ndata)
}
