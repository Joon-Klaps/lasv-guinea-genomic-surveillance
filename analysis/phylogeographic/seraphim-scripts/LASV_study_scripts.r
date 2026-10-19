library(diagram)
library(lubridate)
library(rgdal)
library(seraphim)

study_areas = c("MRU")
segments = c("L","S")
clades = c("MRU","NIG1","NIG2","NIG3")
e_NGA = extent(2, 15.5, 3.3, 14.4)
e_MRU = extent(-15, -4, 4, 13.04)
nberOfExtractionFiles = 1000

# 1. Extraction of the spatio-temporal information contained in posterior trees

options(scipen=9); burnIns = c(100000000, 314000000)
mostRecentSamplingDatum = decimal_date(ymd("2021-09-18"))

	# 1.1. Selecting 1,000 trees from the post-burn-in distribution

for (i in 1:length(segments))
	{
		allTrees = scan(file=paste0("RRW_phylogeography/LASV_",segments[i],"_alignment_4.trees"), what="", sep="\n", quiet=T, blank.lines.skip=F)
		lineIndex1 = which(grepl(paste0("tree STATE_0 "),allTrees))
		lineIndex2 = which(grepl(paste0("tree STATE_",burnIns[i]," "),allTrees))
		lineIndex3 = which(grepl(paste0("End;"),allTrees))[2]
		samplingFrequency = floor(((lineIndex3-1)-(lineIndex2))/nberOfExtractionFiles)
		treesToSample = seq(((lineIndex3-1)-((nberOfExtractionFiles-1)*samplingFrequency)),(lineIndex3-1),samplingFrequency)
		selectedTrees = c(allTrees[1:(lineIndex1-1)],allTrees[treesToSample],"End;")
		write(selectedTrees, paste0("LASV_seg",segments[i],"_",nberOfExtractionFiles,".trees"))
		log = scan(file=paste0("RRW_phylogeography/LASV_",segments[i],"_alignment_4.log"), what="", sep="\n", quiet=T, blank.lines.skip=F)
		logStates = rep(NA, length(log))
		for (j in 1:length(logStates))
			{
				logStates[j] = unlist(strsplit(log[j],"\t"))[1]
			}
		selectedTrees = scan(file=paste0("LASV_seg",segments[i],"_",nberOfExtractionFiles,".trees"), what="", sep="\n", quiet=T, blank.lines.skip=F)
		lines = which(grepl("tree STATE_",selectedTrees)); selectedTreeIDs = rep(NA, length(lines))
		for (j in 1:length(selectedTreeIDs))
			{
				selectedTreeIDs[j] = gsub("STATE_","",unlist(strsplit(selectedTrees[lines[j]]," "))[2])
			}
		selectedLogStates = logStates[which(logStates%in%selectedTreeIDs)]
		indices = c(1:5,which(logStates%in%selectedLogStates))
		write(log[indices], paste0("LASV_seg",segments[i],"_",nberOfExtractionFiles,".log"))
	}

	# 1.2. Generating the maximum clade credibility (MCC) trees with TreeAnnotator

for (i in 1:length(segments))
	{
		input = paste0("LASV_seg",segments[i],"_",nberOfExtractionFiles,".trees")
		output = paste0("LASV_seg",segments[i],"_",nberOfExtractionFiles,".tree")
		system(paste0("BEAST_version_1104/bin/treeannotator -burninTrees 0 -heights keep ",input," ",output))
	}

	# 1.3. Extractions of the spatio-temporal information contained in posterior trees

for (i in 1:length(segments))
	{
		localTreesDirectory = paste0("Tree_extractions_L&S/LASV_seg",segments[i],"_",nberOfExtractionFiles)
		allTrees = scan(file=paste0("LASV_seg",segments[i],"_",nberOfExtractionFiles,".trees"), what="", sep="\n", quiet=T, blank.lines.skip=F)
		burnIn = 0; randomSampling = FALSE
		nberOfTreesToSample = nberOfExtractionFiles
		coordinateAttributeName = "location"; nberOfCores = 5
		treeExtractions(localTreesDirectory, allTrees, burnIn, randomSampling, nberOfTreesToSample, mostRecentSamplingDatum, coordinateAttributeName, nberOfCores)
	}

	# 1.4. Extractions of the spatio-temporal information contained in MCC trees

source("mccTreeExtractions.r")
for (i in 1:length(segments))
	{
		mcc_tre = readAnnotatedNexus(paste0("LASV_seg",segments[i],"_",nberOfExtractionFiles,".tree"))
		mcc_tab = mccTreeExtractions(mcc_tre, mostRecentSamplingDatum)
		write.csv(mcc_tab, paste0("LASV_seg",segments[i],"_",nberOfExtractionFiles,".csv"), row.names=F, quote=F)
	}

	# 1.5. Dividing extraction files per clade

source("isolatingTargetClade.r")
for (i in 1:length(segments))
	{
		localTreesDirectory1 = paste0("Tree_extractions_L&S/LASV_seg",segments[i],"_",nberOfExtractionFiles)
		localTreesDirectory2 = paste0("Tree_extractions_L&S/LASV_mru",segments[i],"_",nberOfExtractionFiles)
		dir.create(file.path(localTreesDirectory2), showWarnings=F)
		metadata = read.csv(paste0("RRW_phylogeography/LASV_",segments[i],"_alignment_4.csv"))
		mcc_tre = readAnnotatedNexus(paste0("LASV_seg",segments[i],"_",nberOfExtractionFiles,".tree"))
		tipNodes = mcc_tre$edge[which(!mcc_tre$edge[,2]%in%mcc_tre$edge[,1]),2]
		tipClades = rep(NA, length(tipNodes)); samplingCoordinates = rep(NA, length(tipNodes))
		mcc_tre$tip.label = gsub(" ","",gsub("'","",mcc_tre$tip.label))
		for (j in 1:length(tipNodes))
			{
				tipClades[j] = metadata[which(gsub(" ","",metadata[,1])==mcc_tre$tip.label[tipNodes[j]]),"CO"]
				lon = round(mcc_tre$annotations[[which(mcc_tre$edge[,2]==tipNodes[j])]]$location2,7)
				lat = round(mcc_tre$annotations[[which(mcc_tre$edge[,2]==tipNodes[j])]]$location1,7)
				samplingCoordinates[j] = paste0(lon,"_",lat)
			}
		selectedTips = metadata[which(metadata[,"CO"]!="NGA"),1]
		metadata = cbind(mcc_tre$tip.label, samplingCoordinates, tipClades)
		for (j in 0:nberOfExtractionFiles)
			{
				if (j == 0)
					{
						tab1 = read.csv(paste0("LASV_seg",segments[i],"_",nberOfExtractionFiles,".csv"), header=T)
					}	else	{
						tab1 = read.csv(paste0(localTreesDirectory1,"/TreeExtractions_",j,".csv"), header=T)
					}
				samplingCoordinates = rep(NA, dim(tab1)[1]); node2s = rep(NA, dim(tab1)[1])
				for (k in 1:dim(tab1)[1])
					{
						lon = round(tab1[k,"endLon"],7); lat = round(tab1[k,"endLat"],7)
						samplingCoordinates[k] = paste0(lon,"_",lat); node2s[k] = tab1[k,"node2"]
					}
				indices = which(!tab1[,"node2"]%in%tab1[,"node1"])
				samplingCoordinates = samplingCoordinates[indices]; node2s = node2s[indices]
				indices = which(metadata[,"tipClades"]!="NGA"); tipNodes = rep(NA, length(indices))
				for (k in 1:length(indices))
					{
						index = which(samplingCoordinates == metadata[indices[k],"samplingCoordinates"])
						if (length(index) != 1) cat(i,j,k,"\n")
						tipNodes[k] = node2s[index]
					}
				tab2 = isolatingTargetClade(tab1, tipNodes)
				if (j == 0)
					{
						write.csv(tab2, paste0("LASV_mru",segments[i],"_",nberOfExtractionFiles,".csv"), row.names=F, quote=F)
					}	else	{
						write.csv(tab2, paste0(localTreesDirectory2,"/TreeExtractions_",j,".csv"), row.names=F, quote=F)
					}
			}
	}

	# 1.6. Annotating new tip nodes (newly sequenced genomes)

newSamples = c("G0405","G0274","M00362","M00363","M00364","M00008","M00009","M00539","M00541","M00542","G0683","G0671")
for (i in 1:length(segments))
	{
		metadata = read.csv(paste0("RRW_phylogeography/LASV_",segments[i],"_alignment_4.csv"), head=T)
		newCoordinates = rep(NA, length(newSamples))
		for (j in 1:length(newSamples))
			{
				indices = c()
				for (k in 1:dim(metadata)[1])
					{
						if (grepl(gsub("\\|","_",newSamples[j]),gsub("\\|","_",metadata[k,1]))) indices = c(indices, k)
					}
				if (length(indices) != 1)
					{
						cat(i,j,length(indices),"\n")
					}	else	{
						lon = round(metadata[indices,"longitude"], 3)
						lat = round(metadata[indices,"latitude"], 3)
						newCoordinates[j] = paste0(lon,"_",lat)
					}
			}
		tab1 = read.csv(paste0("LASV_mru",segments[i],"_",nberOfExtractionFiles,".csv"), header=T)
		samplingCoordinates = rep(NA, dim(tab1)[1]); node2s = rep(NA, dim(tab1)[1])
		for (j in 1:dim(tab1)[1])
			{
				lon = round(tab1[j,"endLon"],3); lat = round(tab1[j,"endLat"],3)
				samplingCoordinates[j] = paste0(lon,"_",lat); node2s[j] = tab1[j,"node2"]
			}
		newColumn = matrix(0, nrow=dim(tab1)[1], ncol=1); colnames(newColumn) = "newSample"
		for (j in 1:length(newSamples))
			{
				index = which(samplingCoordinates == newCoordinates[j])
				if (length(index) != 1) cat(i,j,"\n")
				newColumn[index,1] = 1
			}
		if ("newSample"%in%colnames(tab1))
			{
				tab2 = tab1; tab2[,"newSample"] = newColumn
			}	else		{
				tab2 = cbind(tab1, newColumn)
			}
		write.csv(tab2, paste0("LASV_mru",segments[i],"_",nberOfExtractionFiles,"_TO_CHECK.csv"), row.names=F, quote=F)
	}

# 2. Plotting the dispersal history of LASV lineages in MRU

percentage = 80; precision = 20
for (i in 1:length(segments))
	{
		localTreesDirectory = paste0("Tree_extractions_L&S/LASV_mru",segments[i],"_",nberOfExtractionFiles)
		mcc = read.csv(paste0("LASV_mru",segments[i],"_",nberOfExtractionFiles,".csv"))
		startDatum = min(mcc[,"startYear"]); prob = percentage/100
		polygons = suppressWarnings(spreadGraphic2(localTreesDirectory, nberOfExtractionFiles, prob, startDatum, precision))
		for (j in 1:length(polygons))
			{
				year = as.numeric(names(polygons[[j]]))
				layerName = paste0("LASV_mru",segments[i],"_",percentage,"_",precision,"_",year)
				writeOGR(polygons[[j]], dsn=paste0("LASV_80HPD_20yrs"), layer=layerName, driver="ESRI Shapefile")
			}
	}

background = raster("Environmental_files/LASV_rasters/Elevation_MRU_0.008.asc")
background = crop(background, e_MRU); background[!is.na(background[])] = 1; r = background
background_cols = colorRampPalette(c("grey","white"),bias=1)(max(r[],na.rm=T)-min(r[],na.rm=T))[1:(max(r[],na.rm=T)-min(r[],na.rm=T))]
borders = crop(readOGR(dsn="Environmental_files/Natural_Earth/", layer="International_borders"), e_MRU)
lakes = crop(readOGR(dsn="Environmental_files/Natural_Earth/", layer="Natural_Earth_lakes"), e_MRU)
rivers = crop(readOGR(dsn="Environmental_files/LASV_rivers/", layer="FAO_GeoNetwork"), e_MRU)
rivers@data[which(rivers@data[,"Strahler"]<3),"Strahler"] = 0
riverCols = rep("#AFBED5",dim(rivers@data)[1]); riverCols[which(rivers@data[,"Strahler"]<3)] = NA

pdf(paste0("Dispersal_history_NEW1.pdf"), width=12, height=10.2)
par(mfrow=c(2,2), mar=c(0,0,0,0), oma=c(1.5,2.7,1.0,0), mgp=c(0,0.4,0), lwd=0.2, bty="o")
titles = c("Segment L","Segment S"); plottingStrahlerLegend = FALSE
cols = colorRampPalette(brewer.pal(11,'RdYlGn'))(131)[21:121]
cols = gsub("FF","",rev(viridis::viridis(161)[31:161]))
for (i in 1:length(segments))
	{
		mcc = read.csv(paste0("LASV_mru",segments[i],"_",nberOfExtractionFiles,".csv"))
		mcc = mcc[order(mcc[,"startYear"]),]; mcc1 = mcc[1,]; mcc2 = mcc[c(2:dim(mcc)[1]),]
		mcc2 = mcc2[order(mcc2[,"endYear"]),]; mcc = rbind(mcc1,mcc2)
		mcc = mcc[order(1-mcc[,"newSample"]),]
		minYear = min(mcc[,"startYear"]); maxYear = max(mcc[,"endYear"])
		ancestralC = (((mcc[1,"startYear"]-minYear)/(maxYear-minYear))*100)+1
		startYearsC = (((mcc[,"startYear"]-minYear)/(maxYear-minYear))*100)+1
		endYearsC = (((mcc[,"endYear"]-minYear)/(maxYear-minYear))*100)+1
		col_ancestral = cols[ancestralC]; cols_endYears = cols[endYearsC]
		files = list.files("LASV_80HPD_20yrs"); pols = list(); cols_pol = list()
		files = files[which((grepl(paste0("_mru",segments[i]),files))&(grepl(".shp",files)))]
		for (k in 1:length(files))
			{
				pol = shapefile(paste0("LASV_80HPD_20yrs/",files[k]))
				pol@proj4string = CRS("+init=epsg:4326"); pols[[k]] = pol
				date = gsub(".shp","",unlist(strsplit(files[k],"_"))[length(unlist(strsplit(files[k],"_")))])
				precision = 20; date = as.numeric(date)+1-(precision/2)
				yearM = (((date-minYear)/(maxYear-minYear))*100)+1
				cols_pol[[k]] = cols[yearM]
			}
		for (k in 1:length(cols_pol)) cols_pol[[k]] = paste0(cols_pol[[k]],"70")
		plot(background, main="", cex.main=1, cex.axis=0.7, bty="n", box=F, axes=F, legend=F, axis.args=list(cex.axis=0.7), col="white", colNA="gray90")
		for (k in 1:length(pols)) plot(pols[[k]], axes=F, col=cols_pol[[k]], add=T, border=NA)
		plot(lakes, add=T, lwd=0.7, col=unique(riverCols)[!is.na(unique(riverCols))], border=NA)
		plot(rivers, add=T, lwd=rivers@data[,"Strahler"]/3.5, col=riverCols)
		mtext(titles[i], side=3, line=-27.3, at=3.7, cex=0.9, font=1, col="gray30")
		plot(borders, add=T, lwd=2, col="gray60", lty=2); counter = 0
		rect(xmin(background), ymin(background), xmax(background), ymax(background), xpd=T, lwd=0.2)
		axis(1,c(ceiling(xmin(background)), floor(xmax(background))), pos=ymin(background), mgp=c(0,0.7,0),
			 cex.axis=0.9, lwd=0, lwd.tick=0.2, padj=-0.8, tck=-0.01, col.axis="gray30")
		axis(2, c(ceiling(ymin(background)), floor(ymax(background))), pos=xmin(background), mgp=c(0,0.8,0),
			 cex.axis=0.9, lwd=0, lwd.tick=0.2, padj=1, tck=-0.01, col.axis="gray30")
		for (k in 1:dim(mcc)[1])
			{
				curvedarrow(cbind(mcc[k,"startLon"],mcc[k,"startLat"]), cbind(mcc[k,"endLon"],mcc[k,"endLat"]), arr.length=0,
				    		arr.width=0, lwd=0.3, lty=1, lcol="gray30", arr.col=NA, arr.pos=FALSE, curve=0.1, dr=NA, endhead=F)
			}
		for (k in dim(mcc)[1]:1)
			{
				if (mcc[k,"newSample"] == 0)
					{
						points(mcc[k,"endLon"], mcc[k,"endLat"], pch=16, col=cols_endYears[k], cex=1.0)
						points(mcc[k,"endLon"], mcc[k,"endLat"], pch=1, col="gray30", cex=1.0, lwd=0.2)
					}	else		{
						points(mcc[k,"endLon"], mcc[k,"endLat"], pch=16, col=rgb(222,67,39,255,maxColorValue=255), cex=1.3)
						points(mcc[k,"endLon"], mcc[k,"endLat"], pch=1, col="gray30", cex=1.3, lwd=0.2)
					}
				if ((mcc[k,"startYear"] == min(mcc[,"startYear"]))&(counter == 0))
					{
						counter = counter+1
						points(mcc[k,"startLon"], mcc[k,"startLat"], pch=16, col=col_ancestral, cex=1.0)
						points(mcc[k,"startLon"], mcc[k,"startLat"], pch=1, col="gray30", cex=1.0, lwd=0.2)
					}
			}
		legendRast = raster(as.matrix(c(minYear,maxYear)))
		plot(legendRast, legend.only=T, add=T, col=cols, legend.width=0.5, legend.shrink=0.3, smallplot=c(0.040,0.055,0.09,0.40), alpha=1.0,
			 legend.args=list(text="", cex=0.9, line=0.5, col="gray30"), axis.args=list(cex.axis=1.0, lwd=0, lwd.tick=0.2, tck=-0.8,
			 col.axis="gray30", line=0, mgp=c(0,0.6,0)))
		if (plottingStrahlerLegend == TRUE)
			{
				legend(x=-14.5, y=6.3, legend=c("S = 3 ","S = 4 ","S = 5 ","S = 6 ","S = 7 "), lwd=c(3:7)/3.5, cex=0.9, lty=1,
					   col=unique(riverCols)[!is.na(unique(riverCols))], x.intersp=0.5, border=NA, bg="white", box.lwd=0.2, box.col="gray30")
			}
	}
dev.off()

i = 1; pdf(paste0("Dispersal_history_NEW",i+1,".pdf"), width=12, height=10.2) # dev.new(width=12, height=10.2)
par(mfrow=c(2,2), mar=c(0,0,0,0), oma=c(1.5,2.7,1.0,0), mgp=c(0,0.4,0), lwd=0.2, bty="o")
timePeriods = c(1850, 1900, 1950, 2022)
titles = c("Segment L","Segment S"); plottingStrahlerLegend = FALSE
cols = colorRampPalette(brewer.pal(11,'RdYlGn'))(131)[21:121]
cols = gsub("FF","",rev(viridis::viridis(161)[31:161]))
mcc = read.csv(paste0("LASV_mru",segments[i],"_",nberOfExtractionFiles,".csv"))
mcc = mcc[order(mcc[,"startYear"]),]; mcc1 = mcc[1,]; mcc2 = mcc[c(2:dim(mcc)[1]),]
mcc2 = mcc2[order(mcc2[,"endYear"]),]; mcc = rbind(mcc1,mcc2)
mcc = mcc[order(1-mcc[,"newSample"]),]
minYear = min(mcc[,"startYear"]); maxYear = max(mcc[,"endYear"])
ancestralC = (((mcc[1,"startYear"]-minYear)/(maxYear-minYear))*100)+1
startYearsC = (((mcc[,"startYear"]-minYear)/(maxYear-minYear))*100)+1
endYearsC = (((mcc[,"endYear"]-minYear)/(maxYear-minYear))*100)+1
col_ancestral = cols[ancestralC]; cols_endYears = cols[endYearsC]
files = list.files("LASV_80HPD_20yrs"); pols = list(); cols_pol = list()
files = files[which((grepl(paste0("_mru",segments[i]),files))&(grepl(".shp",files)))]
for (k in 1:length(files))
	{
		pol = shapefile(paste0("LASV_80HPD_20yrs/",files[k]))
		pol@proj4string = CRS("+init=epsg:4326"); pols[[k]] = pol
		date = gsub(".shp","",unlist(strsplit(files[k],"_"))[length(unlist(strsplit(files[k],"_")))])
		precision = 20; date = as.numeric(date)+1-(precision/2)
		yearM = (((date-minYear)/(maxYear-minYear))*100)+1
		cols_pol[[k]] = cols[yearM]
	}
for (k in 1:length(cols_pol)) cols_pol[[k]] = paste0(cols_pol[[k]],"70")
for (h in 1:length(timePeriods))
	{
		plot(background, main="", cex.main=1, cex.axis=0.7, bty="n", box=F, axes=F, legend=F, axis.args=list(cex.axis=0.7), col="white", colNA="gray90")
		for (k in 1:length(pols))
			{
				if (as.numeric(gsub("_",".",gsub("X","",names(pols[[k]])))) <= timePeriods[h])
					{
						plot(pols[[k]], axes=F, col=cols_pol[[k]], add=T, border=NA)
					}
			}
		plot(lakes, add=T, lwd=0.7, col=unique(riverCols)[!is.na(unique(riverCols))], border=NA)
		plot(rivers, add=T, lwd=rivers@data[,"Strahler"]/3.5, col=riverCols)
		plot(borders, add=T, lwd=2, col="gray60", lty=2); counter = 0
		rect(xmin(background), ymin(background), xmax(background), ymax(background), xpd=T, lwd=0.2)
		axis(1,c(ceiling(xmin(background)), floor(xmax(background))), pos=ymin(background), mgp=c(0,0.7,0),
			 cex.axis=0.9, lwd=0, lwd.tick=0.2, padj=-0.8, tck=-0.01, col.axis="gray30")
		axis(2, c(ceiling(ymin(background)), floor(ymax(background))), pos=xmin(background), mgp=c(0,0.8,0),
			 cex.axis=0.9, lwd=0, lwd.tick=0.2, padj=1, tck=-0.01, col.axis="gray30")
		selectedBranches = which(mcc[,"endYear"] <= timePeriods[h])
		for (k in 1:dim(mcc)[1])
			{
				if (k %in% selectedBranches)
					{
						curvedarrow(cbind(mcc[k,"startLon"],mcc[k,"startLat"]), cbind(mcc[k,"endLon"],mcc[k,"endLat"]), arr.length=0,
				    				arr.width=0, lwd=0.3, lty=1, lcol="gray30", arr.col=NA, arr.pos=FALSE, curve=0.1, dr=NA, endhead=F)
					}
			}
		for (k in dim(mcc)[1]:1)
			{
				if (k %in% selectedBranches)
					{
						if (mcc[k,"newSample"] == 0)
							{
								points(mcc[k,"endLon"], mcc[k,"endLat"], pch=16, col=cols_endYears[k], cex=1.0)
								points(mcc[k,"endLon"], mcc[k,"endLat"], pch=1, col="gray30", cex=1.0, lwd=0.2)
							}	else		{
								points(mcc[k,"endLon"], mcc[k,"endLat"], pch=16, col=rgb(222,67,39,255,maxColorValue=255), cex=1.3)
								points(mcc[k,"endLon"], mcc[k,"endLat"], pch=1, col="gray30", cex=1.3, lwd=0.2)
							}
						if ((mcc[k,"startYear"] == min(mcc[,"startYear"]))&(counter == 0))
							{
								points(mcc[k,"startLon"], mcc[k,"startLat"], pch=16, col=col_ancestral, cex=1.0)
								points(mcc[k,"startLon"], mcc[k,"startLat"], pch=1, col="gray30", cex=1.0, lwd=0.2)
							}
					}
			}
		if (h == length(timePeriods))
			{
				legendRast = raster(as.matrix(c(minYear,maxYear)))
				plot(legendRast, legend.only=T, add=T, col=cols, legend.width=0.5, legend.shrink=0.3, smallplot=c(0.040,0.055,0.09,0.40), alpha=1.0,
					 legend.args=list(text="", cex=0.9, line=0.5, col="gray30"), axis.args=list(cex.axis=1.0, lwd=0, lwd.tick=0.2, tck=-0.8,
					 col.axis="gray30", line=0, mgp=c(0,0.6,0)))
			}	else	{
				mtext(timePeriods[h], side=3, line=-27.3, at=-14.2, cex=0.9, font=1, col="gray30")
			}
		if (plottingStrahlerLegend == TRUE)
			{
				legend(x=-14.5, y=6.3, legend=c("S = 3 ","S = 4 ","S = 5 ","S = 6 ","S = 7 "), lwd=c(3:7)/3.5, cex=0.9, lty=1,
					   col=unique(riverCols)[!is.na(unique(riverCols))], x.intersp=0.5, border=NA, bg="white", box.lwd=0.2, box.col="gray30")
			}
	}
dev.off()

