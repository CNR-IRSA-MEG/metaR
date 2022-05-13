################
# ESTER Genome #
################

library(rvest)

urlLink <- read.csv("./Desktop/Dati_Ester/Rotifer_data.csv")

metadat <- data.frame(BioSample = NA,
                      SRA = NA,
                      sampleName = NA,
                      collectionDate = NA,
                      broadScaleEnvironmentalContext = NA,
                      localScaleEnvironmentalContext = NA,
                      environmentalMedium = NA,
                      geographicLocation = NA,
                      investigationType = NA,
                      organismCount = NA,
                      projectName = NA,
                      sampleCollectionDeviceOrMethod = NA,
                      sampleMaterialProcessing = NA,
                      host = NA,
                      ENAFirstPublic = NA,
                      ENALastUpdate = NA,
                      ENACHECKLIST = NA,
                      ExternalId = NA,
                      INSDCCenterAlias = NA,
                      INSDCCenterName = NA,
                      INSDCFirstPublic = NA,
                      INSDCLastUpdate = NA,
                      INSDCStatus = NA,
                      SRAAccession = NA,
                      SubmitterId = NA,
                      amountOrSizeSampleCollected = NA ,
                      depth = NA,
                      latitude = NA,
                      longitude = NA,
                      miscellaneousEnvironmentalPackage = NA,
                      pcrPrimers = NA,
                      sequencingMethod = NA,
                      targetGene = NA,
                      targetSubfragment = NA,
                      title = NA)

for(i in 1:nrow(urlLink)){
  metadat.1 <- data.frame(BioSample = NA,
                        SRA = NA,
                        sampleName = NA,
                        collectionDate = NA,
                        broadScaleEnvironmentalContext = NA,
                        localScaleEnvironmentalContext = NA,
                        environmentalMedium = NA,
                        geographicLocation = NA,
                        investigationType = NA,
                        organismCount = NA,
                        projectName = NA,
                        sampleCollectionDeviceOrMethod = NA,
                        sampleMaterialProcessing = NA,
                        host = NA,
                        ENAFirstPublic = NA,
                        ENALastUpdate = NA,
                        ENACHECKLIST = NA,
                        ExternalId = NA,
                        INSDCCenterAlias = NA,
                        INSDCCenterName = NA,
                        INSDCFirstPublic = NA,
                        INSDCLastUpdate = NA,
                        INSDCStatus = NA,
                        SRAAccession = NA,
                        SubmitterId = NA,
                        amountOrSizeSampleCollected = NA ,
                        depth = NA,
                        latitude = NA,
                        longitude = NA,
                        miscellaneousEnvironmentalPackage = NA,
                        pcrPrimers = NA,
                        sequencingMethod = NA,
                        targetGene = NA,
                        targetSubfragment = NA,
                        title = NA)  

metadat.1[,1:2] <- c(urlLink$BioSample[i], urlLink$SRA[i])
  
url <- urlLink$URL[i] # <- URL
tmp <- read_html(url) # Open URL
tmp <- html_nodes(tmp, "body") # GO to the body
# tmp.dl <- html_nodes(tmp, "dl")[1] # Go to the node "dl" to keep the sample name
# tmp.dl.df <- data.frame(row.names = c("BioSample", "SRA"),
#                         values = c(urlLink$BioSample[i], urlLink$SRA[i])
#                         #values = html_text(html_nodes(tmp.dl, "dd"))
#                         )

tmp.th <- html_nodes(tmp, "th") # Name of attibutes 
tmp.dt <- html_nodes(tmp, "td") # Values of attibutes
attrib <- data.frame(row.names = html_text(tmp.th),
                     values = html_text(tmp.dt))

metadat.1$sampleName <- ifelse("sample name" %in% rownames(attrib), attrib$values[which(rownames(attrib) == "sample name")], NA)
metadat.1$collectionDate <- ifelse("collection date" %in% rownames(attrib), attrib$values[which(rownames(attrib) == "collection date")], NA)
metadat.1$broadScaleEnvironmentalContext <- ifelse("broad-scale environmental context" %in% rownames(attrib), attrib$values[which(rownames(attrib) == "broad-scale environmental context")], NA)
metadat.1$localScaleEnvironmentalContext <- ifelse("local-scale environmental context" %in% rownames(attrib), attrib$values[which(rownames(attrib) == "local-scale environmental context")], NA)
metadat.1$environmentalMedium <- ifelse("environmental medium" %in% rownames(attrib), attrib$values[which(rownames(attrib) == "environmental medium")], NA)
metadat.1$geographicLocation <- ifelse("geographic location" %in% rownames(attrib), attrib$values[which(rownames(attrib) == "geographic location")], NA)
metadat.1$investigationType <- ifelse("investigation type" %in% rownames(attrib), attrib$values[which(rownames(attrib) == "investigation type")], NA)
metadat.1$organismCount <- ifelse("organism count" %in% rownames(attrib), attrib$values[which(rownames(attrib) == "organism count")], NA)
metadat.1$projectName <- ifelse("project name" %in% rownames(attrib), attrib$values[which(rownames(attrib) == "project name")], NA)
metadat.1$sampleCollectionDeviceOrMethod <- ifelse("sample collection device or method" %in% rownames(attrib), attrib$values[which(rownames(attrib) == "sample collection device or method")], NA)
metadat.1$sampleMaterialProcessing <- ifelse("sample material processing" %in% rownames(attrib), attrib$values[which(rownames(attrib) == "sample material processing")], NA)
metadat.1$host <- ifelse("host" %in% rownames(attrib), attrib$values[which(rownames(attrib) == "host")], NA)
metadat.1$ENAFirstPublic <- ifelse("ENA first public" %in% rownames(attrib), attrib$values[which(rownames(attrib) == "ENA first public")], NA)
metadat.1$ENALastUpdate <- ifelse("ENA last update" %in% rownames(attrib), attrib$values[which(rownames(attrib) == "ENA last update")], NA)
metadat.1$ENACHECKLIST <- ifelse("ENA-CHECKLIST" %in% rownames(attrib), attrib$values[which(rownames(attrib) == "ENA-CHECKLIST")], NA)
metadat.1$ExternalId <- ifelse("External Id" %in% rownames(attrib), attrib$values[which(rownames(attrib) == "External Id")], NA)
metadat.1$INSDCCenterAlias <- ifelse("INSDC center alias" %in% rownames(attrib), attrib$values[which(rownames(attrib) == "INSDC center alias")], NA)
metadat.1$INSDCCenterName <- ifelse("INSDC center name" %in% rownames(attrib), attrib$values[which(rownames(attrib) == "INSDC center name")], NA)
metadat.1$INSDCFirstPublic <- ifelse("INSDC first public" %in% rownames(attrib), attrib$values[which(rownames(attrib) == "INSDC first public")], NA)
metadat.1$INSDCLastUpdate <- ifelse("INSDC last update" %in% rownames(attrib), attrib$values[which(rownames(attrib) == "INSDC last update")], NA)
metadat.1$INSDCStatus <- ifelse("INSDC status" %in% rownames(attrib), attrib$values[which(rownames(attrib) == "INSDC status")], NA)
metadat.1$SRAAccession <- ifelse("SRA accession" %in% rownames(attrib), attrib$values[which(rownames(attrib) == "SRA accession")], NA)
metadat.1$SubmitterId <- ifelse("Submitter Id" %in% rownames(attrib), attrib$values[which(rownames(attrib) == "Submitter Id")], NA)
metadat.1$amountOrSizeSampleCollected <- ifelse("amount or size of sample collected" %in% rownames(attrib), attrib$values[which(rownames(attrib) == "amount or size of sample collected")], NA)
metadat.1$depth <- ifelse("geographic location (depth)" %in% rownames(attrib), attrib$values[which(rownames(attrib) == "geographic location (depth)")], NA)
metadat.1$latitude <- ifelse("geographic location (latitude)" %in% rownames(attrib), attrib$values[which(rownames(attrib) == "geographic location (latitude)")], NA)
metadat.1$longitude <- ifelse("geographic location (longitude)" %in% rownames(attrib), attrib$values[which(rownames(attrib) == "geographic location (longitude)")], NA)
metadat.1$miscellaneousEnvironmentalPackage <- ifelse("miscellaneous environmental package" %in% rownames(attrib), attrib$values[which(rownames(attrib) == "miscellaneous environmental package")], NA)
metadat.1$pcrPrimers <- ifelse("pcr primers" %in% rownames(attrib), attrib$values[which(rownames(attrib) == "pcr primers")], NA)
metadat.1$sequencingMethod <- ifelse("sequencing method" %in% rownames(attrib), attrib$values[which(rownames(attrib) == "sequencing method")], NA)
metadat.1$targetGene <- ifelse("target gene" %in% rownames(attrib), attrib$values[which(rownames(attrib) == "target gene")], NA)
metadat.1$targetSubfragment <- ifelse("target subfragment" %in% rownames(attrib), attrib$values[which(rownames(attrib) == "target subfragment")], NA)
metadat.1$title <- ifelse("title" %in% rownames(attrib), attrib$values[which(rownames(attrib) == "title")], NA)

metadat <- rbind(metadat, metadat.1)

}; rm(url, tmp, tmp.dt, tmp.th, attrib, metadat.1)

write.csv(metadat, "./Desktop/Dati_Ester/Metadata/metadata.csv")
