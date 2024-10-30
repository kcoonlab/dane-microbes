library(getPass)           
library(httr)              
library(jsonlite)          
library(geojsonio)          
library(geojsonR)         
library(rgdal)            
library(sp)           
library(raster)          
library(rasterVis)       
library(ggplot2)          
library(RColorBrewer)  

API_URL = 'https://lpdaacsvc.cr.usgs.gov/appeears/api/'

user <- getPass(msg = "Enter NASA Earthdata Login Username: ")    

password <- getPass(msg = "Enter NASA Earthdata Login Password: ")  

secret <- jsonlite::base64_enc(paste(user, password, sep = ":")) 

API_URL = 'https://lpdaacsvc.cr.usgs.gov/appeears/api/'  

response <- httr::POST(paste0(API_URL,"login"), add_headers("Authorization" = paste("Basic", gsub("\n", "", secret)),
                                                 "Content-Type" =
                                                   "application/x-www-form-urlencoded;charset=UTF-8"), 
                 body = "grant_type=client_credentials")

response_content <- content(response)                          # Retrieve the content of the request
token_response <- toJSON(response_content, auto_unbox = TRUE)  # Convert the response to the JSON object
remove(user, password, secret, response)                       # Remove the variables that are not needed anymore 
prettify(token_response)      

desired_prods <- c("MYD13Q1.006","MYD11A1.006","MYD17A2HGF.006","MYD17A3HGF.006") 
desired_layers <- c("_250m_16_days_EVI","_250m_16_days_NDVI","_250m_16_days_VI_Quality","LST_Day_1km","LST_Night_1km","QC_Day","QC_Night","Gpp_500m","Npp_500m")     
desired_prods2 <- c("MYD13Q1.006","MYD13Q1.006","MYD13Q1.006","MYD11A1.006","MYD11A1.006","MYD11A1.006","MYD11A1.006","MYD17A2HGF.006","MYD17A3HGF.006") 
layers <- data.frame(product = desired_prods2, layer = desired_layers)           
proj_req <- GET(paste0(API_URL, "spatial/proj"))            
proj_content <- content(proj_req)
proj_response <- toJSON(proj_content, auto_unbox = TRUE)
remove(proj_req, proj_content)
projs <- fromJSON(proj_response)
projection <- projs[projs$Name=="geographic",] 
taskType <- 'area'
projection <- projection$Name
outFormat <- 'geotiff'
recurring <- FALSE
out <- list(projection )
names(out) <- c("projection")
out$format$type <- outFormat
token <- paste("Bearer", fromJSON(token_response)$token) 
task_ids_collated_1 <- data.frame(matrix(ncol = 2, nrow = 0))
colnames(task_ids_collated_1) <- c('task_id', 'task_name')
task_ids_collated_2 <- data.frame(matrix(ncol = 2, nrow = 0))
colnames(task_ids_collated_2) <- c('task_id', 'task_name')
task_ids_collated_3 <- data.frame(matrix(ncol = 2, nrow = 0))
colnames(task_ids_collated_3) <- c('task_id', 'task_name')

for (i in c(4,5,6,13,26,81,83,86,88,158,164,189,218,243,245,247,249,251,253,254,256,319,335,340,369,371,391,397,417,465,513,515,522,523,526,530,552,556,558,559,565,567,583,584,587,590,597,598,607,615,618,619,620,676,700,715,718,719,724,725,900,969,972,973,975,978,1931,1932,2041,3262,3333,3921,3922,4022,4130,8300,8301,8512,9007,9009,9014,9019,9098,9601,9970,9971)) {
	nps <- readOGR(dsn= "mygeodata/dane-sites-master-polygon.shp", layer = "dane-sites-master-polygon")
	nps_sub <- subset(nps, AREANUMB == i) 
	remove(nps)
	nps_json <- geojsonio::geojson_json(nps_sub, geometry = "polygon")
	nps_js <- geojsonR::FROM_GeoJson(nps_json)
	remove(nps_json)
	taskName <- paste('Dane_',i,'_1',sep="")    
	taskType <- 'area'
	outFormat <- 'geotiff'
	startDate <- '05-01-2007'
	endDate <- '09-30-2011'
	date <- data.frame(startDate = startDate, endDate = endDate)
	nps_js$features[[1]]$geometry$coordinates <- list(nps_js$features[[1]]$geometry$coordinates)
	task_info <- list(date, layers, out, nps_js)
	names(task_info) <- c("dates", "layers", "output", "geo")
	task <- list(task_info, taskName, taskType)
	names(task) <- c("params", "task_name", "task_type") 
	task_json <- jsonlite::toJSON(task,auto_unbox = TRUE, digits = 10)
	response <- POST(paste0(API_URL, "task"), body = task_json , encode = "json", 
                 add_headers(Authorization = token, "Content-Type" = "application/json"))
	task_content <- content(response)
	task_response <- toJSON(task_content, auto_unbox = TRUE)
	prettify(task_response)
	task_id <- fromJSON(task_response)[[1]]
	task_ids_collated_1[nrow(task_ids_collated_1) + 1, 1] <- task_id
	task_ids_collated_1[nrow(task_ids_collated_1), 2] <- taskName
}

task_ids_collated_1_mod <- task_ids_collated_1[-c(9, 13, 27, 28, 43, 65, 82), ]
task_ids_collated_1_mod2 <- task_ids_collated_1_mod
task_ids_collated_1_mod2[nrow(task_ids_collated_1_mod2) + 1,] = c("56055e1f-a264-4161-b6d2-b390a4984c52", "Dane_88_1")
task_ids_collated_1_mod2[nrow(task_ids_collated_1_mod2) + 1,] = c("4619eab4-17db-40ec-93d2-ca4a772f3a39", "Dane_218_1")
task_ids_collated_1_mod2[nrow(task_ids_collated_1_mod2) + 1,] = c("c15f4c71-f3a3-4f6a-9b9d-c3dc0b8ead00", "Dane_391_1")
task_ids_collated_1_mod2[nrow(task_ids_collated_1_mod2) + 1,] = c("bf28f790-d6b0-416f-b964-b021786fde6a", "Dane_397_1")
task_ids_collated_1_mod2[nrow(task_ids_collated_1_mod2) + 1,] = c("50cc3b94-bb4a-4078-b43b-28d15db633f8", "Dane_583_1")
task_ids_collated_1_mod2[nrow(task_ids_collated_1_mod2) + 1,] = c("b0cd25b8-03c4-4c81-8067-6355e04ddf1d", "Dane_975_1")
task_ids_collated_1_mod2[nrow(task_ids_collated_1_mod2) + 1,] = c("b5139c20-a8d7-4fe3-a1a2-d304aa188840", "Dane_9019_1")
row.names(task_ids_collated_1_mod2) <- 1:nrow(task_ids_collated_1_mod2)
write.csv(task_ids_collated_1_mod2,"task_ids_collated_1.csv",row.names = FALSE)
task_ids_collated_1_mod3 <- task_ids_collated_1_mod2[-c(84, 86), ]
row.names(task_ids_collated_1_mod3) <- 1:nrow(task_ids_collated_1_mod3)

for(i in 1:nrow(task_ids_collated_1_mod3)){
	task_id <- task_ids_collated_1_mod3[i,1]
	task_name <- task_ids_collated_1_mod3[i,2]
	outDir <- file.path(paste('env-data/MODIS/MODIS_1/',task_name,sep=""))
	suppressWarnings(dir.create(outDir))  
	response <- GET(paste0(API_URL, "bundle/", task_id), add_headers(Authorization = token))
	bundle_response <- prettify(toJSON(content(response), auto_unbox = TRUE))
	bundle <- fromJSON(bundle_response)$files
	bundle <- bundle[which(bundle$file_type=="csv"),]
		for (j in bundle$file_id){
  			filename <- bundle[bundle$file_id == j,]$file_name           
  			filepath <- paste(outDir,filename, sep = "/")
  			suppressWarnings(dir.create(dirname(filepath)))
  			response <- GET(paste0(API_URL, "bundle/", task_id, "/", j), 
                  write_disk(filepath, overwrite = TRUE), progress(),
                  add_headers(Authorization = token))
		}
	}

find . -type f -name '*-Psn-QC-500m-lookup.csv' -exec rm {} +
find . -type f -name '*-Psn-QC-500m-Statistics-QA.csv' -exec rm {} +
find . -type f -name '*-Npp-QC-500m-Statistics-QA.csv' -exec rm {} +

for (i in c(4,5,6,13,26,81,83,86,88,158,164,189,218,243,245,247,249,251,253,254,256,319,335,340,369,371,391,397,417,465,513,515,522,523,526,530,552,556,558,559,565,567,583,584,587,590,597,598,607,615,618,619,620,676,700,715,718,719,724,725,900,969,972,973,975,978,1931,1932,2041,3262,3333,3921,3922,4022,4130,8300,8301,8512,9007,9009,9014,9019,9098,9601,9970,9971)) {
	nps <- readOGR(dsn= "mygeodata/dane-sites-master-polygon.shp", layer = "dane-sites-master-polygon")
	nps_sub <- subset(nps, AREANUMB == i) 
	remove(nps)
	nps_json <- geojsonio::geojson_json(nps_sub, geometry = "polygon")
	nps_js <- geojsonR::FROM_GeoJson(nps_json)
	remove(nps_json)
	taskName <- paste('Dane_',i,'_2',sep="")    
	taskType <- 'area'
	outFormat <- 'geotiff'
	startDate <- '05-01-2012'
	endDate <- '09-30-2016'
	date <- data.frame(startDate = startDate, endDate = endDate)
	nps_js$features[[1]]$geometry$coordinates <- list(nps_js$features[[1]]$geometry$coordinates)
	task_info <- list(date, layers, out, nps_js)
	names(task_info) <- c("dates", "layers", "output", "geo")
	task <- list(task_info, taskName, taskType)
	names(task) <- c("params", "task_name", "task_type") 
	task_json <- jsonlite::toJSON(task,auto_unbox = TRUE, digits = 10)
	response <- POST(paste0(API_URL, "task"), body = task_json , encode = "json", 
                 add_headers(Authorization = token, "Content-Type" = "application/json"))
	task_content <- content(response)
	task_response <- toJSON(task_content, auto_unbox = TRUE)
	prettify(task_response)
	task_id <- fromJSON(task_response)[[1]]
	task_ids_collated_2[nrow(task_ids_collated_2) + 1, 1] <- task_id
	task_ids_collated_2[nrow(task_ids_collated_2), 2] <- taskName
}

task_ids_collated_2_mod <- task_ids_collated_2[-c(9, 13, 27, 28, 43, 65, 82), ]
task_ids_collated_2_mod2 <- task_ids_collated_2_mod
task_ids_collated_2_mod2[nrow(task_ids_collated_2_mod2) + 1,] = c("fa057cfc-3766-42b4-8208-3d89157d8f1e", "Dane_88_2")
task_ids_collated_2_mod2[nrow(task_ids_collated_2_mod2) + 1,] = c("0aa33054-f854-4bc7-a8d8-d174afab11bf", "Dane_218_2")
task_ids_collated_2_mod2[nrow(task_ids_collated_2_mod2) + 1,] = c("25a7737e-a6e0-4c94-bef2-eca0c462628f", "Dane_391_2")
task_ids_collated_2_mod2[nrow(task_ids_collated_2_mod2) + 1,] = c("5941a3d2-882a-492f-8762-ba142a8c565c", "Dane_397_2")
task_ids_collated_2_mod2[nrow(task_ids_collated_2_mod2) + 1,] = c("c49a0b12-84a9-4395-8579-6fbf5bd14884", "Dane_583_2")
task_ids_collated_2_mod2[nrow(task_ids_collated_2_mod2) + 1,] = c("975f28ef-8179-4040-9695-46647d49de8b", "Dane_975_2")
task_ids_collated_2_mod2[nrow(task_ids_collated_2_mod2) + 1,] = c("2e97c96d-af7b-43ec-b2e9-3533aca776b9", "Dane_9019_2")
row.names(task_ids_collated_2_mod2) <- 1:nrow(task_ids_collated_2_mod2)
write.csv(task_ids_collated_2_mod2,"task_ids_collated_2.csv",row.names = FALSE)
task_ids_collated_2_mod3 <- task_ids_collated_2_mod2[-c(84, 86), ]
row.names(task_ids_collated_2_mod3) <- 1:nrow(task_ids_collated_2_mod3)

for(i in 1:nrow(task_ids_collated_2_mod3)){
	task_id <- task_ids_collated_2_mod3[i,1]
	task_name <- task_ids_collated_2_mod3[i,2]
	outDir <- file.path(paste('env-data/MODIS/MODIS_2/',task_name,sep=""))
	suppressWarnings(dir.create(outDir))  
	response <- GET(paste0(API_URL, "bundle/", task_id), add_headers(Authorization = token))
	bundle_response <- prettify(toJSON(content(response), auto_unbox = TRUE))
	bundle <- fromJSON(bundle_response)$files
	bundle <- bundle[which(bundle$file_type=="csv"),]
		for (j in bundle$file_id){
  			filename <- bundle[bundle$file_id == j,]$file_name           
  			filepath <- paste(outDir,filename, sep = "/")
  			suppressWarnings(dir.create(dirname(filepath)))
  			response <- GET(paste0(API_URL, "bundle/", task_id, "/", j), 
                  write_disk(filepath, overwrite = TRUE), progress(),
                  add_headers(Authorization = token))
		}
	}

find . -type f -name '*-Psn-QC-500m-lookup.csv' -exec rm {} +
find . -type f -name '*-Psn-QC-500m-Statistics-QA.csv' -exec rm {} +
find . -type f -name '*-Npp-QC-500m-Statistics-QA.csv' -exec rm {} +

for (i in c(4,5,6,13,26,81,83,86,88,158,164,189,218,243,245,247,249,251,253,254,256,319,335,340,369,371,391,397,417,465,513,515,522,523,526,530,552,556,558,559,565,567,583,584,587,590,597,598,607,615,618,619,620,676,700,715,718,719,724,725,900,969,972,973,975,978,1931,1932,2041,3262,3333,3921,3922,4022,4130,8300,8301,8512,9007,9009,9014,9019,9098,9601,9970,9971)) {
	nps <- readOGR(dsn= "mygeodata/dane-sites-master-polygon.shp", layer = "dane-sites-master-polygon")	
	nps_sub <- subset(nps, AREANUMB == i) 
	remove(nps)
	nps_json <- geojsonio::geojson_json(nps_sub, geometry = "polygon")
	nps_js <- geojsonR::FROM_GeoJson(nps_json)
	remove(nps_json)
	taskName <- paste('Dane_',i,'_3',sep="")    
	taskType <- 'area'
	outFormat <- 'geotiff'
	startDate <- '05-01-2017'
	endDate <- '09-30-2021'
	date <- data.frame(startDate = startDate, endDate = endDate)
	nps_js$features[[1]]$geometry$coordinates <- list(nps_js$features[[1]]$geometry$coordinates)
	task_info <- list(date, layers, out, nps_js)
	names(task_info) <- c("dates", "layers", "output", "geo")
	task <- list(task_info, taskName, taskType)
	names(task) <- c("params", "task_name", "task_type") 
	task_json <- jsonlite::toJSON(task,auto_unbox = TRUE, digits = 10)
	response <- POST(paste0(API_URL, "task"), body = task_json , encode = "json", 
                 add_headers(Authorization = token, "Content-Type" = "application/json"))
	task_content <- content(response)
	task_response <- toJSON(task_content, auto_unbox = TRUE)
	prettify(task_response)
	task_id <- fromJSON(task_response)[[1]]
	task_ids_collated_3[nrow(task_ids_collated_3) + 1, 1] <- task_id
	task_ids_collated_3[nrow(task_ids_collated_3), 2] <- taskName
}

task_ids_collated_3_mod <- task_ids_collated_3[-c(9, 13, 27, 28, 43, 65, 82), ]
task_ids_collated_3_mod2 <- task_ids_collated_3_mod
task_ids_collated_3_mod2[nrow(task_ids_collated_3_mod2) + 1,] = c("53436105-3aed-4ccb-9fa8-287dbd2cfbeb", "Dane_88_3")
task_ids_collated_3_mod2[nrow(task_ids_collated_3_mod2) + 1,] = c("3b72c9db-c548-4f80-b8c8-c7b7a72173ad", "Dane_218_3")
task_ids_collated_3_mod2[nrow(task_ids_collated_3_mod2) + 1,] = c("64e76f76-21bb-458d-a689-30d4bf4df145", "Dane_391_3")
task_ids_collated_3_mod2[nrow(task_ids_collated_3_mod2) + 1,] = c("5bf965ad-de91-49b3-94ac-2cff79051ea6", "Dane_397_3")
task_ids_collated_3_mod2[nrow(task_ids_collated_3_mod2) + 1,] = c("814747dd-f376-4b8e-8f8c-fd54f749840c", "Dane_583_3")
task_ids_collated_3_mod2[nrow(task_ids_collated_3_mod2) + 1,] = c("0a5dbd31-1412-4444-9b46-0ae80eebdd02", "Dane_975_3")
task_ids_collated_3_mod2[nrow(task_ids_collated_3_mod2) + 1,] = c("790e3b1c-a00f-47ff-a346-3402578da592", "Dane_9019_3")
row.names(task_ids_collated_3_mod2) <- 1:nrow(task_ids_collated_3_mod2)
write.csv(task_ids_collated_3_mod2,"task_ids_collated_3.csv",row.names = FALSE)
task_ids_collated_3_mod3 <- task_ids_collated_3_mod2[-c(84, 86), ]
row.names(task_ids_collated_3_mod3) <- 1:nrow(task_ids_collated_3_mod3)

for(i in 1:nrow(task_ids_collated_3_mod2)){
	task_id <- task_ids_collated_3_mod2[i,1]
	task_name <- task_ids_collated_3_mod2[i,2]
	outDir <- file.path(paste('env-data/MODIS/MODIS_3/',task_name,sep=""))
	suppressWarnings(dir.create(outDir))  
	response <- GET(paste0(API_URL, "bundle/", task_id), add_headers(Authorization = token))
	bundle_response <- prettify(toJSON(content(response), auto_unbox = TRUE))
	bundle <- fromJSON(bundle_response)$files
	bundle <- bundle[which(bundle$file_type=="csv"),]
		for (j in bundle$file_id){
  			filename <- bundle[bundle$file_id == j,]$file_name           
  			filepath <- paste(outDir,filename, sep = "/")
  			suppressWarnings(dir.create(dirname(filepath)))
  			response <- GET(paste0(API_URL, "bundle/", task_id, "/", j), 
                  write_disk(filepath, overwrite = TRUE), progress(),
                  add_headers(Authorization = token))
		}
	}

find . -type f -name '*-Psn-QC-500m-lookup.csv' -exec rm {} +
find . -type f -name '*-Psn-QC-500m-Statistics-QA.csv' -exec rm {} +
find . -type f -name '*-Npp-QC-500m-Statistics-QA.csv' -exec rm {} +

task_ids_collated_4 <- data.frame(matrix(ncol = 2, nrow = 0))
colnames(task_ids_collated_4) <- c('task_id', 'task_name')
task_ids_collated_4[nrow(task_ids_collated_4) + 1, 1] <- "33644fd7-a6ed-4f55-91c0-fed5249efdab"
task_ids_collated_4[nrow(task_ids_collated_4), 2] <- "Dane_9019_1"
task_ids_collated_4[nrow(task_ids_collated_4) + 1, 1] <- "bc934c89-0b5b-44e8-9734-1d58be225644"
task_ids_collated_4[nrow(task_ids_collated_4), 2] <- "Dane_9019_2"
task_ids_collated_4[nrow(task_ids_collated_4) + 1, 1] <- "7bc34e10-3140-4dd1-81fb-a37d393bcf5b"
task_ids_collated_4[nrow(task_ids_collated_4), 2] <- "Dane_9019_3"
task_ids_collated_4[nrow(task_ids_collated_4) + 1, 1] <- "c61d01b7-ffa5-4c87-b55e-6d41e456f5d0"
task_ids_collated_4[nrow(task_ids_collated_4), 2] <- "Dane_583_1"
task_ids_collated_4[nrow(task_ids_collated_4) + 1, 1] <- "d63e689e-40a2-4d8b-a503-fb12053531b0"
task_ids_collated_4[nrow(task_ids_collated_4), 2] <- "Dane_583_2"
task_ids_collated_4[nrow(task_ids_collated_4) + 1, 1] <- "cb992e37-bb6a-4e38-bef9-5ca85d447863"
task_ids_collated_4[nrow(task_ids_collated_4), 2] <- "Dane_583_3"
task_ids_collated_4[nrow(task_ids_collated_4) + 1, 1] <- "a42fa2ad-e424-4865-8868-4a276c5cc544"
task_ids_collated_4[nrow(task_ids_collated_4), 2] <- "Dane_583_4"
task_ids_collated_4[nrow(task_ids_collated_4) + 1, 1] <- "c4165b09-770e-4182-ae63-32e237486d84"
task_ids_collated_4[nrow(task_ids_collated_4), 2] <- "Dane_583_5"
write.csv(task_ids_collated_4,"task_ids_collated_4.csv",row.names = FALSE)

for(i in 1:nrow(task_ids_collated_4)){
	task_id <- task_ids_collated_4[i,1]
	task_name <- task_ids_collated_4[i,2]
	outDir <- file.path(paste('env-data/MODIS/MODIS_4/',task_name,sep=""))
	suppressWarnings(dir.create(outDir))  
	response <- GET(paste0(API_URL, "bundle/", task_id), add_headers(Authorization = token))
	bundle_response <- prettify(toJSON(content(response), auto_unbox = TRUE))
	bundle <- fromJSON(bundle_response)$files
	bundle <- bundle[which(bundle$file_type=="csv"),]
		for (j in bundle$file_id){
  			filename <- bundle[bundle$file_id == j,]$file_name           
  			filepath <- paste(outDir,filename, sep = "/")
  			suppressWarnings(dir.create(dirname(filepath)))
  			response <- GET(paste0(API_URL, "bundle/", task_id, "/", j), 
                  write_disk(filepath, overwrite = TRUE), progress(),
                  add_headers(Authorization = token))
		}
	}

find . -type f -name '*-Psn-QC-500m-lookup.csv' -exec rm {} +
find . -type f -name '*-Psn-QC-500m-Statistics-QA.csv' -exec rm {} +
find . -type f -name '*-Npp-QC-500m-Statistics-QA.csv' -exec rm {} +

library(dplyr)
library(fasstr)
library(lubridate)
library(purrr)
library(zoo)

MYD11A1.collated <- data.frame(matrix(ncol = 5, nrow = 0))
colnames(MYD11A1.collated) <- c('AreaNumb', 'Date', 'Avg_LST_Day7', 'Avg_LST_Night7')

for (i in c(4,5,6,13,26,81,83,86,88,158,164,189,218,243,245,247,249,251,253,254,256,319,335,340,369,371,391,397,417,465,513,515,522,523,526,530,552,556,558,559,565,567,584,587,590,597,598,607,615,618,619,620,676,700,715,718,719,724,725,900,969,972,973,975,978,1931,1932,2041,3262,3333,3921,3922,4022,4130,8300,8301,8512,9007,9009,9014,9098,9601,9970,9971)) {
	for (j in c(1,2,3)){
		FileName <- paste('env-data/MODIS/MODIS_',j,'/Dane_',i,'_',j,'/MYD11A1-006-Statistics.csv',sep="")
		MYD11A1 <- read.csv(FileName)
		MYD11A1.Day <- MYD11A1[which(MYD11A1$Dataset=="LST_Day_1km"),]
		MYD11A1.Night <- MYD11A1[which(MYD11A1$Dataset=="LST_Night_1km"),]
		data.Day <- MYD11A1.Day %>%
 			mutate(Date = ymd(Date)) %>%
  			fill_missing_dates(dates = Date, fill_end_years = FALSE) %>%
  			mutate(Avg_LST_Day7 = map_dbl(Date, ~mean(Mean[(Date > . - 7) & (Date < .)], na.rm=TRUE)))
  		data.Day <- data.Day[-c(1, 2, 3, 4, 5, 6, 7), ] 		
		data.Day <- na.locf(data.Day, na.rm = TRUE)
		data.Night <- MYD11A1.Night %>%
 			mutate(Date = ymd(Date)) %>%
  			fill_missing_dates(dates = Date, fill_end_years = FALSE) %>%
  			mutate(Avg_LST_Night7 = map_dbl(Date, ~mean(Mean[(Date > . - 7) & (Date < .)], na.rm=TRUE)))
  		data.Night <- data.Night[-c(1, 2, 3, 4, 5, 6, 7), ] 		
		data.Night <- na.locf(data.Night, na.rm = TRUE)
		MYD11A1.final <- merge(data.Day, data.Night, by = 'Date', all=TRUE)
		MYD11A1.final$AreaNumb <- paste('S',i,sep="")
		MYD11A1.final2=subset(MYD11A1.final,select=c(AreaNumb,Date,Avg_LST_Day7,Avg_LST_Night7))
		MYD11A1.collated <- rbind(MYD11A1.final2,MYD11A1.collated)
		}
	}
	
for (i in c(583,9019)) {
	for (j in c(1,2,3)){
		FileName <- paste('env-data/MODIS/MODIS_4','/Dane_',i,'_',j,'/MYD11A1-006-Statistics.csv',sep="")
		MYD11A1 <- read.csv(FileName)
		MYD11A1.Day <- MYD11A1[which(MYD11A1$Dataset=="LST_Day_1km"),]
		MYD11A1.Night <- MYD11A1[which(MYD11A1$Dataset=="LST_Night_1km"),]
		data.Day <- MYD11A1.Day %>%
 			mutate(Date = ymd(Date)) %>%
  			fill_missing_dates(dates = Date, fill_end_years = FALSE) %>%
  			mutate(Avg_LST_Day7 = map_dbl(Date, ~mean(Mean[(Date > . - 7) & (Date < .)], na.rm=TRUE)))
  		data.Day <- data.Day[-c(1, 2, 3, 4, 5, 6, 7), ] 		
		data.Day <- na.locf(data.Day, na.rm = TRUE)
		data.Night <- MYD11A1.Night %>%
 			mutate(Date = ymd(Date)) %>%
  			fill_missing_dates(dates = Date, fill_end_years = FALSE) %>%
  			mutate(Avg_LST_Night7 = map_dbl(Date, ~mean(Mean[(Date > . - 7) & (Date < .)], na.rm=TRUE)))
  		data.Night <- data.Night[-c(1, 2, 3, 4, 5, 6, 7), ] 		
		data.Night <- na.locf(data.Night, na.rm = TRUE)
		MYD11A1.final <- merge(data.Day, data.Night, by = 'Date', all=TRUE)
		MYD11A1.final$AreaNumb <- paste('S',i,sep="")
		MYD11A1.final2=subset(MYD11A1.final,select=c(AreaNumb,Date,Avg_LST_Day7,Avg_LST_Night7))
		MYD11A1.collated <- rbind(MYD11A1.final2,MYD11A1.collated)
		}
	}
	
write.csv(MYD11A1.collated,"env-data/MYD11A1-collated.csv",row.names = FALSE)

MYD13Q1.collated <- data.frame(matrix(ncol = 5, nrow = 0))
colnames(MYD13Q1.collated) <- c('AreaNumb', 'Date', 'Mean_EVI', 'Mean_NDVI')

for (i in c(4,5,6,13,26,81,83,86,88,158,164,189,218,243,245,247,249,251,253,254,256,319,335,340,369,371,391,397,417,465,513,515,522,523,526,530,552,556,558,559,565,567,584,587,590,597,598,607,615,618,619,620,676,700,715,718,719,724,725,900,969,972,973,975,978,1931,1932,2041,3262,3333,3921,3922,4022,4130,8300,8301,8512,9007,9009,9014,9098,9601,9970,9971)) {
	for (j in c(1,2,3)){
		FileName <- paste('env-data/MODIS/MODIS_',j,'/Dane_',i,'_',j,'/MYD13Q1-006-Statistics.csv',sep="")
		MYD13Q1 <- read.csv(FileName)
		MYD13Q1.EVI <- MYD13Q1[which(MYD13Q1$Dataset=="_250m_16_days_EVI"),]
		MYD13Q1.NDVI <- MYD13Q1[which(MYD13Q1$Dataset=="_250m_16_days_NDVI"),]
		data.EVI <- MYD13Q1.EVI %>%
 			mutate(Date = ymd(Date)) %>%
  			fill_missing_dates(dates = Date, fill_end_years = FALSE)
  		data.EVI <- data.EVI[-c(1:16), ] 		
		data.EVI <- na.locf(data.EVI, na.rm = TRUE)
		data.EVI <- rename(data.EVI, Mean_EVI = Mean)
		data.NDVI <- MYD13Q1.NDVI %>%
  			mutate(Date = ymd(Date)) %>%
  			fill_missing_dates(dates = Date, fill_end_years = FALSE)
  		data.NDVI <- data.NDVI[-c(1:16), ] 		
		data.NDVI <- na.locf(data.NDVI, na.rm = TRUE)
		data.NDVI <- rename(data.NDVI, Mean_NDVI = Mean)
		MYD13Q1.final <- merge(data.EVI, data.NDVI, by = 'Date', all=TRUE)
		MYD13Q1.final$AreaNumb <- paste('S',i,sep="")
		MYD13Q1.final2=subset(MYD13Q1.final,select=c(AreaNumb,Date,Mean_EVI,Mean_NDVI))
		MYD13Q1.collated <- rbind(MYD13Q1.final2,MYD13Q1.collated)
		}
	}

for (i in c(583,9019)) {
	for (j in c(1,2,3)){
		FileName <- paste('env-data/MODIS/MODIS_4','/Dane_',i,'_',j,'/MYD13Q1-006-Statistics.csv',sep="")
		MYD13Q1 <- read.csv(FileName)
		MYD13Q1.EVI <- MYD13Q1[which(MYD13Q1$Dataset=="_250m_16_days_EVI"),]
		MYD13Q1.NDVI <- MYD13Q1[which(MYD13Q1$Dataset=="_250m_16_days_NDVI"),]
		data.EVI <- MYD13Q1.EVI %>%
 			mutate(Date = ymd(Date)) %>%
  			fill_missing_dates(dates = Date, fill_end_years = FALSE)
  		data.EVI <- data.EVI[-c(1:16), ] 		
		data.EVI <- na.locf(data.EVI, na.rm = TRUE)
		data.EVI <- rename(data.EVI, Mean_EVI = Mean)
		data.NDVI <- MYD13Q1.NDVI %>%
  			mutate(Date = ymd(Date)) %>%
  			fill_missing_dates(dates = Date, fill_end_years = FALSE)
  		data.NDVI <- data.NDVI[-c(1:16), ] 		
		data.NDVI <- na.locf(data.NDVI, na.rm = TRUE)
		data.NDVI <- rename(data.NDVI, Mean_NDVI = Mean)
		MYD13Q1.final <- merge(data.EVI, data.NDVI, by = 'Date', all=TRUE)
		MYD13Q1.final$AreaNumb <- paste('S',i,sep="")
		MYD13Q1.final2=subset(MYD13Q1.final,select=c(AreaNumb,Date,Mean_EVI,Mean_NDVI))
		MYD13Q1.collated <- rbind(MYD13Q1.final2,MYD13Q1.collated)
		}
	}

write.csv(MYD13Q1.collated,"env-data/MYD13Q1-collated.csv",row.names = FALSE)

MYD17A2HGF.collated <- data.frame(matrix(ncol = 5, nrow = 0))
colnames(MYD17A2HGF.collated) <- c('AreaNumb', 'Date', 'Mean_Gpp')

for (i in c(4,5,6,13,26,81,83,86,88,158,164,189,218,243,245,247,249,251,253,254,256,319,335,340,369,371,391,397,417,465,513,515,522,523,526,530,552,556,558,559,565,567,584,587,590,597,598,607,615,618,619,620,676,700,715,718,719,724,725,900,969,972,973,975,978,1931,1932,2041,3262,3333,3921,3922,4022,4130,8300,8301,8512,9007,9009,9014,9098,9601,9970,9971)) {
	for (j in c(1,2,3)){
		FileName <- paste('env-data/MODIS/MODIS_',j,'/Dane_',i,'_',j,'/MYD17A2HGF-006-Statistics.csv',sep="")
		if (!file.exists(FileName)) next
		MYD17A2HGF <- read.csv(FileName)
		data <- MYD17A2HGF %>%
  			mutate(Date = ymd(Date)) %>%
  			fill_missing_dates(dates = Date, fill_end_years = FALSE)
  		data <- data[-c(1:8), ] 		
		data <- na.locf(data, na.rm = TRUE)
		data <- rename(data, Mean_Gpp = Mean)
		data$AreaNumb <- paste('S',i,sep="")
		data2=subset(data,select=c(AreaNumb,Date,Mean_Gpp))
		MYD17A2HGF.collated <- rbind(data2,MYD17A2HGF.collated)
		}
	}

for (i in c(583,9019)) {
	for (j in c(1,2,3)){
		FileName <- paste('env-data/MODIS/MODIS_4','/Dane_',i,'_',j,'/MYD17A2HGF-006-Statistics.csv',sep="")
		if (!file.exists(FileName)) next
		MYD17A2HGF <- read.csv(FileName)		
		data <- MYD17A2HGF %>%
  			mutate(Date = ymd(Date)) %>%
  			fill_missing_dates(dates = Date, fill_end_years = FALSE)
  		data <- data[-c(1:8), ] 		
		data <- na.locf(data, na.rm = TRUE)
		data <- rename(data, Mean_Gpp = Mean)
		data$AreaNumb <- paste('S',i,sep="")
		data2=subset(data,select=c(AreaNumb,Date,Mean_Gpp))
		MYD17A2HGF.collated <- rbind(data2,MYD17A2HGF.collated)
		}
	}

write.csv(MYD17A2HGF.collated,"env-data/MYD17A2HGF-collated.csv",row.names = FALSE)

MYD17A3HGF.collated <- data.frame(matrix(ncol = 5, nrow = 0))
colnames(MYD17A3HGF.collated) <- c('AreaNumb', 'Date', 'Mean_Npp')

for (i in c(4,5,6,13,26,81,83,86,88,158,164,189,218,243,245,247,249,251,253,254,256,319,335,340,369,371,391,397,417,465,513,515,522,523,526,530,552,556,558,559,565,567,584,587,590,597,598,607,615,618,619,620,676,700,715,718,719,724,725,900,969,972,973,975,978,1931,1932,2041,3262,3333,3921,3922,4022,4130,8300,8301,8512,9007,9009,9014,9098,9601,9970,9971)) {
	for (j in c(1,2,3)){
		FileName <- paste('env-data/MODIS/MODIS_',j,'/Dane_',i,'_',j,'/MYD17A3HGF-006-Statistics.csv',sep="")
		if (!file.exists(FileName)) next
		MYD17A3HGF <- read.csv(FileName)		
		data <- MYD17A3HGF %>%
  			mutate(Date = ymd(Date)) %>%
  			fill_missing_dates(dates = Date, fill_end_years = FALSE)
		data <- na.locf(data, na.rm = TRUE)
		data <- rename(data, Mean_Npp = Mean)
		data$AreaNumb <- paste('S',i,sep="")
		data2=subset(data,select=c(AreaNumb,Date,Mean_Npp))
		MYD17A3HGF.collated <- rbind(data2,MYD17A3HGF.collated)
		}
	}

for (i in c(583,9019)) {
	for (j in c(1,2,3)){
		FileName <- paste('env-data/MODIS/MODIS_4','/Dane_',i,'_',j,'/MYD17A3HGF-006-Statistics.csv',sep="")
		if (!file.exists(FileName)) next
		MYD17A3HGF <- read.csv(FileName)		
		data <- MYD17A3HGF %>%
  			mutate(Date = ymd(Date)) %>%
  			fill_missing_dates(dates = Date, fill_end_years = FALSE)
		data <- na.locf(data, na.rm = TRUE)
		data <- rename(data, Mean_Npp = Mean)
		data$AreaNumb <- paste('S',i,sep="")
		data2=subset(data,select=c(AreaNumb,Date,Mean_Npp))
		MYD17A3HGF.collated <- rbind(data2,MYD17A3HGF.collated)
		}
	}

write.csv(MYD17A3HGF.collated,"env-data/MYD17A3HGF-collated.csv",row.names = FALSE)

data <- read.csv("env-data/MYD11A1-collated.csv", header = TRUE)
data <- data %>% mutate(Date = ymd(Date))

data2 <- read.csv("env-data/MYD13Q1-collated.csv", header = TRUE)
data2 <- data2 %>% mutate(Date = ymd(Date))

data3 <- read.csv("env-data/MYD17A2HGF-collated.csv", header = TRUE)
data3 <- data3 %>% mutate(Date = ymd(Date))

data4 <- read.csv("env-data/MYD17A3HGF-collated.csv", header = TRUE)
data4 <- data4 %>% mutate(Date = ymd(Date))

dane=read.csv("Combined-Master.csv",header=TRUE)
dane$AreaNumb <- paste('S',dane$AreaNumb,sep="")
dane <- dane %>%
  mutate(SampDate = mdy(SampDate)) %>%  
  rename(
    Date = SampDate
    )

lookup <- unique(dane)
dane2 <- merge(data, lookup, by = c('AreaNumb','Date'), all = TRUE)
dane3 <- merge(data2, lookup, by = c('AreaNumb','Date'), all = TRUE)
dane3 <- subset(dane3,select=c(AreaNumb,Date,Mean_EVI,Mean_NDVI))
dane4 <- merge(data3, lookup, by = c('AreaNumb','Date'), all = TRUE)
dane4 <- subset(dane4,select=c(AreaNumb,Date,Mean_Gpp))
dane5 <- merge(data4, lookup, by = c('AreaNumb','Date'), all = TRUE)
dane5 <- subset(dane5,select=c(AreaNumb,Date,Mean_Npp))
dane.modis.collated <- merge(dane2, dane3, by = c('AreaNumb', 'Date'), all = TRUE)
dane.modis.collated <- merge(dane.modis.collated, dane4, by = c('AreaNumb', 'Date'), all = TRUE)
dane.modis.collated <- merge(dane.modis.collated, dane5, by = c('AreaNumb', 'Date'), all = TRUE)
dane.modis.collated <- dane.modis.collated[, c(1, 6, 2, 5, 7:19, 3, 4, 20:23)]
dane.modis.collated <- dane.modis.collated[order(dane.modis.collated[,1], dane.modis.collated[,3] ),]
dane.modis.collated <- dane.modis.collated[which(!is.na(dane.modis.collated$SiteVis)),]

write.csv(dane.modis.collated,"dane-modis-collated.csv",row.names = FALSE)

data=read.csv("dane-modis-collated.csv",header=TRUE)
data.hist=data[which(data$Year==2007|data$Year==2008|data$Year==2009|data$Year==2010|data$Year==2011|data$Year==2012|data$Year==2013|data$Year==2014|data$Year==2015|data$Year==2016|data$Year==2017|data$Year==2018),]
data.hist=subset(data.hist,select=c(AreaNumb,Avg_LST_Day7,Avg_LST_Night7,Mean_EVI,Mean_NDVI,Mean_Gpp,Mean_Npp))
avg.by.site=aggregate(.~AreaNumb,data.hist,function(x) mean(x,na.rm=TRUE),na.action=na.pass)
AvgLSTDay=avg.by.site$Avg_LST_Day7
AvgLSTNight=avg.by.site$Avg_LST_Night7
AvgEVI=avg.by.site$Mean_EVI
AvgNDVI=avg.by.site$Mean_NDVI
AvgGpp=avg.by.site$Mean_Gpp
AvgNpp=avg.by.site$Mean_Npp
AreaNumb=avg.by.site$AreaNumb

library(gtools)
df<-data.frame(AreaNumb,AvgLSTDay,AvgLSTNight,AvgEVI,AvgNDVI,AvgGpp,AvgNpp)
clust.data=read.csv("alpha-div.csv",header=TRUE)
clust.data$AreaNumb <- paste('S',clust.data$AreaNumb,sep="")
data.total=merge(clust.data,df,by="AreaNumb", all = TRUE)
data.total=subset(data.total,select=-c(faith_pd,observed_features,inv_simpson,cfus_per_ml))
data.total=data.total[order(nchar(data.total$AreaNumb), data.total$AreaNumb),]
row.names(data.total) <- 1:nrow(data.total)
data.total$AvgGpp[is.nan(data.total$AvgGpp)]<-NA
data.total$AvgNpp[is.nan(data.total$AvgNpp)]<-NA

write.csv(data.total,"hist-env-data.csv",row.names = FALSE)

data=read.csv("dane-modis-collated.csv",header=TRUE)
data.curr=data[which(data$Year==2019|data$Year==2020|data$Year==2021),]
data.curr=subset(data.curr,select=c(AreaNumb,Year,Avg_LST_Day7,Avg_LST_Night7,Mean_EVI,Mean_NDVI,Mean_Gpp,Mean_Npp))
avg.by.site.year=aggregate(.~AreaNumb+Year,data.curr,function(x) mean(x,na.rm=TRUE),na.action=na.pass)
AvgLSTDay=avg.by.site.year$Avg_LST_Day7
AvgLSTNight=avg.by.site.year$Avg_LST_Night7
AvgEVI=avg.by.site.year$Mean_EVI
AvgNDVI=avg.by.site.year$Mean_NDVI
AvgGpp=avg.by.site.year$Mean_Gpp
AvgNpp=avg.by.site.year$Mean_Npp
AreaNumb=avg.by.site.year$AreaNumb
Year=avg.by.site.year$Year

library(gtools)
df<-data.frame(AreaNumb,Year,AvgLSTDay,AvgLSTNight,AvgEVI,AvgNDVI,AvgGpp,AvgNpp)
df=df[order(nchar(df$AreaNumb), df$AreaNumb),]
row.names(df) <- 1:nrow(df)
df$AvgGpp[is.nan(df$AvgGpp)]<-NA
df$AvgNpp[is.nan(df$AvgNpp)]<-NA

write.csv(df,"env-data/MODIS/curr-env-data-by-site-year.csv",row.names = FALSE)











































































