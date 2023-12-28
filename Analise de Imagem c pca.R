library(EBImage)
library(tiff)
library(dplyr)
library(openxlsx)
rm(list = ls())

setwd("~/Documents/Image analysis /Teste/test_images")
study_path <- getwd() #This pulls your working directory out so we don't have to keep pulling that info
image_h_path <- paste0(study_path, "/image_h")
image_d_path <- paste0(study_path, "/image_d")
segment_h_path <- paste0(study_path, "/segment_h")
segment_d_path <- paste0(study_path, "/segment_d")
features_h_path <- paste0(study_path, "/features_h")
features_d_path <- paste0(study_path, "/features_d")
overlay_h_path <- paste0(study_path, "/overlay_h")
overlay_d_path <- paste0(study_path, "/overlay_d")
test_images_path <- paste0(study_path, "/test_images")

#Set Paths
my.list <- list.files(path = test_images_path, pattern = ".tif")
my_files <- my.list
my_intensity <- c('Low', 'Low', 'Low', 'Low', 'Low', 'Low', 'Low', 'Low', 'Low', 'Low',
                  'Mid', 'Mid', 'Mid', 'Mid', 'Mid', 'Mid', 'Mid', 'Mid', 'Mid', 'Mid',
                  'No', 'No', 'No', 'No', 'No', 'No', 'No', 'No', 'No', 'No',
                  'High', 'High', 'High', 'High', 'High', 'High', 'High', 'High', 'High')
meu_list_pca <- data.frame()
meu_list_matrix <- list()
new_df <- data.frame()

for(i in 1:length(my_files)) {#i = 12
  im <- readTIFF(my_files[i])
  im <- rgbImage(red = im[, ,1] , green = im[,,2] , blue = im[,,3])
  im_col <- im
  
  colorMode(im) <- Grayscale
  im1 <- 0.299 * im[,,1] + 0.587 * im[,,2] + 0.114 *im[,,3]
  
  im2 <- 1-im1
  test_median <-   gblur(im2, sigma = 2)
  kern <- makeBrush(5, shape='diamond') 
  test_median2 <-  erode(test_median, kern)
  test_median <- test_median>otsu(test_median)
  test_median2 <- test_median2>otsu(test_median2)
  test_median <- fillHull(test_median)
  test_median2 <- fillHull(test_median2)
  display(test_median2)
  
  
  #dab1_w <- watershed(distmap(test_median))
  dab1_w2 <- watershed(distmap(test_median2))
  
  
  b <- paintObjects(x = dab1_w2, tgt = rgbImage("green"= im2), col = "white", thick = TRUE)  
  
  display(b)
  dab1_w <- bwlabel(dab1_w2)
  my_compute_r <- computeFeatures(x = dab1_w2 , ref = im[, ,1])
  colnames(my_compute_r)<- paste("red", colnames(my_compute_r), sep = "_")
  my_compute_g <- computeFeatures(x = dab1_w2 , ref = im[, ,2])
  colnames(my_compute_g)<- paste("green", colnames(my_compute_g), sep = "_")
  
  my_compute_b <- computeFeatures(x = dab1_w2 , ref = im[, ,3])
  colnames(my_compute_b)<- paste("blue", colnames(my_compute_b), sep = "_")
  
  my_compute <- cbind(my_compute_r, my_compute_g, my_compute_b)
  pdf(file = paste("numbplot_", gsub(pattern = ".tif", replacement = "",x =  my_files[i] ), ".pdf", sep =   "") )
  
  plot(im_col)
  text(my_compute[, "red_x.0.m.cx"],
       my_compute[, "red_x.0.m.cy"], 
       labels = seq_len(nrow(my_compute )),
       col = "green4", cex = .2
  )
  dev.off()
  
  my_compute <- as.data.frame(my_compute)
  #my_compute$image <- rep(gsub(pattern = ".tif", replacement = "",x =  my_files[i] ))
  meu_list_pca <- prcomp(my_compute, scale = TRUE)
  #names(meu_list_pca)[5] <- my_files[i]
  meu_list_matrix<- my_compute[1:20,] %>% as.data.frame
  names(meu_list_matrix) <- names(my_compute)
  meu_list_matrix$slide<- my_files[i]
  meu_list_matrix$intensity<- my_intensity[i]

  new_df <- rbind(new_df, meu_list_matrix)
  #saveRDS(meu_list_matrix, file = "meu_list_matrix.R")
  print(i)
}

write.xlsx(new_df, "rf_model_pca.xlsx")

dim(new_df)

library(randomForest)
library(caret)
library(caTools)
setwd("~/Documents/Image analysis ")
seg_features <- read.xlsx("rf_model/rf_model_pca.xlsx")
seg_features$Type <- as.factor(seg_features$intensity)

#RF modeling from labeled nuclei from positive and negative nuclei
set.seed(105)
sample <- sample.split(seg_features$Type, SplitRatio = 0.7) 
train = subset(seg_features, sample==T)
test = subset(seg_features, sample==F)
rf1 <- randomForest(Type~ B.a.b.mean + B.a.b.q05 + B.a.b.q095 + B.a.b.q099, data = train, ntree=2000, na.action = na.exclude)

#Prediction
predict<-predict(rf1,test)
output <- confusionMatrix(as.factor(predict), as.factor(test$Type))
output

#Load model
#rf1 <- readRDS("rf_model/RF_model.RDS")

#Read in DAB features
d_features_list <- list.files(path = "features_d")

#Make new updated features folder
dir.create("updated_features_d/")

for (i in 1:length(d_features_list)){#i=1
  print(i)
  d_features <- read.csv(paste0("features_d/", d_features_list[i]))
  intensity <- predict(rf1, d_features)
  d_features$Intensity <- intensity
  write.csv(d_features, paste0("updated_features_d/", d_features_list[i]))
}

#Stain Segments
d_features_list <- list.files(path = "updated_features_d")
#Nuclei Segments
h_features_list <- list.files(path = "features_h")

#Remove control segmentations
d_features_list <- d_features_list [!grepl("pc", d_features_list)]
d_features_list <- d_features_list [!grepl("nc", d_features_list)]

h_features_list <- h_features_list [!grepl("pc", h_features_list)]
h_features_list <- h_features_list [!grepl("nc", h_features_list)]

d_features_list
h_features_list

ratio_frame <- NULL
frame_d <- NULL
frame_h <- NULL
i=1
for (i in 1:length(d_features_list)){ 
  print(i)  #Keep track of where you are
  ff <- read.csv(paste0("updated_features_d/", d_features_list[i])) #Read in D features
  pull <- ff[!ff$Intensity == "No",] #Filter No Segmentations
  pull2 <- pull[!pull$Intensity == "Low",] #Filter No/Low Segmentations
  pull$Image <- d_features_list[i] #Add image label for reference
  frame_d <- rbind(frame_d, pull) #Store segmentations (we don't need to store the low filter data since the intensity is stored anyways)
  
  hf <- read.csv(paste0("features_h/", h_features_list[i])) #Read in H features
  frame_h <- rbind(frame_h, hf) #Store in data table
  
  #Caclulate Ratios
  rat <- sum(pull$B.0.s.area)/sum(hf$B.0.s.area) #Caculate D:H Ration
  rat2 <- sum(pull2$B.0.s.area)/sum(hf$B.0.s.area) #Caculate D:H Ration (no low as well)
  
  store <- data.frame("Ratio_N" = rat, "Ration_NL" = rat2, "Image_D" = d_features_list[i], "Image_H" = h_features_list[i]) #I like to store file name to make sure I didnt compare the wrong filtes together.
  
  #Add other sums
  store$Area_D_N <- sum(pull$B.0.s.area)
  store$Area_D_NL <- sum(pull2$B.0.s.area)
  store$Area_H <- sum(hf$B.0.s.area)
  ratio_frame <- rbind(ratio_frame, store)
}

write.csv(ratio_frame, "Ratio_teste.csv")
