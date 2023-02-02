
library(Morpho)
library(Rvcg)
library(rgl)
library(caret)
library(mixOmics)
library(fsbrain)
library(RVAideMemoire)
library(doParallel)

source("VTK_ops.R")



########## set values and load data ########## 

# freesurfer DIR
subjects_dir = "~/data4/FTD_Lund/ftd_subjects_FS6/"

# folder which corresponds to where the processing will occur
base_folder<-"~/data/deformetrica" 

# file name to prefix all the analysis
test_name<-paste("prefix_")

# parameters for number of folds and repeats for analysis
Nrepeat.outer=10
Nfold.outer<-5
Nrepeat.inner=10
Nfold.inner<-5

#### load the CSV file, assumes there's headings called "Group", "Age" and "ICV" but these can be changed and added to in the code below
csvfile<-read.csv(file.path(base_folder,"data.csv")) 

test_var<-"Group_variable" # from CSV file that you want to test, or change to correlation variable and change the value below
run_group=1 # change to 0 for correlation analysis

Cort<-as.matrix(csvfile[,c(1:68)]) # from csv file select which columns to include in cortical calculation.
DTI<-as.matrix(csvfile[c(1:8)]) # from csv file select which columns to include in DTI calculation.

# this tells us how many steps to loop through when assessing sparsity
test.keepX = list (Thal = c(seq(20, 800, 10)),
                   Cort = c(seq(10, 68, 2)),
                   DTI=seq(1,8,1))

# need to change this to point to deformetrica Qt5 install location
QT_PATH<-"~/USER/anaconda3/envs/deformetrica/lib/python3.8/site-packages/PyQt5/Qt5/lib"  


########## rest below shouldn't need changing  ########## 

#### load the momenta

side="rh"
mom_file=file.path(paste("output_",side,"_",struct,sep=""),"DeterministicAtlas__EstimatedParameters__Momenta.txt")
mom_points_rh<-as.integer(read.table(mom_file,nrows = 1,header = FALSE )[2])
mom_points<-mom_points_rh
mom_read<-read.table(mom_file,skip=2,header = FALSE )
mom<-array(0,dim=c(nrow(csvfile),mom_points*3))
for (i in 0:nrow(csvfile)-1) { mom[i+1,]<-unlist(mom_read[(i*mom_points+1):(i*mom_points+mom_points),]) } #inelegant, but it reshapes the data as needed. All x variables first, then y, then z
rh_mom<-mom

side="lh"
mom_file=file.path(paste("output_",side,"_",struct,sep=""),"DeterministicAtlas__EstimatedParameters__Momenta.txt")
mom_points_lh<-as.integer(read.table(mom_file,nrows = 1,header = FALSE )[2])
mom_points<-mom_points_lh
mom_read<-read.table(mom_file,skip=2,header = FALSE )
mom<-array(0,dim=c(nrow(csvfile),mom_points*3))
for (i in 0:nrow(csvfile)-1) { mom[i+1,]<-unlist(mom_read[(i*mom_points+1):(i*mom_points+mom_points),]) } #inelegant, but it reshapes the data as needed. All x variables first, then y, then z
lh_mom<-mom


########## ANALYSIS ########## 

dataset.comb<-cbind(lh_mom,rh_mom)
Thal<-residuals(lm(dataset.comb~scale(csvfile$Age)+scale(csvfile$ICV)))   ###### if needed can change covariates here
Cort<-residuals(lm(Cort~scale(csvfile$Age)))
DTI<-residuals(lm(DTI~scale(csvfile$Age)))

rownames(Thal)<-NULL
rownames(Cort)<-NULL
rownames(DTI)<-NULL

X<-list(Thal=Thal,Cort=Cort,DTI=DTI)

##### Group Analysis ####
if (run_group == 1 ) {
  
  Nrepeat<-Nrepeat.outer
  Nfold<-Nfold.outer
  Mrepeat<-Nrepeat.inner
  Mfold<-Nfold.inner
  
  Y<-csvfile[[test_var]]
  design = matrix(1, ncol = length(X)+1, nrow = length(X)+1,dimnames = list(c(names(X),"Y"), c(names(X),"Y")))
  diag(design) = 0
  
  nest_repeat<-matrix(NA,nrow = Nrepeat,ncol = 4)
  nest_fold<-matrix(NA,nrow = Nfold,ncol = 4)
  
  folds <- cut(seq(1,nrow(X[[1]])),breaks=Nfold,labels=FALSE) #create folds
  
  start_time <- Sys.time()
  pb <- txtProgressBar(min = 0, max = Nrepeat, style = 3)
  
  for (i_repeat in 1:Nrepeat){
    setTxtProgressBar(pb, i_repeat)
    
    folds <- createFolds(factor(Y), k = 5, list = FALSE)
    indx <- which(folds==1,arr.ind=TRUE)
    X.test<-lapply(X, "[",indx,)
    X.train<-lapply(X, "[",-indx,)
    Y.test<-(Y[indx])
    Y.train<-(Y[-indx])
    
    for (i_fold in 1:Nfold) {
      indx <- which(folds==i_fold,arr.ind=TRUE)
      X.test<-lapply(X, "[",indx,)
      X.train<-lapply(X, "[",-indx,)
      Y.test<-(Y[indx])
      Y.train<-(Y[-indx])
      
      tune.diablo = tune.block.splsda(X = X.train, Y = Y.train, ncomp = 1,
                                      test.keepX = test.keepX, design = design,
                                      validation = 'Mfold', folds = Mfold, nrepeat = Mrepeat,
                                      dist = "centroids.dist")
      
      mod.diablo = block.splsda(X = X.train, Y = Y.train, ncomp = 1,
                                keepX = tune.diablo$choice.keepX, design = design)
      test<-predict(mod.diablo,X.test, method="centroids.dist")
      
      a<-tryCatch(table(Y.test,test$WeightedPredict.class$max.dist)[1,1], error = function(e) 0)
      b<-tryCatch(table(Y.test,test$WeightedPredict.class$max.dist)[1,2], error = function(e) 0)
      c<-tryCatch(table(Y.test,test$WeightedPredict.class$max.dist)[2,1], error = function(e) 0)
      d<-tryCatch(table(Y.test,test$WeightedPredict.class$max.dist)[2,2], error = function(e) 0)
      
      BER<-0.5*(b/(a+b) + c/(c+d))
      
      nest_fold[i_fold,]<-c(unlist(tune.diablo$choice.keepX),min(tune.diablo$error.rate))
    }
    
    nest_repeat[i_repeat,]<-nest_fold[which(nest_fold[,4]==min(nest_fold[,4],na.rm = T)), ]
  }
  
  close(pb)
  
  keepX <-  list (Thal = nest_repeat[which(nest_repeat[,4]==min(nest_repeat[,4],na.rm = T)),1],
                  Cort = nest_repeat[which(nest_repeat[,4]==min(nest_repeat[,4],na.rm = T)),2],
                  DTI  = nest_repeat[which(nest_repeat[,4]==min(nest_repeat[,4],na.rm = T)),3])  
  
  mod.diablo = block.splsda(X = X, Y = Y, ncomp = 1, keepX = keepX, design = design)
  
  mod.diablo.perm<-DIABLO.test(mod.diablo,nperm=100,cpus=1)
  mod<-mod.diablo
  
}

#### Correlation analysis ####
if (run_group == 0) {
  
  Y<-as.matrix(csvfile[[test_var]])
  design = matrix(1, ncol = length(X)+1, nrow = length(X)+1,dimnames = list(c(names(X),"Y"), c(names(X),"Y")))
  diag(design) = 0
  
  MSE_array<-cbind(expand.grid(test.keepX),0)
  
  MSE.inner<-cbind(expand.grid(c(test.keepX,1:Nrepeat.inner)),0)
  MSE.inner.temp<-numeric(length=Nfold.inner)
  MSE.outer<-matrix(nrow=Nfold.outer*Nrepeat.outer,ncol = 4)
  
  for (i.Nrepeat.outer in 1:Nrepeat.outer) {
    
    folds.outer <- caret::createFolds(Y, k = Nfold.outer, list = FALSE) # define index of outer folds
    
    for (i.fold.outer in 1:Nfold.outer) { # for each fold in the outer loop
      
      indx.outer <- which(folds.outer==i.fold.outer,arr.ind=TRUE)  # use the index to make test/training outer set
      X.test.outer<-lapply(X, "[",indx.outer,)
      X.train.outer<-lapply(X, "[",-indx.outer,)
      Y.test.outer<-as.matrix(Y[indx.outer])
      Y.train.outer<-as.matrix(Y[-indx.outer])
      
      # loop through KeepX options
      for (i.MSE in 1:nrow(MSE_array)){           
        
        i.keepX<-MSE_array[i.MSE,1:3]
        
        
        for (i.Nrepeat.inner in 1:Nrepeat.inner) {  # repeats of inner loop
          
          folds.inner <- caret::createFolds(Y.train.outer, k = Nfold.inner, list = FALSE)  # again create inner index variable
          
          for (i.fold.inner in 1:Nfold.inner) { # for each fold in the inner loop
            
            indx.inner <- which(folds.inner==i.fold.inner,arr.ind=TRUE) # split into test and training inner set
            X.test.inner<-lapply(X.train.outer, "[",indx.inner,)
            X.train.inner<-lapply(X.train.outer, "[",-indx.inner,)
            Y.test.inner<-as.matrix(Y.train.outer[indx.inner])
            Y.train.inner<-as.matrix(Y.train.outer[-indx.inner])
            
            train.spls<-mixOmics::block.spls(X = X.train.inner, Y = Y.train.inner, keepX=i.keepX, ncomp = 1, design = design, tol=0.1, all.outputs = F)
            test<-predict(train.spls,X.test.inner)
            MSE.inner.temp[i.fold.inner] <- mean( (matrix(Y.test.inner) - matrix(test$WeightedPredict))^2 ) #MSE
            
          }
          MSE.inner[i.MSE,3+i.Nrepeat.inner]<-mean(MSE.inner.temp)
          
        }
        print(c(i.Nrepeat.outer,i.fold.outer,i.MSE))  #print it
      }
      MSE.inner[,ncol(MSE.inner)]<-rowMeans(MSE.inner[,4:(3+Nrepeat.inner)])
      
      
      
      keepX.outer = list(Thal=MSE.inner[which(MSE.inner[,ncol(MSE.inner)]==min(MSE.inner[,ncol(MSE.inner)]))[1],1],   # select the lowest MSE of the inner folds & repeats
                         Cort=MSE.inner[which(MSE.inner[,ncol(MSE.inner)]==min(MSE.inner[,ncol(MSE.inner)]))[1],2],
                         DTI=MSE.inner[which(MSE.inner[,ncol(MSE.inner)]==min(MSE.inner[,ncol(MSE.inner)]))[1],3])  
      
      train.spls<-mixOmics::block.spls(X = X.train.outer, Y = Y.train.outer, keepX=keepX.outer, ncomp = 1, design = design, tol=0.1)  # test it against the outer testing sample
      test<-predict(train.spls,X.test.outer)
      
      MSE.outer[(i.Nrepeat.outer-1)*Nfold.outer+i.fold.outer,]<-c(unlist(keepX.outer) ,
                                                                  mean((matrix(Y.test.outer) - matrix(test$WeightedPredict))^2) )   # save the outer tested MSE and keepX tested for all the outer folds and repeats 
      
    }
  }
  MSE.repeat <- MSE.outer
  
  MSE.repeat[order(MSE.repeat[,4]),]
  MSE.repeat[which(MSE.repeat[,4]==min(MSE.repeat[,4]))[1],]
  
  
  keepX<-list(Thal=MSE.repeat[which(MSE.repeat[,4]==min(MSE.repeat[,4]))[1],1],   # select the lowest MSE and keepx for all the inner folds and repeats
              Cort=MSE.repeat[which(MSE.repeat[,4]==min(MSE.repeat[,4]))[1],2],
              DTI=MSE.repeat[which(MSE.repeat[,4]==min(MSE.repeat[,4]))[1],3])           
  
  mod.spls<-block.spls(X = X, Y = Y, keepX=keepX, ncomp = 1, design = design)
  #mod.spls.mean<-mean((matrix(Y)-matrix(predict(mod.spls,X)$WeightedPredict))^2)
  mod.spls.mean<-as.vector(cor(matrix(Y),matrix(predict(mod.spls,X)$WeightedPredict))^2) #Rsquared is better
  
  nrep=1000
  spls.perm<-seq(1,nrep)
  pb <- txtProgressBar(min = 0, max = nrep, style = 3)
  for (rep in 1:nrep) {
    setTxtProgressBar(pb, rep)
    Y.perm<-matrix(sample(Y),dimnames=list(rownames(Y),"Y"))
    perm.mod.spls<-block.spls(X = X, Y = Y.perm, keepX=keepX, ncomp = 1, design = design)
    #spls.perm[rep]<-mean((matrix(Y.perm)-matrix(predict(perm.mod.spls,X)$WeightedPredict))^2) #MSE, but gives very asymmetrical values
    spls.perm[rep]<-cor(matrix(Y.perm),matrix(predict(perm.mod.spls,X)$WeightedPredict))^2 #Rsquared is better
  }
  close(pb)
  
  #hist(spls.perm,breaks=nrep/5)
  #abline(v = mod.spls.mean,col="red")
  p_value<-length(which(c(spls.perm,mod.spls.mean)>=mod.spls.mean))/nrep # for MSE is <, for R2 is >
  p_value
  
  spls.predict<-as.vector(cor(matrix(Y),matrix(predict(mod.spls,X)$WeightedPredict))) #r - predictive value is pretty good
  spls.r2<-as.vector(cor(matrix(Y),matrix(predict(mod.spls,X)$WeightedPredict))^2) #R squared
  
  mod<-mod.spls
  
}


########### END OF ANALYSIS ########### 
########### WRITE VALUES ##############

#### write shape ####

lh_mom_coef<-mod$loadings$Thal[1:(mom_points_lh*3)]
rh_mom_coef<-mod$loadings$Thal[(((mom_points_lh*3)+1):(length(mod$loadings$Thal)))]

for (side in side_list) {
  deform_dir<-file.path(paste("shoot_",side,"_",struct,sep=""))
  unlink(deform_dir,recursive = TRUE)
  mom_coef<-matrix(as.vector(get(paste(side,"_mom_coef",sep=""))),ncol=3,byrow=FALSE)  ## *-1 this was there before and I'm not sure why?
  
  #write momenta file
  dir.create(file.path(deform_dir))
  mom_file_new=file.path(deform_dir,"new_momenta.txt")
  write(paste("1 ",get(paste("mom_points_",side,sep=""))," 3\n",sep=""),file=mom_file_new)
  write.table(mom_coef,col.names=FALSE,row.names=FALSE,append=TRUE,file=mom_file_new,quote=FALSE,na="")*3
  
  #copy/add in the other needed files
  file.copy(file.path(paste("output_",side,"_",struct,sep=""),"DeterministicAtlas__EstimatedParameters__Template_mesh.vtk"),
            file.path(paste("shoot_",side,"_",struct,sep="")))
  file.copy(file.path(paste("output_",side,"_",struct,sep=""),"DeterministicAtlas__EstimatedParameters__ControlPoints.txt"),
            file.path(paste("shoot_",side,"_",struct,sep="")))
  sink(file.path(paste("shoot_",side,"_",struct,sep=""),paste("model.xml",sep="")))
  cat("<?xml version=\"1.0\"?>",sep="","\n")
  cat("<model>",sep="","\n")
  cat("<model-type>Shooting</model-type>",sep="","\n")
  cat("<dimension>3</dimension>",sep="","\n")
  cat("<template>",sep="","\n")
  cat("<object id=\"mesh\">",sep="","\n")
  cat("<deformable-object-type>SurfaceMesh</deformable-object-type>",sep="","\n")
  cat("<filename>DeterministicAtlas__EstimatedParameters__Template_mesh.vtk</filename>",sep="","\n")
  cat("</object>",sep="","\n")
  cat("</template>",sep="","\n")
  cat("<initial-control-points>DeterministicAtlas__EstimatedParameters__ControlPoints.txt</initial-control-points>",sep="","\n")
  cat("<initial-momenta>new_momenta.txt</initial-momenta>",sep="","\n")
  cat("<deformation-parameters>",sep="","\n")
  cat("<kernel-width>",kernel,"</kernel-width>",sep="","\n")
  cat("<kernel-type>torch</kernel-type>",sep="","\n")
  cat("</deformation-parameters>",sep="","\n")
  cat("</model>",sep="","\n")
  sink()
  
  system(paste("export LD_LIBRARY_PATH=",QT_PATH," && cd shoot_",side,"_",struct," && ~/anaconda3/envs/deformetrica/bin/deformetrica compute model.xml -v DEBUG",sep=""))
  
  avg<-read_VTK(paste(deform_dir,"/DeterministicAtlas__EstimatedParameters__Template_mesh.vtk",sep=""))
  mesh<-read_VTK(paste(deform_dir,"/output/Shooting__GeodesicFlow__","mesh","__tp_10__age_1.00.vtk",sep=""))
  
  diff<-mesh$vb[1:3,]-avg$vb[1:3,]
  dot<-(diff[1,]*mesh$normals[1,]) + (diff[2,]*mesh$normals[2,]) + (diff[3,]*mesh$normals[3,])
  mag<-sqrt((mesh$normals[1,]*dot)^2+(mesh$normals[2,]*dot)^2+(mesh$normals[3,]*dot)^2)
  mag<-ifelse(dot < 0, mag*-1, mag)
  
  coef_name<-paste("PLS_",formatC(1, width = 2, format = "d", flag = "0"),sep="")
  
  output_file<-paste(base_folder,"/",test_name,"_",side,"_",struct,".vtk",sep="")
  write_VTK(output_file,read_VTK(paste(deform_dir,"/DeterministicAtlas__EstimatedParameters__Template_mesh.vtk",sep="")))
  write_VTK_scalar(output_file,coef_name,mag)
  
  assign(paste0(side,"_avg"),read_VTK(paste(deform_dir,"/DeterministicAtlas__EstimatedParameters__Template_mesh.vtk",sep="")))
  assign(paste0(side,"_mag"),mag)
}

#### write cortical ####


r3dDefaults$windowRect=c(100, 100, 800, 800);
atlas = 'aparc'
template_subject = 'fsaverage'


atlas_region_names = get.atlas.region.names(atlas, template_subjects_dir=subjects_dir, template_subject=template_subject)

region_value_list<-t(mod$loadings$Cort)

lh_region_value_list=c(0,
                       as.numeric(t(mod$loadings$Cort)[1:3]),
                       0,
                       as.numeric(t(mod$loadings$Cort)[4:34])) #need an extra leading zero to cover "unamed" area and at position 5 to cover corpuscallosum for fsbrain
rh_region_value_list=c(0,
                       as.numeric(t(mod$loadings$Cort)[35:37]),
                       0,
                       as.numeric(t(mod$loadings$Cort)[38:68])) #need an extra leading zero to cover "unamed" area and at position 5 to cover corpuscallosum for fsbrain

lh_region_value_list[which(lh_region_value_list==0)]<-NA #better to have zeros as NA, renders them as white later on
rh_region_value_list[which(rh_region_value_list==0)]<-NA #better to have zeros as NA, renders them as white later on
#rh_region_value_list<-rep(NA,36)
colFn_diverging = function(n) { grDevices::hcl.colors(n, palette = "Blue-Red 3"); }
makecmap_options = list('colFn'=colFn_diverging, 'symm'=TRUE);

vis.region.values.on.subject(subjects_dir, template_subject, atlas, lh_region_value_list, rh_region_value_list,
                             surf='pial',bg='sulc_light', views = 't9',
                             makecmap_options=makecmap_options,draw_colorbar=TRUE, rglactions=list('trans_fun'=clip.data));

##### write output ####
if (run_group == 1) {
  sink(paste(base_folder,"/",test_name,"_",struct,"_diagnostics.txt",sep=""))
  print("Tested keepX")
  print(test.keepX)
  print("Used keepX")
  print(keepX)
  print("Explained variance")
  mod.diablo$prop_expl_var
  print("AUC by component")
  auroc(mod.diablo)
  pred.mod.diablo<-predict(mod.diablo,X)
  print("Expected confusion matrix")
  print(ftable(as.vector(Y),as.vector(Y)))
  print("Predicted confusion matrix")
  pred_matrix<-ftable(as.vector(pred.mod.diablo$MajorityVote$max.dist),as.vector(Y))
  print(pred_matrix)
  print("sensitivity")
  sens<-pred_matrix[1,1]/sum(pred_matrix[1:2,1])*100
  print(sens)
  print("specificity")
  spec<-pred_matrix[2,2]/sum(pred_matrix[1:2,2])*100
  print(spec)
  print("DTI loadings")
  print(mod$loadings$DTI)
  print(mod.diablo.perm)
  sink()
}
if (run_group == 0) {
  sink(paste(base_folder,"/",test_name,"_diagnostics.txt",sep=""))
  print("Tested keepX")
  print(test.keepX)
  print("Used keepX")
  print(keepX)
  print("Explained variance - r")
  print(spls.predict) #r - predictive value is pretty good
  print("Explained variance - r^2")
  print(spls.r2)
  print("DTI loadings")
  print(mod$loadings$DTI)
  print(p_value)
  sink()
}
save.image(file=paste(base_folder,"/",test_name,"_mod_diablo.Rdata",sep=""))
