library(ANTsR)
library(Morpho)
library(Rvcg)
library(rgl)

source("VTK_ops.R")

# assumes there's a CSV file with an ID variable which corresponds to the subject ID from freesurfer
csvfile<-read.csv(file.path(base_folder,"data.csv")) 

#for example
csvfile<-data.frame(matrix(ncol=1,nrow=1, c("bert"), dimnames=list(NULL, c("ID"))))

# folder which corresponds to where the processing will occur
base_folder<-"~/data/deformetrica" 

# change this to the freesurfer subject directory if wishing to use freesurfer aseg instead of manual segmentations
# if using manual segmentations, they should be in a folder named "source/" in the format "lh_thal_ID" such as "lh_thal_CON001" or "lh_thal_bert", etc.
subjects_dir<-"~/freesurfer/subjects/"

# need to change this to point to deformetrica Qt5 install location
QT_PATH<-"~/USER/anaconda3/envs/deformetrica/lib/python3.8/site-packages/PyQt5/Qt5/lib"  

# these are deformetrica control parameters
kernel=3
noise=3
kernel_deform=5


# ###### optionally extract from freesurfer aseg if no ####
# if (file.exists(file.path(base_folder,"source")) == FALSE) {
#   library(freesurferformats)
#   dir.create(file.path(base_folder,"source"))
#   for(x in csvfile$ID){
#     side="lh"
#     system(paste("mri_binarize",
#                  "--i",file.path("subjects_dir",x,"mri","aparc+aseg.mgz"),
#                  "--o",file.path(base_folder,"source",paste(side,"_",struct,"_",x,".nii.gz",sep="")),
#                  "--match",10,"--binval 1"))
#     side="rh"
#     system(paste("mri_binarize",
#                  "--i",file.path("subjects_dir",x,"mri","aparc+aseg.mgz"),
#                  "--o",file.path(base_folder,"source",paste(side,"_",struct,"_",x,".nii.gz",sep="")),
#                  "--match",49,"--binval 1"))
#   }
# }


###### make template for each structure for each side ####
if (file.exists(file.path(base_folder,"template")) == FALSE) {
  dir.create(file.path(base_folder,"template"))
  for (side in side_list){
    loop_reg=3
    for (y in 1:loop_reg) {
      cat("loop",y,"\n")
      if (file.exists(file.path(base_folder,"template",paste(side,"_",struct,"_average_rigid.nii.gz",sep=""))) == TRUE) {  # if there's a template, use it
        fi<-antsImageRead(file.path(base_folder,"template",paste(side,"_",struct,"_average_rigid.nii.gz",sep="")))
      } else { 
        fi<-antsImageRead(file.path(base_folder,"source",paste(side,"_",struct,"_",csvfile$ID[1],".nii.gz",sep=""))) 
      }
      for(x in csvfile$ID){
        cat(x," ")
        mi<-antsImageRead(file.path(base_folder,"source",paste(side,"_",struct,"_",x,".nii.gz",sep="")))
        mytx<-antsRegistration(fixed=fi , moving=mi, typeofTransform = c('QuickRigid'))
        antsImageWrite(antsApplyTransforms( fixed=fi, moving=mi, transformlist=mytx$fwdtransforms,interpolator="nearestNeighbor" ),
                       file.path(base_folder,"template",paste(side,"_",struct,"_",x,".nii.gz",sep="")))
      }
      file.remove(file.path(base_folder,"template",paste(side,"_",struct,"_average_rigid.nii.gz",sep=""))) #delete the avg file because it's going to be re-written
      cat("\n","make average","\n")
      avg<-antsAverageImages(list.files(file.path(base_folder,"template"),pattern = paste(side,"_",struct,sep=""), full.names = TRUE))
      avg_thresh<-thresholdImage(avg,0.2,1)  # inverse relationship - thresh 0.2 means more included
      antsImageWrite(avg_thresh, file.path(base_folder,"template",paste(side,"_",struct,"_average_rigid.nii.gz",sep="")))
      cat("\n")
    }
  }
}

#### make VTK ####
if (file.exists(file.path(base_folder,"mesh")) == FALSE) {
  dir.create(file.path(base_folder,"mesh"))
  for (side in side_list){
    avg_thresh<-antsImageRead(file.path(base_folder,"template",paste(side,"_",struct,"_average_rigid.nii.gz",sep="")))
    avg_thresh_mesh<-vcgSmooth(vcgIsosurface(as.array(avg_thresh),threshold=1),"HC",iteration = 10)
    avg_thresh_mesh<-vcgClean(vcgQEdecim(avg_thresh_mesh,tarface=2000),sel=c(1:7,1))
    write_VTK(file.path(base_folder,"mesh",paste(side,"_",struct,"_average_rigid.vtk",sep="")), avg_thresh_mesh)
    for(x in csvfile$ID){
      cat(x," ")
      mi<-antsImageRead(file.path(base_folder,"template",paste(side,"_",struct,"_",x,".nii.gz",sep="")))
      avg_thresh_mesh<-vcgSmooth(vcgIsosurface(as.array(mi),threshold=1),"HC",iteration = 5)
      avg_thresh_mesh<-vcgClean(vcgQEdecim(avg_thresh_mesh,tarface=2000),sel = c(1:7,1))
      write_VTK(file.path(base_folder,"mesh",paste(side,"_",struct,"_",x,".vtk",sep="")), avg_thresh_mesh)
    }
  }
}

#### run deformetrica atlas/registration  ####
if (run_deform == 1) {
  for (side in side_list) {
    
    unlink(file.path(base_folder,paste0("output_",side,"_",struct)),recursive = TRUE)
    
    sink(file.path(base_folder,paste(side,"_",struct,"_data_set.xml",sep="")))
    cat("<?xml version=\"1.0\" ?>",sep="","\n")
    cat("<data-set>",sep="","\n",append=TRUE)
    for (subject in csvfile$ID) {
      cat("<subject id=\"",subject,"\">","\n",sep="")
      cat("<visit id=\"baseline\">",sep="","\n")
      cat("<filename object_id=\"mesh\">","mesh/",side,"_",struct,"_",subject,".vtk","</filename>",sep="","\n")
      cat("</visit>",sep="","\n")
      cat("</subject>",sep="","\n")
    }
    cat("</data-set>",sep="","\n")
    sink()
    
    sink(file.path(base_folder,paste(side,"_",struct,"_model.xml",sep="")))
    cat("<?xml version=\"1.0\"?>",sep="","\n")
    cat("<model>",sep="","\n")
    cat("<model-type>DeterministicAtlas</model-type>",sep="","\n")
    cat("<dimension>3</dimension>",sep="","\n")
    cat("<template>",sep="","\n")
    cat("<object id=\"mesh\">",sep="","\n")
    cat("<deformable-object-type>SurfaceMesh</deformable-object-type>",sep="","\n")
    cat("<attachment-type>Current</attachment-type>",sep="","\n")
    cat("<noise-std>",noise,"</noise-std>",sep="","\n")
    cat("<kernel-type>torch</kernel-type>",sep="","\n")
    cat("<kernel-width>",kernel,"</kernel-width>",sep="","\n")
    cat("<filename>mesh/",side,"_",struct,"_average_rigid.vtk</filename>",sep="","\n")
    cat("</object>",sep="","\n")
    cat("</template>",sep="","\n")
    cat("<deformation-parameters>",sep="","\n")
    cat("<kernel-width>",kernel_deform,"</kernel-width>",sep="","\n")
    cat("<kernel-type>torch</kernel-type>",sep="","\n")
    cat("</deformation-parameters>",sep="","\n")
    cat("</model>",sep="","\n")
    sink()
    
    system(paste0("export LD_LIBRARY_PATH=",QT_PATH," && ~/anaconda3/envs/deformetrica/bin/deformetrica estimate ",side,"_",struct,"_model.xml ",side,"_",struct,"_data_set.xml -p optimization_parameters.xml -v DEBUG --output=output_",side,"_",struct))
  }
}
