# Load libraries
library(neuprintr)
library(hemibrainr)
library(nat.templatebrains)
library(nat.jrcbrains)
library(elmr)
source("R/startup/functions.R")

# Get registrations if you have not done so already
# nat.jrcbrains::download_saalfeldlab_registrations()
nat.jrcbrains::register_saalfeldlab_registrations()

# Get directories
obj.dir <- file.path("data","obj","hemibrain")
dir.create(obj.dir, showWarnings = FALSE, recursive = TRUE)
swc.dir <- file.path("data","swc","hemibrain")
dir.create(swc.dir, showWarnings = FALSE, recursive = TRUE)

# Select AOTU19
hb <- neuprintr::neuprint_search(search = "AOTU019", field = "type")
hb.ids <- as.character(hb$bodyid)


# Read mesh data
# hb.meshes <- hemibrainr::hemibrain_neuron_meshes(hb.ids,cloudvolume.url=cloudvolume.url) # for some reason does not retain normals
download_hemibrain_obj(hb.ids, save.obj = obj.dir)
objs <- list.files(obj.dir, full.names = TRUE)
hb.meshes <- nat::neuronlist()
for(obj in objs){
  message("Working on: ", obj)
  mesh <- Morpho::obj2mesh(obj)
  mesh <- nat::as.neuronlist(list(mesh))
  names(mesh) <- gsub("\\.obj","",basename(obj))
  hb.meshes <- c(hb.meshes, mesh)
}

# Transform meshes into um
hb.meshes.t <- xform_brain(hb.meshes, reference = "JRC2018F", sample = "JRCFIB2018Fraw")

# Save transformed obj file
for(m in names(hb.meshes.t)){
  mesh <- hb.meshes.t[[m]]
  filename <- file.path(obj.dir,m)
  Rvcg::vcgObjWrite(mesh = mesh, 
                    filename = filename, 
                    writeNormals = TRUE)
}

# Skeletonise with radius information.
## see ?skeletor for arguments
## At this point you can try using instead: https://neuromorpho.org/xyz2swc/ui/
## May need to tweak these parameters
hb.meshes.skels <- nat::neuronlist()
for(obj in objs){
  message("Working on: ", obj)
  hb.meshes.skel<- skeletor(
    segments = NULL,
    obj = obj,
    mesh3d = FALSE,
    save.obj = NULL,
    cloudvolume.url = getOption("fafbseg.cloudvolume.url"),
    operator = c("umbrella", "contangent"),
    clean = TRUE,
    remove_disconnected = 10,
    theta = 0.01,
    radius = TRUE,
    ratio = 0.1,
    SL = 10,
    WH0 = 2,
    iter_lim = 4,
    epsilon = 0.05,
    precision = 1e-06,
    validate = TRUE,
    method.radii = "knn", # c("knn", "ray"),
    method = c("wavefront"), # "vertex_clusters", "edge_collapse", "teasar", "tangent_ball"),
    heal = TRUE,
    heal.k = 10L,
    heal.threshold = Inf,
    reroot = TRUE,
    k.soma.search = 10,
    radius.soma.search = 2500,
    reroot_method = "density",#c("direction", "density"),
    brain = NULL,
    n = 5,
    n_rays = 20,
    projection = "sphere", #c("sphere", "tangents"),
    fallback = "knn",
    waves = 2,
    step_size = 1,
    sampling_dist = 500,
    cluster_pos = "median", #c("median", "center"),
    shape_weight = 1,
    sample_weight = 0.1,
    inv_dist = 100,
    cpu = Inf,
    elapsed = Inf,
  )
  hb.meshes.skels <- c(hb.meshes.skels,hb.meshes.skel)
}

# save swc
write.neurons(nl = hb.meshes.skels,
              dir = swc.dir,
              files = names(hb.meshes.skels),
              Force = TRUE,
              format = "swc") 

# check it looks okay
nopen3d()
plot3d(JRC2018F, alpha = 0.1)
plot3d(hb.meshes.skels, lwd = 0.1)
for(n in 1:length(hb.meshes.skels)){
  swc <- hb.meshes.skels[[n]]$d
  p <- nat::xyzmatrix(swc)
  spheres3d(p, radius = swc$W, color=rainbow(2)[n])
}
plot3d(hb.meshes.t, alpha = 0.3)

