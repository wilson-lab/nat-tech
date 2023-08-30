# Load libraries
library(fafbseg)
library(nat.templatebrains)
library(nat.jrcbrains)
library(elmr)

# Get registrations if you have not done so already
# nat.jrcbrains::download_saalfeldlab_registrations()
nat.jrcbrains::register_saalfeldlab_registrations()

# Get directories
obj.dir <- file.path("data","obj","flywire")
dir.create(obj.dir, showWarnings = FALSE, recursive = TRUE)
swc.dir <- file.path("data","swc", "flywire")
dir.create(swc.dir, showWarnings = FALSE, recursive = TRUE)

# Get meta data, if this does not work for you obtain data here: https://github.com/flyconnectome/flywire_annotations
#ft <- fafbseg::flytable_query("select _id, root_id, root_630, supervoxel_id, proofread, status, pos_x, pos_y, pos_z, nucleus_id, soma_x, soma_y, soma_z, side, ito_lee_hemilineage, hartenstein_hemilineage, top_nt, flow, super_class, cell_class, cell_type, hemibrain_type, root_duplicated from info")
#ft <- as.data.frame(ft)
ft <-  read.table("data/Supplemental_file1_annotations.tsv", header=TRUE, sep="\t", quote = "")

# Select AOTU19
fw.select <- subset(ft, hemibrain_type=="AOTU019")
fw.ids <- fw.select$root_id

# Read mesh data
# fw.meshes <- fafbseg::read_cloudvolume_meshes(fw.ids) # for some reason does not retain normals
fafbseg::download_neuron_obj(fw.ids, save.obj = obj.dir)
objs <- list.files(obj.dir, full.names = TRUE)
fw.meshes <- nat::neuronlist()
for(obj in objs){
  message("Working on: ", obj)
  mesh <- Morpho::obj2mesh(obj)
  mesh <- nat::as.neuronlist(list(mesh))
  names(mesh) <- gsub("\\.obj","",basename(obj))
  fw.meshes <- c(fw.meshes, mesh)
}

# Transform meshes into um
fw.meshes.t <- xform_brain(fw.meshes, reference = "JRC2018F", sample = "FlyWire")

# Save transformed obj file
for(m in names(fw.meshes.t)){
  mesh <- fw.meshes.t[[m]]
  filename <- file.path(obj.dir,m)
  Rvcg::vcgObjWrite(mesh = mesh, 
                    filename = filename, 
                    writeNormals = TRUE)
}

# Skeletonise with radius information.
## see ?skeletor for arguments
## At this point you can try using instead: https://neuromorpho.org/xyz2swc/ui/
fw.meshes.skels <- nat::neuronlist()
for(obj in objs){
  message("Working on: ", obj)
  fw.meshes.skel<- skeletor(
    segments = NULL,
    obj = obj,
    mesh3d = FALSE,
    save.obj = NULL,
    cloudvolume.url = getOption("fafbseg.cloudvolume.url"),
    operator = c("umbrella", "contangent"),
    clean = FALSE,
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
  fw.meshes.skels <- c(fw.meshes.skels,fw.meshes.skel)
}

# save swc
write.neurons(nl = fw.meshes.skels,
              dir = swc.dir,
              files = names(fw.meshes.skels),
              Force = TRUE,
              format = "swc") 

# check it looks okay
nopen3d(userMatrix = structure(c(0.998503506183624, 0.028934620320797, 
                                 -0.0464092344045639, 0, 0.0174863673746586, -0.972944557666779, 
                                 -0.230376064777374, 0, -0.0518193989992142, 0.229219749569893, 
                                 -0.971994340419769, 0, 0, 0, 0, 1), dim = c(4L, 4L)), zoom = 0.281241029500961, 
        windowRect = c(38L, 47L, 1182L, 921L))
plot3d(JRC2018F, alpha = 0.1)
plot3d(fw.meshes.skels, lwd = 0.1)
for(n in 1:length(fw.meshes.skels)){
  swc <- fw.meshes.skels[[n]]$d
  p <- nat::xyzmatrix(swc)
  spheres3d(p, radius = swc$W, color=rainbow(2)[n])
}
plot3d(fw.meshes.t, alpha = 0.3)

