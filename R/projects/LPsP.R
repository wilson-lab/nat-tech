# Load libraries
source("R/startup/packages.R")
source("R/startup/functions.R")

###########
# flywire #
###########

# Get data
ft <- fafbseg::flytable_query("select _id, root_id, root_630, supervoxel_id, proofread, status, pos_x, pos_y, pos_z, nucleus_id, side, ito_lee_hemilineage, hartenstein_hemilineage, top_nt, top_nt_conf, flow, super_class, cell_class, cell_sub_class, cell_type, hemibrain_type, root_duplicated from info")
fw.lpsp <- subset(ft, grepl("LPsP",hemibrain_type))
fw.lpsp.ids <- unique(fw.lpsp$root_id)

#  Get plots
p.right <- flywire_ntplot(fw.lpsp.ids[1])
p.right <- p.right +
  ggplot2::theme_minimal()+
  ggplot2::theme(
    plot.background = ggplot2::element_rect(fill = "black"),
    panel.background = ggplot2::element_rect(fill = "darkgrey"),
    axis.text = ggplot2::element_text(color = "white"),
    axis.title = ggplot2::element_text(color = "white"),
    plot.title = ggplot2::element_text(color = "white"),
    legend.background = ggplot2::element_rect(fill = "black"),
    legend.text = ggplot2::element_text(color = "white"),
    legend.title = ggplot2::element_text(color = "white")
  )
p.left <- flywire_ntplot(fw.lpsp.ids[2]) 
p.left <- p.left +
  ggplot2::theme_minimal()+
  ggplot2::theme(
    plot.background = ggplot2::element_rect(fill = "black"),
    panel.background = ggplot2::element_rect(fill = "darkgrey"),
    axis.text = ggplot2::element_text(color = "white"),
    axis.title = ggplot2::element_text(color = "white"),
    plot.title = ggplot2::element_text(color = "white"),
    legend.background = ggplot2::element_rect(fill = "black"),
    legend.text = ggplot2::element_text(color = "white"),
    legend.title = ggplot2::element_text(color = "white")
  )

# Save
ggplot2::ggsave(file="/Users/abates/projects/wilson-lab/nat-tech/images/nt_predictions/flywire_LPsP_right.png",
       plot=p.right,
       units="px",
       width = 3000,
       height = 3000)
ggplot2::ggsave(file="/Users/abates/projects/wilson-lab/nat-tech/images/nt_predictions/flywire_LPsP_left.png",
       plot=p.left,
       units="px",
       width = 3000,
       height = 3000)

#############
# hemibrain #
#############

# Select LPsP
hb <- neuprintr::neuprint_search(search = "LPsP", field = "type")
hb.ids <- as.character(hb$bodyid)

# Get synapses
# remote <- "/Users/abates/Dropbox (HMS)/hemibrainr"
# remote <- "/Users/abates/Library/CloudStorage/GoogleDrive-ab2248@cam.ac.uk/Shared drives/hemibrainr"
# options(Gdrive_hemibrain_data = remote)
# hb.syns <- hemibrain_synapses()
hb.syns <- readr::read_csv("/Users/GD/LMBD/Papers/synister/synister_analysis/data/synapses/hemibrain_nt_presynapses_2023-01-20.csv",
                           col_types = hemibrainr:::sql_col_types)
hb.syns.lpsp <- hb.syns %>%
  dplyr::filter(bodyid %in% hb.ids) %>%
  dplyr::mutate(syn_top_p = as.numeric(syn_top_p))
hb.syns.lpsp.right <- subset(hb.syns.lpsp, bodyid == hb.ids[2])
hb.syns.lpsp.left <- subset(hb.syns.lpsp, bodyid == hb.ids[1])

# Plot right
nts <- c("gaba", "acetylcholine", "glutamate", "octopamine", "serotonin", "dopamine")
x = dplyr::filter(hb.syns.lpsp.right, .data$confidence >= 0.5 & 
                    .data$syn_top_nt %in% nts) %>%
  dplyr::distinct(connector_id, .keep_all = TRUE)
ntcols = c(gaba = "#E6A749", acetylcholine = "#4B506B", glutamate = "#70B657", 
           octopamine = "#7A4F98", serotonin = "#93A3CF", dopamine = "#CF6F6C", 
           neither = "grey70")[nts]
p.hb.right <- ggplot2::ggplot(x, ggplot2::aes(.data$syn_top_p, fill = .data$syn_top_nt)) + 
  ggplot2::geom_histogram() + ggplot2::scale_x_continuous(name = "probability") + 
  ggplot2::scale_fill_manual("nt", values = ntcols, breaks = names(ntcols)) +
  ggplot2::theme_minimal()+
  ggplot2::theme(
    plot.background = ggplot2::element_rect(fill = "black"),
    panel.background = ggplot2::element_rect(fill = "darkgrey"),
    axis.text = ggplot2::element_text(color = "white"),
    axis.title = ggplot2::element_text(color = "white"),
    plot.title = ggplot2::element_text(color = "white"),
    legend.background = ggplot2::element_rect(fill = "black"),
    legend.text = ggplot2::element_text(color = "white"),
    legend.title = ggplot2::element_text(color = "white")
  )

# Plot left
x = dplyr::filter(hb.syns.lpsp.left, .data$confidence >= 0.5 & 
                    .data$syn_top_nt %in% nts) %>%
  dplyr::distinct(connector_id, .keep_all = TRUE)
p.hb.left <- ggplot2::ggplot(x, ggplot2::aes(.data$syn_top_p, fill = .data$syn_top_nt)) + 
  ggplot2::geom_histogram() + ggplot2::scale_x_continuous(name = "probability") + 
  ggplot2::scale_fill_manual("nt", values = ntcols, breaks = names(ntcols)) +
  ggplot2::theme_minimal()+
  ggplot2::theme(
    plot.background = ggplot2::element_rect(fill = "black"),
    panel.background = ggplot2::element_rect(fill = "darkgrey"),
    axis.text = ggplot2::element_text(color = "white"),
    axis.title = ggplot2::element_text(color = "white"),
    plot.title = ggplot2::element_text(color = "white"),
    legend.background = ggplot2::element_rect(fill = "black"),
    legend.text = ggplot2::element_text(color = "white"),
    legend.title = ggplot2::element_text(color = "white")
  )

# Save
ggplot2::ggsave(file="/Users/abates/projects/wilson-lab/nat-tech/images/nt_predictions/hemibrain_LPsP_right.png",
                plot=p.hb.right,
                units="px",
                width = 3000,
                height = 3000)
ggplot2::ggsave(file="/Users/abates/projects/wilson-lab/nat-tech/images/nt_predictions/hemibain_LPsP_left.png",
                plot=p.hb.left,
                units="px",
                width = 3000,
                height = 3000)


