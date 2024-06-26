## From spatialLIBD/global.R
library('ggplot2')
library('cowplot')
library('Polychrome')
library('viridisLite')

## Should be part of the spatialLIBD package... when I make it :P


geom_spatial <-  function(mapping = NULL,
    data = NULL,
    stat = "identity",
    position = "identity",
    na.rm = FALSE,
    show.legend = NA,
    inherit.aes = FALSE,
    ...) {
    GeomCustom <- ggproto(
        "GeomCustom",
        Geom,
        setup_data = function(self, data, params) {
            data <- ggproto_parent(Geom, self)$setup_data(data, params)
            data
        },
        
        draw_group = function(data, panel_scales, coord) {
            vp <- grid::viewport(x = data$x, y = data$y)
            g <- grid::editGrob(data$grob[[1]], vp = vp)
            ggplot2:::ggname("geom_spatial", g)
        },
        
        required_aes = c("grob", "x", "y")
        
    )
    
    layer(
        geom = GeomCustom,
        mapping = mapping,
        data = data,
        stat = stat,
        position = position,
        show.legend = show.legend,
        inherit.aes = inherit.aes,
        params = list(na.rm = na.rm, ...)
    )
}


get_colors <- function(colors = NULL, clusters) {
    n_clus <- length(unique(clusters))
    
    if (is.null(colors) | n_clus > length(colors)) {
        ## Original ones
        # colors <- c("#b2df8a","#e41a1c","#377eb8","#4daf4a","#ff7f00","gold",
        # "#a65628", "#999999", "black", "grey", "white", "purple")
        
        ## From https://medialab.github.io/iwanthue/
        ## which I found the link to from
        ## https://stackoverflow.com/questions/15282580/how-to-generate-a-number-of-most-distinctive-colors-in-r
        ## Used the colorblind friendly and default palette
        
        ## From https://developer.r-project.org/Blog/public/2019/11/21/a-new-palette-for-r/index.html
        ## they point to https://cran.r-project.org/web/packages/Polychrome/vignettes/polychrome.html
        
        
        colors <-
            if (n_clus > 12)
                Polychrome::palette36.colors(n_clus)
        else
            c(
                "#b2df8a",
                "#e41a1c",
                "#377eb8",
                "#4daf4a",
                "#ff7f00",
                "gold",
                "#a65628",
                "#999999",
                "black",
                "grey",
                "white",
                "purple"
            )
        names(colors) <- seq_len(length(colors))
        
    }
    return(colors)
}



sce_image_clus <- function(sce,
    sampleid,
    clustervar,
    colors = c(
        "#b2df8a",
        "#e41a1c",
        "#377eb8",
        "#4daf4a",
        "#ff7f00",
        "gold",
        "#a65628",
        "#999999",
        "black",
        "grey",
        "white",
        "purple"
    ),
    spatial = TRUE,
    ...) {
    sce_sub <- sce[, sce$sample_name == sampleid]
    d <- as.data.frame(colData(sce_sub))
    sce_image_clus_p(
        sce = sce_sub,
        d = d,
        clustervar = clustervar,
        sampleid = sampleid,
        spatial = spatial,
        title = paste0(sampleid, ...),
        colors = get_colors(colors, d[, clustervar])
    )
}

sce_image_clus_p <-
    function(sce,
        d,
        clustervar,
        sampleid,
        colors,
        spatial,
        title) {
        p <- ggplot(d,
            aes(
                x = imagecol,
                y = imagerow,
                fill = factor(!!sym(clustervar)),
                
            ))
        if (spatial) {
            p <-
                p + geom_spatial(
                    data = subset(metadata(sce)$image, sample == sampleid),
                    aes(grob = grob),
                    x = 0.5,
                    y = 0.5
                )
        }
        p <- p +
            geom_point(shape = 21,
                size = 1.25,
                stroke = 0.25) +
            coord_cartesian(expand = FALSE) +
            scale_fill_manual(values = colors) +
            xlim(0, max(sce$width)) +
            ylim(max(sce$height), 0) +
            xlab("") + ylab("") +
            labs(fill = clustervar) +
            guides(fill = guide_legend(override.aes = list(size = 3))) +
            ggtitle(title) +
            theme_set(theme_bw(base_size = 10)) +
            theme(
                panel.grid.major = element_blank(),
                panel.grid.minor = element_blank(),
                panel.background = element_blank(),
                axis.line = element_line(colour = "black"),
                axis.text = element_blank(),
                axis.ticks = element_blank()
            )
        return(p)
    }

sce_image_clus_gene <-
    function(sce,
        sampleid,
        geneid = genes[17856],
        spatial = TRUE,
        assayname = 'logcounts',
        minCount = 0,
        ...) {
        sce_sub <- sce[, sce$sample_name == sampleid]
        d <- as.data.frame(colData(sce_sub))
        
        if(geneid %in% colnames(colData(sce_sub))) {
            d$COUNT <- colData(sce_sub)[[geneid]]
        } else {
            d$COUNT <- assays(sce_sub)[[assayname]][which(genes == geneid),]
        }
        d$COUNT[d$COUNT <= minCount] <- NA
        sce_image_clus_gene_p(
            sce = sce_sub,
            d = d,
            sampleid = sampleid,
            spatial = spatial,
            title = paste(
                sampleid,
                geneid,
                if(!geneid %in% colnames(colData(sce_sub))) assayname,
                paste0('min Count: >', minCount),
                ...,
                sep = " - "
            )
        )
    }


sce_image_clus_gene_p <-
    function(sce, d, sampleid, spatial, title) {
        p <-
            ggplot(d,
                aes(
                    x = imagecol,
                    y = imagerow,
                    fill = COUNT,
                    color = COUNT,
                    key =  key
                ))
        
        if (spatial) {
            p <-
                p + geom_spatial(
                    data = subset(metadata(sce)$image, sample == sampleid),
                    aes(grob = grob),
                    x = 0.5,
                    y = 0.5
                )
        }
        p <- p +
            geom_point(shape = 21,
                size = 1.25,
                stroke = 0.25) +
            coord_cartesian(expand = FALSE) +
            scale_fill_gradientn(
                colors = viridis(21)
            ) +
            scale_color_gradientn(
                colors = viridis(21)
            ) +
            xlim(0, max(sce$width)) +
            ylim(max(sce$height), 0) +
            xlab("") + ylab("") +
            labs(fill = "COUNT") +
            ggtitle(title) +
            theme_set(theme_bw(base_size = 10)) +
            theme(
                panel.grid.major = element_blank(),
                panel.grid.minor = element_blank(),
                panel.background = element_blank(),
                axis.line = element_line(colour = "black"),
                axis.text = element_blank(),
                axis.ticks = element_blank()
            )
        return(p)
    }

# sce_image_clus_gene(sce, '151507')
# sce_image_clus_gene(sce, '151507', minCount = 3)
# sce_image_clus_gene(sce, '151507', minCount = 3, assayname = 'counts')


sort_clusters <- function(clusters, map_subset = NULL) {
    if (is.null(map_subset)) {
        map_subset <- rep(TRUE, length(clusters))
    }
    map <-
        rank(length(clusters[map_subset]) - table(clusters[map_subset]), ties.method = 'first')
    res <- map[clusters]
    factor(res)
}

sce_image_grid <-
    function(sce,
        clusters,
        pdf_file,
        sort_clust = TRUE,
        colors = NULL,
        return_plots = FALSE,
        spatial = TRUE,
        ...) {
        colors <- get_colors(colors, clusters)
        sce$Clus <- if (sort_clust)
            sort_clusters(clusters)
        else
            clusters
        plots <-
            lapply(unique(sce$sample_name), function(sampleid) {
                sce_image_clus(sce,
                    sampleid,
                    'Clus',
                    colors = colors,
                    spatial = spatial,
                    ...)
            })
        names(plots) <- unique(sce$sample_name)
        if (!return_plots) {
            pdf(pdf_file, height = 24, width = 36)
            print(cowplot::plot_grid(plotlist = plots))
            dev.off()
            return(pdf_file)
        }
        else {
            return(plots)
        }
    }

sce_image_grid_gene <-
    function(sce,
        geneid = genes[17856],
        pdf_file,
        assayname = 'logcounts',
        minCount = 0,
        return_plots = FALSE,
        spatial = TRUE,
        ...) {
        plots <- lapply(unique(sce$sample_name), function(sampleid) {
            sce_image_clus_gene(sce, sampleid, geneid, spatial, assayname, minCount, ...)
        })
        names(plots) <- unique(sce$sample_name)
        if (!return_plots) {
            pdf(pdf_file, height = 24, width = 36)
            print(cowplot::plot_grid(plotlist = plots))
            dev.off()
            return(pdf_file)
        }
        else {
            return(plots)
        }
    }

# plots <- sce_image_grid_gene(sce, minCount = 3, return_plots = TRUE)
# print(cowplot::plot_grid(plotlist = plots))

# sce_image_grid_comp <-
#     function(sce,
#         clus1,
#         clus2,
#         pdf_file,
#         map_subset = NULL,
#         colors = NULL,
#         return_plots = FALSE,
#         ...) {
#         clus1 <- sort_clusters(clus1, map_subset)
#         clus2 <- sort_clusters(clus2, map_subset)
#         if (is.null(colors)) {
#             colors <- c('FALSE' = 'red', 'TRUE' = 'grey80')
#         }
#         clus_agree <-
#             factor(clus1 == clus2, levels = c('FALSE', 'TRUE'))
#         sce_image_grid(
#             sce,
#             clus_agree,
#             pdf_file,
#             sort_clust = FALSE,
#             colors = colors,
#             return_plots = return_plots,
#             ...
#         )
#     }


sce_image_grid_by_clus <-
    function(sce,
        clusters,
        pdf_file,
        colors = NULL,
        spatial = TRUE,
        ...) {
        if (is.null(colors)) {
            colors <- c('FALSE' = 'transparent', 'TRUE' = 'red')
        }
        clusters_uni <- sort(unique(clusters))
        pdf(pdf_file, height = 24, width = 36)
        lapply(clusters_uni, function(clus) {
            curr_clus <- factor(clusters == clus, levels = c('FALSE', 'TRUE'))
            plots <-
                sce_image_grid(
                    sce,
                    curr_clus,
                    sort_clust = FALSE,
                    colors = colors,
                    return_plots = TRUE,
                    spatial = spatial,
                    ... = paste(..., '- cluster', clus)
                )
            print(cowplot::plot_grid(plotlist = plots))
            return(clus)
        })
        dev.off()
        return(pdf_file)
    }
