function(input, output, session) {
  
  observeEvent(input$slide, {
    updateSelectInput(session, "tumor_only_ROI",
                      choices = gsub(" Unspecified", "",
                                     unique(paste(ROI[ROI$SlideName == input$slide, "ROINumber"],
                                                  ROI[ROI$SlideName == input$slide, "ROILabel"]))))
    updateNumericInput(session, "clusters",
                       value = 2,
                       max = min(c(15, sum(ROI$SlideName == input$slide))))
    hide("outline_deconv")
    hide("outline_cluster")
    hide("heatmap")
    hide("table")
  })
  
  # SpatialDecon calculations
  spatialdecon_results <- eventReactive(input$run, {
    specific_expression = expression[expression$SlideName == input$slide,]
    specific_expression$ROIName <- paste(specific_expression$ROINumber, specific_expression$ROILabel, sep = "_")
    specific_expression <- t(specific_expression)
    colnames(specific_expression)  <- as.vector(specific_expression["ROIName", ])
    specific_expression.names.remove <- c("SlideName", "ROINumber", "ROILabel", "ROIName")
    specific_expression <- specific_expression[!(row.names(specific_expression) %in% specific_expression.names.remove), ]
    specific_expression <- type.convert(specific_expression)
    background_scan <- derive_GeoMx_background(
      norm = specific_expression,
      probepool = rep(1, nrow(specific_expression)),
      negnames = "NegProbe"
    )
    list_nuclei <- ROI[ROI$SlideName == input$slide, ]$ROINucleiCount
    if (input$immune == TRUE) {
      results <- spatialdecon(
        norm = specific_expression,
        bg = background_scan,
        X = safeTME,
        cell_counts = list_nuclei,
        cellmerges = safeTME.matches
      )
    } else {
      list_pure <- unlist(lapply(colnames(specific_expression),
                                 function(x) gsub("_", " ", gsub("_Unspecified", "", x)) %in% input$tumor_only_ROI))
      results <- spatialdecon(
        norm = specific_expression,
        bg = background_scan,
        X = safeTME,
        is_pure_tumor = list_pure,
        cell_counts = list_nuclei,
        cellmerges = safeTME.matches
      )
    }
    results$SlideName <- input$slide
    results
  })
  
  observeEvent(spatialdecon_results(), {
    show("outline_deconv")
    show("outline_cluster")
    show("heatmap")
    show("table")
    })
 
  # Florets plots:
  output$fluorescent_deconv <- renderPlot(
    {slide_fluo <- load.image(paste("data/jpg_fluorescent/", input$slide, ".jpg", sep = ""))
    par(mar=c(0,0,0.1,0.1))
    plot(slide_fluo, axes = FALSE, frame=FALSE, xaxs="i", yaxs="i", rescale=FALSE)
    })
  
  outline_deconv_plot <- eventReactive(c(spatialdecon_results(), input$slide), {
    slide_outline <- load.image(paste("data/jpg_outline/", isolate(input$slide), "_outline.jpg", sep = ""))
    par(mar=c(0.1,0,0.1,0))
    plot(slide_outline, axes = FALSE, frame=FALSE, xaxs="i", yaxs="i", rescale=FALSE)
    if (spatialdecon_results()$SlideName == input$slide) {
      ROI_specific <- ROI[ROI$SlideName == isolate(input$slide), ]
      ROI_specific_stroma <- ROI_specific$ROILabel %in% c("Stroma", "Unspecified")
      ROICoordinateX_stroma <- ROI_specific[ROI_specific_stroma, ]$ROICoordinateX
      ROICoordinateY_stroma <- ROI_specific[ROI_specific_stroma, ]$ROICoordinateY
      florets(x = ROICoordinateX_stroma, y = ROICoordinateY_stroma, b = spatialdecon_results()$beta[, ROI_specific_stroma], cex = 2.5, add = TRUE,
              col = c("#4cac48", "#b9e398", "#3a87bd", "#b6d9e9", "#fbd9eb", "#bae27e", "#fbb86d", "#97c1dc", "#fb958a", "#c4c3df", "#fdfebd", "#90d8cc", "#fbc478", "#e73d3f"))
    }
  })
  
  output$outline_deconv <- renderPlot({
    outline_deconv_plot()
    })
  
  # Cluster plots
  output$fluorescent_cluster <- renderPlot(
    {slide_fluo <- load.image(paste("data/jpg_fluorescent/", input$slide, ".jpg", sep = ""))
    par(mar=c(0,0,0.1,0.1))
    plot(slide_fluo, axes = FALSE, frame=FALSE, xaxs="i", yaxs="i", rescale=FALSE)})
  
  outline_cluster_plot <- eventReactive(c(spatialdecon_results(), input$slide, input$clusters), {
    slide_outline <- load.image(paste("data/jpg_outline/", input$slide, "_outline.jpg", sep = ""))
    par(mar=c(0.1,0,0.1,0))
    plot(slide_outline, axes = FALSE, frame=FALSE, xaxs="i", yaxs="i", rescale=FALSE)
    if (spatialdecon_results()$SlideName == input$slide) {
      cluster_celltype <- as.data.frame(cutree(hclust(dist(t(spatialdecon_results()$beta)),
                                                      method="complete"), input$clusters))
      colnames(cluster_celltype) <- c("cluster_num")
      cluster_celltype$ROIName <- rownames(cluster_celltype)
      cluster_celltype <- merge(cluster_celltype, color_cluster, all.x = TRUE)
      ROI_specific <- ROI[ROI$SlideName == input$slide, ]
      ROI_specific$ROIName <- paste(ROI_specific$ROINumber, ROI_specific$ROILabel, sep = "_")
      ROI_specific <- merge(ROI_specific, cluster_celltype, all.x = TRUE)
      ROI_specific_stroma <- ROI_specific$ROILabel %in% c("Stroma", "Unspecified")
      ROICoordinateX_stroma <- ROI_specific[ROI_specific_stroma, ]$ROICoordinateX
      ROICoordinateY_stroma <- ROI_specific[ROI_specific_stroma, ]$ROICoordinateY
      ROI_color_stroma <- ROI_specific[ROI_specific_stroma, ]$color
      for (i in 1:sum(ROI_specific_stroma)) {
        color <- col2rgb(ROI_color_stroma[i])
        draw.circle(x = ROICoordinateX_stroma[i], y = ROICoordinateY_stroma[i], radius = 10,
                    col = rgb(color[1], color[2], color[3], max = 255, alpha = 200),
                    border = rgb(color[1], color[2], color[3], max = 255, alpha = 200))
      }}
  })
    
  output$outline_cluster <- renderPlot({
    outline_cluster_plot()
  })
  
  heatmap_plot <- eventReactive(c(spatialdecon_results(), input$slide, input$clusters), {
    if (spatialdecon_results()$SlideName == input$slide) {
      cluster_celltype <- hclust(dist(t(spatialdecon_results()$beta)), method="complete")
      heatmap.2(spatialdecon_results()$beta, Rowv = FALSE, col = rev(brewer.pal(11,"RdBu")),
                Colv=as.dendrogram(cluster_celltype) %>% 
                  set("branches_k_color", value = as.vector(color_cluster$color[1:input$clusters]), k = input$clusters) %>% 
                  set("branches_lwd", 3), dendrogram = "column",
                scale="row", density.info="none", trace="none")
    } else {
      plot.new()
    }
  })
  
  output$heatmap <- renderPlot({
    heatmap_plot()
  })
  
  # Table
  
  results_table <- eventReactive({
    spatialdecon_results()
    input$slide}, {
    if (spatialdecon_results()$SlideName == input$slide) {
      beta_matrix <- spatialdecon_results()$beta
      colnames(beta_matrix) <- gsub("_", " ", gsub("_Unspecified", "", colnames(beta_matrix)))
      beta_matrix
    } else {
      # print("empty df")
      data.frame()
    }
  })
  
  output$table <- DT::renderDataTable({
    results_table()},
    options = list(scrollX=TRUE, scrollCollapse=TRUE, scrollY=TRUE, paging = FALSE), width = "100%")
}