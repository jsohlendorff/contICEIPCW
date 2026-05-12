### plot_treatment_breaks.R --- 
#----------------------------------------------------------------------
## Author: Johan Sebastian Ohlendorff
## Created: Apr 29 2026 (15:09) 
## Version: 
## Last-Updated: Apr 29 2026 (16:51) 
##           By: Johan Sebastian Ohlendorff
##     Update #: 70
#----------------------------------------------------------------------
## 
### Commentary: 
## 
### Change Log:
#----------------------------------------------------------------------
## 
### Code:
plot_treatment_breaks <- function(dt_list, treatments = "A", baseline_cols = "A_0", subset_size = NULL, type = "facet") {
    base_dt <- dt_list[[1]][, c("id", baseline_cols), with = FALSE]
    setnames(base_dt, old = baseline_cols, new = treatments)
    base_dt[, c("time", "event") := .(0, "A")]
    time_dt <- dt_list[[2]][event != "L", c("id", "time", "event", treatments), with = FALSE]
    time_dt <- rbindlist(list(base_dt, time_dt), use.names = TRUE, fill = TRUE)
    setkey(time_dt, id, time)
    time_dt[, time_next := shift(time, type = "lead"), by = id]

    terminal_events_dt <- dt_list[[2]][event %in% c("Y", "D", "tauend", "C"), .(id, event, time)]
    setnames(terminal_events_dt, c("time", "event"),  c("terminal_time", "terminal_event"))
    ## Add to data table for terminal events
    time_dt <- merge(time_dt, terminal_events_dt, by = "id", all.x = TRUE)
    time_dt <- time_dt[event == "A"]
    setnames(time_dt, c("time", "time_next"), c("treatment_start", "treatment_end"))

    out <- list()
    for (i in seq_along(treatments)) {
        df_temp <- time_dt[time_dt[[treatments[i]]] == 1, c("id", "treatment_start", "event", treatments[i], "treatment_end", "terminal_time", "terminal_event"), with = FALSE]
        set(
            df_temp, j = c(treatments[i], "event"), value = NULL
        )
        df_temp[,treatment := treatments[i]]
        out[[i]] <- df_temp
    }
    out_dt <- rbindlist(out, use.names = TRUE, fill = TRUE)
    if (type == "segment") {
        out_dt[, offset := as.numeric(factor(treatment)) * 0.15]
    }
    if (!is.null(subset_size)) {
        subset <- unique(out_dt$id)[1:subset_size]
        out_dt <- out_dt[id %in% subset]
    }

    n_ids <- length(unique(out_dt$id))

    # Define adaptive sizes (tune constants as needed)
    lw <- max(0.10, 5 / sqrt(n_ids))   # segment thickness
    pt <- max(0.2, 15 / sqrt(n_ids))     # point size
    

    
    if (type == "segment"){
        p<-ggplot(out_dt, aes(y = as.numeric(factor(id))+ offset)) +
            geom_segment(aes(x = treatment_start,
                                  xend = treatment_end,
                                  color = treatment),
                              linewidth = lw)
    } else {
        p <- ggplot(out_dt, aes(y = factor(id)))+
            geom_segment(aes(x = treatment_start,
                                  xend = treatment_end,
                                  yend = factor(id)),
                              linewidth = lw, color = "steelblue")
    }
    p <- p +
        
        # Terminal time as a dot
        geom_point(aes(x = terminal_time, color = terminal_event),
                   size = pt)
    if (type == "facet") {
        p <- p +
            facet_grid(~ treatment, scales = "free") 
    }
    p <- p + labs(x = "Time",
         y = "ID")
    if (n_ids > 50 || type == "segment") {
        p <- p + theme(axis.text.y = element_blank(),
                       axis.ticks.y = element_blank())
    } 
    p
}


######################################################################
### plot_treatment_breaks.R ends here
