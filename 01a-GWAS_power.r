library(data.table)
library(ggplot2)
library(ggpubr)

N <- 5976
alpha <- 5e-8
H2 <- 0.005

threshold <- qchisq(alpha, df = 1, lower.tail = FALSE)
power <- pchisq(threshold, df = 1, lower.tail = FALSE, ncp = N * H2)
power + 0.93*power

power_func <- function(freq, h2) {
    f = freq
    y.explained = h2
    b.alt = sqrt(y.explained / (2*f*(1-f)) ) #this is beta that explains 0.5%
    sigma = sqrt(1 - y.explained) #error sd after SNP effect is accounted for
    ns = seq(1000,10000,10) #candidate n
    ses = sigma / sqrt( ns*2*f*(1-f) ) #SE corresponding to each n
    q.thresh = qchisq(5e-8, df = 1, ncp = 0, lower = F) #threshold corresp alpha=5e-8
    pwr = pchisq(q.thresh, df = 1, ncp = (b.alt/ses)^2, lower = F) #power at alpha=5e-8
    #plot(ns,pwr, col = "darkgreen", xlab = "n", ylab = "power", 
    #     main = paste0("QT sd=1; MAF=",f,"; beta=",b.alt), t = "l", lwd = 1.5)
    #abline( h = 0.9, lty = 2 )
    pwr_df <- data.table(ns=ns, ses=ses, pwr=pwr, y.explained=y.explained)
    return(pwr_df)
}

pwr_df_1 <- power_func(0.25, 0.01)
pwr_df_2 <- power_func(0.25, 0.0075)
pwr_df_3 <- power_func(0.25, 0.005)
pwr_df_4 <- power_func(0.25, 0.0025)
pwr_df <- rbind(pwr_df_1, pwr_df_2, pwr_df_3, pwr_df_4)
pwr_df$y.explained <- factor(pwr_df$y.explained, levels = c(0.0100, 0.0075, 0.0050, 0.0025))

power_plot <- ggplot(data=pwr_df, aes(x=ns, y=pwr, colour = y.explained)) + 
    geom_line(size=1.3) +
    scale_x_continuous(breaks = pretty(pwr_df$ns, n = 10)) +
    scale_y_continuous(breaks = pretty(pwr_df$pwr, n = 10)) +
    ylab("Power (%)") + xlab("Sample-size") +
    guides(color=guide_legend(title=bquote('SNP'~h[2]))) +
    geom_hline(yintercept = 0.8, linetype="dashed", color = "black", size=1) + 
    geom_vline(xintercept = 5976, linetype="dashed", color = "firebrick", size=1) + theme_pubr() + 
    theme(axis.title = element_text(size = 16), legend.text = element_text(size = 16), legend.title = element_text(size = 16)) 

ggsave(plot = power_plot, 
       filename = "#.tiff", 
        device="tiff", width = 280, height = 190, dpi=300, compression = "lzw", units="mm")


