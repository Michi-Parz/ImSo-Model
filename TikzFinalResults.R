library(tikzDevice)

tikzi <- theme(
  legend.background = element_rect(fill = "transparent"),
  legend.spacing.y = unit(-1.5, "cm"),
  axis.text = element_text(size = 15),
  axis.text.x = element_text(hjust = 0.7),
  axis.title = element_text(size = 15),
  plot.title = element_text(size = 15),
  legend.text = element_text(size = 15),
  legend.title = element_text(size = 15, margin = margin(b = 5)),
  panel.grid.major = element_line(linewidth = 1))


xlag_nam <- "ISoX"
ylag_nam <- "ISoY"

plot_madj <- 4.2


# Daten -------------------------


gg_tikz_save <- gg_basic(ntp_o_long, aes(x = fpHz, y = NTP, group = ID),
                         ymin = 0, ymax = 70, isl = T, lab = F)+
  geom_line(alpha = 0.5)+
  tikzi+
  labs(x = xlag_nam, y = ylag_nam)


gg_tikz_save
tikz("Tikz/Data.tex", width = plot_madj, height = plot_madj)
gg_tikz_save
dev.off()



# Data adjust--------------

gg_o <- gg_basic(ntp_o_long, aes(x = fpHz, y = NTP, group = ID),
                 ymin = -20, ymax = 70, isl = T, lab = F)+
  geom_line(alpha = 0.5)+
  labs(title = "Original", x = xlag_nam, y = ylag_nam)+
  theme(axis.text =  element_text(size = 10),
        axis.title = element_text(size = 12))

gg_a <- gg_basic(ntp_long, aes(x = fpHz, y = NTP, group = ID),
                 ymin = -20, ymax = 70, isl = T, lab = F)+
  geom_line(alpha = 0.5)+
  labs(title = "Adjusted", x = xlag_nam, y = ylag_nam)+
  theme(axis.text =  element_text(size = 10),
        axis.title = element_text(size = 12))


ggarrange(gg_o, gg_a)


ggarrange(gg_o, gg_a)
tikz("Tikz/DataAdjusted.tex", width = plot_madj+2, height = plot_madj)
ggarrange(gg_o, gg_a)
dev.off()



# Einfluss Daten Anpassung auf Fit-------------


gg_tikz_save <- gg_basic(vgl_data, aes(x = fpHz, y = Fit, color = Data, group = Data),
         ymin = 0, ymax = 60, isl = T)+
  geom_line()+
  labs(color = "Data basis", x = xlag_nam, y = ylag_nam)+
  theme(
    legend.position = c(0.4, 0.2),
    legend.justification = c(0.5, 0.5),
    legend.background = element_rect(fill = "transparent"),
    legend.spacing.y = unit(-0.7, "cm"),
    panel.grid.major = element_line(linewidth = 1))+
  tikzi


gg_tikz_save
tikz("Tikz/AdjustmentInfluence.tex", width = plot_madj, height = plot_madj)
gg_tikz_save
dev.off()



# Segmentierte Regression wie physikalsiches Modell-------------


gg_tikz_save <- gg_basic(
  lnTestDF,
  aes(y = IMso, x = fpHz, color = "L'n Model"),
  ymin = 0, ymax = 65,
  isl = T, lab = F
)+geom_line() + geom_point()+
  geom_line(aes(y = func_1b(fpHz,
                            segRegCoeff[1],
                            segRegCoeff[2],
                            segRegCoeff[3],
                            segRegCoeff[4]),
                color = "Regression"))+
  labs(color = "", x = xlag_nam, y = ylag_nam)+
  theme(
    legend.position = c(0.4, 0.2),
    legend.justification = c(0.5, 0.5),
    legend.background = element_rect(fill = "transparent"),
    legend.spacing.y = unit(-0.7, "cm"),
    panel.grid.major = element_line(linewidth = 1))+
  tikzi


gg_tikz_save
tikz("Tikz/SegRegLikeStandard.tex", width = plot_madj, height = plot_madj)
gg_tikz_save
dev.off()






# Modelkurve-------------



gg_tikz_save <- gg_basic(ntp_long,
                         ymin = -20, ymax = 70,
                         isl = T, lab = F)+
  geom_line(aes(x = fpHz, y = NTP, group = ID), alpha = 0.1)+
  geom_line(
    mapping = aes(x = fpHz, y = NTP),
    data = ntp_fit, linewidth = 0.8, color = "red"
  )+
  labs(x = xlag_nam, y = ylag_nam)+
  tikzi


gg_tikz_save
tikz("Tikz/FixEff.tex", width = plot_madj, height = plot_madj)
gg_tikz_save
dev.off()


# Residuen Analyse-----------

gg_tikz_save <- resfreq+
  labs(x = xlag_nam)+
  tikzi+
  theme(axis.text.x = element_text(size = 12))

gg_tikz_save
tikz("Tikz/ResiduenAnalyse.tex", width = plot_madj+3, height = plot_madj)
gg_tikz_save
dev.off()



# Posterior beta samples-------------




ggarrange(ggb0,ggb1, ggb2, ggb3)
tikz("Tikz/PosterBetaSamples.tex", width = plot_madj, height = plot_madj)
ggarrange(ggb0,ggb1, ggb2, ggb3)
dev.off()

# Posterior tau samples-------------

ggarrange(ggt1,ggt2)
tikz("Tikz/PosterTauSamples.tex", width = plot_madj, height = plot_madj)
ggarrange(ggt1,ggt2)
dev.off()




# Marginale Standardabweichung-------



gg_tikz_save <- ggplot(stdev_df)+
  aes(x = fpHz)+
  geom_line(aes(y = q50))+
  geom_ribbon(aes(ymin = q2.5, ymax = q97.5), alpha = 1/4)+
  scale_x_continuous(trans = log2_trans(),
                     minor_breaks = c(freq,6300,8000,10000),
                     breaks = c(63,125,250,500,1000,2000,4000,8000),
                     labels = c(63,125,250,500,1000,2000,4000,8000),
  )+
  labs(y = "Standard deviation / dB",
       x = xlag_nam)+
  theme(legend.position = "bottom")+
  tikzi


gg_tikz_save
tikz("Tikz/StDevMarginal.tex", width = plot_madj, height = plot_madj)
gg_tikz_save
dev.off()


# Corrmatrix--------------


corrplot::corrplot(mean_corrmat, "color", tl.col = "black")
tikz("Tikz/Correlation.tex", width = plot_madj, height = plot_madj)
corrplot::corrplot(mean_corrmat, "color", tl.col = "black")
dev.off()





# RE corrmat-------------


corrplot::corrplot(ReMeanCorr,"color",addCoef.col = 'black',tl.col = "black")
tikz("Tikz/RandoEffCorr.tex", width = plot_madj, height = plot_madj)
corrplot::corrplot(ReMeanCorr,"color",addCoef.col = 'black',tl.col = "black")
dev.off()


# RE stdev----------



gg_tikz_save <- ggplot(reSdDflong)+geom_boxplot(aes(x = Coeff, y = RE))+
  labs(y = "Standard deviation random effect",
       x = "Coefficient")+
  tikzi

gg_tikz_save
tikz("Tikz/RandoEffStdev.tex", width = plot_madj+2, height = plot_madj)
gg_tikz_save
dev.off()


# FitREeff------------

gg_tikz_save <- gg_basic(fit_re, aes(x = fpHz, group = ID),
                         ymin = -20, ymax = 70, isl = T, lab = F)+
  geom_line(aes(y = NTP, color = "Data"), alpha = 0.5)+
  geom_line(aes(y = Fit, color = "Fit"), alpha = 0.5)+
  labs(color = "", x = xlag_nam, y = ylag_nam)+
  scale_color_manual(values = c("black", "red"))+
  tikzi+
  theme(
    legend.position = c(0.25, 0.3),
    legend.justification = c(0.5, 0.5),
    legend.background = element_rect(fill = "transparent"),
    legend.spacing.y = unit(-0.7, "cm"),
    panel.grid.major = element_line(linewidth = 1))

gg_tikz_save
tikz("Tikz/FitREeff.tex", width = plot_madj, height = plot_madj)
gg_tikz_save
dev.off()


# Softplus---------------


gg_tikz_save <- ggplot(df_all, aes(x = x, y = y, color = curve, linetype = curve)) +
  geom_line(size = 1) +
  scale_color_manual(
    values = c(
      "softplus 1"   = "#d95f02",
      "softplus 2"   = "#7570b3",
      "softplus 10"  = "#e7298a",
      "softplus 20"  = "#66a61e",
      "ReLU" = "black"
    )
  ) +
  scale_linetype_manual(
    values = c(rep("solid", length(alphas)), "solid"),
    guide = "none" # falls du Linientyp nicht in der Legende willst
  ) +
  labs(
    title = expression(paste("softplus vs. ReLU")),
    x = "x",
    y = "Function value"
  ) +
  theme_minimal(base_size = 12)+
  theme(
    legend.position = c(0.2, 0.75),
    legend.justification = c(0.5, 0.5),
    legend.background = element_rect(fill = "transparent"),
    legend.spacing.y = unit(-0.7, "cm"),
    panel.grid.major = element_line(linewidth = 1))

gg_tikz_save
tikz("Tikz/Softplus.tex", width = plot_madj, height = plot_madj)
gg_tikz_save
dev.off()



# Smallest whitned residual---------


gg_tikz_save <- whiteResOutlier +
  tikzi +
  labs(x = xlag_nam, y = ylag_nam)

gg_tikz_save
tikz("Tikz/SmalestWhitnedResidalCurve.tex", width = plot_madj, height = plot_madj)
gg_tikz_save
dev.off()


# PPC Freq---------




gg_tikz_save <- ppcFreq +
  tikzi +
  theme(
    legend.position = c(0.4, 0.2),
    legend.justification = c(0.5, 0.5),
    legend.background = element_rect(fill = "transparent"),
    legend.spacing.y = unit(-0.7, "cm"),
    legend.title = element_text(size = 12),
    panel.grid.major = element_line(linewidth = 1))+
  labs(x = xlag_nam, y = ylag_nam)


gg_tikz_save
tikz("Tikz/PPCfreq.tex", width = plot_madj+1, height = plot_madj+1)
gg_tikz_save
dev.off()


# PPC Freq Biomonal---------




gg_tikz_save <- ppcFreqBinom +
  tikzi +
  labs(x = xlag_nam, y = "P(.Binomial - Expected. .= .OOR - Expected.)")


gg_tikz_save
tikz("Tikz/PPCfreqBinom.tex", width = plot_madj+3, height = plot_madj+1)
gg_tikz_save
dev.off()


# PPC Freq StDev---------




gg_tikz_save <- ppcFreqSt +
  tikzi +
  labs(x = xlag_nam)


gg_tikz_save
tikz("Tikz/PPCfreqSD.tex", width = plot_madj+3, height = plot_madj+1)
gg_tikz_save
dev.off()


# PPC Freq and StDev------

gg_tikz_save1 <- ppcFreq+
  tikzi+ labs(x = xlag_nam, y = ylag_nam)

gg_tikz_save2 <- ppcFreqSt+
  tikzi+ labs(x = xlag_nam)+
  theme(
        axis.text.x = element_text(angle =40, hjust = 1))


ggarrange(gg_tikz_save1, gg_tikz_save2, legend = "bottom", common.legend = T)
tikz("Tikz/PPCfreqBoth.tex", width = plot_madj+3, height = plot_madj+1)
ggarrange(gg_tikz_save1, gg_tikz_save2, legend = "bottom", common.legend = T)
dev.off()




# PPC White---------




gg_tikz_save <- PPCwhite +
  tikzi 

gg_tikz_save
tikz("Tikz/PPCwhite.tex", width = plot_madj+3, height = plot_madj+1)
gg_tikz_save
dev.off()



# SMD outlier---------


gg_tikz_save1 <- mahalanobisOutlier1+
  tikzi+
  labs(x = xlag_nam, y = ylag_nam)

gg_tikz_save2 <- mahalanobisOutlier2+
  tikzi+
  labs(x = xlag_nam, y = ylag_nam)

ggarrange(
  gg_tikz_save1, gg_tikz_save2,
  legend = "bottom",
  legend.grob = get_legend(gghelp+tikzi)
)
tikz("Tikz/PPCsmdCurves.tex", width = plot_madj+3, height = plot_madj+1)
ggarrange(
  gg_tikz_save1, gg_tikz_save2,
  legend = "bottom",
  legend.grob = get_legend(gghelp+tikzi)
)
dev.off()

