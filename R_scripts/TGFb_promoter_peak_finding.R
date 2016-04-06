#finding peaks in TGFb promoter
tp1 <- turnpoints(subsetByOverlaps(bA.cov.h2az.emt_markers.tgfb, gr.plot[2])$mean[1500:3500])
summary(tp1)
plot(tp1)

pdf("TGFb-peak_distances.pdf")
plot(1:2001, subsetByOverlaps(bA.cov.h2az.emt_markers.tgfb, gr.plot[2])$mean[1500:3500], type = "l", axes = F, xlab = "Relative distance [bp]", ylab = "[RPM]")
lines(tp1)
abline(v = 456)
abline(v = 1037)
abline(v = 1363)
axis(1, at = c(456, 1037, 1363), tick = T)
axis(2)
dev.off()


