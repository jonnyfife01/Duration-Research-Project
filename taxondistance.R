library(chronosphere)
# library(divDyn)
library(geosphere)




# Get the data from pbdb
# updated may 17, 2026
dat <- chronosphere::fetch("pbdb", ser="occs4", ver="20260517")

# Filter records not identified at least to genus
dat <- dat[dat$accepted_rank %in% c("genus", "subgenus", "species", "subspecies"),]

# Omit non informative genus entries
dat <- dat[dat$genus != "",]

# Filter phyla
marineNoPlant <- c("",
                   "Agmata",
                   "Annelida",
                   "Bilateralomorpha",
                   "Brachiopoda",
                   "Bryozoa",
                   "Calcispongea",
                   "Chaetognatha",
                   "Cnidaria",
                   "Ctenophora",
                   "Echinodermata",
                   "Entoprocta",
                   "Foraminifera",
                   "Hemichordata",
                   "Hyolitha",
                   "Mollusca",
                   "Nematoda",
                   "Nematomorpha",
                   "Nemertina",
                   "Onychophora",
                   "Petalonamae",
                   "Phoronida",
                   "Platyhelminthes",
                   "Porifera",
                   "Problematica",
                   "Rhizopodea",
                   "Rotifera",
                   "Sarcomastigophora",
                   "Sipuncula",
                   "Uncertain",
                   "Vetulicolia",
                   ""
)

# which rows?
bByPhyla <- dat$phylum%in% marineNoPlant

#B. classes
needClass <- c(
  "Artiopoda",
  "Branchiopoda",
  "Cephalocarida",
  "Copepoda",
  "Malacostraca",
  "Maxillopoda",
  "Megacheira",
  "Merostomoidea",
  "Ostracoda",
  "Paratrilobita",
  "Pycnogonida",
  "Remipedia",
  "Thylacocephala",
  "Trilobita",
  "Xiphosura"
)

# which rows?
bNeedClass <- dat$class %in% needClass

# subset the data
dat <- dat[bByPhyla | bNeedClass , ]

# remove potential homonymy problems
dat$clgen <- paste(dat$class, dat$genus)

# number of remaining rows
nrow(dat)

# Check to make sure only marine invertebrates are left
unique(dat$phylum)

# remove non-marine environments
omitEnv <- c(
  "\"floodplain\"",
  "alluvial fan",
  "cave",
  "\"channel\"",
  "channel lag" ,
  "coarse channel fill",
  "crater lake",
  "crevasse splay",
  "dry floodplain",
  "delta plain",
  "dune",
  "eolian indet.",
  "fine channel fill",
  "fissure fill",
  "fluvial indet.",
  "fluvial-lacustrine indet.",
  "fluvial-deltaic indet.",
  "glacial",
  "interdune",
  "karst indet.",
  "lacustrine - large",
  "lacustrine - small",
  "lacustrine delta front",
  "lacustrine delta plain",
  "lacustrine deltaic indet.",
  "lacustrine indet.",
  "lacustrine interdistributary bay",
  "lacustrine prodelta",
  "levee",
  "loess",
  "mire/swamp",
  "pond",
  "sinkhole",
  "spring",
  "tar",
  "terrestrial indet.",
  "wet floodplain"
)

# actual omission
dat<-dat[!dat$environment%in%omitEnv, ]

# Filter genera with only one entry
gencount <- as.data.frame(table(dat$genus))

multgen <- gencount$Var1[gencount$Freq > 1]
multgen <- as.character(multgen)

dat <- dat[dat$genus %in% multgen, ]

# Duration per genus

duration <- c()

for (i in unique(dat$genus)) {
  # use upper max age and lower min age to be generous
  maxage <- max(dat$max_ma[dat$genus == i])
  minage <- min(dat$min_ma[dat$genus == i])
  duration[i] <- maxage - minage
}

head(duration)

# Histogram of durations
hist(duration)

# Find the distance between oldest and youngest occurrence

distance <- c()

for (i in unique(dat$genus)) {
  oldlat <- dat$paleolat[dat$genus == i & (dat$max_ma[dat$genus == i] == max(dat$max_ma[dat$genus == i]))]
  oldlong <- dat$paleolng[dat$genus == i & (dat$max_ma[dat$genus == i] == max(dat$max_ma[dat$genus == i]))]
  younglat <- dat$paleolat[dat$genus == i & (dat$min_ma[dat$genus == i] == min(dat$min_ma[dat$genus == i]))]
  younglong <- dat$paleolng[dat$genus == i & (dat$min_ma[dat$genus == i] == max(dat$min_ma[dat$genus == i]))]
  # Find haversine distance
  # Sometimes there are multiple occurrences of same age so only include first entry
  # So that the distance can actually be calculated and not return an error
  point1 <- c(oldlong[1], oldlat[1])
  point2 <- c(younglong[1], younglat[1])
  
  distance[i] <- distHaversine(point1, point2)
  
}

head(distance)

# Histogram of the distances
hist(distance, main = "Genus Distance from Origin at Extinction",
     xlab = "Distance(m)")

results <- cbind(duration, distance)
results <- as.data.frame(results)

# log results
logresults <- as.data.frame(log10(results))
results$genus <- unique(dat$genus)
logresults$genus <- unique(dat$genus)

# Clean data to remove
# Kinda problematic because it is possible that there are enough 0/infinities to affect results
# The dataset is so big this hopefully won't happen
logresults <- logresults[!is.na(logresults$distance),]
logresults <- logresults[logresults$distance != -Inf, ]

hist(logresults$duration)
hist(logresults$distance)

# Statistics of Results
resultstat <- lm(distance ~ duration, results)
summary(resultstat)

logresultstat <- lm(distance ~ duration, logresults)
summary(logresultstat)

# Plot Distance and Duration
plot(duration, distance,
     xlab = "Genus Lifespan (Ma)", ylab = "Dispersal Distance (m)",
     main = "Relationship between Genus Lifespan and Dispersion",
     col = "#1A611E50")
abline(resultstat)
legend("topright",legend = bquote(r^2 == .(summary(resultstat)$r.squared)))

plot(logresults$duration, logresults$distance,
     xlab = expression("log(Genus Lifespan (Ma))"), ylab = "log(Dispersal Distance (m))",
     main = "Relationship between Genus Lifespan and Dispersion (log)",
     col = "#7C7CD950")
abline(logresultstat)
legend("topleft",legend = bquote(r^2 == .(summary(logresultstat)$r.squared)))

# Test for correlation between duration and distance

cor.test(duration, distance)

cor.test(logresults$duration, logresults$distance)




####################
# Compare different groups

# Comparing taxa which originate in tropics to those which do not
troplat <- 23.4 # establishing tropics at 23.4 degrees latitude

# use a for loop to determine which genera originate in the tropics
tropicgenera <- c()
nontropicgenera <- c()

for(i in unique(dat$genus)){
  if(min(dat$paleolat[dat$genus == i & (dat$max_ma[dat$genus == i] == max(dat$max_ma[dat$genus == i]))] < troplat, na.rm = TRUE)) {
    tropicgenera[i] <- i
  } else{
    nontropicgenera[i] <- i
  }
}

length(unique(dat$genus))

length(tropicgenera) + length(nontropicgenera)

# Compare duration/distances for tropical/nontropical genera
tropicresults <- results[results$genus %in% tropicgenera, ]

nontropicresults <- results[results$genus %in% nontropicgenera, ]

logtropicresults <- tropicresults
logtropicresults$duration <- log10(logtropicresults$duration)
logtropicresults$distance <- log10(logtropicresults$distance)
logtropicresults <- logtropicresults[!is.na(logtropicresults),]
logtropicresults <- logtropicresults[logtropicresults$distance != -Inf,]

lognontropicresults <- nontropicresults
lognontropicresults$duration <- log10(lognontropicresults$duration)
lognontropicresults$distance <- log10(lognontropicresults$distance)
lognontropicresults <- lognontropicresults[!is.na(lognontropicresults),]
lognontropicresults <- lognontropicresults[lognontropicresults$distance != -Inf,]

# Histograms of duration and distance for tropic and nontropic genera
hist(tropicresults$duration)
hist(tropicresults$distance)

hist(logtropicresults$duration)
hist(logtropicresults$distance)

hist(nontropicresults$duration)
hist(nontropicresults$distance)

hist(lognontropicresults$duration)
hist(lognontropicresults$distance)

# Statistics of tropic results
tropicstat <- lm(distance ~ duration, logtropicresults)
summary(tropicstat)

nontropicstat <- lm(distance ~ duration, lognontropicresults)
summary(nontropicstat)

# Scatterplots for results
plot(logtropicresults$duration, logtropicresults$distance,
     main = "Dispersion of Tropical Genera",
     xlab = "log(Genus Duration(Ma))",
     ylab = "Log(Dispersal Distance (m))",
     col = "#34eb3480", pch = 16)
abline(tropicstat)
legend("topleft",legend = bquote(r^2 == .(summary(tropicstat)$r.squared)))

plot(lognontropicresults$duration, lognontropicresults$distance,
     main = "Dispersion of Non-Tropical Genera",
     xlab = "log(Genus Duration(Ma))",
     ylab = "Log(Dispersal Distance(m))",
     col = "#38e8e850", pch = 16)
abline(nontropicstat)
legend("topleft",legend = bquote(r^2 == .(summary(nontropicstat)$r.squared)))

################################
# Compare results for taxa which are regarded as disperers vs nondispersers
# Using Bivalves vs Brachiopods
bivalves <- dat[dat$class == "Bivalvia", ]
bivalvegenus <- unique(bivalves$genus)

brachiopods <- dat[dat$phylum == "Brachiopoda", ]
brachgenus <- unique(brachiopods$genus)


# Get logged stats for both  groups
bivalveresults <- logresults[bivalvegenus %in% logresults$genus, ]
bivalvestats <- lm(data = bivalveresults, distance ~ duration)
summary(bivalvestats)

brachresults <- logresults[brachgenus %in% logresults$genus, ]
brachstats <- lm(distance ~ duration, brachresults)
summary(brachstats)

# Plot results
plot(bivalveresults$duration, bivalveresults$distance,
     main = "Dispersal of Fossil Bivalves",
     ylab = "log(Dispersal Distance(m))",
     xlab = "log(Duration of Genus(Ma))",
     col = "#0004fc50", pch = 16)
abline(bivalvestats)
legend("topleft",legend = bquote(r^2 == .(summary(bivalvestats)$r.squared)))

plot(brachresults$duration, brachresults$distance,
     main = "Dispersal of Fossil Brachiopods",
     ylab = "log(Dispersal Distance(m))",
     xlab = "log(Duration of Genus(Ma))",
     col = "#fc500070", pch = 16)
abline(brachstats)
legend("topleft",legend = bquote(r^2 == .(summary(brachstats)$r.squared)))

###########################################
# Put all the results into a pdf
# Make sure to change path if editing for personal use
pdf("Results/graphs.pdf")

hist(duration, main = "Distribution of Genus Durations",
     xlab = "Genus Duration (Ma)")

hist(distance, main = "Genus Distance from Origin at Extinction",
     xlab = "Distance(m)")

plot(duration, distance,
     xlab = "Genus Lifespan (Ma)", ylab = "Dispersal Distance (m)",
     main = "Relationship between Genus Lifespan and Dispersion",
     col = "#1A611E50")
abline(resultstat)
legend("topright",legend = bquote(r^2 == .(summary(resultstat)$r.squared)))
print(summary(resultstat))

plot(logresults$duration, logresults$distance,
     xlab = expression("log(Genus Lifespan (Ma))"), ylab = "log(Dispersal Distance (m))",
     main = "Relationship between Genus Lifespan and Dispersion (log)",
     col = "#7C7CD950")
abline(logresultstat)
legend("topleft",legend = bquote(r^2 == .(summary(logresultstat)$r.squared)))


plot(logtropicresults$duration, logtropicresults$distance,
     main = "Dispersion of Tropical Genera",
     xlab = "log(Genus Duration(Ma))",
     ylab = "Log(Dispersal Distance (m))",
     col = "#34eb3480", pch = 16)
abline(tropicstat)
legend("topleft",legend = bquote(r^2 == .(summary(tropicstat)$r.squared)))

plot(lognontropicresults$duration, lognontropicresults$distance,
     main = "Dispersion of Non-Tropical Genera",
     xlab = "log(Genus Duration(Ma))",
     ylab = "Log(Dispersal Distance(m))",
     col = "#38e8e850", pch = 16)
abline(nontropicstat)
legend("topleft",legend = bquote(r^2 == .(summary(nontropicstat)$r.squared)))

plot(bivalveresults$duration, bivalveresults$distance,
     main = "Dispersal of Fossil Bivalves",
     ylab = "log(Dispersal Distance(m))",
     xlab = "log(Duration of Genus(Ma))",
     col = "#0004fc50", pch = 16)
abline(bivalvestats)
legend("topleft",legend = bquote(r^2 == .(summary(bivalvestats)$r.squared)))

plot(brachresults$duration, brachresults$distance,
     main = "Dispersal of Fossil Brachiopods",
     ylab = "log(Dispersal Distance(m))",
     xlab = "log(Duration of Genus(Ma))",
     col = "#fc500070", pch = 16)
abline(brachstats)
legend("topleft",legend = bquote(r^2 == .(summary(brachstats)$r.squared)))

dev.off()


###############
# Making a bunch of pngs of the graphs
png("Results/primary.png")
plot(duration, distance,
     xlab = "Genus Lifespan (Ma)", ylab = "Dispersal Distance (m)",
     main = "Relationship between Genus Lifespan and Dispersion",
     col = "#1A611E50")
abline(resultstat)
legend("topright",legend = bquote(r^2 == .(summary(resultstat)$r.squared)))
dev.off()

png("Results/log primary.png")
plot(logresults$duration, logresults$distance,
     xlab = expression("log(Genus Lifespan (Ma))"), ylab = "log(Dispersal Distance (m))",
     main = "Relationship between Genus Lifespan and Dispersion (log)",
     col = "#7C7CD950")
abline(logresultstat)
legend("topleft",legend = bquote(r^2 == .(summary(logresultstat)$r.squared)))
dev.off()

png("Results/tropic.png")
plot(logtropicresults$duration, logtropicresults$distance,
     main = "Dispersion of Tropical Genera",
     xlab = "log(Genus Duration(Ma))",
     ylab = "Log(Dispersal Distance (m))",
     col = "#34eb3480", pch = 16)
abline(tropicstat)
legend("topleft",legend = bquote(r^2 == .(summary(tropicstat)$r.squared)))
dev.off()


png("Results/nontropic.png")
plot(lognontropicresults$duration, lognontropicresults$distance,
     main = "Dispersion of Non-Tropical Genera",
     xlab = "log(Genus Duration(Ma))",
     ylab = "Log(Dispersal Distance(m))",
     col = "#38e8e850", pch = 16)
abline(nontropicstat)
legend("topleft",legend = bquote(r^2 == .(summary(nontropicstat)$r.squared)))
dev.off()

png("Results/bivalves.png")
plot(bivalveresults$duration, bivalveresults$distance,
     main = "Dispersal of Fossil Bivalves",
     ylab = "log(Dispersal Distance(m))",
     xlab = "log(Duration of Genus(Ma))",
     col = "#0004fc50", pch = 16)
abline(bivalvestats)
legend("topleft",legend = bquote(r^2 == .(summary(bivalvestats)$r.squared)))
dev.off()

png("Results/brachiopods.png")
plot(brachresults$duration, brachresults$distance,
     main = "Dispersal of Fossil Brachiopods",
     ylab = "log(Dispersal Distance(m))",
     xlab = "log(Duration of Genus(Ma))",
     col = "#fc500070", pch = 16)
abline(brachstats)
legend("topleft",legend = bquote(r^2 == .(summary(brachstats)$r.squared)))
dev.off()
