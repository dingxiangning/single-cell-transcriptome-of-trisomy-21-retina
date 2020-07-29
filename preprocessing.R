## Load and normalize

library(data.table)

load_other = function(){
    other = fread(file = "Adult_QC_Pass_Samples_TPM.txt")
    .genes = other$Gene.ID
    other = other[,-c(1,2)]
    other = data.frame(t(other), stringsAsFactors = F)
    colnames(other) = .genes
    return(other)
}

load_mixed = function(){
    mixed = fread(file = "rawUMI.txt")
    .genes = mixed$Gene.ID
    mixed = mixed[,-1]
    mixed = data.frame(t(mixed), stringsAsFactors = F)
    colnames(mixed) = .genes
    return(mixed)
}

mixed = load_mixed()
mixed = scale(mixed)
mixed = t(mixed)
mixed = mixed[!is.nan(mixed[,1]),]

other = load_other()
other = scale(other)
other = t(other)
other = other[!is.nan(other[,1]),]

## Remove batch effect

int.genes = intersect(rownames(mixed), rownames(other))
mixed = mixed[int.genes,]
other = other[int.genes,]

merged = cbind(mixed, other)
meta = data.frame(
    name = c(colnames(mixed), colnames(other)),
    batch = c(rep(1, ncol(mixed)), rep(2, ncol(other))),
    stringsAsFactors = F
)

library(sva)

model = model.matrix(~1, data = meta)
batch = meta$batch
merged = ComBat(dat = merged, batch = batch, mod = model, par.prior = T, prior.plots = F)

## silico

merged = t(merged)
merged = scale(merged)
merged = t(merged)
range(merged)

silico = rnorm(nrow(merged) * as.integer(ncol(merged) * 0.3))
silico = matrix(silico, nrow = nrow(merged))
rownames(silico) = rownames(merged)
colnames(silico) = paste0("silico_", 1:ncol(silico))

## load variable gene2 of chr21

library(readxl)

vgs = read_excel("Table_13_Integrated Quantitative Transcriptome Maps of Human Trisomy 21 Tissues and Cells.XLSX", skip = 2)
vg.names = intersect(vgs$Gene_name, rownames(merged))

vg.sel = sample(vg.names, length(vg.names)/2)
ce.sel = sample(1:ncol(silico), ncol(silico)/2)
silico[vg.sel, ce.sel] = silico[vg.sel, ce.sel] * 10

combined = cbind(merged, silico)
group = c(batch, rep(3, ncol(silico)))
names(group) = rep(c('mixed', 'normal', 'silico'), c(length(which(group==1)), length(which(group==2)), ncol(silico)))

saveRDS(combined, file = "combined.rds")
saveRDS(group, file = "group.rds")

rm(list = ls())
combined = readRDS("combined.rds")
group = readRDS("group.rds")

library(randomForest)
set.seed(10)

idx = which(names(group)!='mixed')
X = as.data.frame(t(combined[,idx]))
Y = names(group)[idx]
TX = as.data.frame(t(combined[,-idx]))

idx = sample(1:length(Y), length(Y))
X = X[idx,]
Y = Y[idx]
label = rep(1, length(Y))
label[Y == 'normal'] = 0
X$label = as.factor(label)

idx = sample(1:length(label), length(label)*0.8)
x.train = X[idx,]
y.train = label[idx]
x.test = X[-idx,]
y.test = label[-idx]

colnames(x.train) = gsub('-', '_', colnames(x.train))
colnames(x.test) = gsub('-', '_', colnames(x.test))

rf = randomForest(label ~ ., data = x.train, importance = T, proximity = T)

predictions = predict(rf, x.test)

plot(rf)

colnames(TX) = gsub('-', '_', colnames(TX))
predictions = predict(rf, TX)
table(predictions)

saveRDS(predictions, file = 'predictions.rds')
