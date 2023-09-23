setwd('/home/zcqsdhj/TranscriptionalSignatures/Revisions/Data')

load('TCGA_BRCA.ExprDeconvolution_050823_p0.79_training.Rdata')
z2 <- Z.tumor_training[,apply(Z.tumor_training, 2, sum) > 0]
z3 <- log2(z2+1)
z4 <- apply(z3, 2, scale)
rownames(z4) <- rownames(z3)

Z.tumor_training <- z4
save(Z.tumor_training, file = 'TCGA_BRCA.ExprDeconvolution_270723_p0.79_training.Rdata')
