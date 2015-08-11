print(GeneFC_cntrl_p45ko$Direction)
plot(GeneFC_cntrl_p45ko$PostFC, EBOut$PPDE, xlim =c(0,5), ylim=c(0,1), main="Control/Experimental FC vs. PPDE", 
     sub=GeneFC_cntrl_p45ko$Direction, xlab="EBSeq Posterior Fold Change", ylab="EBSeq posterior prob of DE")
abline(h=0.95)
#identify(GeneFC_cntrl_p45ko$PostFC, EBOut$PPDE, labels=names(GeneFC_cntrl_p45ko$PostFC))
