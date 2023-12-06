rm(list=ls())

################################################################################
# Nonparametric Result
################################################################################

RRR <- read.csv("NP/Result_Merge_NP.csv")
Eff <- apply(RRR[,2:4],2,median)
SE  <- apply((RRR[,2:4] - matrix(Eff,dim(RRR)[1],3,byrow=T))^2 + RRR[,5:7]^2,
             2,median)^(0.5)
round(Eff,3); round(SE,3)

################################################################################
# Parametric Result
################################################################################

load("Result_Parametric.RData")

################################################################################
# Table
################################################################################

RRR <- rbind(
  sprintf("\\multirow{3}{*}{Semiparametric COCA} & Estimate & \\multicolumn{1}{c|}{$%0.3f$} & \\multicolumn{1}{c|}{$%0.3f$} & $%0.3f$ \\\\ \\cline{2-5}",
          Eff[2], Eff[1], Eff[3]),
  sprintf(" & SE & \\multicolumn{1}{c|}{$%0.3f$} & \\multicolumn{1}{c|}{$%0.3f$} & $%0.3f$ \\\\ \\cline{2-5}",
          SE[2], SE[1], SE[3]),
  sprintf(" & 95\\%% CI & \\multicolumn{1}{c|}{($%0.3f$,$%0.3f$)} & \\multicolumn{1}{c|}{($%0.3f$,$%0.3f$)} & ($%0.3f$,$%0.3f$) \\\\ \\hline",
          Eff[2]-qnorm(0.975)*SE[2], Eff[2]+qnorm(0.975)*SE[2],
          Eff[1]-qnorm(0.975)*SE[1], Eff[1]+qnorm(0.975)*SE[1],
          Eff[3]-qnorm(0.975)*SE[3], Eff[3]+qnorm(0.975)*SE[3]),
  sprintf("\\multirow{3}{*}{EPS parametric COCA} & Estimate & \\multicolumn{1}{c|}{$%0.3f$} & \\multicolumn{1}{c|}{$%0.3f$} & $%0.3f$ \\\\ \\cline{2-5}",
          Result.COCA.Para[[2]][1,1], 
          Result.COCA.Para[[1]][1,1], 
          Result.COCA.Para[[3]][1,1]),
  sprintf(" & SE & \\multicolumn{1}{c|}{$%0.3f$} & \\multicolumn{1}{c|}{$%0.3f$} & $%0.3f$ \\\\ \\cline{2-5}",
          Result.COCA.Para[[2]][1,2], 
          Result.COCA.Para[[1]][1,2], 
          Result.COCA.Para[[3]][1,2]),
  sprintf(" & 95\\%% CI & \\multicolumn{1}{c|}{($%0.3f$,$%0.3f$)} & \\multicolumn{1}{c|}{($%0.3f$,$%0.3f$)} & ($%0.3f$,$%0.3f$) \\\\ \\hline",
          Result.COCA.Para[[2]][1,3], Result.COCA.Para[[2]][1,4],  
          Result.COCA.Para[[1]][1,3], Result.COCA.Para[[1]][1,4],  
          Result.COCA.Para[[3]][1,3], Result.COCA.Para[[3]][1,4]),
  sprintf("\\multirow{3}{*}{Outcome bridge function parametric COCA} & Estimate & \\multicolumn{1}{c|}{$%0.3f$} & \\multicolumn{1}{c|}{$%0.3f$} & $%0.3f$ \\\\ \\cline{2-5}",
          Result.COCA.Para[[2]][2,1], 
          Result.COCA.Para[[1]][2,1], 
          Result.COCA.Para[[3]][2,1]),
  sprintf(" & SE & \\multicolumn{1}{c|}{$%0.3f$} & \\multicolumn{1}{c|}{$%0.3f$} & $%0.3f$ \\\\ \\cline{2-5}",
          Result.COCA.Para[[2]][2,2], 
          Result.COCA.Para[[1]][2,2], 
          Result.COCA.Para[[3]][2,2]),
  sprintf(" & 95\\%% CI & \\multicolumn{1}{c|}{($%0.3f$,$%0.3f$)} & \\multicolumn{1}{c|}{($%0.3f$,$%0.3f$)} & ($%0.3f$,$%0.3f$) \\\\ \\hline",
          Result.COCA.Para[[2]][2,3], Result.COCA.Para[[2]][2,4],  
          Result.COCA.Para[[1]][2,3], Result.COCA.Para[[1]][2,4],  
          Result.COCA.Para[[3]][2,3], Result.COCA.Para[[3]][2,4]),
  sprintf("\\multirow{3}{*}{Doubly-robust parametric COCA} & Estimate & \\multicolumn{1}{c|}{$%0.3f$} & \\multicolumn{1}{c|}{$%0.3f$} & $%0.3f$ \\\\ \\cline{2-5}",
          Result.COCA.Para[[2]][3,1], 
          Result.COCA.Para[[1]][3,1], 
          Result.COCA.Para[[3]][3,1]),
  sprintf(" & SE & \\multicolumn{1}{c|}{$%0.3f$} & \\multicolumn{1}{c|}{$%0.3f$} & $%0.3f$ \\\\ \\cline{2-5}",
          Result.COCA.Para[[2]][3,2], 
          Result.COCA.Para[[1]][3,2], 
          Result.COCA.Para[[3]][3,2]),
  sprintf(" & 95\\%% CI & \\multicolumn{1}{c|}{($%0.3f$,$%0.3f$)} & \\multicolumn{1}{c|}{($%0.3f$,$%0.3f$)} & ($%0.3f$,$%0.3f$) \\\\ \\hline",
          Result.COCA.Para[[2]][3,3], Result.COCA.Para[[2]][3,4],  
          Result.COCA.Para[[1]][3,3], Result.COCA.Para[[1]][3,4],  
          Result.COCA.Para[[3]][3,3], Result.COCA.Para[[3]][3,4]),
  sprintf("\\multirow{3}{*}{Standard DiD under parallel trends} & Estimate & \\multicolumn{1}{c|}{$%0.3f$} & \\multicolumn{1}{c|}{$%0.3f$} & $%0.3f$ \\\\ \\cline{2-5}",
          Result.PT[2,1], 
          Result.PT[1,1], 
          Result.PT[3,1]),
  sprintf(" & SE & \\multicolumn{1}{c|}{$%0.3f$} & \\multicolumn{1}{c|}{$%0.3f$} & $%0.3f$ \\\\ \\cline{2-5}",
          Result.PT[2,2], 
          Result.PT[1,2], 
          Result.PT[3,2]),
  sprintf(" & 95\\%% CI & \\multicolumn{1}{c|}{($%0.3f$,$%0.3f$)} & \\multicolumn{1}{c|}{($%0.3f$,$%0.3f$)} & ($%0.3f$,$%0.3f$) \\\\ \\hline",
          Result.PT[2,3], Result.PT[2,4],
          Result.PT[1,3], Result.PT[1,4],
          Result.PT[3,3], Result.PT[3,4]) )

print(data.frame(RRR),row.names=F)

RRR <- rbind(
  sprintf("%0.3f & %0.3f & %0.3f",
          Eff[2], Eff[1], Eff[3]),
  sprintf("%0.3f & %0.3f & %0.3f",
          SE[2], SE[1], SE[3]),
  sprintf("(%0.3f,%0.3f) & (%0.3f,%0.3f) & (%0.3f,%0.3f)",
          Eff[2]-qnorm(0.975)*SE[2], Eff[2]+qnorm(0.975)*SE[2],
          Eff[1]-qnorm(0.975)*SE[1], Eff[1]+qnorm(0.975)*SE[1],
          Eff[3]-qnorm(0.975)*SE[3], Eff[3]+qnorm(0.975)*SE[3]),
  sprintf("%0.3f & %0.3f & %0.3f",
          Result.COCA.Para[[2]][3,1], 
          Result.COCA.Para[[1]][3,1], 
          Result.COCA.Para[[3]][3,1]),
  sprintf("%0.3f & %0.3f & %0.3f",
          Result.COCA.Para[[2]][3,2], 
          Result.COCA.Para[[1]][3,2], 
          Result.COCA.Para[[3]][3,2]),
  sprintf("(%0.3f,%0.3f) & (%0.3f,%0.3f) & (%0.3f,%0.3f)",
          Result.COCA.Para[[2]][3,3], Result.COCA.Para[[2]][3,4],  
          Result.COCA.Para[[1]][3,3], Result.COCA.Para[[1]][3,4],  
          Result.COCA.Para[[3]][3,3], Result.COCA.Para[[3]][3,4]),
  sprintf("%0.3f & %0.3f & %0.3f",
          Result.PT[2,1], 
          Result.PT[1,1], 
          Result.PT[3,1]),
  sprintf("%0.3f & %0.3f & %0.3f",
          Result.PT[2,2], 
          Result.PT[1,2], 
          Result.PT[3,2]),
  sprintf("(%0.3f,%0.3f) & (%0.3f,%0.3f) & (%0.3f,%0.3f)",
          Result.PT[2,3], Result.PT[2,4],
          Result.PT[1,3], Result.PT[1,4],
          Result.PT[3,3], Result.PT[3,4]) )

print(data.frame(RRR),row.names=F)













RRR <- rbind(
  sprintf("\\multicolumn{2}{|c|}{\\multirow{3}{*}{Semiparametric}} & Estimator & \\multicolumn{1}{c|}{$%0.3f$} & \\multicolumn{1}{c|}{$%0.3f$} & $%0.3f$ \\\\ \\cline{3-6} ",
          Eff[2], Eff[1], Eff[3]),
  sprintf("\\multicolumn{1}{|c}{} & & SE & \\multicolumn{1}{c|}{$%0.3f$} & \\multicolumn{1}{c|}{$%0.3f$} & $%0.3f$ \\\\ \\cline{3-6} ",
          SE[2], SE[1], SE[3]),
  sprintf("\\multicolumn{1}{|c}{} & & 95\\%% CI & \\multicolumn{1}{c|}{$(%0.3f,%0.3f)$} & \\multicolumn{1}{c|}{$(%0.3f,%0.3f)$} & $(%0.3f,%0.3f)$ \\\\ \\hline ",
          Eff[2]-qnorm(0.975)*SE[2], Eff[2]+qnorm(0.975)*SE[2],
          Eff[1]-qnorm(0.975)*SE[1], Eff[1]+qnorm(0.975)*SE[1],
          Eff[3]-qnorm(0.975)*SE[3], Eff[3]+qnorm(0.975)*SE[3]),
  sprintf("\\multicolumn{1}{|c|}{\\multirow{9}{*}{Parametric}} & \\multirow{3}{*}{EPS} & Estimator & \\multicolumn{1}{c|}{$%0.3f$} & \\multicolumn{1}{c|}{$%0.3f$} & $%0.3f$ \\\\ \\cline{3-6} ",
          Result.COCA.Para[[2]][1,1], 
          Result.COCA.Para[[1]][1,1], 
          Result.COCA.Para[[3]][1,1]),
  sprintf("\\multicolumn{1}{|c|}{} & & SE & \\multicolumn{1}{c|}{$%0.3f$} & \\multicolumn{1}{c|}{$%0.3f$} & $%0.3f$ \\\\ \\cline{3-6} ",
          Result.COCA.Para[[2]][1,2], 
          Result.COCA.Para[[1]][1,2], 
          Result.COCA.Para[[3]][1,2]),
  sprintf("\\multicolumn{1}{|c|}{} & & 95\\%% CI & \\multicolumn{1}{c|}{$(%0.3f,%0.3f)$} & \\multicolumn{1}{c|}{$(%0.3f,%0.3f)$} & $(%0.3f,%0.3f)$ \\\\ \\cline{2-6} ",
          Result.COCA.Para[[2]][1,3], Result.COCA.Para[[2]][1,4],  
          Result.COCA.Para[[1]][1,3], Result.COCA.Para[[1]][1,4],  
          Result.COCA.Para[[3]][1,3], Result.COCA.Para[[3]][1,4]),
  sprintf("\\multicolumn{1}{|c|}{} & \\multirow{3}{*}{\\begin{tabular}[c]{@{}c@{}}COCA \\\\ bridge\\\\ function\\end{tabular}} & Estimator & \\multicolumn{1}{c|}{$%0.3f$} & \\multicolumn{1}{c|}{$%0.3f$} & $%0.3f$ \\\\ \\cline{3-6} ",
          Result.COCA.Para[[2]][2,1], 
          Result.COCA.Para[[1]][2,1], 
          Result.COCA.Para[[3]][2,1]),
  sprintf("\\multicolumn{1}{|c|}{} & & SE & \\multicolumn{1}{c|}{$%0.3f$} & \\multicolumn{1}{c|}{$%0.3f$} & $%0.3f$ \\\\ \\cline{3-6} ",
          Result.COCA.Para[[2]][2,2], 
          Result.COCA.Para[[1]][2,2], 
          Result.COCA.Para[[3]][2,2]),
  sprintf("\\multicolumn{1}{|c|}{} & & 95\\%% CI & \\multicolumn{1}{c|}{$(%0.3f,%0.3f)$} & \\multicolumn{1}{c|}{$(%0.3f,%0.3f)$} & $(%0.3f,%0.3f)$ \\\\ \\cline{2-6} ",
          Result.COCA.Para[[2]][2,3], Result.COCA.Para[[2]][2,4],  
          Result.COCA.Para[[1]][2,3], Result.COCA.Para[[1]][2,4],  
          Result.COCA.Para[[3]][2,3], Result.COCA.Para[[3]][2,4]),
  sprintf("\\multicolumn{1}{|c|}{} & \\multirow{3}{*}{Doubly-robust} & Estimator & \\multicolumn{1}{c|}{$%0.3f$} & \\multicolumn{1}{c|}{$%0.3f$} & $%0.3f$ \\\\ \\cline{3-6} ",
          Result.COCA.Para[[2]][3,1], 
          Result.COCA.Para[[1]][3,1], 
          Result.COCA.Para[[3]][3,1]),
  sprintf("\\multicolumn{1}{|c|}{} & & SE & \\multicolumn{1}{c|}{$%0.3f$} & \\multicolumn{1}{c|}{$%0.3f$} & $%0.3f$ \\\\ \\cline{3-6} ",
          Result.COCA.Para[[2]][3,2], 
          Result.COCA.Para[[1]][3,2], 
          Result.COCA.Para[[3]][3,2]),
  sprintf("\\multicolumn{1}{|c|}{} & & 95\\%% CI & \\multicolumn{1}{c|}{$(%0.3f,%0.3f)$} & \\multicolumn{1}{c|}{$(%0.3f,%0.3f)$} & $(%0.3f,%0.3f)$ \\\\ \\hline",
          Result.COCA.Para[[2]][3,3], Result.COCA.Para[[2]][3,4],  
          Result.COCA.Para[[1]][3,3], Result.COCA.Para[[1]][3,4],  
          Result.COCA.Para[[3]][3,3], Result.COCA.Para[[3]][3,4]),
  sprintf("\\multicolumn{2}{|c|}{\\multirow{3}{*}{DiD under parallel trends}} & Estimator & \\multicolumn{1}{c|}{$%0.3f$} & \\multicolumn{1}{c|}{$%0.3f$} & $%0.3f$ \\\\ \\cline{3-6} ",
          Result.PT[2,1], 
          Result.PT[1,1], 
          Result.PT[3,1]),
  sprintf("\\multicolumn{1}{|c}{} & & SE & \\multicolumn{1}{c|}{$%0.3f$} & \\multicolumn{1}{c|}{$%0.3f$} & $%0.3f$ \\\\ \\cline{3-6} ",
          Result.PT[2,2], 
          Result.PT[1,2], 
          Result.PT[3,2]),
  sprintf("\\multicolumn{1}{|c}{} & & 95\\%% CI & \\multicolumn{1}{c|}{$(%0.3f,%0.3f)$} & \\multicolumn{1}{c|}{$(%0.3f,%0.3f)$} & $(%0.3f,%0.3f)$ \\\\ \\hline ",
          Result.PT[2,3], Result.PT[2,4],
          Result.PT[1,3], Result.PT[1,4],
          Result.PT[3,3], Result.PT[3,4]) )

print(data.frame(RRR),row.names=F)

################################################################################
# Over-identification
################################################################################

RRR <- read.csv("NP/Result_Merge_NP.csv")

# H0: 2013 is a valid NCO = [2014] - [2013,2014]
VV <- RRR[,2]-RRR[,4]
SS <- (RRR$VAR33 + RRR$VAR11 - 2*RRR$VAR13)^(0.5)
Eff <- median(RRR[,2])-median(RRR[,4])
SE  <- median((VV - Eff)^2 + SS^2)^(0.5)
round(Eff,3); round(SE,3); round(Eff/SE,3)

# H0: 2014 is a valid NCO = [2013] - [2013,2014]
VV <- RRR[,3]-RRR[,4]
SS <- (RRR$VAR33 + RRR$VAR22 - 2*RRR$VAR23)^(0.5)
Eff <- median(RRR[,3])-median(RRR[,4])
SE  <- median((VV - Eff)^2 + SS^2)^(0.5)
round(Eff,3); round(SE,3); round(Eff/SE,3)

PP <- function(V){
  c(round(V$Diff,3),
    round(V$SE,3),
    round(V$Diff/V$SE,3))
}

cbind( PP(Result.COCA.Para.PS.Compare.23),
       PP(Result.COCA.Para.PS.Compare.13))

cbind( PP(Result.COCA.Para.OR.Compare.23),
       PP(Result.COCA.Para.OR.Compare.13))

cbind( PP(Result.COCA.Para.DR.Compare.23),
       PP(Result.COCA.Para.DR.Compare.13))


