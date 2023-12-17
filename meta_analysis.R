library(netmeta)
library(dplyr)
library(officer)

# make a dataframe
# study 1 to 2: RGH-MD-16 placebo(-13.29(1.82), n=148), 3mg(-21.45(1.74), n=140), risperidone(-29.27(1.74), n=138)
# study 3 to 4: RGH-MD-04 placebo(-14.3(1.5), n=149), 3mg(-20.2(1.5), n=151), 6mg(-23.0(1.5), n=154)
# study 5 to 7: A002-A4 placebo(-9.5(2.1), n=125), 3mg(-10.4(2.0), n=126), 6mg(-13.8(2.0), n=131), risperidone(-20.2(3.0), n=55)

meta_df <- data.frame(studlab = c("Study 1", "Study 2", "Study 3", "Study 4", "Study 5", "Study 6", "Study 7"),
                      t1 = c("3mg", "Risperidone", "3mg", "6mg", "3mg", "6mg", "Risperidone"),
                      y1 = c(-21.45, -29.27, -20.2, -23.0, -10.4, -13.8, -20.2),
                      sd1 = c(1.74, 1.74, 1.5, 1.5, 2.0, 2.0, 3.0),
                      n1 = c(140, 138, 151, 154, 126, 131, 55),
                      t2 = rep("Placebo", 7),
                      y2 = c(-13.29, -13.29, -14.3, -14.3, -9.5, -9.5, -9.5),
                      sd2 = c(1.82, 1.82, 1.5, 1.5, 2.1, 2.1, 2.1),
                      n2= c(148, 148, 149, 149, 125, 125, 125))


pw1 <- pairwise(treat = list(t1, t2), n = list(n1, n2), TE = list(y1, y2), seTE = list(sd1, sd2), studlab = studlab, data = meta_df)

net1 <- netmeta(TE, seTE, treat1, treat2, studlab,
                data = pw1, sm = "SMD", reference = "Placebo", baseline.reference = FALSE,
                comb.random = FALSE)

results <- summary(net1, digits = 3)
forexport <- results$x # result
forexport

forest(net1) # plotting
