#mammal 
mtree2 <- read.nexus('/Users/sunny/Desktop/project/1.mammal/tree.with.length/mtreebl.tre')
mtree <- read.nexus('project/1.mammal/tree.no.length/mtree.tre')

m_stats <- read.csv('/Users/sunny/Desktop/project/1.mammal/stats/mammalcstats.csv')
mtree2 <- multi2di(mtree2) #change the tree in to binary tree,resolve polytomy
mchar <- ReadCharacters('/Users/sunny/Desktop/project/1.mammal/analysisdata/mammal-character.nex')
unique(as.vector(mchar))
l <- subtrees(mtree2,wait = F)
#1] "0"    "?"    "1"    "(01)" "2"    "-"    "3"    "(12)" "4"    "(23)" "5"    "(34)" "(02)" "(03)" "6"   
# "7"    "(15)"

#seal
setwd('/Users/Sunny/Desktop')
stree <- read.nexus('project/seal/tress/stree.tre')
schar <- ReadCharacters('/Users/sunny/Desktop/project/seal/character.matrix/sealchar.nex')
unique(as.vector(schar)) #"1"    "?"    "0"    "2"    "-"    "(01)" "(1)"  "3"    "(12)"
stree <- multi2di(stree)
s_stats <- read.csv('/Users/Sunny/Desktop/project/seal/stats.seal.csv')
s <- subtrees(stree,wait = F)

#bird-dinosaur
btree <-  read.nexus('project/2.bird-dinosaur/trees/birdtreeble.tre')
btree <- multi2di(btree)
bchar <- ReadCharacters ('/Users/sunny/Desktop/project/2.bird-dinosaur/matrix/bird.matrix2.nex')
b_stats <- read.csv('project/2.bird-dinosaur/birdstats.csv')
unique(as.vector(bchar))
b <- subtrees(btree,wait = F)

#tetrapod
ttree <- read.nexus('project/tetrapod/tetrapod.tre')
ttree <- multi2di(ttree)
tchar <- ReadCharacters('project/tetrapod/tetrapod.nex')
t_stats <- read.csv('project/tetrapod/t.stats.csv')
unique(as.vector(tchar))
t <- subtrees(ttree,wait = F)

#ganthostome
gtree <- read.nexus('project/gnathostome/gtree2.tre')
gtree <- multi2di(gtree)
gchar <- ReadCharacters('project/gnathostome/gtree.nex')
unique(as.vector(tchar))
g_stats <- read.csv('project/gnathostome/g.stats.csv')
g <- subtrees(gtree, wait = F)

#convert character states into phydat 
char_phyDat <- phyDat(mchar, type = "USER", contrast=contrast, return.index = TRUE)
char_phdDat_s <- phyDat(schar,type = "USER", contrast=contrast,return.index = TRUE )
char_phyDat_b <- phyDat(bchar, type = "USER", contrast=contrast,return.index = TRUE)
char_phyDat_t <-phyDat(tchar, type = "USER", contrast=contrast,return.index = TRUE)
char_phyDat_g <- phyDat(gchar, type = "USER", contrast=contrast,return.index = TRUE)
#then get the ancestral states data
ances <- ancestral.pars(mtree2,char_phyDat , type = "MPR")
ances_s <- ancestral.pars(stree,char_phdDat_s , type = "MPR" )
ances_b <- ancestral.pars(btree, char_phyDat_b , type = "MPR")
ances_t <- ancestral.pars(ttree, char_phyDat_t, type = "MPR")
ances_g <- ancestral.pars(gtree, char_phyDat_g, type = "MPR")
# Number of characters
num_characters <- ncol(mchar)  # Adjust to 583 characters
num_characters_seal <- ncol(schar)
num_characters_bird <- ncol(bchar)
num_char_t <- ncol(tchar)
num_char_g <- ncol(gchar)
####finding synapomorphy:

#char_changes <- c(char_changes, sprintf('%s -> %s', parent_label, child_label))
test <- process_ancestral_data(mtree2, ances, num_characters,char_phyDat)
test_bird <- process_ancestral_data(btree, ances_b, num_characters_bird,char_phyDat_b)
test_seal <- process_ancestral_data(stree, ances_s, num_characters_seal,char_phdDat_s)
test_t <- process_ancestral_data(ttree, ances_t, num_char_t,char_phyDat_t)
test_g <- process_ancestral_data(gtree, ances_g, num_char_g,char_phyDat_g)
#for bird, node label might need to be changed


#######getting the goodness of fit at branch level
#1. get the subtrees info from list of subtrees, then convert into dataframe
expanded_m <- combined_function(l)
expanded_seal <- combined_function(s)
expanded_b <- combined_function(b)
expanded_t <- combined_function(t)
expanded_g <- combined_function(g)

# 2. Function to extract nodes, synapomorphy and create a dataframe from 'test'
#datalist is text freviously, move up to test! can be combined
#merge fitstas based on character
result_list <- list_to_list(test)
#result_test <- list_to_list(testm)
test_list <- list_to_list(testm)
result_s <- list_to_list(test_seal)
result_b <- list_to_list(test_bird)
result_t <- list_to_list(test_t)
result_g <- list_to_list(test_g)

#test_df <- stack(result_test)%>% as.data.frame() %>%setNames(c("Character", "node")) %>%merge(m_stats, by = "Character")

result_m <- stack(result_list) %>%
  as.data.frame() %>%
  setNames(c("Character", "node")) %>%
  merge(m_stats, by = "Character")

result_s <- stack(result_s) %>%
  as.data.frame() %>%
  setNames(c("Character", "node")) %>%
  merge(s_stats, by = "Character")

result_b <- stack(result_b) %>%
  as.data.frame() %>%
  setNames(c("Character", "node")) %>%
  merge(b_stats, by = "Character")

result_t <- stack(result_t) %>%
  as.data.frame() %>%
  setNames(c("Character", "node")) %>%
  merge(t_stats, by = "Character")

result_g <- stack(result_g) %>%
  as.data.frame() %>%
  setNames(c("Character", "node")) %>%
  merge(g_stats, by = "Character")
# Print the resulting DataFrame
#test_df <- summarise_data(expanded_m,test_df )
summary_df <- summarise_data(expanded_m, result_m)
summary_s <- summarise_data(expanded_seal, result_s)
summary_b <- summarise_data(expanded_b, result_b)
summary_t <- summarise_data(expanded_t, result_t)
summary_g <- summarise_data(expanded_g,result_g)


#####Linear model
#character shapes the tree, so the imbalance is dependent 
# Standardize the variables using scale()
GI_standardized <- scale(summary_df$avg_Gfit)
RI_standardized <- scale(summary_df$avg_RI)
imbalance <- scale(summary_df$imbalance)
#perform correlation test:
cor_m <- cor.test(summary_df$avg_Gfit, summary_df$sum_RI, method=c('spearman'))
print(cor_m)#0.72

# Perform linear regression with standardized variables
filtered_df <- summary_df %>% filter(imbalance != 0 & imbalance != 1)

modelm <- lm( GI_standardized ~ imbalance, data = summary_df) #both not
modeltest <- lm(avg_Gfit ~ imbalance, data = test_df)
modelm2 <- lm( avg_RI~ imbalance , data = summary_df )
summary(modelm) #0.19
summary(modelm2) #0.13

plot(summary_df$imbalance,summary_df$avg_Gfit , 
     pch = 19, 
     col = "blue", 
     xlab = "imbalance", 
     ylab = "G fit", 
     main = "evolution of mammal ear")
abline(modelm, col = "darkolivegreen3", lwd = 2)

plot(summary_df$imbalance, summary_df$avg_RI, 
     pch = 19, 
     col = "darkred", 
     xlab = "imbalance", 
     ylab = "RI", 
     main = "evolution of mammal ear")
# Add the regression line
abline(modelm2, col = "darkolivegreen3", lwd = 2)
par(mfrow = c(2, 2))  # Arrange the plots in a 2x2 layout
plot(modelm)
# Obtain residuals and fitted values
residuals <- residuals(modelm)
fitted_values <- fitted(modelm)


GI_standardized <- scale(summary_s$avg_Gfit)
RI_standardized <- scale(summary_s$avg_RI)
cor_s <- cor.test(summary_s$avg_Gfit, summary_s$sum_RI, method=c('spearman'))
print(cor_s) #0.434
models <- lm(avg_GI ~ imbalance , data = summary_s) #not 
summary(models) #0.78
models2 <- lm( avg_RI ~ imbalance ,data = summary_s) #sig
summary(models2) #0.178

plot(summary_s$imbalance, summary_s$avg_Gfit, 
     pch = 19, 
     col = "blue", 
     xlab = "imbalance", 
     ylab = "G fit", 
     main = "evolution of seal")
abline(models, col = "darkolivegreen3", lwd = 2)

plot(summary_s$imbalance, summary_s$avg_RI, 
     pch = 19, 
     col = "darkred", 
     xlab = "imbalance", 
     ylab = "RI", 
     main = "evolution of seal")
# Add the regression line
abline(models2, col = "darkolivegreen3", lwd = 2)

###bird dino
cor_s <- cor.test(summary_s$avg_Gfit, summary_s$sum_RI, method=c('spearman'))
print(cor_s) #0,435
modelb <- lm(avg_GI ~ imbalance , data = summary_b)
summary(modelb) #0.445
plot(summary_b$imbalance, summary_b$avg_Gfit, 
     pch = 19, 
     col = "blue", 
     xlab = "imbalance", 
     ylab = "G fit", 
     main = "dinosaur to bird transition")
abline(modelb, col = "darkolivegreen3", lwd = 2)
modelb2 <- lm( avg_RI ~ imbalance ,data = summary_b)
summary(modelb2) #0.04
plot(summary_b$imbalance, summary_b$avg_RI, 
     pch = 19, 
     col = "darkred", 
     xlab = "imbalance", 
     ylab = "RI", 
     main = "dinosaur to bird transition")
abline(modelb2, col = "darkolivegreen3", lwd = 2)

####tetrapod 
cor_t <- cor.test(summary_t$avg_Gfit, summary_t$sum_RI, method = "spearman")
print(cor_t)#0.73

modelt <- lm(avg_Gfit ~ imbalance , data = summary_t)
summary(modelt)  #0.0169
plot(summary_t$imbalance, summary_t$avg_Gfit, 
     pch = 19, 
     col = "blue", 
     xlab = "imbalance", 
     ylab = "G fit", 
     main = "tetrpaod evolution")
abline(modelt, col = "darkolivegreen3", lwd = 2)
modelt2 <- lm(avg_RI ~ imbalance , data = summary_t)
summary(modelt2) #0.042
plot(summary_t$imbalance, summary_t$avg_RI, 
     pch = 19, 
     col = "darkred", 
     xlab = "imbalance", 
     ylab = "RI", 
     main = "tetrpaod evolution")
abline(modelt2, col = "darkolivegreen3", lwd = 2)

#####gnathostome
cor_g <- cor.test(summary_g$avg_Gfit, summary_g$sum_RI, method=c('kendall'))
print(cor_g) #0.400
modelg <- lm(avg_Gfit ~ imbalance , data = summary_g)
summary(modelg)  ##0.0382
plot(summary_g$imbalance, summary_g$avg_Gfit, 
     pch = 19, 
     col = "blue", 
     xlab = "imbalance", 
     ylab = "G fit", 
     main = "gnathosome evolution")
abline(modelg, col = "darkolivegreen3", lwd = 2)
modelg2 <- lm(avg_RI ~ imbalance , data = summary_g)
summary(modelg2) #0.0275
plot(summary_b$imbalance, summary_b$avg_RI, 
     pch = 19, 
     col = "darkred", 
     xlab = "imbalance", 
     ylab = "RI", 
     main = "gnathosome evolution")
abline(modelg2, col = "darkolivegreen3", lwd = 2)


####creat plot of Gfit and ri on tree
#1.mammal get base node lable of each subtree
node_m <- expanded_m %>%
  group_by(subtrees) %>%
  summarise(nodes = first(node_labels)) %>% mutate(parent = as.integer(nodes))

node_m$Gfit <-  round(summary_df$avg_Gfit, 2)
node_m$RI <- round(summary_df$avg_RI, 2)

tree_data <- fortify(mtree2)
mtree2_plot <- merge(tree_data, node_m, by = "parent", all.x = TRUE)
mtree2_plot <- as.treedata(mtree2_plot)

ggtree(mtree2_plot, aes(color = Gfit),size = 1) +
  scale_color_continuous(low = 'bisque2', high ='tomato4' ) +
  geom_text2(aes(subset = !isTip, label = Gfit), 
             size = 2, hjust = 1, vjust = -2) +
  geom_tiplab(size = 2, align = F) +
  theme(legend.position = "right")

ggtree(mtree2_plot, aes(color = RI),size = 1) +
  scale_color_continuous(low = 'lightskyblue1', high ='darkolivegreen' ) +
  geom_text2(aes(subset = !isTip, label = RI), 
             size = 2, hjust = 1, vjust = -2) +
  geom_tiplab(size = 2, align = F)+
  theme(legend.position = "right")

#2.bird
node_b <- expanded_b%>%
  group_by(subtrees) %>%
  summarise(nodes = first(node_labels)) %>% mutate(parent = as.integer(nodes))
node_b$Gfit <- round(summary_b$avg_Gfit, 2)
node_b$RI <- round(summary_b$avg_Gfit, 2)
btree_data <- fortify(btree)
bplot <- merge(btree_data, node_b, by = "parent", all.x = TRUE)

ggtree(btree, aes(color = Gfit), size = 1) %<+% bplot + 
  scale_color_continuous(low = 'bisque2', high = 'tomato4') +
  geom_tiplab(size = 2, align = F) +
  geom_text2(aes(subset = !isTip, label = Gfit), 
             size = 2, hjust =1, vjust = 1) +
  theme(legend.position = "right")

ggtree(bplot, aes(color = RI),size = 1) +
  scale_color_continuous(low = 'lightskyblue1', high ='darkolivegreen' ) +
  geom_tiplab(size = 2, align = F) +
  geom_text2(aes(subset = !isTip, label = Gfit), 
             size = 2, hjust = 1, vjust = -0.7) +
  theme(legend.position = "right")

#3.seal
node_s <- expanded_seal%>%
  group_by(subtrees) %>%
  summarise(nodes = first(node_labels)) %>% mutate(parent = as.integer(nodes))
node_s$Gfit <- round(summary_s$avg_Gfit, 2)
node_s$RI <- round(summary_s$avg_RI, 2)
stree_data <- fortify(stree)
splot <- merge(stree_data, node_s, by = "parent", all.x = TRUE)

ggtree(splot, aes(color = Gfit), size = 1) + 
  scale_color_continuous(low = 'bisque2', high = 'tomato4') +
  geom_tiplab(size = 2, align = F) +
  geom_text2(aes(subset = !isTip, label = Gfit), 
             size = 2, hjust = -0.3, vjust = 0.3) +
  theme(legend.position = "right")

ggtree(splot, aes(color = RI), size = 1) +
  scale_color_continuous(low = 'lightskyblue1', high = 'darkolivegreen') +
  geom_tiplab(size = 2, align = F) +
  geom_text2(aes(subset = !isTip, label = RI), 
             size = 2, hjust = -0.3, vjust = 0.3) +  # Adjust size and position as needed
  theme(legend.position = "right")

#4.gnathostome
node_g <- expanded_g%>%
  group_by(subtrees) %>%
  summarise(nodes = first(node_labels)) %>% mutate(parent = as.integer(nodes))
node_g$Gfit <- round(summary_g$avg_Gfit, 2)
node_g$RI <- round(summary_g$avg_RI, 2)

gtree_data <- fortify(gtree)
gplot <- merge(gtree_data, node_g, by = "parent", all.x = TRUE)

ggtree(gplot, aes(color = Gfit), size = 1) + 
  scale_color_continuous(low = 'bisque2', high = 'tomato4') +
  geom_tiplab(size = 2, align = F) +
  geom_text2(aes(subset = !isTip, label = Gfit), 
             size = 2.2, hjust = 1, vjust = -1) +
  theme(legend.position = "right")

ggtree(gplot, aes(color = RI),size = 1) +
  scale_color_continuous(low = 'lightskyblue1', high ='darkolivegreen' ) +
  geom_text2(aes(subset = !isTip, label = RI), 
             size = 2.4, hjust = 1, vjust = 1) +
  geom_tiplab(size = 2, align = F) +
  theme(legend.position = "right")

#5.tetrapod
node_t <- expanded_t%>%
  group_by(subtrees) %>%
  summarise(nodes = first(node_labels)) %>% mutate(parent = as.integer(nodes))
node_t$Gfit <- round(summary_t$avg_Gfit, 2)
node_t$RI <- round(summary_t$avg_RI, 2)
ttree_data <- fortify(ttree)
tplot <- merge(ttree_data, node_t, by = "parent", all.x = TRUE)

ggtree(tplot, aes(color = Gfit), size = 1) + 
  scale_color_continuous(low = 'bisque2', high = 'tomato4') +
  geom_text2(aes(subset = !isTip, label = Gfit), 
             size = 2.4, hjust = 1.3, vjust = 0) +
  geom_tiplab(size = 2, align = F) +
  theme(legend.position = "right")

ggtree(tplot, aes(color = RI),size = 1) +
  scale_color_continuous(low = 'lightskyblue1', high ='darkolivegreen' ) +
  geom_text2(aes(subset = !isTip, label = RI), 
             size = 2.4, hjust = 1.2, vjust = 0) +
  geom_tiplab(size = 2, align = F) +
  theme(legend.position = "right")

#getting the completeness of the character state data
m.pro <- char_proportion(mchar) #0.870
b.pro <- char_proportion(bchar) #0.142
s.pro <-char_proportion(schar) #0.857
t.pro <- char_proportion(tchar) #0.163
g.char <- char_proportion(gchar) #0.633
