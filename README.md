# microbiome_randomisations
randomisation algorithms for bacterial communities from a known source

The R script randomises bacterial communities isolated from a source using (1) a swapping algorithm and a (2) lottery system whereby the identity of the bacteria is exchanged with one present in the source community. The randomisations are stratified within communities that belong to the same source (constrained by source). 
The swapping algorithm preserves both the matrix fill (the number of zero occurrences) and the column and row abundance totals. The lottery system only introduces bacteria from the source into the isolated communities if those bacteria have at least once have been observed in the isolated communities (i.e. if bacteria can survive in the conditions that the isolated communities are exposed to). The chance that a bacterial is selected from the source is proportional to its abundance in the appropriate source community.
Code for the comparison of the different communities after randomisation is provided. 
Sample data is provided to illustrate the workings of the code. More information on the data is provided in the R script and in the manuscript (submitted)
