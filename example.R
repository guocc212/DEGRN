# exmpale 
# Notes:
# If you want to draw the network for DEGRN, you can do as followed:


# Step 1. load the scripts of functions for DEGRN
source("plot_network_for_DEGRN.R")

# Step 2. choose the interested gene functions from DEGRN.
# for example: response to chitin
# we choose the GO, named as "response to chitin", to show the potential network predicted from the DEGRN. 
# Firstly, we select the mode == "go", and add the go="response to chitin" with the parameter of 30 NOVEL TFs and top 30 of the most enriched target genes.
DEGRN_plot_pipeline(modes = "go", go="response to chitin", tf.number = 30, target.number = 30)

# Step 3. choose the interested TFs for novel gene functions from DEGRN.
# for example: AT3G23250 (ATMYB15)
# we choose the TF ATMYB15 to show the potential network predicted from the DEGRN. 
# Firstly, we select the mode == "go", and add the go="response to chitin" with the parameter of 30 NOVEL TFs and top 30 of the most enriched target genes.
DEGRN_plot_pipeline(modes = "gene", gene="AT3G23250", go.number = 30, target.number = 30)
