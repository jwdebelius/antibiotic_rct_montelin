Code for **Gut microbiome perturbations following intake of nitrofurantoin, pivmecillinam, co-trimoxazole, ceftibuten, and ciprofloxacin; a randomized controlled trial with healthy adults** by Montelin, Debelius et al.

# Overview

In this study, we performed a randomized control trial of the effect of five commonly used antibiotics on the gut microbiome of 86 participants. After filtering and sample lost, we retained 473 samples from 74 individuals. These notebooks outline processing for metagenomic data.

Raw data is avalaible from ENA under accession PRJEB105465.
The origingal data was processed through the StaG pipeline, avalaible [here](https://github.com/ctmrbio/stag-mwc).

# Run enviroment

The conda enviroment can be found under conda-env.txt. Analysis is based on a [qiime2-2025.7 amplicon](https://library.qiime2.org/quickstart/qiime2) enviroment, with the [q2-boots](https://library.qiime2.org/plugins/caporaso-lab/q2-boots) and [gemelli](https://github.com/biocore/gemelli) plugins. Jupyter lab was also installed. The [emmeans](https://cran.r-project.org/web/packages/emmeans/index.html) package was installed in R and is required to run. This should be installed by opening R **inside** the QIIME 2 conda enviroment.

Figure 3 requires teh use of q2-empress, which can be installed in an earlier QIIME 2 release. See the [empress](https://github.com/biocore/empress) plugin for more details. 


# Notebooks

* **A.BuildDataset.ipynb**<br>Takes metadata and tables, imports into QIIME 2, and performs diversity calculations
* **B.TrackSampleLoss.ipynb**<br>Determines the number of samples retained at each step and raesons for exclusion
	- Figure S1. Participant and sample inclusion through sample processing.
		* `output/figures/figure_s1_sample_keep.png`
		* `output/figures/figure_s1_sample_keep.pdf`
* **C.Baseline.ipynb**<br>Plots the taxonomic composition at baseline (before treatment)
	-  Figure S2. Taxonomic microbiome composition at baseline.
		* `output/figures/figure_s2_taxa_baseline.png`
		* `output/figures/figure_s2_taxa_baseline.pdf`
* **D.AMR_baseline.ipynb**<br>Looks at the baseline AMR levels
	- Figure S3. Resistome composition at baseline.
		* `output/figures/rigure_S3_amr_baseline.png`
		* `output/figures/rigure_S3_amr_baseline.pdf`
* **E.IndividualTaxaTraces.ipynb**<br>Calculates the individual taxonomic traces at the phylum and family levels for participants.
	- Figure S4. Stacked bar plot of individual microbiome composition.
		* `output/figures/figure_S4_indv_taxa.png`
		* `output/figures/figure_S4_indv_taxa.pdf`
* **F.WithinIndividualBeta.ipynb**<br>Generates the orbital plot showing within individual beta diversity traces
	- Figure S5. Antibiotic treatment is insufficient to perturb individual microbial identity.
		* `output/figures/figure_s5_intraindividual.png`
        * `output/figures/figure_s5_intraindividual.pdf`
* **G.AlphaFits.ipynb**<br>Fits all three alpha diveristy metrics
	- Figure S6. Change in alpha diversity over time among participants.
         * `output/figures/fig_s6_alpha_trace.png`
         * `output/figures/fig_s6_alpha_trace.pdf`
	- Figure S7. Individual effects on the change in alpha diversity in response to antibiotic treatment
        * `output/figures/fig_s7_alpha_fit.png`
        * `output/figures/fig_s7_alpha_fit.pdf`
	- Table S3. Post-hoc comparisons of alpha diversity recovery after antibiotic treatment
        * `output/tables/table_s3_alpha_posthoc.tsv`
* **H.BetaFits.ipynb**<br>Beta Diveristy traces by time and treatment
	- Figure S8. Average change in beta diversity following antibiotic treatment.
        * `output/figures/figure_s8_beta_traces.png`
        * `output/figures/figure_s8_beta_traces.pdf`
	- Figure S9. Change in beta diversity following antibiotic treatment by period.
        - `output/figures/figure_s9_beta_boxes.png`
        - `output/figures/figure_s9_beta_boxes.pdf`
	- Figure S10. Effect of antibiotic treatment and time since intervention on beta diversity.
        - `output/figures/figure_s10_beta_forest.png`
        - `output/figures/figure_s10_beta_forest.pdf`
	- Table S4. Post-hoc comparisons of beta diversity recovery after antibiotic treatment.
        - `output/tables/table_s4_beta_posthoc.tsv`
* **I.Figure2.ipynb**<br>Main text figure 2 with alpha diveristy, beta diveristy
	- Figure 2
        - `output/figures/figure2.png`
        - `output/figures/figure2.pdf`
* **J.Targeted.ipynb**<br>Targeted taxonomic diversity
	- Table S5. Targeted differential abundance results for organisms of interest.
* **K.L.ANCOM-BC.ipynb**<br>Differential Abunandance
	- Figure 3. (*This requires Empress for tree generation and Illustrator for figure assumby*)
        - `output/figures/Figure_3_ancom_bc_rough.pdf`
        - To accompany this, the notebook generates `final-data/07.ancom-bc/emperor_meta.tsv` which is then processed with the tree through the Empress plugin. Empress screenshots are filled in the boxes on the rough data for the final figure.
	- Figure S11. Forest plot of rate of recovery of species with a statistically significant difference in ANCOM-BC.
        -  
	- Table S6. Untargeted differential abundance of species recovery following antibiotic treatment
    - 
* **L.AMR.ipynb**<br>AMR carraige
	- Figure S12. Portion of reads mapped to antimicrobial resistance genes.
	- Figure S13. Relative risk of carrying AMR pathways following treatment.
    - Table S7. Differentially prevalent antimicrobial resistance pathways across recovery periods

