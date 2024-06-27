# MetimgastDesign

Codes to replicate the results of the article "Design of a phase II trial for digestive cancer therapeutics"

Folders:

- R = old code without cleaning (I let it here if you are curious)
- R script = codes to replicate the article's results:
    - functions.R = functions used in the other scripts
    - metadata.R = data used in other scripts and centralized in this script (for example, the scenarios)
    - thresholds.R = finding the optimal thresholds
    - thresholds_postreview.R = finding the optimal thresholds with about 5% type I error rate as asked by reviewer to have comparable designs
    - resultats_these.R = simulate the results of simulated trials
    - figures.R = generate the figures in the paper
    - simu_corr_data.R = how we simulated the correlated data
- data = folder with the stored results (it takes time to replicate the simulations)
