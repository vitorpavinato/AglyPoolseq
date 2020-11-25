##

### Description:

This analysis examines the abundance and geographic distribution of polymorphisms relative to their frequency.

### Analysis steps
1. Subset VCF down to non-problematic populations. First, pull data from google sheet an reformat for bcftools input.
  `DEST/Analyses/GeographicEndemism/pullProblematicPops.R`

2. Run awk script to extract out SNP-level metrics for polymorphic populations
  `DEST/Analyses/GeographicEndemism/get_geoEndemic.SNAPE.sh`
  `DEST/Analyses/GeographicEndemism/get_geoEndemic.PoolSNP.sh`

3. Parse output
  `DEST/Analyses/GeographicEndemism/geographic_clustering_endemic.sh`

4. Collect and summarize
  `DEST/Analyses/GeographicEndemism/geographic_clustering_endemic.collect.sh`

5. Make object for plotting
  `DEST/Analyses/GeographicEndemism/geographic_clustering_endemic.make_plot_data.R`

6. Make plot
  `DEST/Analyses/GeographicEndemism/geographic_clustering_endemic.plot.R`
