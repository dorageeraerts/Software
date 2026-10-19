# PEDESTAL and 1PHE VALUES

- These values were provided by Gábor, the scripts are in his possession. Contact: Gábor Nyitrai <nyiti28@gmail.com>. The plot for each channels and telescope of this analysis are availabe on the MURAVES Team on Microsoft Teams (https://uclouvain.sharepoint.com/:f:/s/O365G-MURAVESTaskforce/IgD_j0vG1JEgSJjxvK0yWCKWAdzogIsc1Mx18P2Dd-R6_T4?e=miZJjK).
    - The values a retrived by building a high statistic pedestal histogram called banchmark and fit each single channel with a function composed buy a line (noise) + 4 gaussians at *fixed distance* one another, to fit photoelectron peaks.
    NB: for some channels the final 1phe value is 'maunally' tuned so something more appropriate then the result of the fit. These values are not 100% reproduciable.
    - Some of the values are NaN or 0:
        - NaN pedestal/1phe values are replaced with the average of the other
          (non-NaN) values for that column in the same input file
        - A 1phe value of 0 is replaced with 1000 and flagged (flag=1), since it
          usually indicates a dead channel.
- The values provided are not mapped according the real logic (`/user/abiolchi/Software/muraves/tracks_reconstruction/muraves_cfg_files/spiroc-hybrid-map.cfg`). They are just in the same order as saved on the raw file. 
*NB: The script PedestalReader.py do not perfom the mapping, because it is performed directly (and only) before the clustering, in the main reconstruction script*.
- Gábor checked that the provided pedestal and 1phe are reasonable for most of the run ~95%. The scripts, the banchmark dataset and the results of this analysis are available on MURAVES Team on Microsoft Teams (https://uclouvain.sharepoint.com/:f:/s/O365G-MURAVESTaskforce/IgD_j0vG1JEgSJjxvK0yWCKWAdzogIsc1Mx18P2Dd-R6_T4?e=miZJjK)
