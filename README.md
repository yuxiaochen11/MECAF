# MECAF
Differential abundance test
In this study, we propose a novel test for identifying the difference between two high-dimensional microbiome abundance data matrices based on the centred log-ratio transformation of the microbiome compositions. The test p-value can be calculated directly with a closed-form solution from the derived asymptotic null distribution. We also investigate the asymptotic statistical  power against sparse alternatives which are typically encountered in microbiome studies. The proposed Maximum-type test is Equal-Covariance-Assumption-Free, making it widely applicable to studies that compare microbiome compositions between conditions. Our simulation studies demonstrated that the proposed MECAF test achieves desirable power than competing methods while having the type I error rate well controlled under various scenarios.The usefulness of the proposed test is further illustrated with two real microbiome data analyses.


In the file 'Simulation_1.R', we calculate the type I error and power of two compared groups with different covariance matrix.

In the file 'Simulation_2.R', we calculate the type I error and power of two compared groups with same banded covariance matrix.

In the file 'Simulation_3.R', we calculate the type I error and power of two compared groups with same sparse covariance matrix.

In the file 'Real data_1.R', we make differencial abundances analysis of the CDI microbiome study.

In the file 'Real data_2.R', we calculate the type I error of two groups by procanova.

In the file 'Real data_3.R', we make differencial abundances analysis in female and male mice of the murine microbiome study.
