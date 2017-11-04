# git_RBioFS
A simple to use package for implementing machine learning feature selection in biological and biochemical research

To cite in publication
  
    Zhang J, Hadj-Moussa H, Storey KB. 2016. Current progress of high-throughput microRNA differential expression analysis and random forest gene selection for model and non-model systems: an R implementation. J Integr Bioinform. 13: 306.


Installation

  - Install devtools (if not already done)
  
        install.packages("devtools")
        
  - Install bioconductor (if not already done)
        
        source("https://bioconductor.org/biocLite.R")
      
        biocLite()
        
  - Install the package
        
        devtools::install_github("jzhangc/git_RBioFS/RBioFS", repos = BiocInstaller::biocinstallRepos())
        

Change log

    0.4.5 (Nov.3.2017)
        - fs_csv_generator() added
        - Bug fixes

    0.4.4
        - Bug fixes

    0.4.3
        - Web app verion of rbioFS_PCA() added: rbioFS_PCA_app()
        - Code update for rbioFS_PCA() for better input data compatibility
        - Parallel computing functionality added for rbioFS_app()
        - Quantile normalization functionality added for rbioFS_app()
        - Clear screen button added for rbioFS_app()
        - Parallel computing modules updated with the more efficient foreach method for rbioFS(), rbioRF_initialFS() and rbioRF_SFS()
        - Bug fixes

    0.4.2
        - Bug fixes for rbioFS_app()

    0.4.1
        - Progress bar added for rbioFS_app()
        - Plotting buttons added for rbioFS_app() 
        - UI elements re-arranged for a better presentation for rbioFS_app()
        - Small icons added for the buttons for rbioFS_app()
        - "Summary" tabs re-labelled as "Results" tabs for rbioFS_app()
        - NA check added for the initial FS module in rbioFS_app()
        - Other minor bug fixes

    0.4.0
        - Web app version of the main fuction rbioFS() added: rbioFS_app()
        - A bug fixed for the plot subsetting functionality for rbioRF_SFS()
        - Other minor bug fixes

    0.3.3
        - zzz.R file added
    
    0.3.1 - 0.3.2 
        - Bug fixes
  
    0.3.0
        - Principal Component Analysis (PCA) and visualization function rbioFS_PCA() added
        - Citation information added
        - Bug fixes
    
    0.2.1 - 0.2.5
        - Bug fixes
    
    0.2.0
        - All-in-one FS function added
        - Bug fixes
    
    0.1.12
        - Bug fixes
    
    0.1.11
        - Output results as txt files functionality added
    
    0.1.10
        - Random Forest data imputation method added to the data imputation function
        - Round down now used for mtry augment
    
    0.1.9
        - Text fixes
    
    0.1.8
        - Name changed to RBioFS
    
    0.1.7
        - Bug fixes
    
    0.1.6
        - rbioRF_iterOOB() updated
    
    0.1.5
        - Iterative OOB error rate computation function added
    
    0.1.4
        - Initial FS function completed
    
    0.1.3
        - Parallel computing added for rbioRF_vi()
    
    0.1.2
        - rbioRF_vi and rbioRF_viplot functions combined to steramline the workflow
    
    0.1.1
        - Data imputation function added
        - File processing functions added
    
    0.1.0
        - Initial commit
