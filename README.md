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
        
  - Install stable release
        
        devtools::install_github("jzhangc/git_RBioFS/RBioFS", repos = BiocInstaller::biocinstallRepos())
        
  - Install development build
        
        devtools::install_github("jzhangc/git_RBioFS/RBioFS", repos = BiocInstaller::biocinstallRepos(), ref = "beta")
        

Change log

    0.5.2 (Jun. 6. 2018)
        - New PLS-DA functions added (non-Shiny):
          - rbioFS_plsda_predict(): use the plsda model to calcualte predicted values for unknown data.
          - rbioFS_plsda_classification(): use the predicted values (produced by rbioFS_plsda_predict) to classify. Note: current probability method is "softmax". A "Bayesian" method will be added later. 
          
        - Updates to RF-FS functions:
          - Plotting module separated from the functions
          - Boxplot for the initial_FS object now horizontal
          - Plot file suffix for both VI boxplot and OOB plot now ".rffs.ifs.plot.pdf" and ".rffs.sfs.plot.pdf", respectively
          - Classes "rf_ifs" and "rf_sfs" created for the output of rbioRF_initial_FS() and rbioRF_SFS(), respectively
          - Display messages added for the functions
          - plots title and axis title are in bold
          
        - RF-FS analysis plotting module separated from rbioFS() function as functions:
          - rbioRF_initialFS_plot() 
          - rbioRF_SFS_plot()
          
        - Updates to PLS-DA functions (non-Shiny):
          - Smooth functionality added for rbioFS_plsda_roc_auc()
          - Legend position now customizable for multi-plot for the relevant functions
          - rbioFS_plsda_jackknife(): options added for hiding the x-axis tick labels (useful in the case of many variables)
          - rbioFS_plsda_jackknife(): for plotting, x-axis margin adjusted
          - rbioFS_plsda_VIP(): options added for hiding the x-axis tick labels (useful in the case of many variables)
          - rbioFS_plsda_VIP() how outputs a list with VIP values as well as a vector containing features above the threshold
          - Error handling added for all the functions featuring multiplot, excluding the rbioFS_plsda_scoreplot()
          - Plot theme adjusted for functions:
            - rbioFS_plsda_ncomp_select()
            - rbioFS_plsda_tuplot()
          - Output (to R environment) object name suffix adjusted for all the relevant functions with added "_plsda"
          - xLabelSize and yLabelSize added for the functions
          - rbioFS_plsda_scoreplot(): sample labeling functionality added
          - Size option added for ggrepel label for the relavent functions


    0.5.1 (May.24.2018)
        - New PLS-DA functions added (non-Shiny):
          - rbioFS_plsda_VIP(): VIP, or variable importance in projection, is plsda's version of VI. Can be used independently from plsda functions
          - rbioFS_plsda_q2r2(): Q2-R2 calculation and plotting
          - rbioFS_plsda_aoc_auc(): ROC and AUC analysis and plotting
          
        - Updates to PLS-DA functions (non-Shiny):
          - More information added to the manual page for rbioFS_plsda_ncomp_select()
          - Additional arugment checks added to rbioFS_plsda_ncomp_select()
          - A bug fixed where rbioFS_plsda_jackknife() fails if no coefs are > (or <) 0
          - Small changes made to message display pattern in rbioFS_plsda_jackknife()
          - Plot property auguments names unified for function only produce one type of plot
        
        - Updates to RF-FS functions:
          - rbioFS() now accepts R objects, in addition to csv files
          - Group variable now customizable for rbioFS()
          - Arugment check added for rbioFS()
          - rbioFS() output element "SFS_matrix" changed to "SFS_training_data_matrix"
          - Function message display feature added to rbioFS()
        
        - Bug fixes
        

    0.5.0 
        - Data preprocessing functions added for modelling precedures such as PLS-DA, sPLS-DA, PCA, SVM, etc.
          - center_scale()
          - dummy()
          
        - PLS-DA functions added (non-Shiny):
          - rbioFS_plsda()
          - rbioFS_plsda_ncomp_select()
          - rbioFS_plsda_tuplot()
          - rbioFS_plsda_scoreplot()
          - rbioFS_plsda_jackknife()
          
        - rbioFS_PCA() re-written with the follwoing new functionalities:
          - The function now outputs a PCA object to the environment
          - PCA scoreplot now supports single component curve
          - PCA scoreplot now supports paired matrix, i.e. more than two components
          - PCA boxplot y upper limited adjusted for both shiny and non-shiny versions
          - PCA score plot now supports sample names for the samples
          
        - Rightside y-axis now uses a function from RBioplot pakcage, which now is a dependency
        
        - Bug fixes
        

    0.4.6
        - All the settings return to default upon "clear"
        - A bug fixed for rbioFS_app() where plots can't be regenerated on a new dataset upon "clear"
        - Other bug fixes
        

    0.4.4 - 0.4.5
        - fs_csv_generator() added
        - Bug fixes
        

    0.4.3
        - Web app verion of rbioFS_PCA() added: rbioFS_PCA_app()
        - Code update for rbioFS_PCA() for better input data compatibility
        - Parallel computing functionality added for rbioFS_app()
        - Quantile normalization functionality added for rbioFS_app()
        - Clear screen button added for rbioFS_app()
        - Parallel computing modules updated with the more efficient foreach method for rbioFS(), rbioRF_initialFS() and rbioRF_SFS()
        - Bug fixes
        

    0.4.1 - 0.4.2
        - Progress bar added for rbioFS_app()
        - Plotting buttons added for rbioFS_app() 
        - UI elements re-arranged for a better presentation for rbioFS_app()
        - Small icons added for the buttons for rbioFS_app()
        - "Summary" tabs re-labelled as "Results" tabs for rbioFS_app()
        - NA check added for the initial FS module in rbioFS_app()
        - Bug fixes for rbioFS_app()
        - Other minor bug fixes
        

    0.4.0
        - Web app version of the main fuction rbioFS() added: rbioFS_app()
        - A bug fixed for the plot subsetting functionality for rbioRF_SFS()
        - Other minor bug fixes
        

    0.3.3
        - zzz.R file added
        
  
    0.3.0 - 0.3.2
        - Principal Component Analysis (PCA) and visualization function rbioFS_PCA() added
        - Citation information added
        - Bug fixes
        
    
    0.2.0 - 0.2.5
        - All-in-one FS function added
        - Bug fixes
        
    
    0.1.11 - 0.1.12
        - Output results as txt files functionality added
        - Bug fixes
        
    
    0.1.10
        - Random Forest data imputation method added to the data imputation function
        - Round down now used for mtry augment
        
    
    0.1.9
        - Text fixes
        
    
    0.1.7 - 0.1.8
        - Name changed to RBioFS
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
