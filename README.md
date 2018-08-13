# git_RBioFS
A comprehensive yet straightforward machine learning package for biological and biochemical research

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

    Feature blueprint
        - ANN (artificial neural network) functions added (non-Shiny) (tentative):
          - rbioClass_ann()
          - rbioClass_ann_plot()
        
        - KNN (k-nearest neighbors) functions added (non-Shiny) (tentative):
          - rbioClass_knn()
          - rbioClass_knn_plot()
        
        - sPLS-DA functions added (non-Shiny) (tentative):
          - rbioClass_splsda()
          - rbioClass_splsda_plot()
        
        - All shiny apps' interface updated with a new look
        
        - Updates to rbioFS_PCA():
          - 3D score plot option for rbioFS_PCA (non-shiny)
          - 3D socre plot option for rbioFS_PCA_app (shiny)
          - Legend style adjusted for the sample labels
          - PCA function rbioFS_PCA() updated with S3 method

        - ROC-AUC now a seperate function that can be used for other classification/FS methods
        

    0.6.1 (Feature preview)
      (ICEBOX)
        - All-in-one classification solution functions added:
          - rbioClass()
          
        - RF-class (random forest classification) functions added (non-Shiny) (tentative):
          - rbioClass_rf_class()
          - rbioClass_rf_class_plot()

        - New PLS-DA functions added (non-Shiny):
          - rbioFS_plsda_loadingplot(): with y loading as well
        
        - Update to SVM functions:
          - rbioClass_svm_perm() plot output file name fixed
          - S3 print method for relevant functions
        
        - Updates to PLS-DA functions:
          - Multi-category Y modelling now possible with rbioFS_plsda
          - rbioFS_plsda_vip_plot() plot output file name fixed
          - rbioClass_plsda_perm() plot output file name fixed
          - S3 print method for relevant functions
          
        - Updates to RF-FS functions:
          - S3 print method for relevant functions

          
      (ADDED)
        - New SVM functions:
          - rbioClass_svm_ncv_fs(): nested SVM cross-validation function with feature selection functionality
        
        - New PLS-DA functions:
          - rbioFS_plsda_vip_plot(): the function only accepts "rbiomvr_vip" class object
          
        - Updates to SVM functions (non-Shiny):
          - Additional argument checks
          
        - Updates to PLS-DA functions (non-Shiny):
          - rbioFS_plsda_VIP() changed to rbioFS_plsda_vip()
          - Bootstraping option added for rbioFS_plsda_vip() so that VIP can use bootstrap data for SD/SEM errorbars
          - rbioFS_plsda_vip() now outputs a "rbiomvr_vip" class object
          - Plot module removed from rbioFS_plsda_vip() and now a separated function: rbioFS_plsda_vip_plot()
          - rbioClass_plsda_roc_auc() now accepts custom newdata
          - Relevant functions now also output results tst file to the directory
          - Additional argument checks
          
        - Updates to RF-FS functions:
          - Boxplot for the rf_ifs object now has a horizontal line indicating the selection result
          - rf_ifs object now contains: feature_initial_FS, vi_at_threshold, vi_summary, initial_FS_OOB_err_summary, training_initial_FS
          
        - Other update
          - Dependency ggplot2 now requires version 3.0.0
          - Added corresponding print method for the relevant S3 classes
          
        - Bug fixes
          

    0.6.0 (July.20.2018)
        - New generic functions:
          - Generic plot function for permutation test: rbioUtil_perm_plot(). Current supported classes: rbiomvr_perm, rbiosvm_perm
          - Generic plot function for classification: rbioUtil_classplot(). Current supported class: prediction 
        
        - SVM functions added (non-Shiny):
          - rbioClass_svm()
          - rbioClass_svm_roc_auc()
          - rbioClass_svm_perm()
          - rbioClass_svm_predict()
        
        - New PLS-DA functions added (non-Shiny):
          - rbioClass_plsda_perm(): permutation test for plsda models, with two permutation methods
        
        - Updates to SVM functions:
          - Class weight determination functionality added to rbioClass_svm()
          - Additional items added for the modelling settings to the SVM model object (i.e. rbiosvm object)
        
        - Updates to PLS-DA functions (non-Shiny):
          - All PLS-DA function names updated with new prefix: rbioClass_, except for the VIP function, which is a FS function
          - Data centering and scaling argument "center.newdata"" added to rbioClass_plsda_predict()
            - When "center.newdata = TRUE", the function applies training data's col.mean and col.sd to the test data 
          - rbioClass_plsda_predict() new supports data.frame object as newdata
          - Additional argument checking added for rbioClass_plsda_perm()
          - rbioClass_plsda_roc_auc() now correctly uses the centered data from the rbiomvr object for ROC-AUC analysis
          - Argument checking functionality adjusted with correct class checking for all PLS-DA functions
          - Output object of rbioClass_class_perm() is now defined as "rbiomvr_perm" object
          - rbioClass_plsda_perm() updated with plotting capability, using rbioUtil_perm_plot method for class "rbiomvr_perm"
          - rbioClass_plsda_classplot() now changed to a generic function rbioUtil_classplot() applicable to other classifier predictions
        
        - Updates to RF-FS functions:
          - All RF-FS function names updated with new prefix: rbioFS_rf_
          
        - verbose argument added for all the relavent functions so that user can silence the messages
        
        - Overall code base optimization
          
        - Bug fixes
          

    0.5.3
        - Updates to rbioFS_PCA():
          - Legend style adjusted for the sample labels
        
        - Updates to RF-FS functions:
          - For consistency, the dashed line indicators now in red in all the relevant functions

        - Updates to PLS-DA functions (non-Shiny):
          - Bayesian probability option added to rbioFS_plsda_predict()
          - The entire probability calculation and classification module of rbioFS_plsda_classification() merged to rbioFS_plsda_predict()
          - The prediction object will now contain the following sections: predicted.value, probability.summary, probability.method
          - rbioFS_plsda_classification() is now changed to a plotting plot function: rbioFS_plsda_classplot()
          - Options to move probability labels out of the pies adedd for rbioFS_plsda_classplot()
          - rbioFS_plsda_predict(): the legend adjusted to "within threshold" and "outside of threshold"
          - Legend style adjusted for the sample labels for the relevant functions
        
        - Bug fixes
          

    0.5.2
        - New PLS-DA functions added (non-Shiny):
          - rbioFS_plsda_predict(): use the plsda model to calcualte predicted values for unknown data.
          - rbioFS_plsda_classification(): use the predicted values (produced by rbioFS_plsda_predict) to classify. Note: current probability method is "softmax". A "Bayesian" method will be added later. 
          
        - Updates to RF-FS functions:
          - Plotting module separated from the functions
          - Boxplot for the rf_ifs object now horizontal
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
          - Size option added for ggrepel label for the relevant functions
        
        - Bug fixes


    0.5.1
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
