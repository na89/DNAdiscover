# DNAdiscover
Predictor for technical bias in DNA variants in NGS data

Installation
==============

To install the package use *devtools* package:

    devtools::install_github("na89/DNAdiscover")
    
Testing the predictor
==============
A test dataset with annotations for 100 variants from gnomAD genomes could be used for testing the predictor

    DNAdiscover(test_data)
    Test data AUC=0.897370177744595
    
You can check how quality of the predictor changes if the number of features is changed:
    
    DNAdiscover(test_data[, c('lcr', 'BaseQRankSum')])
    Test data AUC=0.630644027385401
