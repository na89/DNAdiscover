# DNAdiscover
DNAdiscover uses GATK variant annotations to predict the variants with discordant allele frequencies between different platforms of high-throughput DNA sequencing. 

Predictor design
==============

Predictor uses the following list of features, extracted from GATK annotations and that could be found in gnomAD(https://gnomad.broadinstitute.org):
- Low-complexity DNA region (TRUE/FALSE)
- Decoy DNA sequence (TRUE/FALSE)
- Variant falls within a segmental duplication region (TRUE/FALSE)
- Variant (on sex chromosome) falls outside a pseudoautosomal region
- Variant type (snv, indel, multi-snv, multi-indel, or mixed)
- Allele type (insertion, deletion, snv)
- Variant type was mixed (TRUE/FALSE)
- Has * allele (TRUE/FALSE)
- MappingQualityRankSumTest (MQRankSum)
- StrandOddsRatio (SOR)
- Likelihood-based test for the consanguinity among samples (InbreedingCoeff)
- RMSMappingQuality (MQ)
- QualByDepth (QD)
- ReadPosRankSumTest (ReadPosRankSum)
- FisherStrand (FS)
- ClippingRankSumTest (ClippingRankSum)
- ReadPosRankSum
- BaseQRankSum
- Variant quality (QUAL)
- Probability to fail gnomAD RandomForest filter (see https://broadinstitute.github.io/gnomad_methods/api_reference/variant_qc/random_forest.html)

Predictor could use any combination or any number of these features that are available in your dataset.

Installation
==============

To install the package use *devtools* package:

    devtools::install_github("na89/DNAdiscover")
    
Testing the predictor
==============

A test dataset with annotations for 100 variants from gnomAD genomes could be used for testing the predictor

    DNAdiscover::test_data[1:3, ]
    
    lcr decoy segdup nonpar variant_type allele_type was_mixed has_star     qd info_MQRankSum info_SOR
    TRUE FALSE   TRUE  FALSE        mixed         ins      TRUE    FALSE 4.6593          0.751    5.511
    FALSE FALSE  FALSE  FALSE    multi-snv         snv     FALSE     TRUE 2.1140         -0.102    8.509
    FALSE FALSE  FALSE  FALSE        mixed         del      TRUE    FALSE 3.0381         -0.033    0.818
    
    info_InbreedingCoeff info_FS info_QD info_MQ info_DP rf_probability was_split     score    qual BaseQRankSum
    0.2358  72.069   31.59   49.05  396170      0.0401890      TRUE 0.0401890 47210.0       -1.380
    0.0529  79.885    1.43   60.00  547367      0.0064855      TRUE 0.0064855  4817.3       -1.849
    -0.0008   1.024    3.99   60.00  635655      0.1921200      TRUE 0.1921200  4470.1        0.705
    
    ClippingRankSum     FS InbreedingCoeff    MQ MQRankSum    QD ReadPosRankSum          P                loc
    0.000 72.069          0.2358 49.05     0.751 31.59          0.529 5.7543e-01 12:21623282:T:TAAA
    0.046 79.885          0.0529 60.00    -0.102  1.43         -0.924 1.2497e-01     17:7752300:G:C
    0.079  1.024         -0.0008 60.00    -0.033  3.99          0.126 3.3226e-11     17:708263:CG:C


DNAdiscover will produce the classification of variants into 'bad' and 'good' classes, with 'bad' implying similarity to variants with discordant allele frequencies between exomes and genomes in gnomAD. You will also obtain the ROC AUC estimation for the model using the test data. 
    
    results <- DNAdiscover(test_data)
    Test data AUC=0.897370177744595
    
    as.character(results[1:3])
    [1] "good" "good" "bad" 
    
You can check how quality of the predictor changes if the number of features is changed:
    
    results <- DNAdiscover(test_data[, c('lcr', 'BaseQRankSum')])
    Test data AUC=0.630644027385401


# Discordance analysis results


For gnomAD and AllofUs data are available [here](https://www.synapse.org/#!Synapse:syn52587531)

Citation
==============
E. Atkinson*, M. Artomov*, AA Loboda, HL Rehm, DG MacArthur, KJ Karczewski, BM Neale, MJ Daly. "Discordant calls across genotype discovery approaches elucidate variants with systematic errors". Genome Research, 2023. doi: 10.1101/gr.277908.123
