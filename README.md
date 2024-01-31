# lifecourse-MR-brain
Scripts written for manuscript employing a two-sample multivariable lifecourse Mendelian Randomization approach in UK Biobank (UKB) to investigate the independent effects of childhood and adult body size on midlife brain structure. Secondary analyses assessing similar brain traits in childhood were conducted in the Adolescent Brain Cognitive Development Study (ABCD).

This project has been my first foray into MR as well as the first time I have used R and Bash to conduct any sort of analyses, so I'm well aware that the way I've approached much of the coding is extremely long-winded and inefficient; this was mainly because I was learning each necessary technique to conduct the analysis as I went before combining everything together to create my final analysis. It does the job though and hopefully will be refined in future projects using these methods.

Paper shortly to be submitted for peer review and page will be updated once it is eventually accpeted and published somewhere. Current abstract posted below though for info.

ABSTRACT

Importance: Obese adults are often reported to have smaller brain volumes than their non-obese peers. Whether this represents evidence of midlife obesity-driven atrophy or is instead a legacy of developmental differences established earlier in the lifespan remains unclear.

Objective: Investigate whether early-life differences in body size explain brain volumetric differences commonly attributed to obesity-related atrophy.

Design, Setting, and Participants: Two-sample lifecourse Mendelian Randomization study in 37,501 adults recruited to UK Biobank (UKB) imaging centers from 2014, with secondary analyses in 6,996 children assessed in the Adolescent Brain Cognitive Development Study (ABCD) recruited from 2018.

Exposures: Genetic instruments for childhood (266 variants) and adult (470 variants) body size derived from a GWAS of 407,741 UKB participants.

Main Outcomes and Measures: Primary outcomes are adult total brain volume; grey matter volume, thickness, and surface area; white matter volume and hyperintensities; and hippocampus, amygdala, and thalamus volumes at mean age 55 in UKB. Secondary outcomes are equivalent childhood measures collected at mean age 11 in ABCD.

Results: In UKB, individuals who were genetically-predicted to have been plumper in childhood were found to have multiple smaller adult brain volumes relative to intracranial volume (e.g. z-score difference in normalised total volume per category increase in body size [95%CI] = - 0.20 [-0.28,-0.12]; p=4x10-6). These effect sizes remained essentially unchanged after accounting for birthweight or current adult size in multivariable models, whereas most adult effects attenuated towards null (e.g. adult z-score [95%CI] for total volume = 0.06 [-0.05,0.17]; p=0.272). Observational analyses in ABCD showed a similar pattern of changes already present in those with a high BMI by age 11 (z-score [95%CI] = - 0.10 [-0.13,-0.07]; p=8x10-13), with follow-up MR analyses providing some evidence for a causal effect already at this early age. Sensitivity analyses revealed that many of these effects were likely due to the persistence of larger head sizes established in those who gained excess weight in childhood, rather than smaller brain sizes per se.

Conclusions and Relevance: Persistence of early-life differences in brain volumes across the lifecourse may underlie numerous neuroimaging traits commonly attributed to obesity-related atrophy in later life.

