# Childhood Adiposity Underlies Numerous Adult Brain Traits Commonly Attributed to Midlife Obesity: A Lifecourse Mendelian Randomization Study
_Scott T. Chiesa, Lydia Rader, Victoria Garfield, Isabelle Foote, Sana Suri, George Davey Smith, Alun D. Hughes, and Tom G. Richardson_

Scripts written for manuscript employing a two-sample multivariable lifecourse Mendelian Randomization approach in UK Biobank (UKB) to investigate the independent effects of childhood and adult adiposity on midlife brain structure. Secondary analyses assessing similar brain traits in childhood were conducted in the Adolescent Brain Cognitive Development Study (ABCD).

This project has been my first foray into MR as well as the first time I have used R and Bash to conduct any sort of analyses, so I'm well aware that the way I've approached much of the coding is _extremely_ long-winded and inefficient. This was mainly because I was learning each necessary technique as I went before combining everything together to create my final analysis. It does the job though and hopefully will be refined in future projects using these methods.

Manuscript was accepted to Brain in June 2024 and DOI to final article will be provided here once published.

## ABSTRACT

Obese adults are often reported to have smaller brain volumes than their non-obese peers. Whether this represents evidence of accelerations in obesity-driven atrophy or is instead a legacy of developmental differences established earlier in the lifespan remains unclear. 

This study aimed to investigate whether early-life differences in adiposity explain differences in numerous adult brain traits commonly attributed to mid-life obesity.

We utilised a two-sample lifecourse Mendelian Randomization study in 37,501 adults recruited to UK Biobank (UKB) imaging centers from 2014, with secondary analyses in 6,996 children assessed in the Adolescent Brain Cognitive Development Study (ABCD) recruited from 2018. Exposures were genetic variants for childhood (266 variants) and adult (470 variants) adiposity derived from a GWAS of 407,741 UKB participants. Primary outcomes were adult total brain volume; grey matter volume, thickness, and surface area; white matter volume and hyperintensities; and hippocampus, amygdala, and thalamus volumes at mean age 55 in UKB. Secondary outcomes were equivalent childhood measures collected at mean age 10 in ABCD.

In UKB, individuals who were genetically-predicted to have had higher levels of adiposity in childhood were found to have multiple smaller adult brain volumes relative to intracranial volume (e.g. z-score difference in normalised brain volume per category increase in adiposity [95%CI] = -0.20 [-0.28,-0.12]; p=4x10-6). These effect sizes remained essentially unchanged after accounting for birthweight or current adult size in multivariable models, whereas most observed adult effects attenuated towards null (e.g. adult z-score [95%CI] for total volume = 0.06 [-0.05,0.17]; p=0.3). Observational analyses in ABCD showed a similar pattern of changes already present in those with a high BMI by age 10 (z-score [95%CI] = -0.10 [-0.13,-0.07]; p=8x10-13), with follow-up genetic risk score analyses providing some evidence for a causal effect already at this early age. Sensitivity analyses revealed that many of these effects were likely due to the persistence of larger head sizes established in those who gained excess weight in childhood (childhood z-score [95%CI] for intracranial volume = 0.14 [0.05,0.23]; p = 0.002), rather than smaller brain sizes per se.

Our data suggest that persistence of early-life developmental differences across the lifecourse may underlie numerous neuroimaging traits commonly attributed to obesity-related atrophy in later life.
