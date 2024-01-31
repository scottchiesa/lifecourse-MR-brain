** SCRIPTS FOR PLINK ANALYSES **


##Loop over phenotype files and run associations between SNPs & outcomes

#MRI Outcomes
for i in UKB_*.txt;
do ~/plink2 \
--ci 0.95 \
--covar /home/rmgpstc/Scratch/UKB_MRBRAIN/covars.txt `#use PCs file in same folder` \
--gen /home/rmgpstc/Scratch/UKB_MRBRAIN/birthweight_snps.gen ref-last `#gen file created in earlier step`\
--glm cols=+a1freq omit-ref `#tell plink to run a glm`\
--pheno ${i} `#tell plink to loop phenos, or if only one pheno use other .sh file`\
--sample /home/rmgpstc/Scratch/UKB_MRBRAIN/71702_ukb_QC_IDlist_F.sample `#use sample ID file in same folder`\
--out /home/rmgpstc/Scratch/UKB_MRBRAIN/${i}_all `#state output directory`
done

##Loop over and keep only additive results
for i in UKB_*.txt_all.PHENO1.glm.linear;
do
awk '($8=="ADD")' ./${i} > ./${i}_ADD; ##select any rows in column 8 containing ADD and then write only these to a new file called *_ADD
done

##Add headers
for i in *_ADD;  do  cat Header_plink ${i} > ${i}_header; done ##concatenate file with headers to new files missing headers

##Keep relevant columns from one of the association results files (TV) I THINK I CAN REMOVE THIS...
##awk 'NR>0{print $1,$2, $3, $4, $5, $6, $7, $9, $10, $11, $12, $13, $15}' TV_ADD_header > TV_MR_input

##For the other results files, retain columns 9-13+15 for MR input
##for j in *_header;
##do
##awk '{print $3, $9, $10, $11, $12, $13, $15}' ${j} > ${j}_MR_input; ##for each file created above, take only necessary columns and print to new file
##done

mv UKB_icv.txt_all.PHENO1.glm.linear_ADD_header bw_icv.txt
mv UKB_tv.txt_all.PHENO1.glm.linear_ADD_header bw_tv.txt
mv UKB_gm.txt_all.PHENO1.glm.linear_ADD_header bw_gm.txt
mv UKB_area.txt_all.PHENO1.glm.linear_ADD_header bw_area.txt
mv UKB_thick.txt_all.PHENO1.glm.linear_ADD_header bw_thick.txt
mv UKB_wm.txt_all.PHENO1.glm.linear_ADD_header bw_wm.txt
mv UKB_wmh.txt_all.PHENO1.glm.linear_ADD_header bw_wmh.txt
mv UKB_logwmh.txt_all.PHENO1.glm.linear_ADD_header bw_logwmh.txt
mv UKB_thalamus.txt_all.PHENO1.glm.linear_ADD_header bw_thalamus.txt
mv UKB_hippocampus.txt_all.PHENO1.glm.linear_ADD_header bw_hippo.txt
mv UKB_amygdala.txt_all.PHENO1.glm.linear_ADD_header bw_amygdala.txt
mv UKB_tvn.txt_all.PHENO1.glm.linear_ADD_header bw_tvn.txt
mv UKB_gmn.txt_all.PHENO1.glm.linear_ADD_header bw_gmn.txt
mv UKB_arean.txt_all.PHENO1.glm.linear_ADD_header bw_arean.txt
mv UKB_thickn.txt_all.PHENO1.glm.linear_ADD_header bw_thickn.txt
mv UKB_wmn.txt_all.PHENO1.glm.linear_ADD_header bw_wmn.txt
mv UKB_wmhn.txt_all.PHENO1.glm.linear_ADD_header bw_wmhn.txt
mv UKB_logwmhn.txt_all.PHENO1.glm.linear_ADD_header bw_logwmhn.txt
mv UKB_thaln.txt_all.PHENO1.glm.linear_ADD_header bw_thaln.txt
mv UKB_hippon.txt_all.PHENO1.glm.linear_ADD_header bw_hippon.txt
mv UKB_amygn.txt_all.PHENO1.glm.linear_ADD_header bw_amygn.txt
mv UKB_gmb.txt_all.PHENO1.glm.linear_ADD_header bw_gmb.txt
mv UKB_areab.txt_all.PHENO1.glm.linear_ADD_header bw_areab.txt
mv UKB_thickb.txt_all.PHENO1.glm.linear_ADD_header bw_thickb.txt
mv UKB_wmb.txt_all.PHENO1.glm.linear_ADD_header bw_wmb.txt
mv UKB_logwmhb.txt_all.PHENO1.glm.linear_ADD_header bw_logwmhb.txt
mv UKB_thalb.txt_all.PHENO1.glm.linear_ADD_header bw_thalb.txt
mv UKB_hippob.txt_all.PHENO1.glm.linear_ADD_header bw_hippob.txt
mv UKB_amygb.txt_all.PHENO1.glm.linear_ADD_header bw_amygb.txt

rm *txt_all*
