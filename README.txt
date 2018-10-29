cd ${youtpath}/PSAP_cluster

Here is an example for trio psap analysis 
bash psap/family_psap_pipeline.sh   test/test.trio.vcf  test/test.trio  test/test.trio.ped

Here is an example for singleton psap analysis
bash  psap/individual_psap_pipeline.sh  test/test.singleton.vcf  test/test.singleton  test/test.singleton.ped

the scripts under script/ is used for cohort PSAP analysis, it requires tabix vcf.gz  and pedigree file for the cohort.
