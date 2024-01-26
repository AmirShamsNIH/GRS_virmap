snakemake --snakefile /data/shamsaddinisha/Development/GRS_virmap/snakemake/pre_process.smk \
--configfile /data/shamsaddinisha/Test_Space/GRS_virmap/grs_virmap.json --keep-incomplete --cores



snakemake --snakefile /data/shamsaddinisha/Development/GRS_virmap/snakemake/alignment.smk \
--configfile /data/shamsaddinisha/Test_Space/GRS_virmap/grs_virmap.json --keep-incomplete --cores



python /data/shamsaddinisha/Development/GRS_virmap/grs_virmap.py -i grs_virmap.json