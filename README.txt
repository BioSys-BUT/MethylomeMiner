
How to run methylome_miner

Run script methylome_miner.py in cmd with arguments:
1) --input_bed_file		path to BED file to be analyzed
2) --input_bed_dir		path to dir with all BED files
3) --input_annot_file		path to file with annotation
4) --output_files_prefix	output file path + prefix of new files
5) optional: --min_coverage		value of min coverage for filtering
6) optional: --min_percent_modified	value for min percent modified

Example for PyCharm Run configuration:
--input_bed_file
input_bed_files\KP825_b53_4mC_5mC_6mA_calls_modifications_whole_run_aligned_sorted_pileup_SUP_qscore.bed
--input_bed_dir
input_bed_files
--input_annot_file
input_roary_output_files_corrected\KP825_genome.gff
--output_files_prefix
methylome_miner_output\KP825
--min_coverage
37


python -m pip install -e .

MethylomeMiner --input_bed_file methylations_data\input_bed_files\KP825_b53_4mC_5mC_6mA_calls_modifications_whole_run_aligned_sorted_pileup_SUP_qscore.bed --input_bed_dir methylations_data\input_bed_files --input_annot_file methylations_data\input_roary_output_files_corrected\KP825_genome.gff --work_dir MethylomeMiner_output --file_name KP825 --min_coverage 37

CoreMethylomeMiner --input_bed_dir methylations_data\input_bed_files --input_annot_dir methylations_data\input_roary_output_files_corrected --roary_output_file methylations_data\roary_output\gene_presence_absence.csv --min_coverage 37 --work_dir MethylomeMiner_output --matrix_values presence