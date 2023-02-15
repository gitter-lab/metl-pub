run_sa: sa_optimize_files.tar.gz
	condor_submit chtc_sa_optimize.sub

sa_optimize_files.tar.gz: design_tools.py simple_inference run_sa.py models_8852736 configs *.py
	tar -zcvf sa_optimize_files.tar.gz design_tools.py simple_inference run_sa.py models_8852736 configs *.py *.yaml

run_hill_climbing: 
	condor_submit chtc_sa_optimize.sub

hill_climbing_optimize_files.tar.gz: design_tools.py simple_inference run_hill_climbing.py saved_models configs
	tar -zcvf hill_climbing_optimize_files.tar.gz design_tools.py simple_inference run_hill_climbing.py saved_models configs

