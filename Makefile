DIR=$(shell head -n 1 params.txt)

write:
	python Write_files.py

copy:
	scp $(DIR)*.i $(DIR)*.h5m $(DIR)*.sh ljjacobson@aci-service-1.chtc.wisc.edu:~/

retrieve:
	scp ljjacobson@aci-service-1.chtc.wisc.edu:"~/*.io ~/*.out ~/*.err" $(DIR)

parse:
	python Parse.py

clean:
	rm -f $(DIR)zCube_* $(DIR)cdp.sh $(DIR)submit_jobs.sh $(DIR)local_runs.sh $(DIR)timing.txt
	rm -f $(DIR)out.sat $(DIR)out* $(DIR)runtp* $(DIR)fcad* $(DIR)lcad* $(DIR)comou* $(DIR)fort*
