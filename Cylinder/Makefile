DIR=$(shell head -n 1 params.txt)

write:
	python Write_files.py

tar:
	tar -czvf tarball.tar.gz $(DIR)zCyl* $(DIR)cont.i $(DIR)submit_jobs.sh

copy:
	scp $(DIR)*.i $(DIR)*.ir $(DIR)*.h5m $(DIR)*.sh ljjacobson@aci-service-1.chtc.wisc.edu:~/

retrieve:
	scp ljjacobson@aci-service-1.chtc.wisc.edu:"~/*.io ~/*.ir ~/*.lcad ~/*.lcad ~/*.out ~/*.err" $(DIR)

parse:
	python Parse.py

clean:
	rm -f $(DIR)zCyl_* $(DIR)cont.i $(DIR)setup_files.sh
	rm -f $(DIR)submit_jobs.sh $(DIR)local_runs.sh
	rm -f $(DIR)cubit* $(DIR)history*
	rm -f $(DIR)out* $(DIR)runtp* $(DIR)comou* 
	rm -f $(DIR)fcad* $(DIR)lcad* $(DIR)fort*
