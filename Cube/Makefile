DIR=$(shell head -n 1 params.txt)

write:
	python Write_files.py

tar:
	tar -czvf tarball.tar.gz $(DIR)zCube* $(DIR)cont.i $(DIR)submit_jobs.sh

copy:
	scp $(DIR)*.i $(DIR)*.ir $(DIR)*.h5m $(DIR)*.sh ljjacobson@aci-service-1.chtc.wisc.edu:~/

retrieve:
	scp ljjacobson@aci-service-1.chtc.wisc.edu:"~/*.io ~/*.ir ~/*.lcad ~/*.fcad ~/*.out ~/*.err" $(DIR)

parse:
	python Parse.py

clean:
	rm -f $(DIR)zCube_* $(DIR)cont.i $(DIR)setup_file*.sh
	rm -f $(DIR)submit_jobs.sh $(DIR)local_runs.sh
	rm -f $(DIR)cubit* $(DIR)history*
	rm -f $(DIR)out* $(DIR)runtp* $(DIR)comou* 
	rm -f $(DIR)fcad* $(DIR)lcad* $(DIR)fort*
