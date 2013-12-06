write:
	python Write_files.py

copy:
	scp *.i *.sh ljjacobson@aci-service-1.chtc.wisc.edu:~/

retrieve:
	scp ljjacobson@aci-service-1.chtc.wisc.edu:"~/*.io ~/*.out ~/*.err" ./
	find . -maxdepth 1 -name \*.io -print | wc -l

parse:
	python Parse.py
	find . -maxdepth 1 -name \*.io -print | wc -l

clean:
	rm -f zCube_* comou* out* fort* submit_jobs.sh