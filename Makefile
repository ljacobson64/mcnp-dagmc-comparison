write:
	python Write_files.py

copyHPC:
	scp *.i *.sh ljjacobson@aci-service-1.chtc.wisc.edu:~/

retrieveHPC:
	scp ljjacobson@aci-service-1.chtc.wisc.edu:"~/*.io ~/*.out ~/*.err" ./
	find . -maxdepth 1 -name \*.io -print | wc -l

copyHTC:
	scp *.i        ljjacobson@submit.chtc.wisc.edu:/squid/ljjacobson/
	scp *.sh *.cmd ljjacobson@submit.chtc.wisc.edu:~/surface_study/

retrieveHTC:
	scp ljjacobson@submit.chtc.wisc.edu:"~/surface_study/*.io ~/surface_study/*.log" ./
	find . -maxdepth 1 -name \*.io -print | wc -l

parse:
	python Parse.py
	find . -maxdepth 1 -name \*.io -print | wc -l

clean:
	rm -f zCube_* comou* out* fort* submit_jobs.sh