write:
	python Write_files.py

copy:
	scp *.i   ljjacobson@submit.chtc.wisc.edu:/squid/ljjacobson/
	scp *.sh  ljjacobson@submit.chtc.wisc.edu:~/surface_study/
	scp *.cmd ljjacobson@submit.chtc.wisc.edu:~/surface_study/

retrieve:
	rsync --ignore-existing ljjacobson@submit.chtc.wisc.edu:"~/surface_study/*.io ~/surface_study/*.log" .
	find . -name \*.io -print | wc -l

parse:
	python Parse.py
	find . -name \*.io -print | wc -l

clean:
	rm -f zCube_* comou* out* fort*