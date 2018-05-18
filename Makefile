init:
	make -C other_software/nupack3.0.6 all
	pip install -r requirements.txt

test:
	nosetests tests
