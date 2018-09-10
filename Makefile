init:
	make -C other_software/nupack3.0.6 all
	pip install -r requirements.txt

test:
        #! /usr/bin/python2.7.15
	nosetests tests

clean:
	rm -rf build
	rm -rf dist
	rm -rf dnastorage.egg-info
	rm -rf generate.egg-info

install:
	python setup.py install
