
init:
	tcsh init.csh
test:

clean:
	rm -rf build
	rm -rf dist
	rm -rf dnastorage.egg-info
	rm -rf generate.egg-info


install:
	python setup.py install --user

develop:
	python setup.py develop 
