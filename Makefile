
init:
	bash init.sh

test:

clean:
	rm -rf build
	rm -rf dist
	rm -rf dnastorage.egg-info
	rm -rf generate.egg-info


install:
	pip install --user .

develop:
	pip install -e .
