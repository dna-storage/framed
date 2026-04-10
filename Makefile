
init:
	bash init.sh

test:

clean:
	rm -rf build
	rm -rf dist
	rm -rf dnastorage.egg-info
	rm -rf generate.egg-info


install:
	pip3 install --user .

develop:
	pip3 install -e .
