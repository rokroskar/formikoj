# Manages building, testing, and cleaning the code as well as running the code
# to generate the results and figures for the paper.

# CONFIGURATION
###############################################################################

# Set the package name
PACKAGE = formikoj

# TARGETS
###############################################################################

help:
	@echo "Commands:"
	@echo ""
	@echo "  all        runs 'clean' and 'build'"
	@echo "  build      build and install the package"
	@echo "  clean      clean up the package build"
	@echo ""

all: clean build

clean:
	find . -name "*.pyc" -exec rm -v {} \;
	rm -rvf build dist MANIFEST *.egg-info __pycache__ .coverage .cache

build:
	python setup.py install
	
dist:
	python setup.py bdist_wheel
