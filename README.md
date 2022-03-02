# diachrscripts

This is a collection of Python (3.7) scripts for the analysis of simple and twisted read pairs and interactions.

## Documentation

To generate the documentation locally, do the following to set up RTD
```
virtualenv p38  ## first time only
source p38/bin/activate
pip install sphinx
pip install sphinx-rtd-thema
```

To generate documentation
```
source p38/bin/activate
make html
```

You should now find documentation here: ``GIT/diachrscripts/docs/_build/html/index.html``



## Jupyter

To set up a kernel with the required packages for the notebook, do the following.

```
pip install -r requirements.txt
pip install jupyter
ipython kernel install --name "<name>" --user
jupyter notebook
```
(set the kernel to <name>)

## Testing

To setup testing, run the setup.py script
```
source p38/bin/activate
python setup.py test
```
Now test everything at once
```
nosetests tests
```
or specific files
```
 nosetests tests/test_binomial.py 
```
