# to create a new virtual env
virtualenv --python=python2.7 myenv

# active the virtual env

source myenv/bin/activate

# install requirements.txt

cat requirements.txt | xargs -n1 pip install
