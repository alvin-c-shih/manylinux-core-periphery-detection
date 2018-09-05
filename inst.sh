#sudo rm -rf cpalgorithm.egg-info/ dist/ build/
#python3 setup.py sdist
#python3 -m twine upload --repository-url https://upload.pypi.org/legacy/ dist/*
python3 -m twine upload --repository-url https://test.pypi.org/legacy/ dist/*

#sudo pip3 uninstall cpalgorithm
#sudo pip3 install cpalgorithm

