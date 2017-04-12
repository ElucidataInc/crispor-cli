
crispor_env: crispor_env/bin/activate

crispor_env/bin/activate: requirements.txt
	test -d crispor_env || virtualenv crispor_env
	crispor_env/bin/pip install -Ur requirements.txt
	touch crispor_env/bin/activate

devbuild: crispor_env
	crispor_env/bin/python setup.py install
