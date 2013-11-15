try:
	from setuptools import setup
except ImportError:
	from distutils.core import setup

config = {
		'description': 'My Project',
		'author': 'Colin A. McFaul',
		'url': 'URL to get it.',
		'download_url': 'Where to download it.',
		'author_email': 'cmcfaul@tulane.edu',
		'version': '0.1',
		'install_requires': ['nose'],
		'packages': ['NAME'],
		'scripts': [],
		'name': 'projectname'
}

setup(**config)
