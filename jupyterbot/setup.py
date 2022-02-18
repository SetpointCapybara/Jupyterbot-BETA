from setuptools import setup

setup(
   name='jupyterbot',
   version='0.1',
   description='A useful module',
   author='Man Foo',
   author_email='foomail@foo.com',
   packages=['foo'],  #same as name
   install_requires=['quadprog', 'colour']
)