from setuptools import setup, Extension

symnmf_module = Extension('symnmf',
                          sources=['symnmfmodule.c', 'symnmf.c'])


setup(name='symnmf',
      version='1.0',
      description='Symnmf Implementation',
      ext_modules=[symnmf_module])