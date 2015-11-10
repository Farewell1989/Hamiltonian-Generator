#from numpy.distutils.core import setup, Extension
def configuration(parent_package='',top_path=None):
    from numpy.distutils.misc_util import Configuration
    config=Configuration('Hamiltonian',parent_package,top_path)
    config.add_subpackage('Core')
    config.make_config_py()
    return config

if __name__ == '__main__':
    from numpy.distutils.core import setup
    setup(configuration=configuration)
