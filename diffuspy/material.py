"""Module that contains material properties in a class

"""


class Material(object):
    """Instanciate a material object

    """
    def __init__(self, **kwargs):
        self.__dict__.update(kwargs)
