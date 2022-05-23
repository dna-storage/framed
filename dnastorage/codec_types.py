
'''
These are useful classes to tag codecs with to provide more information on what function the codec serves, neccessary if using a pipeline based building interface
'''

class BaseConversionType(object):
    def __init__(self):
        pass


class DNAtoDNA(BaseConversionType):
    def __init__(self):
        pass

class CWtoDNA(BaseConversionType):
    def __init__(self):
        pass

class CWtoCW(BaseConversionType):
    def __init__(self):
        pass

class CWConsolidate(BaseConversionType):
    def __init__(self):
        pass

class DNAConsolidate(BaseConversionType):
    def __init__(self):
        pass

class Probe(BaseConversionType):
    def __init__(self):
        pass
    
