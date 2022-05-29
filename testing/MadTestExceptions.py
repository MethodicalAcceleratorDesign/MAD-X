#!/usr/bin/python


class MadException(Exception):
    pass

class BuildException(MadException):
    pass

class TestException(MadException):
    pass
