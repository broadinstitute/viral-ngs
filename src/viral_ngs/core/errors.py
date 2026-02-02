#!/usr/bin/env python

class QCError(RuntimeError):
    '''Indicates a failure at a QC step.'''
    
    def __init__(self, reason):
        super(QCError, self).__init__(reason)

class InvalidBamHeaderError(ValueError):
    '''Indicates a malformed or otherwise unusable bam file header'''
    pass