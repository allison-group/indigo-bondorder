class IndigoException(Exception):
    """Base class for exceptions in Indigo."""
    
class IndigoError(IndigoException):
    """Exception for a serious error in Indigo."""
    
class IndigoWarning(Warning):
    """Base class for Indigo warnings."""
    
class IndigoFileError(IndigoError):
    """Exception for errors caused by handling of files."""
    
class IndigoMissingParameters(IndigoError):
    """Exception for functions not being passed all required parameters when 
    function has multiple modes of operations."""
    
class IndigoPointlessConcept(IndigoError):
    """Exception for if something silly and pointless is attempted."""
    
class IndigoUnfeasibleComputation(IndigoError):
    """Exception for when the input information falls outside the bounds of
    computational ability."""
    
class IndigoSearchError(IndigoError):
    """Exception for when a search cannot find a result."""
    
class IndigoMissingFunctionality(IndigoWarning):
    """Warning for if some non-critical functionality is missing due to 
    user system setup."""

class IndigoExternalProgramError(IndigoError):
    """Error in an external program called by indigo
    """

