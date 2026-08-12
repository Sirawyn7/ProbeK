class ProbekError(Exception):
    """Base class for errors that should abort the CLI with a clean message."""


class MissingToolError(ProbekError):
    pass


class ReferenceDataMissingError(ProbekError):
    pass


class InputDataError(ProbekError):
    pass
