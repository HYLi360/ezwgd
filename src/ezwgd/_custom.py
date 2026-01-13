class BaseException4EZWGD(BaseException):
    def __init__(self, *args: object) -> None:
        super().__init__(*args)


class BaseWarning4EZWGD(UserWarning):
    pass


# ---all custom class (inheritances from BaseException and Warning) we used.--------------------------
class ProgramNotFoundError(BaseException4EZWGD):
    """Called when the external executable file does not exist."""

    def __init__(self, program: str) -> None:
        self.message = (
            f"[ERROR] Can't Find the Executable File of This Program: {program}"
        )

    def __str__(self) -> str:
        return f"{self.message}"

class SeqFileNotFoundError(BaseException4EZWGD):
    """Called when the external executable file does not exist."""

    def __init__(self, seq_file_path: str) -> None:
        self.message = (
            f"[ERROR] Can't Find This Sequence File: {seq_file_path}"
        )

    def __str__(self) -> str:
        return f"{self.message}"

class GFF3FileNotFoundError(BaseException4EZWGD):
    """Called when the external executable file does not exist."""

    def __init__(self, gff3_file_path: str) -> None:
        self.message = (
            f"[ERROR] Can't Find This GFF3 File: {gff3_file_path}"
        )

    def __str__(self) -> str:
        return f"{self.message}"
