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


class NotTheSameRecordsError(BaseException4EZWGD):
    """Called when the ID from CDS and Protein Records are mismatch."""

    # Such as:
    # CDS Records       Sp1g1   Sp2g2
    # Protein Records           Sp2g2  Sp4g4
    def __init__(self, cds_records: tuple, prot_records: tuple) -> None:
        self.message = "[ERROR] Unable to Match the Two Records from CDS and Protein."
        self.cds = f"CDS Records ID: {cds_records[0].id}, {cds_records[1].id}"
        self.prot = f"Protein Records ID: {prot_records[1].id}, {prot_records[0].id}"

    def __str__(self) -> str:
        return f"{self.message}{self.cds}, {self.prot}"


class NotTheSameOrderWarning(BaseWarning4EZWGD):
    """[WARNING] CDS与蛋白质记录的ID顺序不一致。已自动调换蛋白质两个记录的先后顺序。"""
