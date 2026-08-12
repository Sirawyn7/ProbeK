from __future__ import annotations

from dataclasses import dataclass, field
from typing import Any, Literal


@dataclass(frozen=True)
class ProbeOccurrence:
    """One (target, probe name) pairing that a given sequence appeared under."""

    target: str
    name: str
    row: dict[str, Any]


@dataclass
class SequenceGroup:
    """All occurrences of one unique probe sequence, possibly across targets."""

    sequence: str
    occurrences: list[ProbeOccurrence] = field(default_factory=list)


@dataclass(frozen=True)
class BlastHit:
    sseqid: str
    pident: float
    length: int
    qstart: int
    qend: int
    sstart: int
    send: int
    evalue: float
    qlen: int

    @property
    def coverage(self) -> float:
        return self.length / self.qlen

    @property
    def span(self) -> tuple[int, int]:
        """(start, end) with start <= end, normalized for reverse-strand hits."""
        return (min(self.sstart, self.send), max(self.sstart, self.send))

    def to_dict(self) -> dict[str, Any]:
        return {
            "sseqid": self.sseqid,
            "pident": self.pident,
            "length": self.length,
            "qstart": self.qstart,
            "qend": self.qend,
            "sstart": self.sstart,
            "send": self.send,
            "evalue": self.evalue,
            "qlen": self.qlen,
        }

    @classmethod
    def from_dict(cls, d: dict[str, Any]) -> "BlastHit":
        return cls(
            sseqid=d["sseqid"],
            pident=float(d["pident"]),
            length=int(d["length"]),
            qstart=int(d["qstart"]),
            qend=int(d["qend"]),
            sstart=int(d["sstart"]),
            send=int(d["send"]),
            evalue=float(d["evalue"]),
            qlen=int(d["qlen"]),
        )


Tier = Literal["A", "B", "C"]


@dataclass(frozen=True)
class ClassifiedHit:
    hit: BlastHit
    tier: Tier
    locus_label: str | None  # ERVK repName for tier A, gene symbol for tier C, else None


@dataclass
class TargetResult:
    target: str
    ranked_rows: list[dict[str, Any]]
    selected_names: set[str]
    short: bool
