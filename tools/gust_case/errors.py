"""The error a case description is refused with.

Every refusal names the YAML key (``mesh.cell_size``), its value and the rule
it breaks, as GUST does for its own inputs. Nothing is clamped or silently
replaced.
"""


class CaseError(Exception):
    """One or more problems in the case description (one per line)."""

    def __init__(self, problems):
        if isinstance(problems, str):
            problems = [problems]
        self.problems = list(problems)
        super().__init__("\n".join(self.problems))
