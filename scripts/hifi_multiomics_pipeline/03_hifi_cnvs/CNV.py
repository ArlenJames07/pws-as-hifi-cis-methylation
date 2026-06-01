from pathlib import Path
import runpy


SOURCE = Path(__file__).resolve().parents[2] / "SV_calling" / "CNV.py"


if __name__ == "__main__":
    runpy.run_path(str(SOURCE), run_name="__main__")
