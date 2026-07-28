from __future__ import annotations
import base64
import csv
import html
import io
import re
import subprocess
import tempfile
from datetime import datetime
from pathlib import Path

import streamlit as st


st.set_page_config(
    page_title="Sequencing Alignment Pipeline",
    page_icon="🧬",
    layout="wide",
)

# --------------------------------------------------
# Access settings
# --------------------------------------------------

ALLOWED_DOMAIN = "evoralis.com"

# Add manually approved Google-account email addresses here.
ALLOWED_EMAILS = {
    "asha.webb@evoralis.com",
    "valentine.patterson@evoralis.com",
}


def get_user_value(name: str, default=None):
    """Read a claim from st.user safely."""
    try:
        value = getattr(st.user, name)
        if value is not None:
            return value
    except Exception:
        pass

    try:
        return st.user.get(name, default)
    except Exception:
        return default


def claim_is_true(value) -> bool:
    """Convert a Google/OIDC boolean claim safely."""
    if isinstance(value, str):
        return value.strip().lower() == "true"
    return bool(value)

def get_logo_html() -> str:
    """Load the Evoralis logo from the same directory as this app."""
    app_directory = Path(__file__).resolve().parent

    possible_logos = [
        app_directory / "EvoralisLogo.png",
        app_directory / "cropped-cropped-0_Evoralis_logo_for-emails_final_v2.png",
    ]

    for logo_path in possible_logos:
        if logo_path.exists():
            encoded_logo = base64.b64encode(
                logo_path.read_bytes()
            ).decode("ascii")

            return (
                f'<img class="evoralis-logo" '
                f'src="data:image/png;base64,{encoded_logo}" '
                f'alt="Evoralis">'
            )

    return ""


def apply_branding(
    title: str,
    subtitle: str,
    *,
    narrow: bool = False,
) -> None:
    """Apply the Plate QC colour scheme, layout and logo."""
    maximum_width = "900px" if narrow else "1200px"
    logo_html = get_logo_html()

    st.markdown(
        f"""
        <style>
          .stApp {{
            background: #e8f7f5;
          }}

          .block-container {{
            max-width: {maximum_width};
            padding-top: 2rem;
            padding-bottom: 3rem;
          }}

          .hero {{
            display: flex;
            align-items: center;
            gap: 1.4rem;
            background: white;
            border: 1px solid #b9dfd8;
            border-radius: 18px;
            padding: 1.4rem 1.6rem;
            margin-bottom: 1.2rem;
          }}

          .hero-text {{
            flex: 1;
          }}

          .hero h1 {{
            margin: 0;
          }}

          .hero p {{
            margin: .4rem 0 0 0;
          }}

          .evoralis-logo {{
            width: auto;
            height: 70px;
            max-width: 220px;
            object-fit: contain;
          }}

          [data-testid="stSidebar"] {{
            background: #f4fbfa;
            border-right: 1px solid #b9dfd8;
          }}

          div[data-testid="stFileUploader"] {{
            background: white;
            border: 1px solid #b9dfd8;
            border-radius: 14px;
            padding: 0.8rem;
          }}

          div[data-testid="stDataFrame"] {{
            background: white;
            border: 1px solid #b9dfd8;
            border-radius: 14px;
            overflow: hidden;
          }}

          div.stButton > button,
          div.stDownloadButton > button {{
            border-radius: 10px;
          }}
        </style>

        <div class="hero">
          {logo_html}
          <div class="hero-text">
            <h1>{title}</h1>
            <p>{subtitle}</p>
          </div>
        </div>
        """,
        unsafe_allow_html=True,
    )

# --------------------------------------------------
# Authentication
# --------------------------------------------------

if not st.user.is_logged_in:
    apply_branding(
        "Sequencing Alignment Pipeline",
        "This private tool is available to authorised users.",
        narrow=True,
    )

    st.info("Sign in with Google to continue.")

    if st.button(
        "Sign in with Google",
        key="google_login_button",
        type="primary",
        use_container_width=True,
    ):
        st.login()

    st.stop()


email = str(get_user_value("email", "") or "").strip().lower()
email_verified = claim_is_true(get_user_value("email_verified", False))

allowed_emails = {
    address.strip().lower()
    for address in ALLOWED_EMAILS
    if address.strip()
}

is_authorised = (
    email_verified
    and (
        email.endswith(f"@{ALLOWED_DOMAIN}")
        or email in allowed_emails
    )
)

if not is_authorised:
    st.error(
        "Access denied. Your Google account is not authorised to use this application."
    )

    if email:
        st.write(f"Signed-in email: **{email}**")

    if st.button(
        "Sign out",
        key="unauthorised_logout_button",
        use_container_width=True,
    ):
        st.logout()

    st.stop()


with st.sidebar:
    st.success(f"Signed in as {email}")

    if st.button(
        "Sign out",
        key="authorised_logout_button",
        use_container_width=True,
    ):
        st.logout()


apply_branding(
    "Sequencing Alignment Pipeline",
    "Upload nucleotide sequences, analyse matches and identify mutations.",
)


def safe_record_name(filename: str) -> str:
    name = re.sub(r"\s+", "_", Path(filename).stem.strip())
    name = re.sub(r"[^A-Za-z0-9_.:-]", "_", name)
    return name or "unnamed_sequence"


def read_sequence(raw: bytes, filename: str) -> str:
    try:
        text = raw.decode("utf-8-sig")
    except UnicodeDecodeError as exc:
        raise ValueError(f"{filename}: expected a UTF-8 text sequence file.") from exc

    parts = []
    for line in text.splitlines():
        line = line.strip()
        if line and not line.startswith(">"):
            parts.append(re.sub(r"\s+", "", line))

    sequence = "".join(parts).upper()
    if not sequence:
        raise ValueError(f"{filename}: no sequence was found.")

    allowed = set("ACGTUNRYKMSWBDHVX.-")
    invalid = sorted(set(sequence) - allowed)
    if invalid:
        raise ValueError(
            f"{filename}: unexpected sequence characters: {''.join(invalid[:20])}"
        )
    return sequence


def combine_seq_files(files) -> tuple[bytes, list[dict[str, object]]]:
    output = io.StringIO()
    summary = []
    used = set()

    for uploaded in sorted(files, key=lambda item: item.name.lower()):
        base = safe_record_name(uploaded.name)
        record = base
        number = 2
        while record in used:
            record = f"{base}_{number}"
            number += 1
        used.add(record)

        sequence = read_sequence(uploaded.getvalue(), uploaded.name)
        output.write(f">{record}\n")
        for start in range(0, len(sequence), 80):
            output.write(sequence[start:start + 80] + "\n")

        summary.append(
            {"file": uploaded.name, "sequence": record, "length_bp": len(sequence)}
        )

    return output.getvalue().encode("utf-8"), summary


def run_command(command: list[str], timeout: int = 900) -> subprocess.CompletedProcess:
    try:
        result = subprocess.run(
            command,
            text=True,
            capture_output=True,
            check=False,
            timeout=timeout,
        )
    except FileNotFoundError as exc:
        raise RuntimeError(
            f"Required program was not found: {command[0]}. "
            "Check packages.txt and reboot the app."
        ) from exc
    except subprocess.TimeoutExpired as exc:
        raise RuntimeError(f"{command[0]} exceeded the processing time limit.") from exc

    if result.returncode != 0:
        raise RuntimeError(
            f"{command[0]} failed:\n{result.stderr.strip() or 'No error details returned.'}"
        )
    return result


def parse_minimap_matches(paf_text: str) -> list[dict[str, object]]:
    rows = []
    for line in paf_text.splitlines():
        if not line.strip():
            continue
        fields = line.split("\t")
        if len(fields) < 12:
            continue

        query_length = int(fields[1])
        query_start = int(fields[2])
        query_end = int(fields[3])
        matching_bases = int(fields[9])
        block_length = int(fields[10])

        rows.append(
            {
                "query": fields[0],
                "reference": fields[5],
                "strand": fields[4],
                "query_coverage_percent": (
                    100 * (query_end - query_start) / query_length
                    if query_length
                    else 0
                ),
                "identity_percent": (
                    100 * matching_bases / block_length
                    if block_length
                    else 0
                ),
                "mapping_quality": int(fields[11]),
            }
        )
    return rows


def parse_blastx(tsv_text: str) -> list[dict[str, object]]:
    """Keep the single best protein-reference alignment for each query sequence."""
    best_hits: dict[str, dict[str, object]] = {}

    for line_number, line in enumerate(tsv_text.splitlines(), start=1):
        if not line.strip():
            continue

        row = line.split("\t")
        if len(row) != 13:
            raise ValueError(
                f"Unexpected BLASTX output on line {line_number}: "
                f"expected 13 fields, found {len(row)}."
            )

        (
            query_id,
            subject_id,
            percent_identity,
            alignment_length,
            query_start,
            query_end,
            subject_start,
            subject_end,
            query_frame,
            query_alignment,
            subject_alignment,
            evalue,
            bitscore,
        ) = row

        hit = {
            "query_id": query_id,
            "subject_id": subject_id,
            "query_start": int(query_start),
            "query_end": int(query_end),
            "subject_start": int(subject_start),
            "subject_end": int(subject_end),
            "query_frame": int(query_frame),
            "query_alignment": query_alignment.upper(),
            "subject_alignment": subject_alignment.upper(),
            "percent_identity": float(percent_identity),
            "alignment_length": int(alignment_length),
            "evalue": float(evalue),
            "bitscore": float(bitscore),
        }

        previous = best_hits.get(query_id)

        # BLAST bit score is the primary measure of alignment quality.
        # Ties are resolved by E-value, alignment length, and identity.
        hit_rank = (
            hit["bitscore"],
            -hit["evalue"],
            hit["alignment_length"],
            hit["percent_identity"],
        )

        if previous is None:
            best_hits[query_id] = hit
        else:
            previous_rank = (
                previous["bitscore"],
                -previous["evalue"],
                previous["alignment_length"],
                previous["percent_identity"],
            )
            if hit_rank > previous_rank:
                best_hits[query_id] = hit

    return [best_hits[key] for key in sorted(best_hits)]


def call_substitutions(hit: dict[str, object]) -> str:
    query_aln = str(hit["query_alignment"])
    reference_aln = str(hit["subject_alignment"])

    step = 1 if int(hit["subject_end"]) >= int(hit["subject_start"]) else -1
    position = int(hit["subject_start"]) - step
    mutations = []

    for query_aa, reference_aa in zip(query_aln, reference_aln):
        if reference_aa != "-":
            position += step

        if query_aa == reference_aa:
            continue
        if query_aa in {"X", "*"}:
            continue

        if reference_aa != "-" and query_aa != "-":
            mutations.append(f"{reference_aa}{position}{query_aa}")
        elif reference_aa != "-" and query_aa == "-":
            mutations.append(f"{reference_aa}{position}del")

    return ",".join(mutations) if mutations else "WT"


def make_substitution_rows(best_hits: list[dict[str, object]]) -> list[dict[str, object]]:
    rows = []
    for hit in best_hits:
        rows.append(
            {
                "sequence": hit["query_id"],
                "reference": hit["subject_id"],
                "frame": hit["query_frame"],
                "identity": f'{hit["percent_identity"]:.2f}',
                "alignment_length": hit["alignment_length"],
                "evalue": f'{hit["evalue"]:.3g}',
                "bitscore": f'{hit["bitscore"]:.2f}',
                "mutations": call_substitutions(hit),
            }
        )
    return rows


def to_tsv(rows: list[dict[str, object]], columns: list[str]) -> bytes:
    buffer = io.StringIO()
    writer = csv.DictWriter(
        buffer,
        fieldnames=columns,
        delimiter="\t",
        lineterminator="\n",
        extrasaction="ignore",
    )
    writer.writeheader()
    writer.writerows(rows)
    return buffer.getvalue().encode("utf-8")


def rows_to_html_table(
    rows: list[dict[str, object]],
    columns: list[str],
    empty_message: str,
) -> str:
    """Render rows as a self-contained HTML table."""
    if not rows:
        return f'<p class="empty">{html.escape(empty_message)}</p>'

    header = "".join(f"<th>{html.escape(column)}</th>" for column in columns)
    body_rows = []

    for row in rows:
        cells = "".join(
            f"<td>{html.escape(str(row.get(column, '')))}</td>"
            for column in columns
        )
        body_rows.append(f"<tr>{cells}</tr>")

    return (
        '<div class="table-wrap"><table>'
        f"<thead><tr>{header}</tr></thead>"
        f"<tbody>{''.join(body_rows)}</tbody>"
        "</table></div>"
    )


def build_html_report(
    *,
    email: str,
    sequence_count: int,
    minimap_rows: list[dict[str, object]],
    substitution_rows: list[dict[str, object]],
    minimap_preset: str,
    evalue_limit: float,
    max_targets: int,
) -> bytes:
    """Build a standalone HTML report matching the Plate QC report style."""
    minimap_columns = [
        "query",
        "reference",
        "strand",
        "query_coverage_percent",
        "identity_percent",
        "mapping_quality",
    ]
    substitution_columns = [
        "sequence",
        "reference",
        "frame",
        "identity",
        "alignment_length",
        "evalue",
        "bitscore",
        "mutations",
    ]

    minimap_table = rows_to_html_table(
        minimap_rows,
        minimap_columns,
        "minimap2 did not report any alignments.",
    )
    substitutions_table = rows_to_html_table(
        substitution_rows,
        substitution_columns,
        "No BLASTX hits passed the selected settings.",
    )

    logo_html = get_logo_html()
    generated_at = datetime.now().strftime("%Y-%m-%d %H:%M:%S")

    report = f"""<!doctype html>
<html lang="en">
<head>
<meta charset="utf-8">
<meta name="viewport" content="width=device-width, initial-scale=1">
<title>Sequencing Alignment Report</title>
<style>
:root {{
  --bg: #e8f7f5;
  --panel: #f9fffe;
  --text: #1c2434;
  --muted: #667085;
  --border: #b9dfd8;
  --accent: #1b8f84;
}}
* {{ box-sizing: border-box; }}
body {{
  margin: 0;
  background: var(--bg);
  color: var(--text);
  font: 15px/1.5 system-ui, -apple-system, BlinkMacSystemFont, "Segoe UI", sans-serif;
}}
main {{ max-width: 1250px; margin: auto; padding: 90px 20px 60px; }}
.topnav {{
  position: fixed;
  top: 0;
  left: 0;
  right: 0;
  z-index: 20;
  display: flex;
  gap: 8px;
  flex-wrap: wrap;
  padding: 12px 20px;
  background: #f0e8f7;
  border-bottom: 1px solid var(--border);
  backdrop-filter: blur(8px);
}}
.topnav a {{
  color: var(--text);
  text-decoration: none;
  padding: 8px 12px;
  border-radius: 999px;
  background: #9370DB;
  font-weight: 650;
}}
.report-header {{
  display: flex;
  align-items: center;
  justify-content: space-between;
  gap: 30px;
  margin-bottom: 20px;
}}
.report-header .evoralis-logo {{
  width: auto;
  height: 140px;
  max-width: 320px;
  object-fit: contain;
}}
.report-title {{
  text-align: right;
}}
h1 {{ margin: 0 0 5px; font-size: 38px; }}
.subtitle {{ color: var(--muted); margin-bottom: 24px; }}
.cards {{
  display: grid;
  grid-template-columns: repeat(auto-fit, minmax(190px, 1fr));
  gap: 14px;
  margin-top: 18px;
}}
.card, section, details {{
  background: var(--panel);
  border: 1px solid var(--border);
  border-radius: 14px;
  box-shadow: 0 7px 22px rgba(31, 42, 68, .06);
}}
.card {{ padding: 18px; }}
.card .label {{ color: var(--muted); font-size: 13px; }}
.card .value {{ font-size: 27px; font-weight: 750; margin-top: 3px; }}
section {{ margin-top: 18px; padding: 22px; }}
details {{ margin-top: 18px; }}
summary {{
  cursor: pointer;
  padding: 18px 22px;
  font-size: 20px;
  font-weight: 750;
}}
details .content {{ padding: 0 22px 22px; }}
.table-wrap {{ overflow-x: auto; }}
table {{
  border-collapse: collapse;
  width: 100%;
}}
th, td {{
  border-bottom: 1px solid var(--border);
  padding: 9px 11px;
  text-align: right;
  white-space: nowrap;
}}
th {{
  background: #b9dfd8 !important;
  color: #1c2434 !important;
}}
th:first-child, td:first-child {{ text-align: left; }}
.note, .empty {{ color: var(--muted); }}
footer {{ margin-top: 20px; color: var(--muted); font-size: 13px; }}
@media (max-width: 700px) {{
  .report-header {{
    align-items: flex-start;
    flex-direction: column;
  }}
  .report-title {{ text-align: left; }}
  .report-header .evoralis-logo {{ height: 90px; }}
}}
@media print {{
  .topnav {{ display: none; }}
  main {{ padding-top: 20px; }}
  details {{ display: block; }}
  details > * {{ display: block; }}
}}
</style>
</head>
<body>
<nav class="topnav">
  <a href="#summary">Summary</a>
  <a href="#settings">Settings</a>
  <a href="#minimap">Sequence Matches</a>
  <a href="#blastx">Mutations Identified</a>
</nav>

<main>
  <div class="report-header">
    {logo_html}
    <div class="report-title">
      <h1>Sequencing Alignment Report</h1>
      <div class="subtitle">
        <strong>Prepared by:</strong> {html.escape(email or "Unknown user")}<br>
        <strong>Report generated:</strong> {generated_at}
      </div>
    </div>
  </div>

  <details id="Sequence Hits" open>
    <summary>Sequence Hits</summary>
    <div class="content">{minimap_table}</div>
  </details>

  <details id="Mutations" open>
    <summary>Mutations Identified</summary>
    <div class="content">{substitutions_table}</div>
  </details>

</main>
</body>
</html>"""

    return report.encode("utf-8")


seq_files = st.file_uploader(
    "1. Upload the folder containing `.seq` files",
    type=["seq"],
    accept_multiple_files="directory",
    key="seq_directory",
)

minimap_reference = st.file_uploader(
    "2. Upload nucleotide reference FASTA for sequence hit identification",
    type=["fasta", "fa", "fna"],
    key="minimap_reference",
    help="This is the DNA reference used to identify candidates which were Sanger Sequenced.",
)

blast_reference = st.file_uploader(
    "3. Upload protein reference for mutational analysis",
    type=["fasta", "fa", "faa"],
    key="blast_reference",
    help="This is a protein reference used to identify mutations.",
)

with st.sidebar:
    st.header("Pipeline settings")
    minimap_preset = st.selectbox(
        "minimap2 preset",
        options=["map-ont", "lr:hq", "map-pb"],
        index=0,
    )
    evalue_limit = st.number_input(
        "BLASTX E-value cutoff",
        min_value=0.0,
        value=1e-5,
        format="%.1e",
    )
    max_targets = st.number_input(
        "BLASTX maximum target sequences",
        min_value=1,
        max_value=100,
        value=1,
        step=1,
    )


if seq_files:
    st.write(f"**Sequence files selected:** {len(seq_files)}")

if not seq_files or minimap_reference is None or blast_reference is None:
    st.info(
        "Upload the `.seq` folder, the nucleotide minimap2 reference, "
        "and the separate protein BLASTX reference."
    )
    st.stop()

if st.button("Run sequence identification analysis, then identify potential mutations", type="primary", use_container_width=True):
    try:
        with st.spinner("Combining files and running minimap2..."):
            combined_fasta, sequence_summary = combine_seq_files(seq_files)

            with tempfile.TemporaryDirectory() as temp_dir:
                temp = Path(temp_dir)
                query_path = temp / "combined.fasta"
                minimap_ref_path = temp / "minimap_reference.fasta"
                blast_ref_path = temp / "blast_reference.fasta"
                db_prefix = temp / "protein_db"

                query_path.write_bytes(combined_fasta)
                minimap_ref_path.write_bytes(minimap_reference.getvalue())
                blast_ref_path.write_bytes(blast_reference.getvalue())

                minimap_result = run_command(
                    [
                        "minimap2",
                        "-x",
                        minimap_preset,
                        "-t",
                        str(int(threads)),
                        str(minimap_ref_path),
                        str(query_path),
                    ]
                )

                paf_text = minimap_result.stdout
                minimap_rows = parse_minimap_matches(paf_text)
                minimap_simple = [
                    {"query": row["query"], "reference": row["reference"]}
                    for row in minimap_rows
                ]

                with st.spinner("minimap2 complete. Running BLASTX..."):
                    run_command(
                        [
                            "makeblastdb",
                            "-in",
                            str(blast_ref_path),
                            "-dbtype",
                            "prot",
                            "-out",
                            str(db_prefix),
                        ]
                    )

                    outfmt = (
                        "6 qseqid sseqid pident length qstart qend sstart send "
                        "qframe qseq sseq evalue bitscore"
                    )
                    blast_result = run_command(
                        [
                            "blastx",
                            "-query",
                            str(query_path),
                            "-db",
                            str(db_prefix),
                            "-evalue",
                            str(evalue_limit),
                            "-max_target_seqs",
                            str(int(max_targets)),
                            "-max_hsps",
                            "1",
                            "-num_threads",
                            str(int(threads)),
                            "-outfmt",
                            outfmt,
                        ]
                    )

                blastx_text = blast_result.stdout
                best_hits = parse_blastx(blastx_text)
                substitution_rows = make_substitution_rows(best_hits)

        st.session_state["pipeline_results"] = {
            "combined_fasta": combined_fasta,
            "sequence_count": len(sequence_summary),
            "paf": paf_text.encode("utf-8"),
            "minimap_matches": to_tsv(
                minimap_simple,
                ["query", "reference"],
            ),
            "minimap_detailed": to_tsv(
                minimap_rows,
                [
                    "query",
                    "reference",
                    "strand",
                    "query_coverage_percent",
                    "identity_percent",
                    "mapping_quality",
                ],
            ),
            "blastx": blastx_text.encode("utf-8"),
            "substitutions": to_tsv(
                substitution_rows,
                [
                    "sequence",
                    "reference",
                    "frame",
                    "identity",
                    "alignment_length",
                    "evalue",
                    "bitscore",
                    "mutations",
                ],
            ),
            "minimap_rows": minimap_rows,
            "substitution_rows": substitution_rows,
            "html_report": build_html_report(
                email=email,
                sequence_count=len(sequence_summary),
                minimap_rows=minimap_rows,
                substitution_rows=substitution_rows,
                minimap_preset=minimap_preset,
                evalue_limit=float(evalue_limit),
                max_targets=int(max_targets),
            ),
        }

    except Exception as exc:
        st.error(str(exc))


if "pipeline_results" in st.session_state:
    results = st.session_state["pipeline_results"]

    st.success(
        f"Processed {results['sequence_count']} sequences. "
        f"minimap2 reported {len(results['minimap_rows'])} alignments and "
        f"BLASTX retained {len(results['substitution_rows'])} best hits."
    )

    st.subheader("Report")

    st.download_button(
        "Download HTML report",
        data=results["html_report"],
        file_name=f"sequencing_alignment_report_{datetime.now().strftime('%Y-%m-%d')}.html",
        mime="text/html",
        key="download_html_report",
        type="primary",
        use_container_width=True,
    )

    downloads = [
        ("Combined FASTA", results["combined_fasta"], "combined.fasta", "text/plain"),
        ("Minimap2 PAF", results["paf"], "aln.paf", "text/plain"),
        (
            "Minimap2 matches",
            results["minimap_matches"],
            "matches.tsv",
            "text/tab-separated-values",
        ),
        (
            "Detailed minimap2 matches",
            results["minimap_detailed"],
            "minimap_detailed.tsv",
            "text/tab-separated-values",
        ),
        ("Raw BLASTX", results["blastx"], "blastx.tsv", "text/tab-separated-values"),
        (
            "Best hits and substitutions",
            results["substitutions"],
            "amino_acid_substitutions.tsv",
            "text/tab-separated-values",
        ),
    ]

    st.subheader("Other downloads")
    columns = st.columns(3)
    for index, (label, data, filename, mime) in enumerate(downloads):
        with columns[index % 3]:
            st.download_button(
                f"Download {label}",
                data=data,
                file_name=filename,
                mime=mime,
                key=f"download_{index}",
                use_container_width=True,
            )

    st.subheader("minimap2 matches")
    if results["minimap_rows"]:
        st.dataframe(results["minimap_rows"], use_container_width=True)
    else:
        st.warning("minimap2 did not report any alignments.")

    st.subheader("Best BLASTX hits and amino-acid substitutions")
    if results["substitution_rows"]:
        st.dataframe(results["substitution_rows"], use_container_width=True)
    else:
        st.warning("No BLASTX hits passed the selected settings.")
