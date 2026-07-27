from __future__ import annotations

import csv
import io
import re
import subprocess
import tempfile
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


# --------------------------------------------------
# Authentication
# --------------------------------------------------

if not st.user.is_logged_in:
    st.title("Sequencing Alignment Pipeline")
    st.write(
        "Upload a folder of nucleotide `.seq` files, a nucleotide FASTA for minimap2, "
        "and a separate protein FASTA for BLASTX. The app runs minimap2 first, then BLASTX, "
        "and reports the best protein alignment and amino-acid substitutions."
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


st.title("Sequencing Alignment Pipeline")
st.write(
    "Upload a folder of nucleotide `.seq` files, a nucleotide FASTA for minimap2, "
    "and a separate protein FASTA for BLASTX. The app runs minimap2 first, then BLASTX, "
    "and reports the best protein alignment and amino-acid substitutions."
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


seq_files = st.file_uploader(
    "1. Upload the folder containing `.seq` files",
    type=["seq"],
    accept_multiple_files="directory",
    key="seq_directory",
)

minimap_reference = st.file_uploader(
    "2. Upload nucleotide reference FASTA for minimap2",
    type=["fasta", "fa", "fna"],
    key="minimap_reference",
    help="This is the DNA reference used by minimap2.",
)

blast_reference = st.file_uploader(
    "3. Upload protein reference FASTA for BLASTX",
    type=["fasta", "fa", "faa"],
    key="blast_reference",
    help="This is a separate protein reference used by makeblastdb and blastx.",
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
        value=10,
        step=1,
    )
    threads = st.number_input(
        "Threads",
        min_value=1,
        max_value=2,
        value=2,
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

if st.button("Run minimap2, then BLASTX", type="primary", use_container_width=True):
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
