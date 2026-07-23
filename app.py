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

st.title("Sequencing Alignment Pipeline")
st.write(
    "Upload a folder of nucleotide `.seq` files, combine them into FASTA, "
    "optionally align them to a nucleotide reference with minimap2, then run "
    "BLASTX against a protein reference and report the best hit and amino-acid substitutions."
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


def parse_blastx(tsv_text: str) -> list[dict[str, object]]:
    best_hits: dict[str, dict[str, object]] = {}

    for line_number, line in enumerate(tsv_text.splitlines(), start=1):
        if not line.strip():
            continue

        row = line.split("\t")
        if len(row) != 11:
            raise ValueError(
                f"Unexpected BLASTX output on line {line_number}: "
                f"expected 11 fields, found {len(row)}."
            )

        (
            query_id,
            subject_id,
            percent_identity,
            alignment_length,
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
        if previous is None or (
            hit["bitscore"],
            -hit["evalue"],
            hit["percent_identity"],
        ) > (
            previous["bitscore"],
            -previous["evalue"],
            previous["percent_identity"],
        ):
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
    "Upload the folder containing `.seq` files",
    type=["seq"],
    accept_multiple_files="directory",
    key="seq_directory",
)

protein_reference = st.file_uploader(
    "Upload the protein reference FASTA for BLASTX",
    type=["fasta", "fa", "faa"],
    key="protein_reference",
)

with st.expander("Optional minimap2 nucleotide-reference alignment"):
    run_minimap = st.checkbox("Run minimap2 before BLASTX", value=False)
    nucleotide_reference = st.file_uploader(
        "Upload nucleotide reference FASTA",
        type=["fasta", "fa", "fna"],
        disabled=not run_minimap,
        key="nucleotide_reference",
    )

with st.sidebar:
    st.header("BLASTX settings")
    evalue_limit = st.number_input(
        "E-value cutoff",
        min_value=0.0,
        value=1e-5,
        format="%.1e",
    )
    max_targets = st.number_input(
        "Maximum target sequences per query",
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

ready = bool(seq_files) and protein_reference is not None
if run_minimap and nucleotide_reference is None:
    ready = False

if not ready:
    st.info("Upload the required files to run the pipeline.")
    st.stop()

if st.button("Run sequencing pipeline", type="primary", use_container_width=True):
    try:
        with st.spinner("Combining sequences and running alignments..."):
            combined_fasta, sequence_summary = combine_seq_files(seq_files)

            with tempfile.TemporaryDirectory() as temp_dir:
                temp = Path(temp_dir)
                query_path = temp / "combined.fasta"
                protein_path = temp / "protein_reference.fasta"
                db_prefix = temp / "protein_db"

                query_path.write_bytes(combined_fasta)
                protein_path.write_bytes(protein_reference.getvalue())

                paf_bytes = None
                match_bytes = None

                if run_minimap:
                    nucleotide_path = temp / "nucleotide_reference.fasta"
                    nucleotide_path.write_bytes(nucleotide_reference.getvalue())

                    minimap_result = run_command(
                        [
                            "minimap2",
                            "-x",
                            "map-ont",
                            "-t",
                            str(int(threads)),
                            str(nucleotide_path),
                            str(query_path),
                        ]
                    )
                    paf_bytes = minimap_result.stdout.encode("utf-8")

                    matches = []
                    for line in minimap_result.stdout.splitlines():
                        fields = line.split("\t")
                        if len(fields) >= 6:
                            matches.append(
                                {"query": fields[0], "reference": fields[5]}
                            )
                    match_bytes = to_tsv(matches, ["query", "reference"])

                run_command(
                    [
                        "makeblastdb",
                        "-in",
                        str(protein_path),
                        "-dbtype",
                        "prot",
                        "-out",
                        str(db_prefix),
                    ]
                )

                outfmt = (
                    "6 qseqid sseqid pident length sstart send "
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
                        "-num_threads",
                        str(int(threads)),
                        "-outfmt",
                        outfmt,
                    ]
                )

                blastx_bytes = blast_result.stdout.encode("utf-8")
                best_hits = parse_blastx(blast_result.stdout)
                substitution_rows = make_substitution_rows(best_hits)
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
                substitutions_bytes = to_tsv(
                    substitution_rows,
                    substitution_columns,
                )

        st.session_state["pipeline_results"] = {
            "combined_fasta": combined_fasta,
            "sequence_summary": sequence_summary,
            "paf": paf_bytes,
            "minimap_matches": match_bytes,
            "blastx": blastx_bytes,
            "substitutions": substitutions_bytes,
            "substitution_rows": substitution_rows,
        }

    except Exception as exc:
        st.error(str(exc))


if "pipeline_results" in st.session_state:
    results = st.session_state["pipeline_results"]

    st.success(
        f"Processed {len(results['sequence_summary'])} sequences; "
        f"{len(results['substitution_rows'])} had a retained BLASTX hit."
    )

    downloads = [
        ("Combined FASTA", results["combined_fasta"], "combined.fasta", "text/plain"),
        ("Raw BLASTX", results["blastx"], "blastx.tsv", "text/tab-separated-values"),
        (
            "Best hits and substitutions",
            results["substitutions"],
            "amino_acid_substitutions.tsv",
            "text/tab-separated-values",
        ),
    ]
    if results["paf"] is not None:
        downloads.extend(
            [
                ("Minimap2 PAF", results["paf"], "aln.paf", "text/plain"),
                (
                    "Minimap2 matches",
                    results["minimap_matches"],
                    "matches.tsv",
                    "text/tab-separated-values",
                ),
            ]
        )

    columns = st.columns(min(3, len(downloads)))
    for index, (label, data, filename, mime) in enumerate(downloads):
        with columns[index % len(columns)]:
            st.download_button(
                f"Download {label}",
                data=data,
                file_name=filename,
                mime=mime,
                key=f"download_{index}",
                use_container_width=True,
            )

    st.subheader("Best BLASTX hits and substitutions")
    if results["substitution_rows"]:
        st.dataframe(results["substitution_rows"], use_container_width=True)
    else:
        st.warning("No BLASTX hits passed the selected settings.")

    with st.expander("Uploaded sequence summary"):
        st.dataframe(results["sequence_summary"], use_container_width=True)
