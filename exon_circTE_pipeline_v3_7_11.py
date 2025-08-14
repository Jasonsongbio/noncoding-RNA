#!/usr/bin/env python3
# exon_circTE_pipeline_v3_7_7.py: Fix header length for RepeatMasker and safe parsing of .out file

import os
import argparse
import subprocess
import logging
from Bio import SeqIO
from Bio.Seq import Seq
from collections import defaultdict

def setup_logging(logfile):
    logging.basicConfig(
        filename=logfile,
        filemode='w',
        level=logging.INFO,
        format='%(asctime)s - %(levelname)s - %(message)s'
    )
    console = logging.StreamHandler()
    console.setLevel(logging.INFO)
    formatter = logging.Formatter('%(asctime)s - %(levelname)s - %(message)s')
    console.setFormatter(formatter)
    logging.getLogger('').addHandler(console)

def parse_blast_and_define_introns(blast_file):
    exons_by_gene = defaultdict(list)
    with open(blast_file) as f:
        for line in f:
            parts = line.strip().split("\t")
            if len(parts) < 11:
                continue
            qid, sid = parts[0], parts[1]
            sstart, send = int(parts[8]), int(parts[9])
            exons_by_gene[qid].append((sid, sstart, send))
    introns_by_gene = {}
    for gid, exons in exons_by_gene.items():
        exons = sorted(exons, key=lambda x: min(x[1], x[2]))
        introns = []
        for i in range(len(exons) - 1):
            chrom = exons[i][0]
            start = max(exons[i][1], exons[i][2])
            end = min(exons[i + 1][1], exons[i + 1][2])
            if end > start:
                introns.append((chrom, start + 1, end - 1))
        introns_by_gene[gid] = introns
    return introns_by_gene

def extract_introns(intron_coords, genome_fasta, out_file):
    genome = SeqIO.to_dict(SeqIO.parse(genome_fasta, "fasta"))
    mapping_file = out_file + ".map.tsv"
    with open(out_file, "w") as out, open(mapping_file, "w") as map_out:
        map_out.write("ShortHeader\tFullHeader\n")
        for gene, introns in intron_coords.items():
            for i, (chrom, start, end) in enumerate(introns):
                seq = genome[chrom].seq[start - 1:end]
                full_header = f"{gene}_i{i+1}_{chrom}_{start}_{end}"
                short_header = f"{gene}_i{i+1}"[:50]  # Ensure <50 chars
                out.write(f">{short_header}\n{seq}\n")
                map_out.write(f"{short_header}\t{full_header}\n")

def run_blast(query, db, out, evalue="1e-10"):
    subprocess.run([
        "blastn", "-query", query, "-db", db,
        "-outfmt", "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore",
        "-evalue", evalue, "-out", out
    ])

def make_blast_db(genome, db_name):
    subprocess.run(["makeblastdb", "-in", genome, "-dbtype", "nucl", "-out", db_name])

def run_repeatmasker(fasta, species, outdir):
    subprocess.run(["RepeatMasker", "-xsmall", "-species", species, "-dir", outdir, fasta])

def parse_repeatmasker_out(out_file):
    te_info = {}
    with open(out_file) as f:
        for line in f:
            if line.startswith("SW") or line.strip() == "" or line.startswith("  "):
                continue
            parts = line.strip().split()
            if len(parts) < 10:
                continue
            try:
                start, end = int(parts[5]), int(parts[6])
                intron = parts[4]
                te_name = parts[9]
                te_info.setdefault(intron, []).append((start, end, te_name))
            except ValueError:
                continue
    return te_info

def extract_te_fragments(masked_fa, output_fa, map_file, repeat_out):
    te_data = parse_repeatmasker_out(repeat_out)
    fragments = []
    with open(output_fa, "w") as out, open(map_file, "w") as map_out:
        map_out.write("TE_ID\tSourceIntron\tStart\tEnd\tTE_Name\n")
        for record in SeqIO.parse(masked_fa, "fasta"):
            name = record.id
            seq = str(record.seq)
            for idx, (start, end, te_name) in enumerate(te_data.get(name, [])):
                frag_seq = seq[start-1:end]
                te_id = f"{name}_TE{idx+1}"
                fragments.append((te_id, name, start, end, te_name))
                out.write(f">{te_id}\n{frag_seq}\n")
                map_out.write(f"{te_id}\t{name}\t{start}\t{end}\t{te_name}\n")
    return {t[0]: t[4] for t in fragments}

def parse_blast_for_rc_hits(blast_file, id_threshold):
    hits = []
    seen = set()
    with open(blast_file) as f:
        for line in f:
            parts = line.strip().split("\t")
            if len(parts) < 11:
                continue
            qid, sid, ident, length, evalue = parts[0], parts[1], float(parts[2]), int(parts[3]), parts[10]
            if qid == sid or (qid, sid) in seen or (sid, qid) in seen:
                continue
            if ident >= id_threshold * 100:
                hits.append((qid, sid, ident, length, evalue))
                seen.add((qid, sid))
    return hits


def extract_te_fragments_with_names(masked_fa, out_file, te_fa, te_map_tsv):
    from Bio import SeqIO
    te_coords = []
    te_seqs = []
    te_id_to_name = {}

    # parse .out for TE names and positions
    with open(out_file) as f:
        for line in f:
            if line.startswith(' '):
                parts = line.strip().split()
                if len(parts) < 15:
                    continue
                try:
                    score = float(parts[0])
                    seqid = parts[4]
                    start = int(parts[5].replace('(', ''))
                    end = int(parts[6].replace(')', ''))
                    strand = parts[8]
                    te_name = parts[9]
                except Exception:
                    continue
                te_coords.append((seqid, start, end, strand, te_name))

    # parse .masked file
    intron_dict = SeqIO.to_dict(SeqIO.parse(masked_fa, "fasta"))
    te_idx = 1
    with open(te_fa, "w") as te_out, open(te_map_tsv, "w") as map_out:
        map_out.write("TE_ID\tTE_Name\tSourceIntron\tStart\tEnd\tStrand\n")
        for seqid, start, end, strand, te_name in te_coords:
            if seqid not in intron_dict:
                continue
            subseq = intron_dict[seqid].seq[start - 1:end]
            if not subseq:
                continue
            te_id = f"TE_{te_idx:05d}"
            te_out.write(f">{te_id}\n{subseq}\n")
            map_out.write(f"{te_id}\t{te_name}\t{seqid}\t{start}\t{end}\t{strand}\n")
            te_id_to_name[te_id] = te_name
            te_idx += 1

    return te_id_to_name

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--mrna_dir", required=True)
    parser.add_argument("--genome_dir", required=True)
    parser.add_argument("--output_dir", required=True)
    parser.add_argument("--repeat_species", required=True)
    parser.add_argument("--evalue", type=str, default="1e-10")
    parser.add_argument("--identity", type=float, default=0.6)
    args = parser.parse_args()

    os.makedirs(args.output_dir, exist_ok=True)
    setup_logging(os.path.join(args.output_dir, "log.txt"))

    mrna_files = [f for f in os.listdir(args.mrna_dir) if f.endswith(".fa") or f.endswith(".fasta")]
    for mrna_file in mrna_files:
        prefix = mrna_file.split("_")[0]
        mrna_path = os.path.join(args.mrna_dir, mrna_file)
        genome_path = os.path.join(args.genome_dir, f"{prefix}_chr.fa")
        blast_db = os.path.join(args.output_dir, f"{prefix}_db")
        blast_out = os.path.join(args.output_dir, f"{prefix}_blast.tsv")
        intron_fa = os.path.join(args.output_dir, f"{prefix}_introns.fa")

        make_blast_db(genome_path, blast_db)
        run_blast(mrna_path, blast_db, blast_out, args.evalue)

        intron_coords = parse_blast_and_define_introns(blast_out)
        extract_introns(intron_coords, genome_path, intron_fa)

        if os.path.getsize(intron_fa) == 0:
            logging.warning(f"No introns extracted for {prefix}. Skipping RepeatMasker and RC analysis.")
            continue

        te_dir = os.path.join(args.output_dir, f"{prefix}_TE")
        os.makedirs(te_dir, exist_ok=True)

        run_repeatmasker(intron_fa, args.repeat_species, te_dir)

        masked_path = os.path.join(te_dir, f"{prefix}_introns.fa.masked")
        repeat_out = os.path.join(te_dir, f"{prefix}_introns.fa.out")
        if not os.path.exists(masked_path):
            logging.error(f"RepeatMasker failed: {masked_path} not found. Skipping {prefix}.")
            continue

        te_frag_fa = os.path.join(te_dir, "te_fragments.fa")
        te_map = os.path.join(te_dir, "te_fragments_map.tsv")
        te_id_to_name = extract_te_fragments_with_names(masked_path, repeat_out, te_frag_fa, te_map)
        qid_to_intron = {}
        with open(te_map) as te_map_file:
            next(te_map_file)
            for row in te_map_file:
                parts = row.strip().split('\t')
                if len(parts) >= 3:
                    qid_to_intron[parts[0]] = parts[2]

        if os.path.getsize(te_frag_fa) == 0:
            logging.warning(f"No TE fragments found for {prefix}. Skipping RC analysis.")
            continue

        subprocess.run(["makeblastdb", "-in", te_frag_fa, "-dbtype", "nucl"])
        blast_rc_out = os.path.join(te_dir, "rc_blast.tsv")
        run_blast(te_frag_fa, te_frag_fa, blast_rc_out, args.evalue)

        rc_hits = parse_blast_for_rc_hits(blast_rc_out, args.identity)
        report = os.path.join(te_dir, "rc_te_report.tsv")
        with open(report, "w") as out:
            out.write("TE1\tTE1_Name\tIntron1\tTE2\tTE2_Name\tIntron2\tIdentity(%)\tMatchLength\tE-value\n")
            for qid, sid, ident, aln_len, evalue in rc_hits:
                qname = te_id_to_name.get(qid, "Unknown")
                sname = te_id_to_name.get(sid, "Unknown")
                out.write(f"{qid}\t{qname}\t{qid_to_intron.get(qid, 'NA')}\t{sid}\t{sname}\t{qid_to_intron.get(sid, 'NA')}\t{ident:.1f}\t{aln_len}\t{evalue}\n")

        logging.info(f"Finished {prefix}. TE reverse complementary analysis saved to {report}.")

if __name__ == "__main__":
    main()
