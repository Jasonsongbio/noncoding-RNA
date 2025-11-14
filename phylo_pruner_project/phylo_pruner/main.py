"""Core functionality for the phylo_pruner package."""
from __future__ import annotations

import argparse
import logging
from pathlib import Path
from typing import Iterable, List, Optional

from ete3 import NCBITaxa, Tree


LOGGER = logging.getLogger("phylo_pruner")


def configure_logging(log_file: Optional[str] = None, level: int = logging.INFO) -> None:
    """Configure logging for console and optionally a file."""
    LOGGER.setLevel(level)
    formatter = logging.Formatter("%(asctime)s - %(levelname)s - %(message)s")

    # Clear existing handlers
    for handler in list(LOGGER.handlers):
        LOGGER.removeHandler(handler)

    console_handler = logging.StreamHandler()
    console_handler.setFormatter(formatter)
    LOGGER.addHandler(console_handler)

    if log_file:
        log_path = Path(log_file)
        log_path.parent.mkdir(parents=True, exist_ok=True)
        file_handler = logging.FileHandler(log_path)
        file_handler.setFormatter(formatter)
        LOGGER.addHandler(file_handler)


def setup_taxonomy_db(output_file: str) -> None:
    """Download and cache the NCBI taxonomy database using ete3."""
    configure_logging()
    output_path = Path(output_file)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    LOGGER.info("开始构建/下载 NCBI 分类学数据库 ...")
    LOGGER.info("目标数据库文件: %s", output_path)

    try:
        ncbi = NCBITaxa(dbfile=str(output_path))
        # Trigger access to ensure DB is created/downloaded
        ncbi.get_taxid_translator([1])
    except Exception as exc:  # pragma: no cover - external dependency errors
        LOGGER.error("在创建分类学数据库时发生错误: %s", exc)
        raise

    LOGGER.info("数据库已成功创建并保存在: %s", output_path)


def _resolve_taxids_from_names(ncbi: NCBITaxa, names: Iterable[str]) -> dict[str, int]:
    """Resolve scientific names to taxids using ete3 translator."""
    name_list = list(names)
    translator = ncbi.get_name_translator(name_list)
    resolved = {name: ids[0] for name, ids in translator.items() if ids}
    return resolved


def prune_tree(args: argparse.Namespace) -> None:
    """Prune a Newick tree based on taxonomy filtering criteria."""
    configure_logging(getattr(args, "log_file", None))
    LOGGER.info("启动树修剪流程 ...")
    LOGGER.info("输入参数: %s", args)

    taxonomy_db = Path(args.taxonomy_db)
    LOGGER.info("加载分类数据库: %s", taxonomy_db)
    if not taxonomy_db.exists():
        LOGGER.error("分类数据库不存在。请先运行 'setup-taxonomy' 命令。")
        raise FileNotFoundError(
            f"Taxonomy database '{taxonomy_db}' not found. Run setup-taxonomy first."
        )

    try:
        ncbi = NCBITaxa(dbfile=str(taxonomy_db))
    except Exception as exc:  # pragma: no cover - external dependency errors
        LOGGER.error("无法加载分类数据库: %s", exc)
        raise

    LOGGER.info("读取 Newick 树: %s", args.newick_in)
    tree = Tree(args.newick_in, format=1)
    LOGGER.info("在 '%s' 中找到 %d 个叶节点。", args.newick_in, len(tree.get_leaves()))

    leaves = tree.get_leaves()

    if args.leaf_type == "taxid":
        for leaf in leaves:
            try:
                resolved_taxid = int(leaf.name)
                leaf.add_feature("resolved_taxid", resolved_taxid)
            except ValueError:
                LOGGER.warning("叶节点 '%s' 不是有效的整数 TaxID。", leaf.name)
    else:
        unique_names = {leaf.name for leaf in leaves}
        LOGGER.info("解析 %d 个独特的叶节点名称为 TaxID。", len(unique_names))
        resolved_map = _resolve_taxids_from_names(ncbi, unique_names)
        for leaf in leaves:
            taxid = resolved_map.get(leaf.name)
            if taxid is None:
                LOGGER.warning("无法为叶节点 '%s' 解析 TaxID。", leaf.name)
                continue
            leaf.add_feature("resolved_taxid", taxid)

    LOGGER.info("调用 ete3.annotate_ncbi_taxa 进行分类学注释。")
    tree.annotate_ncbi_taxa(taxid_attr="resolved_taxid")
    LOGGER.info("已使用 'resolved_taxid' 特征完成树的分类学注释。")

    target_rank = args.keep_rank.lower()
    LOGGER.info("目标分类级别: %s", target_rank)
    for leaf in leaves:
        if not hasattr(leaf, "lineage") or not leaf.lineage:
            LOGGER.warning("叶节点 '%s' 缺少分类谱系信息。", leaf.name)
            continue
        rank_map = ncbi.get_rank(leaf.lineage)
        name_map = ncbi.get_taxid_translator(leaf.lineage)
        found_name: Optional[str] = None
        for taxid in leaf.lineage:
            rank = rank_map.get(taxid, "").lower()
            if rank == target_rank:
                found_name = name_map.get(taxid)
                break
        if found_name:
            leaf.add_feature("target_rank_name", found_name)
            LOGGER.debug("叶节点 '%s' 在级别 '%s' 上的名称为 '%s'。", leaf.name, target_rank, found_name)
        else:
            LOGGER.info("叶节点 '%s' 在级别 '%s' 上没有匹配名称。", leaf.name, target_rank)

    keep_set = {name.lower() for name in args.keep_names}
    leaves_to_keep: List[str] = []
    for leaf in leaves:
        target_name = getattr(leaf, "target_rank_name", None)
        if target_name and target_name.lower() in keep_set:
            leaves_to_keep.append(leaf.name)
    LOGGER.info("根据规则，将保留 %d 个叶节点。", len(leaves_to_keep))

    if not leaves_to_keep:
        LOGGER.warning("没有叶节点符合保留条件，输出将是一棵空树。")

    tree.prune(leaves_to_keep, preserve_branch_length=True)

    output_path = Path(args.newick_out)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    tree.write(outfile=str(output_path), format=1)
    LOGGER.info("修剪后的树已保存到: %s", output_path)

    if getattr(args, "itol_out", None):
        _write_itol_annotations(tree, Path(args.itol_out))

    LOGGER.info("流程完成。")


def _write_itol_annotations(tree: Tree, output_path: Path) -> None:
    """Write a simple iTOL annotation dataset based on target rank names."""
    output_path.parent.mkdir(parents=True, exist_ok=True)
    LOGGER.info("生成 iTOL 注释文件: %s", output_path)
    with output_path.open("w", encoding="utf-8") as handle:
        handle.write("DATASET_COLORSTRIP\n")
        handle.write("SEPARATOR TAB\n")
        handle.write("DATASET_LABEL\tphylo-pruner\n")
        handle.write("COLOR\t#ff0000\n")
        handle.write("DATA\n")
        for leaf in tree.get_leaves():
            label = getattr(leaf, "target_rank_name", "NA")
            handle.write(f"{leaf.name}\t#1f77b4\t{label}\n")
    LOGGER.info("iTOL 注释文件已保存到: %s", output_path)


def cli_entry() -> None:
    """Entry point for console_scripts."""
    parser = build_parser()
    args = parser.parse_args()
    args.func(args)


def build_parser() -> argparse.ArgumentParser:
    """Build the top-level argument parser for the CLI."""
    parser = argparse.ArgumentParser(
        prog="phylo_pruner",
        description="Annotate and prune phylogenetic trees using the NCBI taxonomy database.",
    )
    subparsers = parser.add_subparsers(dest="command", required=True)

    setup_parser = subparsers.add_parser(
        "setup-taxonomy", help="Download and cache the NCBI taxonomy database."
    )
    setup_parser.add_argument(
        "--output-db",
        required=True,
        help="Path to the SQLite database that will store the NCBI taxonomy cache.",
    )

    def _setup_command(args: argparse.Namespace) -> None:
        setup_taxonomy_db(args.output_db)

    setup_parser.set_defaults(func=_setup_command)

    prune_parser = subparsers.add_parser(
        "prune", help="Annotate a Newick tree and prune leaves by taxonomy."
    )
    prune_parser.add_argument("--newick-in", required=True, help="Input Newick tree file.")
    prune_parser.add_argument(
        "--taxonomy-db", required=True, help="SQLite taxonomy database created by setup-taxonomy."
    )
    prune_parser.add_argument(
        "--leaf-type",
        required=True,
        choices=["name", "taxid"],
        help="Interpretation of leaf labels: scientific names or numeric TaxIDs.",
    )
    prune_parser.add_argument(
        "--keep-rank",
        required=True,
        help="Taxonomic rank to use when filtering leaves (e.g., order, family).",
    )
    prune_parser.add_argument(
        "--keep-names",
        required=True,
        nargs="+",
        help="One or more taxon names to retain at the specified rank.",
    )
    prune_parser.add_argument(
        "--newick-out", required=True, help="Path to write the pruned Newick tree."
    )
    prune_parser.add_argument(
        "--itol-out",
        required=False,
        help="Optional path to write an iTOL color strip annotation dataset.",
    )
    prune_parser.add_argument(
        "--log-file",
        required=False,
        help="Optional path for a detailed execution log file.",
    )
    prune_parser.set_defaults(func=prune_tree)

    return parser


__all__ = [
    "setup_taxonomy_db",
    "prune_tree",
    "build_parser",
    "cli_entry",
]
