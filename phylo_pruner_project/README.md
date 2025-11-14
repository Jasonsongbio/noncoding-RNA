# phylo-pruner

`phylo-pruner` 是一个命令行工具，用于结合本地 NCBI 分类学数据库对 Newick 格式的系统发育树进行注释和修剪。该项目遵循 GitHub 与 Bioconda 分发的最佳实践，提供可重复的分析流程和友好的 CLI 体验。

## 功能概览

- **数据库固定（Database Pinning）**：利用 `ete3.NCBITaxa` 下载并缓存 NCBI taxonomy 数据库，确保后续分析的可重复性。
- **树注释与修剪**：支持根据指定的分类级别和名称过滤 Newick 树的叶节点。
- **iTOL 注释输出**：自动生成 iTOL color strip 数据集，以便可视化修剪结果。

## 安装

### 使用 Conda / Bioconda（推荐）

本项目的目标是发布到 Bioconda，届时可通过以下命令安装：

```bash
conda install -c bioconda phylo-pruner
```

### 从源码安装

克隆项目仓库并使用 `setup.py` 安装：

```bash
git clone https://github.com/your-org/phylo-pruner.git
cd phylo-pruner
python setup.py install
```

## 依赖

- [ETE Toolkit (`ete3`)](http://etetoolkit.org/) — 负责 taxonomy 数据库下载、树结构处理及注释。
- Python 标准库组件：`argparse`、`logging`、`pathlib` 等。

## 快速上手示例

### 步骤 1：下载 / 固定分类学数据库

```bash
phylo_pruner setup-taxonomy --output-db taxa.sqlite
```

### 步骤 2：修剪系统发育树

```bash
phylo_pruner prune \
    --newick-in tree.nwk \
    --taxonomy-db taxa.sqlite \
    --leaf-type name \
    --keep-rank order \
    --keep-names Carnivora Primates \
    --newick-out pruned.nwk \
    --itol-out itol_colors.txt \
    --log-file run.log
```

执行完成后，将得到：

- `pruned.nwk`：修剪后的树文件。
- `itol_colors.txt`：iTOL 注释文件，可用于可视化。
- `run.log`：详细的运行日志。

## CLI 参数详解

### `setup-taxonomy`

| 参数 | 类型 | 说明 |
| --- | --- | --- |
| `--output-db` | `str` (必需) | 指定 SQLite 数据库输出路径。若文件不存在，将自动下载并构建。 |

### `prune`

| 参数 | 类型 | 说明 |
| --- | --- | --- |
| `--newick-in` | `str` (必需) | 输入的 Newick 树文件。 |
| `--taxonomy-db` | `str` (必需) | 使用 `setup-taxonomy` 命令生成的 taxonomy 数据库路径。 |
| `--leaf-type` | `str` (必需, 可选值：`name`, `taxid`) | 叶节点名称类型，是科学名还是 NCBI TaxID。 |
| `--keep-rank` | `str` (必需) | 要保留的分类级别（如 `order`, `family` 等）。 |
| `--keep-names` | `str` (必需, 多值) | 在指定分类级别下需要保留的分类单元名称列表。 |
| `--newick-out` | `str` (必需) | 输出修剪后 Newick 树的文件路径。 |
| `--itol-out` | `str` (可选) | 输出 iTOL 注释文件的路径。 |
| `--log-file` | `str` (可选) | 将详细日志写入该文件，并同时在控制台输出。 |

## 示例数据

项目在 `examples/` 目录中提供了一个 `sample_tree.nwk` 树以及使用说明，可帮助用户快速测试流程。

## 许可证

本项目基于 MIT License 发布，详情请见 [LICENSE](LICENSE) 文件。
