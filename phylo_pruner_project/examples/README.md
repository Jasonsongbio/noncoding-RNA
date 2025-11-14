# Examples

This directory contains sample data demonstrating how to use **phylo-pruner**.

## Files

- `sample_tree.nwk` — A toy Newick tree with a handful of mammal and bird species.

## Quickstart

```bash
# 1. Create or download the taxonomy database (replace the path as needed)
phylo_pruner setup-taxonomy --output-db taxa.sqlite

# 2. Prune the sample tree by keeping specific orders
phylo_pruner prune \
    --newick-in examples/sample_tree.nwk \
    --taxonomy-db taxa.sqlite \
    --leaf-type name \
    --keep-rank order \
    --keep-names Primates Rodentia \
    --newick-out examples/pruned_sample.nwk \
    --itol-out examples/pruned_sample_itol.txt \
    --log-file examples/run.log
```

After running the commands above you will obtain a pruned tree containing only the
specified orders plus corresponding iTOL annotations and a detailed execution log.
