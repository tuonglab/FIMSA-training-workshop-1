# Clonotype Cheat Sheet

Copy and paste each chunk below into your notebook, one at a time.

---

## Part 1: Simple looping through a row at a time

```python
################# simple looping through a row at a time
# Loop through a dataframe one cell at a time
for i, row in df.iterrows():
    for column in df.columns:
        value = row[column]

        print("row:", i)
        print("column:", column)
        print("value:", value)
        break  # Remove this break if you want to iterate through all cells
```

---

## Part 2: Clonotype dictionary using TRB only, row by row

```python
################ Getting to a clonotype dictionary using TRB only row by row
from collections import defaultdict

clone_dict = defaultdict(list)

for i, row in df.iterrows():

    # Access individual cells
    cell_barcode = row["barcode"]
    chain = row["chain"]
    v_gene = row["v_gene"]
    j_gene = row["j_gene"]
    cdr3 = row["cdr3"]

    # Define clonotype using TRB cells
    if chain == "TRB":
        clonotype = v_gene + "_" + j_gene + "_" + cdr3
        clone_dict[clonotype].append(cell_barcode)

for clonotype, barcodes in list(clone_dict.items())[:3]:
    print(clonotype, ":", barcodes)
```

---

## Part 3: Clonotype dictionary using TRB, by subsetting the dataframe first

```python
############### Getting to a clonotype dictionary using TRB by subsetting the dataframe first
from collections import defaultdict

clone_dict = defaultdict(list)

for i in df.index:

    # Identify the cell
    cell_barcode = df.loc[i, "barcode"]

    # Get all contigs belonging to this cell
    cell_contigs = df[df["barcode"] == cell_barcode]

    print("CELL:", cell_barcode)

    # Work through each contig in this cell
    for j in cell_contigs.index:

        chain = df.loc[j, "chain"]
        v_gene = df.loc[j, "v_gene"]
        j_gene = df.loc[j, "j_gene"]
        cdr3 = df.loc[j, "cdr3"]

        print("  CONTIG:", chain, v_gene, j_gene, cdr3)

        # Define clonotype using TRB cells
        if chain == "TRB":
            clonotype = v_gene + "_" + j_gene + "_" + cdr3
            clone_dict[clonotype].append(cell_barcode)

for clonotype, barcodes in list(clone_dict.items())[:3]:
    print(clonotype, ":", barcodes)
```

---

## Part 4: Clonotype dictionary using both TRA + TRB, with filtering rules

```python
############### Getting to a clonotype dictionary using both TRA and TRB
# rules for defining the clonotype
    # 1) same V gene
    # 2) same J gene
    # 3) same CDR3
    # filtering rules:
    # 1) only 1 beta and 1 alpha
    # 2) break ties by selecting the highest UMI

from collections import defaultdict

clone_dict = defaultdict(list)

for barcode in df["barcode"].unique():

    # Get all contigs belonging to this cell
    cell_contigs = df[df["barcode"] == barcode]

    # Split into alpha and beta contigs
    alpha_contigs = cell_contigs[cell_contigs["chain"] == "TRA"]
    beta_contigs = cell_contigs[cell_contigs["chain"] == "TRB"]

    # Filtering rule: only keep cells that have at least one of each chain
    if len(alpha_contigs) == 0 or len(beta_contigs) == 0:
        continue

    # Filtering rule: break ties by keeping the contig with the highest UMI count
    alpha_contig = alpha_contigs.sort_values("umis", ascending=False).iloc[0]
    beta_contig = beta_contigs.sort_values("umis", ascending=False).iloc[0]

    # Define clonotype using V gene + J gene + CDR3 from both chains
    clonotype = (
        alpha_contig["v_gene"] + "_" + alpha_contig["j_gene"] + "_" + alpha_contig["cdr3"]
        + "__"
        + beta_contig["v_gene"] + "_" + beta_contig["j_gene"] + "_" + beta_contig["cdr3"]
    )

    clone_dict[clonotype].append(barcode)

for clonotype, barcodes in list(clone_dict.items())[:3]:
    print(clonotype, ":", barcodes)
```

**Note:** this assumes a UMI count column named `umis`. If your dataframe uses a different name (e.g. `umi_count`, `reads`), swap it into the `.sort_values(...)` calls.
