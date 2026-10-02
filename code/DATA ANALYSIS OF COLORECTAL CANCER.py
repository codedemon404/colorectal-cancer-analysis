#         ============================================================
#                  Computational Analysis of Colorectal Cancer
#                        Author: DHIRAJ KUMAR BISWAL
#                        github: codedemon404
#                               Year: 2026
#
#                     Copyright © 2026 DHIRAJ KUMAR BISWAL.
#
#               This code is part of an academic research project.
#         ============================================================
#                                  Workflow:
#                   1. Load GSE74602 expression data
#                   2. Identify matched tumor-normal samples
#                   3. Map microarray probes to gene names
#                   4. Perform paired differential expression analysis
#                   5. Apply FDR correction
#                   6. Generate differential-expression results
#                   7. Generate volcano plot
#            ============================================================
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import re
from pathlib import Path

from scipy.stats import ttest_rel
from statsmodels.stats.multitest import multipletests


# ============================================================
# 1. FILE PATHS
# ============================================================

# Project folder = folder containing this Python script
PROJECT_PATH = Path(__file__).resolve().parent.parent

# Input and output folders
DATA_PATH = PROJECT_PATH / "data"
RESULTS_PATH = PROJECT_PATH / "results"
FIGURES_PATH = PROJECT_PATH / "figures"

# Create output folders if they don't already exist
RESULTS_PATH.mkdir(exist_ok=True)
FIGURES_PATH.mkdir(exist_ok=True)

# Input files
EXPRESSION_FILE = DATA_PATH / "GSE74602_series_matrix.txt"
ANNOTATION_FILE = DATA_PATH / "GPL6104-11576.txt"

# Output files
MAPPED_FILE = RESULTS_PATH / "GSE74602_mapped_expression.csv"
RESULTS_FILE = RESULTS_PATH / "GSE74602_differential_expression.csv"
SIGNIFICANT_FILE = RESULTS_PATH / "GSE74602_significant_genes.csv"

VOLCANO_FILE = FIGURES_PATH / "GSE74602_volcano_plot.png"


# ============================================================
# 2. READ EXPRESSION MATRIX
# ============================================================

print("\n" + "=" * 70)
print("STEP 1 — READING EXPRESSION MATRIX")
print("=" * 70)

data = pd.read_csv(
    EXPRESSION_FILE,
    sep="\t",
    comment="!",
    index_col=0
)


# Clean column names

data.columns = [
    str(x).strip().strip('"')
    for x in data.columns
]


# Clean probe IDs

data.index = (
    data.index
    .astype(str)
    .str.strip()
    .str.strip('"')
)


print("Expression matrix shape:", data.shape)

print("\nFirst 5 probes:")
print(data.index[:5].tolist())

print("\nFirst 5 samples:")
print(data.columns[:5].tolist())


# ============================================================
# 3. READ GEO SAMPLE INFORMATION
# ============================================================

print("\n" + "=" * 70)
print("STEP 2 — READING SAMPLE INFORMATION")
print("=" * 70)


samples = []
source_names = []
titles = []


with open(
    EXPRESSION_FILE,
    "r",
    encoding="utf-8",
    errors="replace"
) as f:

    for line in f:

        if line.startswith("!Sample_geo_accession"):

            samples = [
                x.strip().strip('"')
                for x in line.rstrip("\n").split("\t")[1:]
            ]

        elif line.startswith("!Sample_source_name_ch1"):

            source_names = [
                x.strip().strip('"')
                for x in line.rstrip("\n").split("\t")[1:]
            ]

        elif line.startswith("!Sample_title"):

            titles = [
                x.strip().strip('"')
                for x in line.rstrip("\n").split("\t")[1:]
            ]


print("Number of samples:", len(samples))
print("Number of source names:", len(source_names))
print("Number of titles:", len(titles))


# ============================================================
# 4. CHECK SAMPLE INFORMATION
# ============================================================

print("\nFirst 10 samples:")

for sample, source in zip(
    samples[:10],
    source_names[:10]
):

    print(
        sample,
        "-->",
        source
    )


# ============================================================
# 5. IDENTIFY TUMOR AND NORMAL SAMPLES
# ============================================================

print("\n" + "=" * 70)
print("STEP 3 — IDENTIFYING TUMOR AND NORMAL SAMPLES")
print("=" * 70)


tumor_samples = {}
normal_samples = {}


for sample, source in zip(
    samples,
    source_names
):

    # Remove quotation marks
    source = source.strip().strip('"')


    # Look for Primary Tumor T<number>

    tumor_match = re.search(
        r"Primary Tumor\s+T(\d+)",
        source,
        re.IGNORECASE
    )


    # Look for Normal N<number>

    normal_match = re.search(
        r"Normal\s+N(\d+)",
        source,
        re.IGNORECASE
    )


    if tumor_match:

        patient_number = int(
            tumor_match.group(1)
        )

        tumor_samples[patient_number] = sample


    elif normal_match:

        patient_number = int(
            normal_match.group(1)
        )

        normal_samples[patient_number] = sample


print(
    "Tumor samples found:",
    len(tumor_samples)
)

print(
    "Normal samples found:",
    len(normal_samples)
)


# ============================================================
# 6. CREATE VERIFIED TUMOR-NORMAL PAIRS
# ============================================================

print("\n" + "=" * 70)
print("STEP 4 — CREATING MATCHED PAIRS")
print("=" * 70)


pairs = []


for patient in sorted(tumor_samples):

    if patient in normal_samples:

        pairs.append({
            "Patient": patient,
            "Tumor": tumor_samples[patient],
            "Normal": normal_samples[patient]
        })


pairs_df = pd.DataFrame(pairs)


print(
    "Number of verified tumor-normal pairs:",
    len(pairs_df)
)


print("\nTumor-normal pairs:")

print(
    pairs_df.to_string(index=False)
)


# Make sure we have the expected number

if len(pairs_df) == 0:

    raise ValueError(
        "No tumor-normal pairs were found. "
        "Check the sample source names."
    )


# ============================================================
# 7. CHECK THAT ALL PAIR SAMPLE NAMES EXIST
# ============================================================

print("\n" + "=" * 70)
print("STEP 5 — CHECKING SAMPLE NAMES")
print("=" * 70)


expression_columns = set(
    data.columns
)


missing_samples = []


for _, pair in pairs_df.iterrows():

    if pair["Tumor"] not in expression_columns:

        missing_samples.append(
            pair["Tumor"]
        )

    if pair["Normal"] not in expression_columns:

        missing_samples.append(
            pair["Normal"]
        )


if len(missing_samples) > 0:

    print("Missing sample names:")

    print(missing_samples)

    raise ValueError(
        "Some paired samples are not present "
        "in the expression matrix."
    )


print("All paired sample names found.")

print(
    "Number of expression columns:",
    len(data.columns)
)


# ============================================================
# 8. READ PLATFORM ANNOTATION
# ============================================================

print("\n" + "=" * 70)
print("STEP 6 — READING PLATFORM ANNOTATION")
print("=" * 70)


with open(
    ANNOTATION_FILE,
    "r",
    encoding="utf-8",
    errors="replace"
) as f:

    annotation_lines = f.readlines()


print(
    "Total annotation lines:",
    len(annotation_lines)
)


# ============================================================
# 9. FIND FIRST ILMN PROBE
# ============================================================

first_data_line = None


for i, line in enumerate(annotation_lines):

    parts = line.rstrip("\n").split("\t")

    if len(parts) == 0:
        continue


    first_field = (
        parts[0]
        .strip()
        .strip('"')
    )


    if re.match(
        r"^ILMN_\d+$",
        first_field
    ):

        first_data_line = i
        break


if first_data_line is None:

    raise ValueError(
        "Could not find ILMN probe IDs "
        "in the annotation file."
    )


print(
    "First ILMN probe found at line:",
    first_data_line + 1
)


# ============================================================
# 10. FIND ANNOTATION HEADER
# ============================================================

header_line = first_data_line - 1


while header_line >= 0:

    parts = (
        annotation_lines[header_line]
        .rstrip("\n")
        .split("\t")
    )


    if len(parts) > 1:

        break


    header_line -= 1


if header_line < 0:

    raise ValueError(
        "Could not find annotation header."
    )


header = [
    x.strip().strip('"')
    for x in annotation_lines[
        header_line
    ].rstrip("\n").split("\t")
]


print(
    "Annotation header found at line:",
    header_line + 1
)


print("\nAnnotation columns:")

for i, column in enumerate(header):

    print(
        i,
        "-->",
        column
    )


# ============================================================
# 11. READ ANNOTATION DATA
# ============================================================

annotation_rows = []


for line in annotation_lines[first_data_line:]:

    parts = (
        line.rstrip("\n")
        .split("\t")
    )


    if len(parts) == 0:
        continue


    probe = (
        parts[0]
        .strip()
        .strip('"')
    )


    # Keep only ILMN probe rows

    if not re.match(
        r"^ILMN_\d+$",
        probe
    ):

        continue


    # Make every row same length as header

    if len(parts) < len(header):

        parts += [""] * (
            len(header) - len(parts)
        )


    elif len(parts) > len(header):

        parts = parts[:len(header)]


    annotation_rows.append(parts)


annotation = pd.DataFrame(
    annotation_rows,
    columns=header
)


print(
    "\nAnnotation table shape:",
    annotation.shape
)


# ============================================================
# 12. FIND PROBE COLUMN
# ============================================================

probe_column = None


for column in annotation.columns:

    name = (
        str(column)
        .strip()
        .lower()
        .replace("_", " ")
    )


    if name in [
        "id",
        "probe",
        "probe id"
    ]:

        probe_column = column
        break


if probe_column is None:

    # First column should contain ILMN IDs
    probe_column = annotation.columns[0]


print(
    "\nProbe column:",
    probe_column
)


# ============================================================
# 13. FIND GENE SYMBOL COLUMN
# ============================================================

gene_column = None


for column in annotation.columns:

    name = (
        str(column)
        .strip()
        .lower()
        .replace("_", " ")
    )


    if (
        name == "symbol"
        or
        "gene symbol" in name
    ):

        gene_column = column
        break


if gene_column is None:

    print("\nAvailable annotation columns:")

    for i, column in enumerate(
        annotation.columns
    ):

        print(
            i,
            "-->",
            column
        )


    raise ValueError(
        "Could not find the gene symbol column."
    )


print(
    "Gene symbol column:",
    gene_column
)


# ============================================================
# 14. CREATE PROBE → GENE MAPPING
# ============================================================

probe_to_gene = annotation[
    [
        probe_column,
        gene_column
    ]
].copy()


probe_to_gene.columns = [
    "Probe",
    "Gene"
]


# Clean values

probe_to_gene["Probe"] = (
    probe_to_gene["Probe"]
    .astype(str)
    .str.strip()
    .str.strip('"')
)


probe_to_gene["Gene"] = (
    probe_to_gene["Gene"]
    .astype(str)
    .str.strip()
    .str.strip('"')
)


# Remove empty/missing gene names

bad_gene_names = [
    "",
    "NA",
    "NaN",
    "nan",
    "---",
    "None"
]


probe_to_gene = probe_to_gene[
    ~probe_to_gene["Gene"].isin(
        bad_gene_names
    )
]


# Remove duplicate probe IDs

probe_to_gene = (
    probe_to_gene
    .drop_duplicates(
        subset="Probe"
    )
)


print("\nFirst 10 probe-gene mappings:")

print(
    probe_to_gene
    .head(10)
    .to_string(index=False)
)


print(
    "\nNumber of mapped probes:",
    len(probe_to_gene)
)


# ============================================================
# 15. MERGE EXPRESSION WITH GENE ANNOTATION
# ============================================================

print("\n" + "=" * 70)
print("STEP 7 — MAPPING PROBES TO GENES")
print("=" * 70)


expression_data = data.copy()


expression_data.index.name = "Probe"


expression_data = (
    expression_data
    .reset_index()
)


merged_data = expression_data.merge(
    probe_to_gene,
    on="Probe",
    how="left"
)


print(
    "Expression rows:",
    len(expression_data)
)

print(
    "Merged rows:",
    len(merged_data)
)

print(
    "Rows with gene names:",
    merged_data["Gene"].notna().sum()
)


# ============================================================
# 16. SAVE MAPPED EXPRESSION DATA
# ============================================================

merged_data.to_csv(
    MAPPED_FILE,
    index=False
)


print(
    "\nMapped data saved to:"
)

print(MAPPED_FILE)


# ============================================================
# 17. CHECK EXPRESSION VALUES
# ============================================================

print("\n" + "=" * 70)
print("STEP 8 — CHECKING EXPRESSION VALUES")
print("=" * 70)


# Expression columns are the GSM columns

expression_columns = [
    column
    for column in merged_data.columns
    if str(column).startswith("GSM")
]


expression_values = (
    merged_data[
        expression_columns
    ]
    .apply(
        pd.to_numeric,
        errors="coerce"
    )
)


print(
    "Number of expression columns:",
    len(expression_columns)
)

print(
    "Minimum:",
    expression_values.min().min()
)

print(
    "Maximum:",
    expression_values.max().max()
)

print(
    "Mean:",
    expression_values.mean().mean()
)


# ============================================================
# 18. PAIRED DIFFERENTIAL EXPRESSION
# ============================================================

print("\n" + "=" * 70)
print("STEP 9 — PAIRED DIFFERENTIAL EXPRESSION")
print("=" * 70)


results = []


for index, row in merged_data.iterrows():

    tumor_values = []
    normal_values = []


    # --------------------------------------------------------
    # Get the 30 tumor-normal measurements
    # --------------------------------------------------------

    for _, pair in pairs_df.iterrows():

        tumor_sample = pair["Tumor"]
        normal_sample = pair["Normal"]


        tumor_value = pd.to_numeric(
            row[tumor_sample],
            errors="coerce"
        )


        normal_value = pd.to_numeric(
            row[normal_sample],
            errors="coerce"
        )


        tumor_values.append(
            tumor_value
        )

        normal_values.append(
            normal_value
        )


    tumor_values = np.array(
        tumor_values,
        dtype=float
    )


    normal_values = np.array(
        normal_values,
        dtype=float
    )


    # --------------------------------------------------------
    # Remove missing pairs
    # --------------------------------------------------------

    valid = (
        ~np.isnan(tumor_values)
        &
        ~np.isnan(normal_values)
    )


    tumor_clean = tumor_values[valid]

    normal_clean = normal_values[valid]


    # --------------------------------------------------------
    # Need at least 3 paired observations
    # --------------------------------------------------------

    if len(tumor_clean) < 3:

        results.append({

            "Probe":
                row["Probe"],

            "Gene":
                row["Gene"],

            "log2FC":
                np.nan,

            "p_value":
                np.nan,

            "n_pairs":
                len(tumor_clean)
        })

        continue


    # --------------------------------------------------------
    # Paired t-test
    # --------------------------------------------------------

    statistic, p_value = ttest_rel(
        tumor_clean,
        normal_clean
    )


    # --------------------------------------------------------
    # log2 fold change
    #
    # The expression values are already on a log2 scale.
    #
    # Therefore:
    #
    # log2FC =
    # mean(tumor expression - normal expression)
    # --------------------------------------------------------

    log2FC = np.mean(
        tumor_clean - normal_clean
    )


    results.append({

        "Probe":
            row["Probe"],

        "Gene":
            row["Gene"],

        "log2FC":
            log2FC,

        "p_value":
            p_value,

        "n_pairs":
            len(tumor_clean)
    })


# ============================================================
# 19. CREATE RESULTS TABLE
# ============================================================

results_df = pd.DataFrame(
    results
)


print(
    "\nTotal probes tested:",
    len(results_df)
)


# ============================================================
# 20. MULTIPLE TESTING CORRECTION
# ============================================================

print("\nApplying Benjamini-Hochberg FDR...")


valid = results_df[
    "p_value"
].notna()


results_df[
    "adjusted_p_value"
] = np.nan


results_df.loc[
    valid,
    "adjusted_p_value"
] = multipletests(
    results_df.loc[
        valid,
        "p_value"
    ],
    method="fdr_bh"
)[1]


# ============================================================
# 21. DETERMINE DIRECTION
# ============================================================

results_df["Direction"] = np.where(

    results_df["log2FC"] > 0,

    "Up in Tumor",

    "Down in Tumor"
)


# ============================================================
# 22. SORT BY FDR
# ============================================================

results_df = results_df.sort_values(
    "adjusted_p_value"
)


# ============================================================
# 23. SIGNIFICANT RESULTS
# ============================================================

significant = results_df[
    (
        results_df["adjusted_p_value"]
        < 0.05
    )
]


print("\n" + "=" * 70)
print("DIFFERENTIAL EXPRESSION SUMMARY")
print("=" * 70)


print(
    "Total probes tested:",
    len(results_df)
)


print(
    "Significant probes (FDR < 0.05):",
    len(significant)
)


print(
    "Up in tumor:",
    (
        significant["log2FC"] > 0
    ).sum()
)


print(
    "Down in tumor:",
    (
        significant["log2FC"] < 0
    ).sum()
)


# ============================================================
# 24. SHOW TOP 20 RESULTS
# ============================================================

print("\n" + "=" * 70)
print("TOP 20 DIFFERENTIALLY EXPRESSED PROBES")
print("=" * 70)


print(
    results_df[
        [
            "Probe",
            "Gene",
            "log2FC",
            "p_value",
            "adjusted_p_value",
            "Direction"
        ]
    ]
    .head(20)
    .to_string(index=False)
)


# ============================================================
# 25. SAVE ALL DIFFERENTIAL EXPRESSION RESULTS
# ============================================================

results_df.to_csv(
    RESULTS_FILE,
    index=False
)


print("\nAll results saved to:")

print(RESULTS_FILE)


# ============================================================
# 26. SAVE SIGNIFICANT RESULTS
# ============================================================

significant.to_csv(
    SIGNIFICANT_FILE,
    index=False
)


print("\nSignificant results saved to:")

print(SIGNIFICANT_FILE)


# ============================================================
# 27. VOLCANO PLOT
# ============================================================

print("\n" + "=" * 70)
print("STEP 10 — CREATING VOLCANO PLOT")
print("=" * 70)


plot_data = results_df.copy()


# Avoid log10(0)

plot_data["plot_p"] = (
    plot_data["adjusted_p_value"]
    .clip(lower=1e-300)
)


plot_data["neg_log10_FDR"] = (
    -np.log10(
        plot_data["plot_p"]
    )
)


# Significant points

is_significant = (
    plot_data["adjusted_p_value"] < 0.05
)


# Create plot

plt.figure(
    figsize=(10, 7)
)


plt.scatter(
    plot_data.loc[
        ~is_significant,
        "log2FC"
    ],
    plot_data.loc[
        ~is_significant,
        "neg_log10_FDR"
    ],
    s=8,
    alpha=0.5
)


plt.scatter(
    plot_data.loc[
        is_significant,
        "log2FC"
    ],
    plot_data.loc[
        is_significant,
        "neg_log10_FDR"
    ],
    s=10,
    alpha=0.7
)


plt.axhline(
    -np.log10(0.05),
    linestyle="--"
)


plt.axvline(
    0,
    linestyle="--"
)


plt.xlabel(
    "log2 Fold Change"
)


plt.ylabel(
    "-log10 Adjusted P-value"
)


plt.title(
    "GSE74602: Tumor vs Matched Normal"
)


plt.tight_layout()


plt.savefig(
    VOLCANO_FILE,
    dpi=300
)


plt.show()


print(
    "\nVolcano plot saved to:"
)

print(VOLCANO_FILE)


# ============================================================
# 28. FINISHED
# ============================================================

print("\n" + "=" * 70)
print("✓ ANALYSIS COMPLETE")
print("=" * 70)

print("\nFiles created:")

print("1.", MAPPED_FILE)
print("2.", RESULTS_FILE)
print("3.", SIGNIFICANT_FILE)
print("4.", VOLCANO_FILE)

print("\n✓ Done.")