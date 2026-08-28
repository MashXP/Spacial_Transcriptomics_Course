# Module 1c: R Graphics and Statistics

> 🚀 **Interactive Google Colab Notebook:** [01c_R_Graphics_and_Statistics.ipynb](./01c_R_Graphics_and_Statistics.ipynb)

Imagine you have successfully cleaned your spatial transcriptomics dataset. Now, you need to present your findings to your team. Staring at columns of numbers doesn't tell a story. More importantly, how do you mathematically prove that the differences in gene expression between tumor and healthy cells are real and not just random noise?

This module covers **Graphics and Statistics in R**. You will learn how to turn raw numbers into beautiful, publication-ready plots (using `ggplot2`) and how to run rigorous statistical tests to validate your biological conclusions.

---

## 1. Data Reshaping (Wide vs. Long Format)

Most scientific plotting packages (like `ggplot2`) require data in **long** (tidy) format.
*   **Wide Format**: Each variable/sample has its own column.
*   **Long Format**: One column for variable types/labels, and one column for values.

Use `reshape2::melt()` to reshape data frames from wide to long:

```R
# Wide data frame
wide_df <- data.frame(
  gene = c("GAPDH", "ACTB"),
  Sample_01 = c(12.5, 44.1),
  Sample_02 = c(14.0, 42.0)
)

# Convert to Long format
library(reshape2)
long_df <- melt(wide_df, id.vars = "gene", variable.name = "sample", value.name = "tpm")
print(long_df)
# Output columns: gene, sample, tpm
```

---

## 2. R Graphics & Plotting

Data visualization follows a 7-step workflow: ingestion -> variable typing -> filtering -> stats check -> plot selection -> aesthetic refinement -> export.

### A. Base R Plotting
Quick, built-in functions for exploratory plots:

```R
# Scatterplot
plot(x = 1:5, y = c(2, 4, 6, 8, 10), main = "Scatterplot", xlab = "X", ylab = "Y", pch = 16, col = "blue")

# Line Chart (type: "s" for steps, "b" for both points and lines)
plot(1:5, c(2, 4, 6, 8, 10), type = "b", lty = 1, lwd = 2, col = "red", main = "Line Chart")

# Dot Chart (visualizing individual data points)
dotchart(c(12, 15, 20), labels = c("S1", "S2", "S3"), main = "Dot Plot")

# Grouped Boxplot using Formula syntax (y ~ x) (notch = TRUE adds median confidence intervals)
boxplot(tpm ~ gene, data = long_df, col = c("salmon", "skyblue"), notch = TRUE, main = "Expression Profiles")

# Histograms with Density Curve Overlay
hist(long_df$tpm, prob = TRUE, breaks = 5, main = "TPM Density", xlab = "TPM", col = "lightgrey")
lines(density(long_df$tpm), col = "blue", lwd = 2) # Overlays a smooth probability density curve
```

### B. Graphics Devices (Exporting Plots)
To write plots directly to disk, open a graphics device before running the plot code, and close it using `dev.off()`:

```R
# 1. Open device
pdf("gene_boxplot.pdf", width = 6, height = 4)

# 2. Run plot code
boxplot(tpm ~ gene, data = long_df, main = "Boxplot", col = "orange")

# 3. Close device to save file
dev.off()
```

### C. Modern Visualizations with ggplot2
`ggplot2` utilizes the **Grammar of Graphics**: layered specifications of Data, Aesthetic Mappings (`aes()` mapping columns to X, Y, Color, or Shape), and Geometric layers (`geom_point`, `geom_boxplot`, etc.).

```R
library(ggplot2)

# 1. Scatterplot with Trendline and Alpha for Transparency
ggplot(data = long_df, aes(x = gene, y = tpm, color = sample)) +
  geom_point(size = 3, alpha = 0.6) + # alpha prevents overplotting issues
  geom_smooth(method = "lm", se = FALSE) + # Linear model trendline
  labs(title = "Gene Expression Comparison", x = "Gene", y = "TPM") +
  theme_classic()

# 2. Faceting
ggplot(long_df, aes(x = sample, y = tpm)) +
  geom_bar(stat = "identity") +
  facet_wrap(~ gene) # Split into separate panels per gene

# 3. Volcano Plot Code with FDR (Benjamini-Hochberg) Adjustment
# Generate dummy dataset for demonstration
de_results <- data.frame(
  gene = paste0("Gene_", 1:100),
  log2FC = rnorm(100, mean = 0, sd = 2),
  pvalue = runif(100, min = 0, max = 0.05)
)
# Reconcile false discoveries by adjusting p-values
de_results$padj <- p.adjust(de_results$pvalue, method = "BH")

# Volcano plot mapping
ggplot(de_results, aes(x = log2FC, y = -log10(padj))) +
  geom_point(aes(color = padj < 0.05 & abs(log2FC) > 1)) +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed") +
  geom_vline(xintercept = c(-1, 1), linetype = "dashed") +
  theme_minimal()

# 4. Heatmaps & matrix generation with zero offsets
# Handle zeros by adding a small constant offset to avoid log(0)
expr_matrix <- matrix(log2(c(12, 0.1, 45, 0.1, 80, 100) + 0.1), nrow = 2)
pheatmap::pheatmap(expr_matrix, cluster_rows = TRUE, cluster_cols = TRUE)
```
---

## 3. R Basic & Descriptive Statistics

Descriptive metrics summarize the central values and dispersion of data distributions.

*   **Central Tendency**: `mean()`, `median()`, and custom Mode helper functions.
*   **Dispersion**: `var()` (variance), `sd()` (standard deviation), and `range()`.

```R
values <- c(2, 5, 8, 12, 12, 25)

mean(values)   # 10.66667
median(values) # 10.0
sd(values)     # 7.865537
```

---

## 4. Assumptions Auditing (Normality Checking)

> **The Problem:** You want to run a standard t-test to check if gene expression differs significantly between two cell types. However, your data is heavily skewed and contains major outliers, meaning a standard t-test will output incorrect, invalid p-values.
>
> **The Solution:** **Normality checking** (Shapiro-Wilk and visual Q-Q plots) allows you to audit the shape of your data first, so you can choose the correct parametric or non-parametric test.

Before choosing a statistical test, you **MUST** verify the distribution assumptions of the dataset.

*   **Visual Check**: Histogram shapes or Q-Q plots (`qqnorm()`, `qqline()`).
*   **Computational Check**: Shapiro-Wilk Test (`shapiro.test()`).
    *   $H_0$: The data is normally distributed.
    *   If $p$-value $> 0.05$: Accept normality (use **parametric** tests).
    *   If $p$-value $\le 0.05$: Reject normality (use **non-parametric** tests).

```R
# Shapiro-Wilk normality test
shapiro.test(values)
# If p-value = 0.449 (> 0.05), data is normally distributed
```

---

## 5. Statistical Test Decision Framework

Choose the correct test based on experimental design:

```
                  +-----------------------------------+
                  |      What is your Data Type?      |
                  +-----------------------------------+
                                    |
         +--------------------------+--------------------------+
         |                                                     |
   [ Numerical ]                                         [ Categorical ]
         |                                                     |
  Normality Audit?                                       Group Counts?
         |                                                     |
    +----+----+                                           +----+----+
    |         |                                           |         |
 [ Yes ]    [ No ]                                    [ 2x2 ]    [ Large ]
(Parametric)(Non-Parametric)                          (Fisher's) (Chi-Sq)
```

### A. Group Comparisons (Differences)

| Numerical Groups | Parametric Test (Normal) | Non-Parametric Test (Non-normal) |
| :--- | :--- | :--- |
| **2 Groups (Independent)** | Student's t-test: `t.test(y ~ x)` | Wilcoxon Rank-Sum: `wilcox.test(y ~ x)` |
| **2 Groups (Paired / Before-After)**| Paired t-test: `t.test(y1, y2, paired=TRUE)` | Wilcoxon Signed-Rank: `wilcox.test(y1, y2, paired=TRUE)`|
| **> 2 Groups (1-way)** | One-way ANOVA: `aov(y ~ x)` | Kruskal-Wallis: `kruskal.test(y ~ x)` |
| **> 2 Groups (2-way)** | Two-way ANOVA: `aov(y ~ x1 + x2)` | Scheirer-Ray-Hare (advanced extension) |

### B. Categorical Frequency Comparisons
*   **Chi-Square Test**: `chisq.test(table)` (compares frequencies in large groups).
*   **Fisher's Exact Test**: `fisher.test(table)` (preferred for small sample sizes).

---

## 6. Correlation & Regression Modeling

### A. Correlation
Correlation measures the strength and direction of the linear relationship between two variables.
*   **Pearson**: Parametric correlation (expects normal distribution).
*   **Spearman / Kendall**: Non-parametric rank-based correlation.

```R
x <- c(1, 2, 3, 4, 5)
y <- c(2, 4, 5, 4, 5)

# Calculate correlation coefficient
cor(x, y, method = "pearson")  # Output: 0.88
cor.test(x, y, method = "pearson") # Provides p-value verification
```

### B. Linear Regression
Regression models the relationship between dependent variable $Y$ and one or more independent variables $X$.

```R
# Simple Linear Regression (Y ~ X)
fit <- lm(y ~ x)
summary(fit) # Yields coefficients, R-squared value, and p-values

# Multiple Linear Regression (Y ~ X1 + X2)
multi_fit <- lm(y ~ x + c(12, 10, 8, 9, 7))
summary(multi_fit)
```

---

## 7. Unsupervised Learning & Clustering

> **The Problem:** You have single-cell data with 20,000 genes. Trying to visualize this in a multi-dimensional scatter plot creates a massive, unreadable cloud of points. You need to identify cell groupings but have no pre-existing cell labels.
>
> **The Solution:** **Dimensionality reduction (PCA)** reduces data complexity to 2D/3D projections, and **unsupervised clustering** (K-Means/Hierarchical) groups cells based on distance metrics to discover cell types.

Bioinformatics datasets with thousands of genes require dimensionality reduction and grouping.

### A. Principal Component Analysis (PCA)
PCA simplifies complex high-dimensional datasets. Crucial: always scale the data before performing PCA so large-value columns do not dominate variance.

```R
# Setup dummy gene expression matrix (columns are genes, rows are cells)
data_matrix <- matrix(rnorm(50), nrow = 5, ncol = 10)

# Run PCA (scale = TRUE is essential)
pca_res <- prcomp(data_matrix, scale = TRUE)

# View variance explained
summary(pca_res)

# Plot PCA results
plot(pca_res$x[, 1], pca_res$x[, 2], xlab = "PC1", ylab = "PC2", main = "PCA Cell Projections")
```

### B. Unsupervised Clustering
*   **Hierarchical Clustering**: Groups elements based on distance metrics (e.g., Euclidean distance).
*   **K-Means Clustering**: Iteratively partitions cells into K clusters based on centroids.

```R
# 1. Hierarchical Clustering
dist_matrix <- dist(data_matrix) # Euclidean distance matrix
h_clust <- hclust(dist_matrix)
plot(h_clust, main = "Hierarchical Clustering Dendrogram")

# 2. K-Means Clustering
km_clust <- kmeans(data_matrix, centers = 2)
km_clust$cluster # View cell assignments
```

---

## 8. Proactive vs. Reactive Programming

In bioinformatics pipelines, write defensive (proactive) code rather than correcting errors later (reactive):
*   **Proactive**: Define column names, check dimensions, and handle `NA` values immediately upon data ingestion (e.g., using `fread` and filtering).
*   **Reactive**: Fixing mismatching shapes or missing values after multiple downstream pipelines have already failed.

---

## 🧑‍💻 Practice Exercises

### Exercise 1: Reshaping & Plotting
A wide dataset of gene expression values:
```R
wide_expr <- data.frame(
  gene = c("SOX2", "POU5F1"),
  Control = c(4.2, 5.0),
  Treated = c(12.8, 14.5)
)
```
1. Reshape the dataset into **long** format using `melt()` (name variable column "condition" and value column "expression").
2. Write a script to export a bar plot of expression values grouped by condition using a graphics device saving to "expression_bars.pdf".

### Exercise 2: Normality Audit & Hypothesis Testing
You are testing if the gene expression of `SOX2` differs significantly between two groups:
*   Group A: `c(12, 15, 14, 11, 13)`
*   Group B: `c(22, 28, 25, 29, 24)`
1. Check the normality assumption of Group A using `shapiro.test()`.
2. Based on the $p$-value result, choose and execute the correct comparison test to evaluate if Group B's expression is significantly higher than Group A's.

---

### Solutions

<details>
<summary>Click to view solutions</summary>

#### Solution 1:
```R
library(reshape2)
# 1. Reshape
long_expr <- melt(wide_expr, id.vars = "gene", variable.name = "condition", value.name = "expression")

# 2. Export Plot
pdf("expression_bars.pdf", width = 5, height = 4)
barplot(expression ~ condition + gene, data = long_expr, 
        beside = TRUE, col = c("blue", "red"),
        main = "Gene Expression by Condition")
dev.off()
```

#### Solution 2:
```R
group_a <- c(12, 15, 14, 11, 13)
group_b <- c(22, 28, 25, 29, 24)

# 1. Normality test
shapiro.test(group_a) # p-value = 0.941 (> 0.05, Normally distributed)
shapiro.test(group_b) # p-value = 0.793 (> 0.05, Normally distributed)

# 2. Execute parametric Student's independent t-test
t.test(group_b, group_a, alternative = "greater")
# Output will confirm highly significant difference (p < 0.01)
```
</details>

---

### Next Modules
- Predecessor: [Module 1a: R Programming Fundamentals](./01a_R_Syntax_and_Variables.md)
- Predecessor: [Module 1b: R Data Manipulation](./01b_R_Data_Manipulation.md)
