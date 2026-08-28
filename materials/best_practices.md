# **Spatial Transcriptomics** Best Practices and Prerequisites

Spatial transcriptomics best practice centers on **matching the platform to the biological question**, protecting sample quality from collection through sectioning, and planning the computational workflow before data generation  (Grases & Porta-Pardo, 2025; Williams et al., 2022; Kleino et al., 2022). The main prerequisites are suitable tissue and preservation chemistry, clear experimental objectives, histology and pathology support, and analysis resources that can handle expression, image, and spatial coordinate data  (Lim et al., 2025; Khan et al., 2025; Kleino et al., 2022).

## Experimental Design

A first prerequisite is to define whether **spatial context is essential** for the biological question, rather than using simpler bulk or single-cell RNA-seq alone  (Wang et al., 2023). Platform choice should then be driven by spatial resolution, tissue area, sensitivity, gene coverage, sample number, and compatible preservation mode  (Williams et al., 2022; Lim et al., 2025; Robles-Remacho et al., 2023).

- Define the **biological question** before choosing a platform  (Lim et al., 2025; Wang et al., 2023).
- Record species, tissue size, cell targets, RNA integrity, and budget during planning  (Lim et al., 2025).
- Consider matched **scRNA-seq reference data** when the tissue can be dissociated  (Williams et al., 2022; Zhang et al., 2022).

## Tissue And Sample Preparation

Sample handling is a decisive determinant of data quality, especially in clinical material where harvesting, preservation, transport, storage, and section processing all affect success  (Liu et al., 2022). Tissue compatibility also constrains chemistry: fresh frozen samples often support oligo-dT whole-transcriptome capture, whereas FFPE samples usually require probe-based targeted panels because of RNA degradation  (Khan et al., 2025; Williams et al., 2022).

- Use **strict quality control** on handling and preservation before committing to expensive assays  (Mirzazadeh et al., 2022; Liu et al., 2022).
- For human lung samples, temporary preservation within **30 minutes** and in dry status improved sample quality; direct OCT cryopreservation was recommended  (Liu et al., 2022).
- Proper tissue orientation, sectioning expertise, and appropriate thickness matter because overlapping cells can obscure spatial relationships  (Khan et al., 2025).

## Data Processing And Analysis

Best practice is to plan preprocessing and downstream analysis before the experiment, because ST produces platform-specific raw data that must be converted into count matrices linked to spatial coordinates and often histology images  (Williams et al., 2022; Kleino et al., 2022). Standard steps include filtering low-quality spots, normalizing counts, reducing dimensionality, and then applying ST-aware methods for spatial domains, cell composition, alignment, or cell-cell communication  (Kleino et al., 2022; Guo et al., 2026).

- Do not assume **scRNA-seq tools** transfer directly, because ST has different data distributions and spatial structure  (Fang et al., 2022; Kleino et al., 2022).
- Use histology images and neighbor information when possible because they improve normalization and spatial inference  (Fang et al., 2022; Kleino et al., 2022; Dries et al., 2021).
- For clustering and integration, choose tools by tissue and platform rather than one universal method  (Kleino et al., 2022; Hu et al., 2024; Guo et al., 2026).

## Reproducibility And Special Cases

Reproducible ST workflows need standardized processing, transparent reporting, and benchmarking with shared reference tissues and harmonized preprocessing pipelines  (Fang et al., 2022; Ge et al., 2025; You et al., 2023). This remains a field-wide weakness: nearly half of tumor ST studies lacked comprehensive processing protocols, which limits reproducibility  (Maciejewski & Czerwińska, 2024).

- Use **reference tissues** and pre-experiments such as permeabilization optimization when the protocol requires them  (You et al., 2023).
- For degraded fresh frozen tissue, specialized rescue workflows can improve mRNA recovery and broaden usable specimens  (Mirzazadeh et al., 2022).
- Multi-slice alignment and 3D reconstruction require additional preprocessing and substantial computational optimization  (Khan et al., 2025).

Spatial transcriptomics prerequisites are therefore both **biological and operational**: a question that truly needs spatial information, tissue prepared for the chosen chemistry, and an analysis plan matched to the platform and study design. Best practice is not one protocol, but a coordinated workflow from sample acquisition to validated, reproducible spatial analysis.
 
_These search results were found and analyzed using Consensus, an AI-powered search engine for research. Try it at https://consensus.app. © 2026 Consensus NLP, Inc. Personal, non-commercial use only; redistribution requires copyright holders’ consent._
 
## References
 
Dries, R., Chen, J. G., Del Rossi, N., Khan, M. M., Sistig, A., & Yuan, G. (2021). Advances in spatial transcriptomic data analysis. *Genome Research, 31*, 1706 - 1718. https://doi.org/10.1101/gr.275224.121
 
Fang, S., Chen, B., Zhang, Y., Sun, H., Liu, L., Liu, S., Li, Y., & Xu, X. (2022). Computational Approaches and Challenges in Spatial Transcriptomics. *Genomics, Proteomics & Bioinformatics, 21*, 24 - 47. https://doi.org/10.1016/j.gpb.2022.10.001
 
Ge, Q., Sheng, Y., Shan, Y., Yang, Y., Jiang, H., & Wang, R. (2025). Enhancing RNA Capture Efficiency in Spatial Transcriptomics: A Review of Innovative Technologies and Strategies. *International Journal of Molecular Sciences, 26*. https://doi.org/10.3390/ijms262211076
 
Grases, D., & Porta-Pardo, E. (2025). A practical guide to spatial transcriptomics: lessons from over 1000 samples.. *Trends in biotechnology*. https://doi.org/10.1016/j.tibtech.2025.08.020
 
Guo, Z., Wu, R., Li, W., Yang, K., Ying, X., Alinejad-Rokny, H., & Ye, Y. (2026). Mapping biology in space: from spatial transcriptomics platforms to analytical tools and databases.. *Science bulletin*. https://doi.org/10.1016/j.scib.2026.01.034
 
Hu, Y., Li, Y., Xie, M., Rao, M., Shen, W., Luo, C., Qin, H., Baek, J., & Zhou, X. (2024). Benchmarking clustering, alignment, and integration methods for spatial transcriptomics. *Genome Biology, 25*. https://doi.org/10.1186/s13059-024-03361-0
 
Khan, M., Arslanturk, S., & Draghici, S. (2025). A comprehensive review of spatial transcriptomics data alignment and integration. *Nucleic Acids Research, 53*. https://doi.org/10.1093/nar/gkaf536
 
Kleino, I., Frolovaitė, P., Suomi, T., & Elo, L. (2022). Computational solutions for spatial transcriptomics. *Computational and Structural Biotechnology Journal, 20*, 4870 - 4884. https://doi.org/10.1016/j.csbj.2022.08.043
 
Lim, H., Wang, Y., Buzdin, A., & Li, X. (2025). A practical guide for choosing an optimal spatial transcriptomics technology from seven major commercially available options. *BMC Genomics, 26*. https://doi.org/10.1186/s12864-025-11235-3
 
Liu, X., Jiang, Y., Song, D., Zhang, L., Xu, G., Hou, R., Zhang, Y., Chen, J., Cheng, Y., Liu, L., Xu, X., Chen, G., Wu, D., Chen, T., Chen, A., & Wang, X. (2022). Clinical challenges of tissue preparation for spatial transcriptome. *Clinical and Translational Medicine, 12*. https://doi.org/10.1002/ctm2.669
 
Maciejewski, K., & Czerwińska, P. (2024). Scoping Review: Methods and Applications of Spatial Transcriptomics in Tumor Research. *Cancers, 16*. https://doi.org/10.3390/cancers16173100
 
Mirzazadeh, R., Andrusivová, Ž., Larsson, L., Newton, P. T., Galicia, L. A., Abalo, X. M., Avijgan, M., Kvastad, L., Denadai-Souza, A., Stakenborg, N., Firsova, A., Shamikh, A., Jurek, A., Schultz, N., Nistér, M., Samakovlis, C., Boeckxstaens, G., & Lundeberg, J. (2022). Spatially resolved transcriptomic profiling of degraded and challenging fresh frozen samples. *Nature Communications, 14*. https://doi.org/10.1038/s41467-023-36071-5
 
Robles-Remacho, A., Sanchez-Martin, R. M., & Díaz-Mochón, J. (2023). Spatial Transcriptomics: Emerging Technologies in Tissue Gene Expression Profiling. *Analytical Chemistry, 95*, 15450 - 15460. https://doi.org/10.1021/acs.analchem.3c02029
 
Wang, Y., Liu, B., Zhao, G., Lee, Y., Buzdin, A., Mu, X., Zhao, J. J., Chen, H., & Li, X. (2023). Spatial transcriptomics: Technologies, applications and experimental considerations. *Genomics, 115*, 110671 - 110671. https://doi.org/10.1016/j.ygeno.2023.110671
 
Williams, C. G., Lee, H. J., Asatsuma, T., Vento-Tormo, R., & Haque, A. (2022). An introduction to spatial transcriptomics for biomedical research. *Genome Medicine, 14*. https://doi.org/10.1186/s13073-022-01075-1
 
You, Y., Fu, Y., Li, L., Zhang, Z., Jia, S., Lu, S., Ren, W., Liu, Y., Xu, Y., Liu, X., Jiang, F., Peng, G., Kumar, A. S., Ritchie, M. E., Liu, X., & Tian, L. (2023). Systematic comparison of sequencing-based spatial transcriptomic methods. *Nature Methods, 21*, 1743 - 1754. https://doi.org/10.1038/s41592-024-02325-3
 
Zhang, L., Chen, D., Song, D., Liu, X., Zhang, Y., Xu, X., & Wang, X. (2022). Clinical and translational values of spatial transcriptomics. *Signal Transduction and Targeted Therapy, 7*. https://doi.org/10.1038/s41392-022-00960-w
 
