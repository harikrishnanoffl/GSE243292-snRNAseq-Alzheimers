# Biological Interpretation — GSE243292 Alzheimer's Disease snRNA-seq Analysis

## Overview

The GSE243292 single-nucleus RNA-seq analysis provides a multi-level view of cellular and molecular changes associated with Alzheimer's disease. The workflow examined cellular composition, cell-type-specific transcriptional changes, functional enrichment, an OPC-to-oligodendrocyte trajectory, and predicted cell-cell communication.

The interpretations below are intentionally cautious. Inferred trajectories and CellChat interactions are not equivalent to experimentally demonstrated lineage or communication events.

## 1. Cellular Composition

Eight populations were identified: Astrocytes, Endothelial cells, Excitatory neurons, Inhibitory neurons, Microglia, Oligodendrocytes, OPCs, and Unresolved cells.

The analyzed dataset contained 122,606 nuclei, including 61,552 Alzheimer's Disease nuclei, 19,892 Normal nuclei, and 41,162 Pathological Aging nuclei. These nuclei come from 15 donors distributed unevenly across groups — 8 Alzheimer's Disease, 5 Pathological Aging, and only 2 Normal donors. All statistical comparisons below are performed at the donor (pseudobulk) level rather than the nucleus level, so this donor count — not the much larger nucleus count — is the effective sample size, and the 2-donor Normal group in particular limits statistical power (see Interpretation Limitations, point 7).

Cell-type composition differed between disease groups. Oligodendrocytes represented approximately 39.50% of Alzheimer's Disease nuclei and 52.31% of Normal nuclei, while excitatory neurons represented approximately 24.67% and 19.76%, respectively.

These are relative proportions of captured nuclei and should not automatically be interpreted as absolute changes in cell numbers.

## 2. Disease-Associated Transcriptional Changes

The primary disease comparison was Alzheimer's Disease versus Normal using sample-level pseudobulk differential expression.

Significant genes were defined using FDR ≤ 0.05 and |log2FC| ≥ 0.25.

| Cell type | Significant genes | Up in AD | Down in AD |
|---|---:|---:|---:|
| Excitatory neurons | 24 | 4 | 20 |
| Inhibitory neurons | 10 | 0 | 10 |
| Astrocytes | 0 | 0 | 0 |
| OPC | 4 | 0 | 4 |
| Microglia | 3 | 0 | 3 |
| Oligodendrocytes | 3 | 0 | 3 |
| Endothelial | 0 | 0 | 0 |

A notable pattern is the predominance of downregulated genes among the significant results.

### Excitatory neurons

Twenty-four significant genes were identified in excitatory neurons. Among the strongest results were CRH, BX547991.1, GPC6-AS1, DND1, AC138956.1, and LINC02015. Several genes had negative log2 fold changes, while ADAMTS2, U2AF1, and HSD11B2 showed positive log2 fold changes.

### Inhibitory neurons

Ten significant genes were identified. CRH, VGF, and SST were among the strongest results and showed negative log2 fold changes.

These results indicate disease-associated transcriptional changes in both excitatory and inhibitory neuronal populations.

## 3. Microglia

Three significant genes were identified in microglia: CCL2, CD83, and CIRBP. All showed negative log2 fold changes in the Alzheimer's Disease versus Normal comparison.

CCL2 showed the strongest statistical signal among the microglial results, with an FDR of approximately 6.21 × 10^-7.

This identifies CCL2 as a prominent microglial disease-associated transcript in this analysis. The result should be interpreted specifically as lower expression in AD samples relative to Normal under the applied analysis, rather than as evidence of generalized microglial activation or suppression.

## 4. Oligodendrocyte Lineage

Disease-associated changes were observed in both OPCs and mature oligodendrocytes.

Four significant genes were detected in OPCs: OPALIN, TP53TG5, STMN1, and CIRBP. All showed negative log2 fold changes.

Three significant genes were detected in oligodendrocytes: CRH, CIRBP, and AC138956.1. These also showed negative log2 fold changes.

The presence of disease-associated changes in both precursor and mature oligodendrocyte populations motivated the trajectory analysis.

## 5. OPC-to-Oligodendrocyte Trajectory

Monocle3 was used to investigate an inferred transcriptional continuum between OPCs and oligodendrocytes.

The trajectory contained 6,794 OPCs and 43,121 oligodendrocytes. OPCs were selected as the root population.

| Cell type | Median pseudotime |
|---|---:|
| OPC | 0.027 |
| Oligodendrocytes | 2.251 |

The ordering followed the expected direction: **OPC → Oligodendrocyte**.

This supports an inferred transcriptional progression from the OPC state toward the oligodendrocyte state within the analyzed dataset.

Pseudotime is an inferred ordering of transcriptional states and should not be interpreted as direct experimental evidence that individual OPCs were observed differentiating into oligodendrocytes.

## 6. Functional Enrichment

GO Biological Process enrichment was performed on significant disease-associated genes.

The analysis produced:

- 79 AD-downregulated GO terms in excitatory neurons
- 98 AD-downregulated GO terms in inhibitory neurons
- 89 AD-downregulated GO terms in microglia

Several other populations had too few significant genes for enrichment under the selected criteria.

The absence of enrichment in a cell type should therefore not be interpreted as proof that the cell type is unaffected; it may reflect the limited number of significant genes available.

## 7. Cell-Cell Communication

CellChat was used to investigate predicted ligand-receptor communication between major cell populations in Alzheimer's Disease and Normal samples.

Seven annotated populations were analyzed, with Unresolved cells excluded.

The workflow included overexpressed signaling gene identification, ligand-receptor interaction identification, communication probability calculation, filtering, network aggregation, and AD versus Normal signaling comparison.

This provides a systems-level assessment of predicted changes in communication between neural and glial populations.

CellChat results should be described as predicted or inferred ligand-receptor communication. They do not establish experimentally validated physical interactions between cells.

## 8. Integrated Interpretation

The analysis suggests that Alzheimer's disease in GSE243292 is associated with changes at multiple biological levels.

**Cellular level:** Relative cellular composition differs between disease groups, with notable differences involving oligodendrocytes and excitatory neurons.

**Transcriptional level:** Excitatory and inhibitory neurons showed the largest numbers of significant pseudobulk disease-associated genes. Significant changes were also detected in microglia, OPCs, and oligodendrocytes.

**Glial level:** The detection of disease-associated genes in microglia and oligodendrocyte-lineage populations indicates that observed molecular changes are not restricted to neuronal populations.

**Lineage/state level:** Monocle3 identified an inferred OPC-to-oligodendrocyte transcriptional continuum, providing a framework for studying oligodendrocyte lineage states.

**Functional level:** GO enrichment provided functional context for disease-associated genes in several populations, particularly excitatory neurons, inhibitory neurons, and microglia.

**Intercellular level:** CellChat extended the analysis to predicted intercellular signaling, allowing AD-associated changes in cell-cell communication to be investigated.

## Overall Biological Conclusion

The GSE243292 analysis indicates that Alzheimer's disease is associated with coordinated changes in cellular composition, cell-type-specific transcriptional states, oligodendrocyte-lineage progression, functional biological processes, and predicted intercellular communication.

The strongest pseudobulk transcriptional changes were observed in neuronal populations, while microglial and oligodendrocyte-lineage populations also showed disease-associated alterations. The inferred OPC-to-oligodendrocyte trajectory provides an additional perspective on oligodendrocyte lineage states, while CellChat provides a framework for investigating changes in predicted intercellular signaling.

Together, these results demonstrate how single-nucleus transcriptomics can connect cellular composition, transcriptional regulation, cellular state transitions, functional processes, and intercellular communication in Alzheimer's disease.

## Interpretation Limitations

1. Relative nucleus proportions do not directly represent absolute cell abundance.
2. Marker-based cell annotation can contain uncertainty.
3. Pseudotime represents an inferred transcriptional ordering rather than experimental lineage tracing.
4. CellChat predicts ligand-receptor communication rather than experimentally proving communication.
5. Pseudobulk analysis is more appropriate for biological replication than treating individual nuclei as independent samples.
6. Some cell populations contained too few significant genes for robust enrichment analysis.
7. **Sample size is limited and unbalanced: the Alzheimer's Disease vs Normal pseudobulk comparison rests on 8 AD samples against only 2 Normal samples.** This constrains statistical power, particularly for the Normal group, and means the analysis is better powered to detect large, consistent effects than subtle ones. A lack of significant genes in a given cell type may reflect this limited power rather than a true absence of disease-associated change.
8. **Donor-level covariates — age, sex, and postmortem interval (PMI) — were not modeled** in the pseudobulk differential expression analysis. These are established sources of variation in postmortem brain transcriptomics, and some genes identified as disease-associated here could partly reflect unmodeled donor differences rather than Alzheimer's disease status specifically.
9. **No batch correction or dataset integration was applied across the 15 samples before clustering.** Cell-type clusters and their disease-status composition could therefore be influenced by sample-of-origin technical variation in addition to genuine biological differences.
10. **No doublet detection or removal step was performed** prior to clustering and annotation; standard QC filtering on feature counts and mitochondrial percentage does not reliably exclude doublets, which could contribute to ambiguous or "Unresolved" cell populations.
11. The observed associations require biological and experimental validation before being interpreted as causal mechanisms.
