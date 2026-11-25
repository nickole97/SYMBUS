# Auditoría de árboles filogenómicos (Snodgrassella / Gilliamella) — 2026-11-25

Tabla por tip: `mag_tree_annotation.tsv` (247 filas = 239 tips en 5 árboles + 8 MAGs en ningún árbol; `in_prelim_bac120` como bandera; `label_in_plot` reproducida con la lógica de los scripts de ploteo). Copias de `.iqtree/.treefile/.contree` del clúster en `cluster/`. **[inf.]** = inferencia; el resto está verificado en archivos/logs.

## (i) Árboles existentes

Todos: IQ-TREE 2.4.0, proteínas, `-m TEST` (ModelFinder por partición; Q.plant mayoritario), supermatriz de OGs single-copy presentes en **todos** los taxa (`03_filter_core_OGs.py`). Treefiles locales = clúster (md5 verificado).

| Árbol | Archivo local = clúster | Fecha | Tips | Genes / sitios (`.iqtree`) | Partición · soporte | Raíz | Scripts · commit |
|---|---|---|---|---|---|---|---|
| Snod-46 (junio) | `New/trees/snodgrassella_phylo.treefile` = `work_meta/09_orthofinder/04_iqtree/snodgrassella/` | 2026-06-08 | 46 MAGs | 415 / 111 019 aa | `-p` · `-B 1000` (solo UFboot) | ninguna | `/dfs10/…/Ch2/scripts/09_orthofinder/` — **fuera de git**; proteínas de `08_annotation/dbcan_all/*/uniInput.faa` |
| Gilli-33 (junio) | `New/trees/gilliamella_phylo.treefile` = `…/gilliamella/` | 2026-06-08 | 33 MAGs | 733 / 219 011 aa | ídem | ninguna | ídem |
| snod_v2 | `Desktop/Ch2/snod_phylo_v2/07_iqtree/` = `work_meta/09_phylogeny/snod_phylo_v2/07_iqtree/` (job 58568003) | 2026-10-28 | 54 = 46 + 8 refs | 411 / 116 873 aa | `-q` (`-p` dio segfault, 07834c5) · `-bb 1000 -alrt 1000` | ninguna (el plot no llama `root()`) | `SYMBUS/…/snod_phylo_v2/scripts/` (ca9750b→07834c5). Treefile y `plot_snod_phylo_v2.R` solo locales |
| snod_v2 corrida previa | sobrescrita (`--redo`) | 2026-10-21 | 56 (incl. ESL0721/0727) | 216 / 61 600 aa (logs `03b_refilter`, `concat_58328407`) | `-p` | — | origen del subtítulo “216 genes · 56 taxa” y del pdf del 21-oct |
| snod_v3 | `SYMBUS/…/snod_phylo_v3/07_iqtree/` = clúster `06_iqtree/` (job 58873520) | 2026-11-04 | 58 = 46 + 8 refs + **2 Gilliamella no declaradas** + 2 outgroups | 170 / 46 896 aa | `-q` · `-bb/-alrt 1000` | `root(getMRCA(Simonsiella, Rivihabitans))` | `snod_phylo_v3/scripts/` (ced0455→073c24c); plot dfcf969. **Treefile local sin trackear** |
| gilli_v2 | `SYMBUS/…/gilli_phylo_v2/06_iqtree/` (job 58873542) | 2026-11-04 | 48 = 33 + 10 refs + 5 outgroups | 484 / 149 136 aa | `-q` · `-bb/-alrt 1000` | `Pasteurella_dagmatis` (rama 1.815 = 54 % del largo total) | `gilli_phylo_v2/scripts/`; plot a7356fa/1fe1a53; treefile en git |
| prelim (bac120) | `Desktop/Ch2/prelim.treefile` = `…/results/04_iqtree_bac120/` | 2026-01-28 | **120 árboles de genes** (115–124 tips; 124 genomas) | nucleótidos, `-S -m MFP -B 1000 --bnni` | — | no es árbol de especies; scripts no localizados |
| snod_v4 | **nunca ejecutado** (no existe `work_meta/09_phylogeny/snod_phylo_v4`) | — | — | — | — | e8d835d. Bug: el plot busca outgroups por accession, pero `01_add_outgroup_seqs.py` los nombra `Simonsiella_muelleri`, etc. |

## (ii) Verificación de problemas sospechados

1. **Outgroups v3 no monofiléticos → raíz en nodo con Gilliamella. Confirmado; consecuencia benigna.** `GCA_029433345.1`+`GCA_029433365.1` (ESL0727/ESL0721) se intercalan: `((ESL,ESL):1.516, Rivihabitans):0.304, Simonsiella):0.127`. El `getMRCA` es el clado de los 4 no-Snodgrassella; enraizando ahí, Snodgrassella (54) es monofilético (100/100) y la primera división interna es colombianos C2+C4 (10 tips, 100/100) vs resto (44, 100/100). La raíz es correcta de facto; el error es que dos Gilliamella aparecen como MAGs `unk_unk_01/02`.
2. **Conteos en subtítulos. Confirmado.** v2 “216 · 56” = corrida 20–21 oct; el treefile (28 oct) tiene 54 taxa / 411 genes. v3 “411 · 56” vs real 170 / 58. phylogenomics.Rmd (415/46, 733/33) es correcto para junio.
3. **216 / 411 / 415. Explicado:** 216 = v2 con 56 taxa (los 2 Gilliamella encogen el core); 411 = v2 final (54); 415 = Snod-46; además 170 = v3 (outgroups), 484 = gilli_v2, 733 = Gilli-33.
4. **ESL0721/0727 en v3 sin metadata. Confirmado; causa:** `snod_phylo_v3/02_prodigal/` contiene los 56 `.faa` de `snod_phylo_v2/03_prodigal/` (md5 idénticos; verificado NV006__54_sub y GCA_029433345.1), incluidos los dos Gilliamella que v2 descargó y luego excluyó (1b9a71f), más los 2 outgroups. No hay log de `01_prodigal` en v3 y su `MAG_DIR` (`06_dereplication/derep_genomes`) no existe → copia manual **[inf.]**. En gilli_v2 el mismo path falló (`01_prodigal_58865179.err`) y se corrigió (073c24c).
5. **Sin `.iqtree/.log/alineamientos` locales. Confirmado**; están en el clúster (`0X_alignments/supermatrix.faa`, `partitions.txt`, `0X_iqtree/`, `logs/`).
6. **Ramas largas. Confirmado.** v2: tallo C2+C4 0.381; ESL0693 0.331 y hermano de C2+C4 (100/100). v3: ESL0693 0.206 y **no** hermano de C2+C4; tallo C2+C4 0.176; outgroups 0.496/0.381; rama interna del par Gilliamella 1.516. NV018__13_sub 0.199/0.228/0.151 (46/v2/v3), siempre hermano del bloque colombiano C2+C4(+ESL0693). La posición de ESL0693 depende del enraizamiento/atracción de ramas largas **[inf.]**.
7. **UFboot < 70. Confirmado:** v3 6/55 (nodo 61 = backbone de 55 tips, 36.3/38; 62; 64 69.7/26; 65 65.2/27; 110; 112). v2 13/51 (incl. nodo 81 = colombianos + ESL0693 + NV018__13_sub, 31/68); Snod-46 9/43; gilli_v2 4/45; Gilli-33 1/30. C1, C2, C4, C2+C4 y S. communis (+tipo) son 100/100 en los tres árboles Snod; S. alvi_D (22) monofilético; los 5 refs de Apis no.
8. **6/52 MAGs Snod en ningún árbol, incl. todo C3. Confirmado** (iii).
9. **NV018__21_sub (inStrain) vs NV018__13_sub (árboles). Confirmado** (iv).
10. **phylogenomics.html (3 nov) = árboles de junio. Confirmado** (“415 core genes”, “46 MAGs”, “733 single-copy”; el Rmd lee `trees/*_phylo.treefile`).

## (iii) MAGs ausentes y por qué

`all_snod_mags.txt` (46) y `all_gilli_mags.txt` (33) de v2/v3/gilli_v2 son idénticas a `snod_mags.txt`/`gill_mags.txt` de junio, generadas por `00_setup.sh` como GTDB ∩ `08_annotation/kofamscan_mag_list.txt` (121 MAGs, 2026-05-12). Los 8 ausentes no están en esa lista (ni en `dbcan_all/`), aunque sí en `06_MAGs_selected/` (140) y en GTDB-Tk (139). Ningún filtro de OrthoFinder los eliminó.

| MAG | Linaje | CheckM compl/cont | En MagsToUse_SD / prelim |
|---|---|---|---|
| NV015__NV015_metabat.8 | **C3** | 83.3 / 1.1 | no / no |
| NV059__NV059_metabat.4 | **C3** | 78.8 / 1.7 | no / no |
| NV018__21_sub | C2 | 98.1 / 4.0 | **sí / sí** — no explicado por calidad ni por lista SD |
| NV034__NV034_maxbin.005 | C1 | 95.5 / 9.3 | no / no |
| 202008210927__41_sub | S. communis | 94.8 / 9.5 | no / no |
| 202008270624__39_sub | S. communis | 98.3 / 9.3 | no / no |
| 202006231018__…_metabat.2 | Gilliamella sp. | 73.7 / 0.3 | no / no |
| 202206071041G__…_metabat.2_sub | G. bombi | 95.5 / 7.3 | no / no |

## (iv) Inconsistencias de nombres

- MAGs: `<Sample>__<bin>` en todos los árboles (OK). Refs: v2/v3 por accession (`GCF_000600005.1`); gilli_v2 y outgroups de v3 por nombre (`Gris1-4`, `Simonsiella_muelleri`).
- `genome_metadata.tsv` v2/v3: NV018__13_sub = “Novel sp. E” (obsoleto: GTDB `s__Snodgrassella sp024103515`, ANI 99 %) y NV015__NV015_metabat.8 = “Novel sp. B” aunque el árbol contiene NV015__NV015_metabat.3 (C2). `instrain_strain_diversity.Rmd` y `New/snodgrassella_functional.Rmd` tratan 13_sub como conocido y 21_sub como C2; como los árboles solo tienen 13_sub, ninguna figura filogenética muestra el MAG novel de NV018.
- Etiquetas `sp_país_NN` dependen del orden de tips: el mismo MAG cambia de número entre figuras (NV015__NV015_metabat.3 = `rubi_COL_03` en Snod-46/v2, `rubi_COL_01` en v3). La tabla da el mapeo.
- Los plots leen `American_bumblebe_md_2025.tsv` y no la csv 2026 autoritativa; para las 55 muestras con MAGs, especie/subgénero/país coinciden.
- `instrain_strain_diversity.Rmd` usa `read.tree("prelim.treefile")[[1]]` como “best ML tree”: es el árbol del **primer marcador bac120**.
- TaxID de ESL0721/0727 intercambiados entre README v2 y `00_download_refs.sh` (accesiones consistentes).

## (v) Recomendación

1. **Gilliamella: gilli_v2 canónico** tal cual (raíz Pasteurella; equivalente: raíz en ESL0721+ESL0727, que es la primera división interna, y `drop.tip(Pasteurella)` para la figura). Documentar las 2 exclusiones (calidad).
2. **Snodgrassella: ninguno es presentable sin cambios.** Mínimo sin re-correr: v3 enraizado como está (la raíz es correcta) con `drop.tip(c("GCA_029433345.1","GCA_029433365.1"))` o añadiéndolas al metadata como Gilliamella; subtítulo “170 genes · 56/58 taxa”.
3. **Re-correr (recomendado, 4 jobs):** v3 sin los 2 `.faa` de Gilliamella en `03_orthofinder/faa_input` **y** con los 6 MAGs faltantes (prodigal desde `06_MAGs_selected/`, corrigiendo `MAG_DIR` en `01_prodigal.sh`): 52 + 8 + 2 = 62 taxa. Como los outgroups reducen el core (411→170), usar `--min_taxa` ≈ 90 % o el esquema v4 (OGs de ingroup + BLAST de outgroups) tras corregir los nombres de outgroup en `plot_snod_phylo_v4.R`. Mantener `-q -m TEST -bb 1000 -alrt 1000`.
4. Actualizar `genome_metadata.tsv` (13_sub → S. sp024103515; NV015 → metabat.3 C2 y metabat.8 C3), trackear `snod_phylo_v3.treefile` y `plot_snod_phylo_v2.R`, corregir subtítulos, sustituir `prelim.treefile[[1]]` en el Rmd de inStrain por el árbol canónico y re-renderizar phylogenomics.Rmd.

## (vi) Qué falta para el cophylo

- **Hospedero:** `TreeFile.nex.con_BombusPhylogeny07` (229 tips, se lee sin raíz; incluye `Apis_meli`, `Apis_dors`). Podar (`keep.tip`) a las 12 especies con MAGs: `atratus, fervidus, funebris, hortulanu, huntii, morrisoni, nevadensi, pensylv1|pensylv2, robustus, rubicundu, sonorusAz|sonorSAC, vosnesens` (decidir qué tip de *pensylvanicus* y *sonorus*), más hospederos de refs: `griseocol` (S. communis tipo, Gris1-4), `pascuorum` (S. gandavensis, ESL0250/0254), `impatiens`, `occidenta`, `terrestri`, `pratorum`, `Apis_meli` (wkB2, K-FP26). No hay tips para *A. cerana/florea/andreniformis* ni abejas sin aguijón (ESL0693/0721/0727): quedan fuera.
- **Homologar nombres:** `species_id` (metadata) → tips truncados a 9 caracteres del nexus; diccionario explícito (arriba).
- **Archivo de enlaces** (`host_tip`, `symbiont_tip`): derivable de `mag_tree_annotation.tsv` (`sample_id → host_species`), una fila por MAG; decidir si colapsar a un tip por cluster ANI/especie GTDB.
- **Árbol simbionte definitivo y enraizado** (v.3) con los 52 MAGs incl. C3; sin eso el cophylo omite dos linajes colombianos.
