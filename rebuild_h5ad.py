import scanpy as sc
import pandas as pd
import numpy as np

# =====================================================
# HUMAN
# =====================================================

human = sc.read_10x_h5(
    r"D:\PhD_work\SNCI_2026\AD_spatial\filtered_feature_bc_matrix.h5"
)

human.var_names_make_unique()

# -----------------------------------------------------
# LOAD SPATIAL COORDINATES
# -----------------------------------------------------

coords_h = pd.read_csv(
    r"D:\PhD_work\SNCI_2026\AD_spatial\spatial\tissue_positions.csv",
    header=None
)

coords_h.columns = [

    "barcode",
    "in_tissue",
    "array_row",
    "array_col",
    "pxl_row_in_fullres",
    "pxl_col_in_fullres"

]

coords_h = coords_h.set_index("barcode")

# -----------------------------------------------------
# MATCH BARCODES
# -----------------------------------------------------

common_h = human.obs_names.intersection(
    coords_h.index
)

human = human[common_h].copy()

coords_h = coords_h.loc[common_h]

# -----------------------------------------------------
# STORE SPATIAL
# -----------------------------------------------------

human.obsm["spatial"] = np.array(

    coords_h[[
        "pxl_col_in_fullres",
        "pxl_row_in_fullres"
    ]]

)

# -----------------------------------------------------
# SAVE
# -----------------------------------------------------

human.write(
    "assets/human/reference/human_adata.h5ad"
)

print("Human rebuilt!")

# =====================================================
# MOUSE
# =====================================================

mouse = sc.read_10x_h5(
    r"D:\PhD_work\SNCI_2026\HD_mouse_spatial\filtered_feature_bc_matrix.h5"
)

mouse.var_names_make_unique()

# -----------------------------------------------------
# LOAD SPATIAL COORDINATES
# -----------------------------------------------------

coords_m = pd.read_csv(
    r"D:\PhD_work\SNCI_2026\HD_mouse_spatial\spatial\tissue_positions_list.csv",
    header=None
)

coords_m.columns = [

    "barcode",
    "in_tissue",
    "array_row",
    "array_col",
    "pxl_row_in_fullres",
    "pxl_col_in_fullres"

]

coords_m = coords_m.set_index("barcode")

# -----------------------------------------------------
# MATCH BARCODES
# -----------------------------------------------------

common_m = mouse.obs_names.intersection(
    coords_m.index
)

mouse = mouse[common_m].copy()

coords_m = coords_m.loc[common_m]

# -----------------------------------------------------
# STORE SPATIAL
# -----------------------------------------------------

mouse.obsm["spatial"] = np.array(

    coords_m[[
        "pxl_col_in_fullres",
        "pxl_row_in_fullres"
    ]]

)

# -----------------------------------------------------
# SAVE
# -----------------------------------------------------

mouse.write(
    "assets/mouse/reference/mouse_adata.h5ad"
)

print("Mouse rebuilt!")