#!/usr/bin/env python3
"""
WM JHU voxel overlay generator

Reads Brain_IDPs_Parsed_for_Mapping.csv, canonicalizes tract names to JHU atlas labels,
builds per-tract voxel overlays (stat-map with intensity = max |z| per tract),
and composes an overview mosaic. Handles hemisphere unification and midline structures
replication for bilateral/NA hemi.

Outputs into Brain_MRI_Maps/WhiteMatterVoxel:
 - WM_JHU_Overlay_<tract>_<hemi>.png
 - WM_JHU_Overlay_Montage.png
 - WM_JHU_Voxel_Coverage_Summary.csv

Dependencies: nilearn, nibabel, numpy, pandas, matplotlib, pillow
"""

import os
import re
import sys
import shutil
import tempfile
import zipfile
import urllib.request
import http.cookiejar
from typing import Dict, List, Tuple, Optional
import xml.etree.ElementTree as ET

import numpy as np
import pandas as pd
try:
    from nilearn import datasets, plotting, image, surface
    import nibabel as nib
    import matplotlib.pyplot as plt
    from PIL import Image
except Exception as e:
    print("Missing dependencies. Please install: nilearn nibabel numpy pandas matplotlib pillow", file=sys.stderr)
    raise


ROOT = os.path.abspath(os.getcwd())
CSV_DEFAULT_PARSED = os.path.join(ROOT, 'Brain_MRI_Maps', 'Brain_IDPs_Parsed_for_Mapping.csv')
OUT_DIR = os.path.join(ROOT, 'Brain_MRI_Maps', 'WhiteMatterVoxel')
OUT_DIR_NETWORK = os.path.join(ROOT, 'Brain_MRI_Maps', 'Network')
OUT_DIR_CORTICAL = os.path.join(ROOT, 'Brain_MRI_Maps', 'Cortical')
OUT_DIR_SUBCORTICAL = os.path.join(ROOT, 'Brain_MRI_Maps', 'Subcortical')


def canonicalize_tract(s: str) -> str:
    v = re.sub(r"[_.]+", " ", str(s or "").strip().lower())
    v = re.sub(r"\s+", " ", v)
    # Mirrors R-side canonicalization (extended)
    if re.search(r"\bexternal capsule\b", v):
        return "External capsule"
    if re.search(r"\banterior limb of internal capsule\b|\balic\b", v):
        return "Anterior limb of internal capsule"
    if re.search(r"\bposterior limb of internal capsule\b|\blic\b|\bplic\b", v):
        return "Posterior limb of internal capsule"
    if re.search(r"\bretrolenticular part of internal capsule\b|\brlic\b", v):
        return "Retrolenticular part of internal capsule"
    if re.search(r"anterior thalamic|\batr\b", v):
        return "Anterior thalamic radiation"
    if re.search(r"posterior thalamic|\bptr\b", v):
        return "Posterior thalamic radiation"
    if re.search(r"\bsagittal stratum\b", v):
        return "Sagittal stratum"
    if "corticospinal" in v:
        return "Corticospinal tract"
    if re.search(r"forceps minor|genu", v):
        return "Genu of corpus callosum"
    if re.search(r"forceps major|splenium", v):
        return "Splenium of corpus callosum"
    if re.search(r"\bcorpus callosum\b", v):
        return "Body of corpus callosum"
    if re.search(r"uncinate|\buf\b", v):
        return "Uncinate fasciculus"
    if re.search(r"inferior fronto[- ]?occipital|\bifo\b", v):
        return "Sagittal stratum"
    if re.search(r"superior longitudinal|\bslf\b", v):
        return "Superior longitudinal fasciculus"
    if re.search(r"inferior longitudinal|\bilf\b", v):
        return "Sagittal stratum"
    if "cingulum" in v and "hippocamp" in v:
        return "Cingulum (hippocampus)"
    if re.search(r"cingulum.*cingulate|\bcg\b", v):
        return "Cingulum (cingulate gyrus)"
    if re.search(r"middle longitudinal|\bmdlf\b", v):
        return "Sagittal stratum"
    if re.search(r"superior fronto[- ]?occipital|\bsfof\b", v):
        return "Superior fronto-occipital fasciculus"
    if re.search(r"middle cerebellar peduncle|\bmcp\b", v):
        return "Middle cerebellar peduncle"
    if re.search(r"superior cerebellar peduncle|\bscp\b", v):
        return "Superior cerebellar peduncle"
    if re.search(r"inferior cerebellar peduncle|\bicp\b", v):
        return "Inferior cerebellar peduncle"
    if re.search(r"fornix.*(cres|crus)|fornix\s+cres", v):
        return "Fornix (cres)"
    if re.search(r"fornix.*body|body\s+of\s+fornix", v):
        return "Fornix (column and body of fornix)"
    if re.search(r"medial\s+lemniscus|\bml\b", v):
        return "Medial lemniscus"
    if re.search(r"anterior corona radiata|\bacr\b", v):
        return "Anterior corona radiata"
    if re.search(r"superior corona radiata|\bscr\b", v):
        return "Superior corona radiata"
    if re.search(r"posterior corona radiata|\bpcr\b", v):
        return "Posterior corona radiata"
    if re.search(r"\bcerebral peduncle\b", v):
        return "Cerebral peduncle"
    if re.search(r"\btapetum\b", v):
        return "Tapetum"
    return s.strip()

def strip_hemi_tokens(name: str) -> str:
    """Remove hemisphere tokens (left/right/lh/rh and standalone L/R/(L)/(R))
    from a region/tract name without touching embedded acronyms like SLF.
    """
    if not name:
        return ''
    n = str(name)
    # Remove common prefixes like 'Left ' / 'Right '
    n = re.sub(r'^(\s*)(left|right)\s+', '', n, flags=re.IGNORECASE)
    # Remove common suffixes ' left' / ' right'
    n = re.sub(r'\s+(left|right)\s*$', '', n, flags=re.IGNORECASE)
    # Remove hemisphere descriptors
    n = re.sub(r'\b(lh|rh|left\s*hemisphere|right\s*hemisphere)\b', '', n, flags=re.IGNORECASE)
    # Remove trailing standalone (L)/(R) or L/R tokens
    n = re.sub(r'\s*[(_-]?\b(l|r)\b[)_-]?\s*$', '', n, flags=re.IGNORECASE)
    # Collapse multiple spaces and punctuation leftovers
    n = re.sub(r'[_\-]+', ' ', n)
    n = re.sub(r'\s{2,}', ' ', n)
    return n.strip()


MIDLINE_TRACTS = {
    "Middle cerebellar peduncle",
    "Forceps minor (anterior)",
    "Forceps major (posterior)",
}


def unify_hemi(hemi_raw: str) -> str:
    if hemi_raw is None:
        return "unknown"
    h = str(hemi_raw).strip().lower()
    if h in {"l", "lh", "left"}:
        return "left"
    if h in {"r", "rh", "right"}:
        return "right"
    if h in {"bilateral", "both", "na", "", "midline", "middle"}:
        return "bilateral"
    return h


def _coerce_bool_series(s: pd.Series) -> pd.Series:
    if s is None:
        return pd.Series([], dtype=bool)
    if pd.api.types.is_bool_dtype(s):
        return s.fillna(False)
    v = s.fillna('').astype(str).str.strip().str.lower()
    return v.isin(['true', 't', '1', 'yes', 'y'])


def load_sig_wm_tract_effects(sig_csv: str) -> Dict[str, Dict[str, float]]:
    if not sig_csv or not os.path.exists(sig_csv):
        return {}
    df = pd.read_csv(sig_csv)
    if 'model_comparison' not in df.columns:
        return {}

    z_col = 'z_change' if 'z_change' in df.columns else ('z_statistic_without_non_imaging' if 'z_statistic_without_non_imaging' in df.columns else None)
    if z_col is None:
        return {}

    is_wm = pd.Series([False] * len(df))
    if 'structure_type' in df.columns:
        is_wm = df['structure_type'].fillna('').astype(str).str.lower().eq('white matter')
    if 'wm_tract' not in df.columns:
        return {}

    sig_mask = pd.Series([True] * len(df))
    if 'selected_for_mapping' in df.columns:
        sig_mask = _coerce_bool_series(df['selected_for_mapping'])
    else:
        ie = _coerce_bool_series(df['independent_effect']) if 'independent_effect' in df.columns else pd.Series([True] * len(df))
        sfdr = _coerce_bool_series(df['significant_fdr']) if 'significant_fdr' in df.columns else pd.Series([False] * len(df))
        sfwe = _coerce_bool_series(df['significant_fwe']) if 'significant_fwe' in df.columns else pd.Series([False] * len(df))
        pfdr = pd.to_numeric(df['p_value_fdr'], errors='coerce') if 'p_value_fdr' in df.columns else pd.Series([np.nan] * len(df))
        pfwe = pd.to_numeric(df['p_value_fwe'], errors='coerce') if 'p_value_fwe' in df.columns else pd.Series([np.nan] * len(df))
        sig_mask = ie & (sfdr | sfwe | (pfdr < 0.05) | (pfwe < 0.05))

    df_wm = df[is_wm & sig_mask].copy()
    if df_wm.empty:
        return {}

    df_wm['wm_tract'] = df_wm['wm_tract'].fillna('').astype(str).str.strip().str.strip(',').str.strip()
    df_wm = df_wm[df_wm['wm_tract'].str.len() > 0]
    if df_wm.empty:
        return {}

    df_wm['_tract_key'] = df_wm['wm_tract'].apply(lambda x: canonicalize_tract(strip_hemi_tokens(x)))
    df_wm['_z'] = pd.to_numeric(df_wm[z_col], errors='coerce')
    df_wm = df_wm[~df_wm['_z'].isna()].copy()
    if df_wm.empty:
        return {}

    out: Dict[str, Dict[str, float]] = {}
    for g, dfg in df_wm.groupby('model_comparison'):
        tract_to_z: Dict[str, float] = {}
        for tract, dft in dfg.groupby('_tract_key'):
            idx = (dft['_z'].abs()).idxmax()
            if pd.isna(idx):
                continue
            tract_to_z[str(tract)] = float(dft.loc[idx, '_z'])
        if tract_to_z:
            out[str(g)] = tract_to_z
    return out


def filter_selected_idps(df: pd.DataFrame) -> pd.DataFrame:
    if df is None or df.empty:
        return df if df is not None else pd.DataFrame()
    if 'selected_for_mapping' in df.columns:
        mask = _coerce_bool_series(df['selected_for_mapping'])
        return df[mask].copy()
    ie = _coerce_bool_series(df['independent_effect']) if 'independent_effect' in df.columns else pd.Series([True] * len(df))
    sfdr = _coerce_bool_series(df['significant_fdr']) if 'significant_fdr' in df.columns else pd.Series([False] * len(df))
    sfwe = _coerce_bool_series(df['significant_fwe']) if 'significant_fwe' in df.columns else pd.Series([False] * len(df))
    p_mask = pd.Series([False] * len(df))
    for pcol in [
        'p_value_fdr',
        'p_value_fwe',
        'p_value_without_non_imaging',
        'p_value_unadjusted',
        'p_value',
    ]:
        if pcol in df.columns:
            p = pd.to_numeric(df[pcol], errors='coerce')
            p_mask = p_mask | (p < 0.05)
    if not (('independent_effect' in df.columns) or ('significant_fdr' in df.columns) or ('significant_fwe' in df.columns) or any(c in df.columns for c in ['p_value_fdr', 'p_value_fwe', 'p_value_without_non_imaging', 'p_value_unadjusted', 'p_value'])):
        raise ValueError("selected-only enabled but no significance/selection columns found in input CSV.")
    mask = ie & (sfdr | sfwe | p_mask)
    return df[mask].copy()


def _infer_z_col(df: pd.DataFrame) -> Optional[str]:
    if df is None or df.empty:
        return None
    if 'z_change' in df.columns:
        return 'z_change'
    if 'z_statistic_without_non_imaging' in df.columns:
        return 'z_statistic_without_non_imaging'
    if 'z_statistic' in df.columns:
        return 'z_statistic'
    return None


def write_uofm_glass_idp_summary(sig_csv: str, group: str, out_dir: str) -> Optional[str]:
    if not sig_csv or not os.path.exists(sig_csv):
        return None
    df = pd.read_csv(sig_csv)
    if df.empty or 'model_comparison' not in df.columns:
        return None
    df = df[df['model_comparison'].astype(str) == str(group)].copy()
    if df.empty:
        return None
    if 'structure_type' in df.columns:
        is_wm = df['structure_type'].fillna('').astype(str).str.lower().eq('white matter')
        df = df[is_wm].copy()
    if df.empty or 'wm_tract' not in df.columns:
        return None
    sig_mask = pd.Series([True] * len(df))
    if 'selected_for_mapping' in df.columns:
        sig_mask = _coerce_bool_series(df['selected_for_mapping'])
    else:
        has_any_sig_cols = (
            ('independent_effect' in df.columns) or ('significant_fdr' in df.columns) or ('significant_fwe' in df.columns)
            or ('p_value_fdr' in df.columns) or ('p_value_fwe' in df.columns)
        )
        if not has_any_sig_cols:
            return None
        ie = _coerce_bool_series(df['independent_effect']) if 'independent_effect' in df.columns else pd.Series([True] * len(df))
        sfdr = _coerce_bool_series(df['significant_fdr']) if 'significant_fdr' in df.columns else pd.Series([False] * len(df))
        sfwe = _coerce_bool_series(df['significant_fwe']) if 'significant_fwe' in df.columns else pd.Series([False] * len(df))
        pfdr = pd.to_numeric(df['p_value_fdr'], errors='coerce') if 'p_value_fdr' in df.columns else pd.Series([np.nan] * len(df))
        pfwe = pd.to_numeric(df['p_value_fwe'], errors='coerce') if 'p_value_fwe' in df.columns else pd.Series([np.nan] * len(df))
        sig_mask = ie & (sfdr | sfwe | (pfdr < 0.05) | (pfwe < 0.05))
    df = df[sig_mask].copy()
    if df.empty:
        return None
    z_col = ('z_change' if 'z_change' in df.columns else ('z_statistic_without_non_imaging' if 'z_statistic_without_non_imaging' in df.columns else ('z_statistic' if 'z_statistic' in df.columns else None)))
    if z_col is None:
        return None
    df['wm_tract'] = df['wm_tract'].fillna('').astype(str).str.strip().str.strip(',').str.strip()
    df = df[df['wm_tract'].str.len() > 0].copy()
    if df.empty:
        return None
    df['_tract_key'] = df['wm_tract'].apply(lambda x: canonicalize_tract(strip_hemi_tokens(x)))
    df['_z'] = pd.to_numeric(df[z_col], errors='coerce')
    df = df[~df['_z'].isna()].copy()
    if df.empty:
        return None
    df['_absz'] = df['_z'].abs()
    picked_idx = df.groupby('_tract_key')['_absz'].idxmax()
    picked = set(int(i) for i in picked_idx.dropna().tolist())
    df['picked_for_uofm_weight'] = [int(i in picked) for i in df.index]
    cols = [
        'model_comparison',
        'variable_name',
        'display_name',
        'hemi',
        'wm_tract',
        '_tract_key',
        z_col,
        '_z',
        '_absz',
        'picked_for_uofm_weight',
    ]
    out = df[[c for c in cols if c in df.columns]].copy()
    out = out.rename(columns={'_tract_key': 'tract_key', '_absz': 'abs_z'})
    out = out.sort_values(by=['picked_for_uofm_weight', 'abs_z'], ascending=[False, False])
    os.makedirs(out_dir, exist_ok=True)
    out_path = os.path.join(out_dir, 'UofM_Glass_IDPs.csv')
    out.to_csv(out_path, index=False)
    return out_path


def _read_jhu_labels_txt(labels_txt: str) -> Tuple[List[str], Dict[str, int]]:
    """Read a JHU labels text file and return (labels, name_to_code).

    The text file usually contains lines like: "1  Posterior thalamic radiation L".
    """
    labels: List[str] = []
    name_to_code: Dict[str, int] = {}
    with open(labels_txt, 'r') as f:
        for ln in f:
            ln = ln.strip()
            if not ln:
                continue
            m = re.match(r"^(\d+)\s+(.+)$", ln)
            if not m:
                continue
            code = int(m.group(1))
            name = m.group(2).strip()
            labels.append(name)
            name_to_code[name] = code
    return labels, name_to_code


def _read_jhu_labels_xml(labels_xml: str) -> Tuple[List[str], Dict[str, int]]:
    """Read JHU labels from an XML file common in FSL installations."""
    labels: List[str] = []
    name_to_code: Dict[str, int] = {}
    try:
        tree = ET.parse(labels_xml)
        root = tree.getroot()
        for elem in root.iter():
            tag = elem.tag.lower()
            if tag in {"label", "region"}:
                idx = elem.get("index") or elem.get("val") or elem.get("id")
                name = (elem.text or "").strip()
                if idx and name:
                    code = int(idx)
                    labels.append(name)
                    name_to_code[name] = code
    except Exception:
        pass
    return labels, name_to_code


def discover_local_jhu(atlas_dir: Optional[str] = None,
                       atlas_nii: Optional[str] = None,
                       labels_txt: Optional[str] = None) -> Optional[Tuple[nib.Nifti1Image, List[str], Dict[str, int], str]]:
    """Attempt to discover JHU atlas files locally.

    Search order:
    - Explicit paths via args: atlas_nii, labels_txt
    - A directory via args/env: atlas_dir or $JHU_ATLAS_DIR
    - FSL installation via $FSLDIR or /usr/local/fsl
    - Project-local drop-in: ./atlases/JHU

    Returns (img, labels, name_to_code, source) or None if not found.
    """
    # If direct files provided
    if atlas_nii and labels_txt and os.path.exists(atlas_nii) and os.path.exists(labels_txt):
        img = nib.load(atlas_nii)
        labels, name_to_code = _read_jhu_labels_txt(labels_txt)
        return img, labels, name_to_code, f"files:{atlas_nii}"

    # Directory-based search
    candidate_dirs: List[str] = []
    if atlas_dir:
        candidate_dirs.append(atlas_dir)
    env_dir = os.environ.get('JHU_ATLAS_DIR')
    if env_dir:
        candidate_dirs.append(env_dir)

    fsldir = os.environ.get('FSLDIR') or '/usr/local/fsl'
    candidate_dirs.append(os.path.join(fsldir, 'data', 'atlases', 'JHU'))
    # Project-local drop-in
    script_dir = os.path.dirname(os.path.abspath(__file__))
    candidate_dirs.append(os.path.join(script_dir, 'atlases', 'JHU'))
    candidate_dirs.append(os.path.join(ROOT, 'atlases', 'JHU'))
    candidate_dirs.append(os.path.join(script_dir, 'ischemic', 'atlases', 'JHU'))
    candidate_dirs.append(os.path.join(ROOT, 'ischemic', 'atlases', 'JHU'))

    NII_CANDIDATES = [
        'JHU-ICBM-tracts-maxprob-1mm.nii.gz',
        'JHU-ICBM-tracts-maxprob-thr0-1mm.nii.gz',
        'JHU-ICBM-tracts-maxprob-2mm.nii.gz',
        'JHU-ICBM-tracts-maxprob-thr0-2mm.nii.gz',
    ]

    for d in candidate_dirs:
        if not os.path.isdir(d):
            continue
        nii_path = None
        # Prefer known filenames
        for nm in NII_CANDIDATES:
            p = os.path.join(d, nm)
            if os.path.exists(p):
                nii_path = p
                break
        # Fallback: find any NIfTI with JHU-ICBM-tracts
        if nii_path is None:
            try:
                for fn in os.listdir(d):
                    if re.search(r'JHU-ICBM-tracts.*\.nii(\.gz)?$', fn, re.I):
                        nii_path = os.path.join(d, fn)
                        break
            except Exception:
                pass
        txt_path = os.path.join(d, 'JHU-ICBM-tracts-labels.txt')
        xml_path = os.path.join(d, 'JHU-ICBM-tracts-labels.xml')
        if nii_path and (os.path.exists(txt_path) or os.path.exists(xml_path)):
            img = nib.load(nii_path)
            labels: List[str] = []
            name_to_code: Dict[str, int] = {}
            if os.path.exists(txt_path):
                labels, name_to_code = _read_jhu_labels_txt(txt_path)
            elif os.path.exists(xml_path):
                labels, name_to_code = _read_jhu_labels_xml(xml_path)
            if labels:
                return img, labels, name_to_code, f"dir:{d}"
    return None


def fetch_jhu(atlas_dir: Optional[str] = None,
              atlas_nii: Optional[str] = None,
              labels_txt: Optional[str] = None) -> Tuple[nib.Nifti1Image, List[str], Optional[Dict[str, int]], str]:
    """Fetch JHU ICBM-DTI-81 white matter labels via nilearn or local files.

    Returns (label_img, labels, name_to_code_or_None, source).
    Raises RuntimeError if not available.
    """
    # Try local files first (fastest, explicit)
    found = discover_local_jhu(atlas_dir=atlas_dir, atlas_nii=atlas_nii, labels_txt=labels_txt)
    if found:
        img, labels, name_to_code, src = found
        return img, labels, name_to_code, src

    # Try nilearn (older versions had fetch_atlas_jhu_icbm_tracts)
    if hasattr(datasets, 'fetch_atlas_jhu_icbm_tracts'):
        atlas = datasets.fetch_atlas_jhu_icbm_tracts()
        img = nib.load(atlas['maps'])
        labels = list(atlas['labels'])
        # No codes provided; assume 1..N mapping, handled downstream.
        return img, labels, None, 'nilearn'

    # Fallback: FSL installation
    found = discover_local_jhu()
    if found:
        img, labels, name_to_code, src = found
        return img, labels, name_to_code, src

    raise RuntimeError(
        "JHU atlas not found. Supply --atlas-dir (or --atlas-nii/--labels-txt),\n"
        "or set JHU_ATLAS_DIR, or install FSL and set FSLDIR.\n"
        "Expected files: JHU-ICBM-tracts-maxprob-1mm.nii.gz and JHU-ICBM-tracts-labels.txt"
    )


def label_index(labels: List[str]) -> Dict[str, int]:
    # labels are like 'Anterior thalamic radiation L' or 'Right ...'
    # Normalize to 'Left/Right <name>' via canonical mapping
    name_to_idx = {}
    for idx, lab in enumerate(labels):
        # In JHU ICBM, index starts at 0 or 1 depending; nilearn labels aligns to map values
        name_to_idx[lab] = idx
    return name_to_idx


def find_label(labels: List[str], tract: str, hemi: str) -> str:
    # Try multiple name forms
    base = canonicalize_tract(strip_hemi_tokens(tract))
    candidates = []
    if hemi == 'left':
        candidates += [f"Left {base}", f"{base} L", f"{base} left", base]
    elif hemi == 'right':
        candidates += [f"Right {base}", f"{base} R", f"{base} right", base]
    else:  # midline/bilateral
        candidates += [base]
    # exact match first
    for c in candidates:
        for lab in labels:
            if lab.lower() == c.lower():
                return lab
    # loose normalized contains match
    def norm(x: str) -> str:
        x = re.sub(r'[_\-]', ' ', x.lower())
        return re.sub(r'\s{2,}', ' ', x).strip()
    for lab in labels:
        if norm(base) in norm(lab):
            if hemi in ('left', 'right'):
                if hemi in lab.lower() or (hemi == 'left' and ' l' in lab.lower()) or (hemi == 'right' and ' r' in lab.lower()):
                    return lab
            else:
                return lab
    # Fallback: token overlap (Jaccard) to find closest label when names differ
    def tokens(x: str) -> set:
        x = norm(x)
        drop = {'left', 'right', 'tract', 'fasciculus', 'radiation', 'peduncle', 'gyrus', 'cortex',
                'anterior', 'posterior', 'superior', 'inferior', 'medial', 'lateral'}
        toks = re.split(r"\W+", x)
        return {t for t in toks if t and t not in drop}
    base_toks = tokens(base)
    best_lab = ''
    best_score = 0.0
    for lab in labels:
        # respect hemisphere when possible
        if hemi in ('left', 'right') and hemi not in lab.lower():
            continue
        lab_no_hemi = re.sub(r'^(left|right)\s+', '', lab.lower())
        lab_toks = tokens(lab_no_hemi)
        if not base_toks or not lab_toks:
            continue
        inter = len(base_toks & lab_toks)
        union = len(base_toks | lab_toks)
        score = inter / union if union else 0.0
        if score > best_score:
            best_score = score
            best_lab = lab
    if best_score >= 0.5:
        return best_lab
    return ''


def build_stat_map(base_img: nib.Nifti1Image,
                   labels: List[str],
                   label_name: str,
                   intensity: float,
                   name_to_code: Optional[Dict[str, int]] = None) -> Optional[nib.Nifti1Image]:
    if not label_name:
        return None
    data = base_img.get_fdata()
    # Determine label code value in the atlas image
    if name_to_code and label_name in name_to_code:
        lbl_idx = int(name_to_code[label_name])
    else:
        # Assume 1..N mapping aligned with labels order
        lbl_idx = labels.index(label_name) + 1
    mask = (data == lbl_idx).astype(np.float32)
    stat = mask * float(intensity)
    return nib.Nifti1Image(stat, base_img.affine)


def save_overlay(stat_img: nib.Nifti1Image,
                 out_path: str,
                 title: str,
                 bg_img: Optional[nib.Nifti1Image] = None,
                 threshold: Optional[float] = None,
                 cmap: str = 'RdBu_r',
                 alpha: float = 0.9,
                 black_bg: bool = True,
                 symmetric_cbar: bool = True,
                 dpi: int = 300):
    if stat_img is None:
        return False
    disp = plotting.plot_stat_map(stat_img,
                                  bg_img=bg_img,
                                  display_mode='ortho',
                                  threshold=threshold,
                                  cmap=cmap,
                                  alpha=alpha,
                                  colorbar=True,
                                  title=title,
                                  black_bg=black_bg,
                                  symmetric_cbar=symmetric_cbar)
    disp.savefig(out_path, dpi=dpi)
    try:
        pdf_path = os.path.splitext(out_path)[0] + '.pdf'
        disp.savefig(pdf_path, dpi=dpi)
    except Exception:
        pass
    disp.close()
    return True


def save_glass_brain(stat_img: nib.Nifti1Image,
                     out_path: str,
                     title: str,
                     threshold: Optional[float] = None,
                     cmap: str = 'RdBu_r',
                     plot_abs: bool = False,
                     black_bg: bool = True,
                     dpi: int = 300) -> bool:
    if stat_img is None:
        return False
    try:
        disp = plotting.plot_glass_brain(stat_img,
                                         display_mode='lyrz',
                                         threshold=threshold,
                                         cmap=cmap,
                                         plot_abs=plot_abs,
                                         black_bg=black_bg,
                                         colorbar=True,
                                         title=title)
        disp.savefig(out_path, dpi=dpi)
        try:
            pdf_path = os.path.splitext(out_path)[0] + '.pdf'
            disp.savefig(pdf_path, dpi=dpi)
        except Exception:
            pass
        disp.close()
        return True
    except Exception as e:
        print(f"Glass brain render failed: {e}")
        return False


def compose_montage(image_paths: List[str], out_path: str, cols: int = 4):
    imgs = [Image.open(p) for p in image_paths if os.path.exists(p)]
    if not imgs:
        return
    w, h = imgs[0].size
    rows = int(np.ceil(len(imgs) / cols))
    canvas = Image.new('RGB', (cols * w, rows * h), color=(255, 255, 255))
    for i, im in enumerate(imgs):
        r = i // cols
        c = i % cols
        canvas.paste(im, (c * w, r * h))
    canvas.save(out_path)
    try:
        pdf_path = os.path.splitext(out_path)[0] + '.pdf'
        canvas.save(pdf_path, format='PDF')
    except Exception:
        pass

def ensure_pdf_companions(dir_path: str) -> None:
    """Ensure every .png under dir_path has a .pdf companion saved via Pillow.
    Useful when upstream rendering produced PNGs earlier or third-party functions
    do not emit PDF directly.
    """
    if not os.path.isdir(dir_path):
        return
    for root, _, files in os.walk(dir_path):
        for f in files:
            if f.lower().endswith('.png'):
                png_fp = os.path.join(root, f)
                pdf_fp = os.path.splitext(png_fp)[0] + '.pdf'
                if os.path.exists(pdf_fp):
                    continue
                try:
                    im = Image.open(png_fp)
                    im.save(pdf_fp, format='PDF')
                except Exception:
                    # best effort; skip on failure
                    pass


# ---------- UofM probability map integration ----------

NITRC_UOFM_V1_LINK_ID = 8036
NITRC_UOFM_V2_LINK_ID = 10102

UOFM_TRACT_KEYWORDS: Dict[str, List[str]] = {
    "Superior longitudinal fasciculus": ["arcuate", "slf"],
}

def _build_cookie_opener(cookie_jar_path: str) -> urllib.request.OpenerDirector:
    cj = http.cookiejar.MozillaCookieJar()
    cj.load(cookie_jar_path, ignore_discard=True, ignore_expires=True)
    return urllib.request.build_opener(urllib.request.HTTPCookieProcessor(cj))


def _nitrc_downloadlink_to_zip(downloadlink_id: int, cookie_jar_path: str, out_zip_path: str) -> None:
    url = f"https://www.nitrc.org/frs/downloadlink.php/{int(downloadlink_id)}"
    opener = _build_cookie_opener(cookie_jar_path)
    req = urllib.request.Request(url, headers={'User-Agent': 'Mozilla/5.0'})
    with opener.open(req, timeout=120) as resp:
        final_url = resp.geturl()
        if 'account/login.php' in final_url:
            raise RuntimeError(
                f"NITRC download requires login; cookie jar did not authenticate. final_url={final_url}"
            )
        ctype = (resp.headers.get('Content-Type') or '').lower()
        cd = (resp.headers.get('Content-Disposition') or '').lower()
        is_zip = ('zip' in ctype) or ('.zip' in final_url.lower()) or ('.zip' in cd)
        if is_zip:
            os.makedirs(os.path.dirname(out_zip_path), exist_ok=True)
            with open(out_zip_path, 'wb') as f:
                shutil.copyfileobj(resp, f)
            return
        data = resp.read(2 * 1024 * 1024)
    html = data.decode('utf-8', errors='ignore')
    candidates = re.findall(r'href=[\"\\\'](https?://[^\"\\\']+)[\"\\\']', html, flags=re.I)
    candidates = [c for c in candidates if 'nitrc.org' not in c.lower()]
    if not candidates:
        raw_urls = re.findall(r'https?://[^\\s<>()\\[\\]\"\\\']+', html, flags=re.I)
        candidates = [u for u in raw_urls if 'nitrc.org' not in u.lower()]
    if not candidates:
        raise RuntimeError("Could not resolve external URL from NITRC downloadlink page.")
    _download_to_file(candidates[0], out_zip_path)


def _download_to_file(url: str, out_path: str) -> None:
    os.makedirs(os.path.dirname(out_path), exist_ok=True)
    req = urllib.request.Request(url, headers={'User-Agent': 'Mozilla/5.0'})
    with urllib.request.urlopen(req, timeout=120) as resp, open(out_path, 'wb') as f:
        shutil.copyfileobj(resp, f)


def install_uofm_atlas_from_nitrc(downloadlink_id: int,
                                  cookie_jar_path: str,
                                  out_dir: str) -> str:
    if not cookie_jar_path or not os.path.exists(cookie_jar_path):
        raise FileNotFoundError(f"Cookie jar not found: {cookie_jar_path}")
    os.makedirs(out_dir, exist_ok=True)
    with tempfile.TemporaryDirectory() as td:
        zip_path = os.path.join(td, f"nitrc_{downloadlink_id}.zip")
        _nitrc_downloadlink_to_zip(downloadlink_id, cookie_jar_path, zip_path)
        with zipfile.ZipFile(zip_path) as z:
            z.extractall(out_dir)
    return out_dir


def _norm_text(s: str) -> str:
    v = re.sub(r'[^a-z0-9]+', ' ', str(s or '').lower())
    return re.sub(r'\s+', ' ', v).strip()


def _match_uofm_map_to_tract_effect(fp: str,
                                   title: str,
                                   effects_by_tract: Dict[str, float]) -> Tuple[Optional[str], Optional[float]]:
    if not effects_by_tract:
        return None, None
    hay = _norm_text(f"{title} {os.path.basename(fp)} {fp}")
    best_tract = None
    best_score = 0
    best_absz = -1.0
    for tract, z in effects_by_tract.items():
        kws = UOFM_TRACT_KEYWORDS.get(tract)
        if not kws:
            toks = [t for t in _norm_text(tract).split(' ') if len(t) >= 4]
            kws = toks[:4] if toks else []
        score = 0
        for kw in kws:
            kw_norm = _norm_text(kw)
            if not kw_norm:
                continue
            if f" {kw_norm} " in f" {hay} ":
                score += 1
        if score <= 0:
            continue
        az = abs(float(z))
        if score > best_score or (score == best_score and az > best_absz):
            best_score = score
            best_absz = az
            best_tract = tract
    if best_tract is None:
        return None, None
    return best_tract, float(effects_by_tract[best_tract])


def discover_uofm_prob_maps(uofm_dir: str) -> List[Tuple[str, str]]:
    """Discover UofM network probability maps (.nii) and derive human-readable titles.
    Returns list of (path, title).
    """
    prob_maps = []
    target_root = os.path.join(uofm_dir, 'Network_Fibers(Group_Probability_Maps)')
    if not os.path.isdir(target_root):
        target_root = uofm_dir
    for root, _, files in os.walk(target_root):
        if os.path.basename(root).lower() == 'nifti_images':
            parts = root.split(os.sep)
            network = parts[-2] if len(parts) >= 2 else ''
            for f in files:
                if f.lower().endswith('.nii') or f.lower().endswith('.nii.gz'):
                    fp = os.path.join(root, f)
                    title = f"UofM {network} - {os.path.splitext(f)[0]}"
                    prob_maps.append((fp, title))
    return prob_maps


def load_bg_img(bg_nii: Optional[str], uofm_dir: Optional[str]) -> Optional[nib.Nifti1Image]:
    """Load background T1 image. Prefer provided path, then UofM template, then MNI152."""
    if bg_nii and os.path.exists(bg_nii):
        return nib.load(bg_nii)
    if uofm_dir:
        tmpl = os.path.join(uofm_dir, 'Template_Images', 'spm8_single_subj_T1_1mm_cubic.nii')
        if os.path.exists(tmpl):
            return nib.load(tmpl)
    try:
        return datasets.load_mni152_template()
    except Exception:
        return None


def _scale_img(img: nib.Nifti1Image, w: float) -> nib.Nifti1Image:
    if w == 1.0:
        return img
    data = np.asarray(img.get_fdata(dtype=np.float32)) * float(w)
    return nib.Nifti1Image(data, img.affine, img.header)


def generate_uofm_overlays(uofm_dir: str,
                           out_dir: str,
                           bg_img: Optional[nib.Nifti1Image],
                           threshold: float = 0.1,
                           make_interactive: bool = True,
                           effects_by_tract: Optional[Dict[str, float]] = None,
                           sig_only: bool = False) -> Tuple[List[str], Optional[str]]:
    """Generate overlays for UofM probability maps. Returns (overlay_paths, montage_path)."""
    os.makedirs(out_dir, exist_ok=True)
    maps = discover_uofm_prob_maps(uofm_dir)
    if not maps:
        print(f"No UofM probability maps found under: {uofm_dir}")
        return [], None
    overlay_paths = []
    for fp, title in maps:
        try:
            pm = nib.load(fp)
            pm_img = pm
            if bg_img is not None:
                try:
                    pm_img = image.resample_to_img(pm, bg_img, interpolation='continuous')
                except Exception:
                    pm_img = pm
            tract_hit, z_hit = _match_uofm_map_to_tract_effect(fp, title, effects_by_tract or {})
            if sig_only and tract_hit is None:
                continue
            pm_to_plot = pm_img
            plot_title = title
            if tract_hit is not None and z_hit is not None:
                w = float(abs(z_hit))
                pm_to_plot = _scale_img(pm_img, w)
                plot_title = f"{title} | {tract_hit} z={z_hit:+.3f}"
            safe = re.sub(r"\W+", "_", os.path.splitext(os.path.basename(fp))[0])
            out_png = os.path.join(out_dir, f"{safe}.png")
            disp = plotting.plot_stat_map(pm_to_plot,
                                          bg_img=bg_img,
                                          display_mode='ortho',
                                          threshold=threshold,
                                          cmap='hot',
                                          colorbar=True,
                                          title=plot_title)
            disp.savefig(out_png, dpi=150)
            try:
                pdf_path = os.path.splitext(out_png)[0] + '.pdf'
                disp.savefig(pdf_path, dpi=150)
            except Exception:
                pass
            disp.close()
            overlay_paths.append(out_png)
            if make_interactive:
                try:
                    html_view = plotting.view_img(pm_to_plot, bg_img=bg_img, threshold=threshold, cmap='hot')
                    out_html = os.path.join(out_dir, f"{safe}.html")
                    html_view.save_as_html(out_html)
                except Exception:
                    pass
        except Exception as e:
            print(f"Failed to render {fp}: {e}")
            continue
    montage_path = None
    if overlay_paths:
        montage_path = os.path.join(out_dir, 'UofM_Network_Montage.png')
        compose_montage(overlay_paths, montage_path, cols=4)
    return overlay_paths, montage_path

def generate_uofm_glass_overlays(uofm_dir: str,
                                 out_dir: str,
                                 threshold: float = 0.1,
                                 make_interactive: bool = True,
                                 effects_by_tract: Optional[Dict[str, float]] = None,
                                 sig_only: bool = False) -> Tuple[List[str], Optional[str]]:
    """Generate glass-brain overlays for UofM probability maps.
    Returns (overlay_paths, montage_path).
    """
    os.makedirs(out_dir, exist_ok=True)
    maps = discover_uofm_prob_maps(uofm_dir)
    if not maps:
        print(f"No UofM probability maps found under: {uofm_dir}")
        return [], None
    overlay_paths = []
    for fp, title in maps:
        try:
            pm = nib.load(fp)
            pm_img = pm
            tract_hit, z_hit = _match_uofm_map_to_tract_effect(fp, title, effects_by_tract or {})
            if sig_only and tract_hit is None:
                continue
            pm_to_plot = pm_img
            plot_title = title
            if tract_hit is not None and z_hit is not None:
                w = float(abs(z_hit))
                pm_to_plot = _scale_img(pm_img, w)
                plot_title = f"{title} | {tract_hit} z={z_hit:+.3f}"
            safe = re.sub(r"\W+", "_", os.path.splitext(os.path.basename(fp))[0])
            out_png = os.path.join(out_dir, f"{safe}.png")
            disp = plotting.plot_glass_brain(pm_to_plot,
                                             threshold=threshold,
                                             cmap='hot',
                                             colorbar=True,
                                             title=plot_title)
            disp.savefig(out_png, dpi=150)
            try:
                pdf_path = os.path.splitext(out_png)[0] + '.pdf'
                disp.savefig(pdf_path, dpi=150)
            except Exception:
                pass
            disp.close()
            overlay_paths.append(out_png)
            if make_interactive:
                try:
                    html_view = plotting.view_img(pm_to_plot, threshold=threshold, cmap='hot')
                    out_html = os.path.join(out_dir, f"{safe}.html")
                    html_view.save_as_html(out_html)
                except Exception:
                    pass
        except Exception as e:
            print(f"Failed to render {fp} (glass): {e}")
            continue
    montage_path = None
    if overlay_paths:
        montage_path = os.path.join(out_dir, 'UofM_Network_Glass_Montage.png')
        compose_montage(overlay_paths, montage_path, cols=4)
    return overlay_paths, montage_path

# ---------- Harvard-Oxford atlas integration (Cortical/Subcortical) ----------

def _discover_harvard_oxford(atlas_kind: str) -> Tuple[nib.Nifti1Image, List[str], Dict[str, int]]:
    """Fetch Harvard-Oxford atlas via nilearn datasets.

    atlas_kind: 'cortical' or 'subcortical'
    Returns (atlas_img, labels, name_to_code)
    """
    if atlas_kind == 'cortical':
        atlas = datasets.fetch_atlas_harvard_oxford('cort-maxprob-thr50-2mm')
    else:
        atlas = datasets.fetch_atlas_harvard_oxford('sub-maxprob-thr50-2mm')
    maps_obj = atlas.maps
    # maps may be a file path (str) or a Niimg-like object depending on nilearn version
    if isinstance(maps_obj, str):
        atlas_img = nib.load(maps_obj)
    else:
        atlas_img = maps_obj  # assume Niimg-like
    # labels list typically includes 'Left XYZ' and 'Right XYZ'
    labels = list(atlas.labels)
    name_to_code: Dict[str, int] = {}
    # The atlas.labels align to integer region codes (starting from 0 for background); map label->index
    for idx, lab in enumerate(labels):
        name_to_code[lab] = idx
    return atlas_img, labels, name_to_code

def _canon_cortical_region(name: str) -> str:
    v = re.sub(r"[_.]+", " ", str(name or "").strip().lower())
    v = re.sub(r"\s+", " ", v)
    # Strip hemisphere tokens from incoming name first
    v = re.sub(r'^(\s*)(left|right)\s+', '', v)
    v = re.sub(r'\s+(left|right)\s*$', '', v)
    # Minimal normalization to align ggseg/common names to Harvard-Oxford labels
    # Occipital/visual
    if re.search(r"lingual", v):
        return "Lingual Gyrus"
    if re.search(r"precuneus", v):
        return "Precuneus"
    if re.search(r"cuneus|pericalcarine|calcarine", v):
        return "Supracalcarine Cortex"
    if re.search(r"lateral\s*occipital|lateraloccipital", v):
        return "Lateral Occipital Cortex, inferior division"
    if re.search(r"occipital\s*pole", v):
        return "Occipital Pole"
    # Central sulcus
    if re.search(r"precentral", v):
        return "Precentral Gyrus"
    if re.search(r"postcentral", v):
        return "Postcentral Gyrus"
    if re.search(r"paracentral", v):
        return "Paracentral Lobule"
    # Temporal
    if re.search(r"heschl|transverse\s*temporal|transversetemporal", v):
        return "Heschl's Gyrus (includes H1 and H2)"
    if re.search(r"planum\s*temporale", v):
        return "Planum temporale"
    if re.search(r"planum\s*polare", v):
        return "Planum polare"
    if re.search(r"fusiform", v):
        return "Temporal Fusiform Cortex, anterior division"
    if re.search(r"middle\s*temporal|middletemporal", v):
        return "Middle Temporal Gyrus, posterior division"
    if re.search(r"inferior\s*temporal|inferiortemporal", v):
        return "Inferior Temporal Gyrus, posterior division"
    if re.search(r"superior\s*temporal|superiortemporal", v):
        return "Superior Temporal Gyrus, posterior division"
    if re.search(r"temporal\s*pole", v):
        return "Temporal Pole"
    if re.search(r"bankssts", v):
        return "Superior Temporal Gyrus, posterior division"
    # Frontal
    if re.search(r"superior\s*frontal|superiorfrontal", v):
        return "Superior Frontal Gyrus"
    if re.search(r"middle\s*frontal|rostral\s*middle\s*frontal|rostralmiddlefrontal", v):
        return "Middle Frontal Gyrus"
    if re.search(r"inferior\s*frontal|inferiorfrontal", v):
        return "Inferior Frontal Gyrus, pars triangularis"
    if re.search(r"pars\s*opercularis", v):
        return "Inferior Frontal Gyrus, pars opercularis"
    if re.search(r"pars\s*triangularis", v):
        return "Inferior Frontal Gyrus, pars triangularis"
    if re.search(r"frontal\s*pole|frontalpole", v):
        return "Frontal Pole"
    if re.search(r"frontal\s*medial|medial\s*orbito\s*frontal|medialorbitofrontal", v):
        return "Frontal Medial Cortex"
    if re.search(r"orbito\s*frontal|lateral\s*orbito\s*frontal|pars\s*orbitalis|parsorbitalis", v):
        return "Frontal Orbital Cortex"
    # Insula/Cingulate
    if re.search(r"insula|insular", v):
        return "Insular Cortex"
    if re.search(r"paracingulate", v):
        return "Paracingulate Gyrus"
    if re.search(r"isthmus.*cingulate|posterior.*cingulate", v):
        return "Cingulate Gyrus, posterior division"
    if re.search(r"caudal\s*anterior\s*cingulate|rostral\s*anterior\s*cingulate|^anterior\s*cingulate", v):
        return "Cingulate Gyrus, anterior division"
    # Parietal
    if re.search(r"supramarginal", v):
        return "Supramarginal Gyrus, anterior division"
    if re.search(r"angular", v):
        return "Angular Gyrus"
    # Subcortical synonyms (Harvard-Oxford subcortical)
    if re.search(r"\bcaudate\b|caudate\s+nucleus", v):
        return "Caudate"
    if re.search(r"\bputamen\b", v):
        return "Putamen"
    if re.search(r"pallid", v):
        return "Pallidum"
    if re.search(r"amygdala", v):
        return "Amygdala"
    if re.search(r"hippocamp", v):
        return "Hippocampus"
    if re.search(r"thalamus", v):
        return "Thalamus"
    if re.search(r"accumbens|nucleus\s+accumbens|accumbens\s+area", v):
        return "Accumbens"
    if re.search(r"entorhinal", v):
        return "Parahippocampal Gyrus, anterior division"
    return name.strip()

def _find_ho_label(labels: List[str], region: str, hemi: str) -> str:
    base = _canon_cortical_region(strip_hemi_tokens(region))
    candidates = []
    if hemi == 'left':
        candidates += [f"Left {base}", f"{base} Left", f"{base} L", base]
    elif hemi == 'right':
        candidates += [f"Right {base}", f"{base} Right", f"{base} R", base]
    else:
        candidates += [base]
    for c in candidates:
        for lab in labels:
            if lab.lower() == c.lower():
                return lab
    # Loose normalized contains match
    def norm(x: str) -> str:
        x = re.sub(r'[_\-]', ' ', x.lower())
        return re.sub(r'\s{2,}', ' ', x).strip()
    for lab in labels:
        if norm(base) in norm(lab):
            if hemi in ('left', 'right'):
                if hemi in lab.lower() or (hemi == 'left' and ' l' in lab.lower()) or (hemi == 'right' and ' r' in lab.lower()):
                    return lab
            else:
                return lab
    # Fallback: token overlap (Jaccard) to find closest label when synonyms differ
    def tokens(x: str) -> set:
        x = norm(x)
        # remove hemisphere words and generic suffixes
        drop = {'left', 'right', 'gyrus', 'cortex', 'lobule', 'lobe', 'area', 'division', 'pole',
                'anterior', 'posterior', 'superior', 'inferior', 'medial', 'lateral', 'pars'}
        toks = re.split(r"\W+", x)
        return {t for t in toks if t and t not in drop}
    base_toks = tokens(base)
    best_lab = ''
    best_score = 0.0
    for lab in labels:
        # respect hemisphere when possible
        if hemi in ('left', 'right') and hemi not in lab.lower():
            continue
        lab_no_hemi = re.sub(r'^(left|right)\s+', '', lab.lower())
        lab_toks = tokens(lab_no_hemi)
        if not base_toks or not lab_toks:
            continue
        inter = len(base_toks & lab_toks)
        union = len(base_toks | lab_toks)
        score = inter / union if union else 0.0
        if score > best_score:
            best_score = score
            best_lab = lab
    if best_score >= 0.5:
        return best_lab
    return ''

def _prepare_cortical_tidy(df: pd.DataFrame) -> pd.DataFrame:
    df = df[df['structure_type'].isin(['Cortical', 'Subcortical/Regional'])]
    # Prefer ggseg_region; fallback to brain_region
    region = df['ggseg_region'].fillna('')
    region = np.where(region != '', region, df.get('brain_region', region))
    df['region_norm'] = [strip_hemi_tokens(str(r).strip()) for r in region]
    hemi_raw = df.get('hemi', pd.Series(['unknown'] * len(df)))
    df['hemi_norm'] = [unify_hemi(h) for h in hemi_raw]
    rows = []
    for _, r in df.iterrows():
        hemi = r['hemi_norm']
        # Infer hemisphere from variable_name/regions if unknown/bilateral
        if hemi in ('unknown', 'bilateral'):
            texts = [str(r.get('variable_name', '')).lower(),
                     str(r.get('ggseg_region', '')).lower(),
                     str(r.get('brain_region', '')).lower(),
                     str(r.get('region_norm', '')).lower()]
            if any(re.search(r'\bleft\b', t) or re.search(r'\blh\b', t) or re.search(r'left\.hemisphere', t) for t in texts):
                hemi = 'left'
            elif any(re.search(r'\bright\b', t) or re.search(r'\brh\b', t) or re.search(r'right\.hemisphere', t) for t in texts):
                hemi = 'right'
        z = float(r.get('z_change', 0.0))
        base_row = dict(region=r['region_norm'], z_value=z)
        # If hemisphere unknown, treat as bilateral
        rows.append({**base_row, 'hemi': hemi if hemi in ('left', 'right') else 'bilateral'})
    tidy = pd.DataFrame(rows)
    tidy = tidy.groupby(['region', 'hemi'], as_index=False)['z_value'].max()
    return tidy

def _render_ho_group(df_group: pd.DataFrame,
                     out_dir_group: str,
                     atlas_kind: str,
                     bg_img: Optional[nib.Nifti1Image]) -> None:
    os.makedirs(out_dir_group, exist_ok=True)
    tidy = _prepare_cortical_tidy(df_group)
    try:
        atlas_img, labels, name_to_code = _discover_harvard_oxford('cortical' if atlas_kind == 'cortical' else 'subcortical')
        image_paths_overlay = []
        glass_paths = []
        cov_rows = []
        match_rows = []
        for _, r in tidy.iterrows():
            hemi_req = r['hemi']
            region = r['region']
            hemi_list = [hemi_req] if hemi_req in ('left','right') else ['left','right']
            for hemi in hemi_list:
                lab = _find_ho_label(labels, region, hemi)
                matched = bool(lab)
                cov_rows.append({'region': region, 'hemi': hemi, 'matched': int(matched)})
                match_rows.append({'atlas': f'Harvard-Oxford-{atlas_kind}', 'region': region, 'hemi_request': hemi_req, 'hemi_assigned': hemi, 'label': lab, 'matched': int(matched)})
                if not matched:
                    continue
                # Harvard-Oxford labels map to integer codes; use name_to_code
                lbl_idx = int(name_to_code.get(lab, -1))
                if lbl_idx <= 0:
                    continue
                data = atlas_img.get_fdata()
                mask = (data == lbl_idx).astype(np.float32)
                stat = mask * float(r['z_value'])
                stat_img = nib.Nifti1Image(stat, atlas_img.affine)
                fname = f"HO_{atlas_kind}_{_canon_cortical_region(region).replace(' ', '_')}_{hemi}.png"
                outp = os.path.join(out_dir_group, fname)
                ok = save_overlay(stat_img, outp,
                                   title=f"{_canon_cortical_region(region)} ({hemi}) | Z={r['z_value']:.2f}",
                                   bg_img=bg_img,
                                   threshold=0.01,
                                   cmap='RdBu_r',
                                   alpha=0.9,
                                   black_bg=True,
                                   symmetric_cbar=True,
                                   dpi=300)
                if ok:
                    image_paths_overlay.append(outp)
                # Glass brain
                glass_png = os.path.join(out_dir_group, f"HO_{atlas_kind}_Glass_{_canon_cortical_region(region).replace(' ', '_')}_{hemi}.png")
                ok_g = save_glass_brain(stat_img, glass_png,
                                        title=f"Glass: {_canon_cortical_region(region)} ({hemi}) | Z={r['z_value']:.2f}",
                                        threshold=0.01,
                                        cmap='RdBu_r',
                                        plot_abs=False,
                                        black_bg=True,
                                        dpi=300)
                if ok_g:
                    glass_paths.append(glass_png)
        cov_df = pd.DataFrame(cov_rows)
        gsum = cov_df.groupby(['hemi'], as_index=False).agg(total=('matched', 'size'), matched=('matched', 'sum'))
        gsum['unmatched'] = gsum['total'] - gsum['matched']
        gsum['coverage_rate'] = np.where(gsum['total'] > 0, gsum['matched'] / gsum['total'], np.nan)
        cov_path = os.path.join(out_dir_group, f'HO_{atlas_kind}_Coverage_Summary.csv')
        gsum.to_csv(cov_path, index=False)
        # Emit label matching diagnostics
        try:
            match_df = pd.DataFrame(match_rows)
            match_df.to_csv(os.path.join(out_dir_group, f'HO_{atlas_kind}_Label_Matches.csv'), index=False)
            umap_df = match_df[match_df['matched'] == 0].copy()
            umap_df['suggest_label'] = ''
            if not umap_df.empty:
                def suggest(row):
                    b = _canon_cortical_region(strip_hemi_tokens(row['region']))
                    def norm(x: str) -> str:
                        x = re.sub(r'[_\-]', ' ', x.lower())
                        return re.sub(r'\s{2,}', ' ', x).strip()
                    options = [lab for lab in labels if norm(b) in norm(lab)]
                    options_hemi = [lab for lab in options if row['hemi_assigned'] in lab.lower()]
                    return options_hemi[0] if options_hemi else (options[0] if options else '')
                umap_df['suggest_label'] = [suggest(r) for _, r in umap_df.iterrows()]
            umap_df.to_csv(os.path.join(out_dir_group, f'HO_{atlas_kind}_Unmapped.csv'), index=False)
        except Exception:
            pass
        montage_path = os.path.join(out_dir_group, f'HO_{atlas_kind}_Overlay_Montage.png')
        compose_montage(image_paths_overlay, montage_path, cols=4)
        print(f"Saved {len(image_paths_overlay)} overlays; montage: {montage_path} | atlas=Harvard-Oxford-{atlas_kind}")

        if glass_paths:
            glass_montage = os.path.join(out_dir_group, f'HO_{atlas_kind}_Glass_Montage.png')
            compose_montage(glass_paths, glass_montage, cols=4)
            print(f"Saved {len(glass_paths)} glass brains; montage: {glass_montage} | atlas=Harvard-Oxford-{atlas_kind}")
    except Exception as e:
        print(f"Harvard-Oxford atlas not available for {atlas_kind}; skipped overlays. Reason: {e}")

def _render_ho_surface_cortical(df_group: pd.DataFrame,
                                out_dir_group: str,
                                bg_img: Optional[nib.Nifti1Image]) -> None:
    """Project HO cortical overlays to fsaverage surfaces and save images & montage."""
    os.makedirs(out_dir_group, exist_ok=True)
    tidy = _prepare_cortical_tidy(df_group)
    try:
        atlas_img, labels, name_to_code = _discover_harvard_oxford('cortical')
        img_paths: List[str] = []
        fsavg = datasets.fetch_surf_fsaverage()
        for _, r in tidy.iterrows():
            hemi_req = r['hemi']
            region = r['region']
            hemi_list = [hemi_req] if hemi_req in ('left','right') else ['left','right']
            for hemi in hemi_list:
                lab = _find_ho_label(labels, region, hemi)
                if not lab:
                    continue
                lbl_idx = int(name_to_code.get(lab, -1))
                if lbl_idx <= 0:
                    continue
                data = atlas_img.get_fdata()
                mask = (data == lbl_idx).astype(np.float32)
                stat_img = nib.Nifti1Image(mask * float(r['z_value']), atlas_img.affine)
                try:
                    if hemi in ('left', 'bilateral'):
                        tex_l = surface.vol_to_surf(stat_img, fsavg.pial_left)
                        fig_l = plotting.plot_surf_stat_map(fsavg.infl_left, tex_l,
                                                             bg_map=fsavg.sulc_left,
                                                             hemi='left', title=f"{_canon_cortical_region(region)} left (surf)",
                                                             cmap='RdBu_r')
                        out_png_l = os.path.join(out_dir_group, f"HO_cortical_Surface_{_canon_cortical_region(region).replace(' ', '_')}_left.png")
                        fig_l.savefig(out_png_l, dpi=300)
                        try:
                            out_pdf_l = os.path.splitext(out_png_l)[0] + '.pdf'
                            fig_l.savefig(out_pdf_l, dpi=300)
                        except Exception:
                            pass
                        plt.close(fig_l)
                        img_paths.append(out_png_l)
                    if hemi in ('right', 'bilateral'):
                        tex_r = surface.vol_to_surf(stat_img, fsavg.pial_right)
                        fig_r = plotting.plot_surf_stat_map(fsavg.infl_right, tex_r,
                                                             bg_map=fsavg.sulc_right,
                                                             hemi='right', title=f"{_canon_cortical_region(region)} right (surf)",
                                                             cmap='RdBu_r')
                        out_png_r = os.path.join(out_dir_group, f"HO_cortical_Surface_{_canon_cortical_region(region).replace(' ', '_')}_right.png")
                        fig_r.savefig(out_png_r, dpi=300)
                        try:
                            out_pdf_r = os.path.splitext(out_png_r)[0] + '.pdf'
                            fig_r.savefig(out_pdf_r, dpi=300)
                        except Exception:
                            pass
                        plt.close(fig_r)
                        img_paths.append(out_png_r)
                except Exception as e:
                    print(f"Surface render failure for {region}: {e}")
                    continue
        if img_paths:
            montage_path = os.path.join(out_dir_group, 'HO_cortical_Surface_Overlay_Montage.png')
            compose_montage(img_paths, montage_path, cols=4)
            print(f"Saved {len(img_paths)} surface overlays; montage: {montage_path}")
    except Exception as e:
        print(f"Harvard-Oxford cortical surface rendering skipped. Reason: {e}")

def _write_hemi_diagnostics(df_group: pd.DataFrame, out_dir_group: str) -> None:
    os.makedirs(out_dir_group, exist_ok=True)
    # White Matter
    df_wm = df_group[df_group['structure_type'] == 'White Matter'].copy()
    df_wm['hemi_norm'] = [unify_hemi(h) for h in df_wm.get('hemi', pd.Series(['unknown'] * len(df_wm)))]
    wm_counts = df_wm['hemi_norm'].value_counts().rename_axis('hemi').reset_index(name='count')
    wm_counts.to_csv(os.path.join(out_dir_group, 'WM_Hemisphere_Distribution.csv'), index=False)
    # Cortical/Subcortical
    df_cs = df_group[df_group['structure_type'].isin(['Cortical', 'Subcortical/Regional'])].copy()
    df_cs['hemi_norm'] = [unify_hemi(h) for h in df_cs.get('hemi', pd.Series(['unknown'] * len(df_cs)))]
    cs_counts = df_cs['hemi_norm'].value_counts().rename_axis('hemi').reset_index(name='count')
    cs_counts.to_csv(os.path.join(out_dir_group, 'Cortical_Subcortical_Hemisphere_Distribution.csv'), index=False)

def _prepare_wm_tidy(df: pd.DataFrame) -> pd.DataFrame:
    """From the unified CSV, keep white matter rows and produce tract/hemi tidy table."""
    df = df[df['structure_type'] == 'White Matter']
    df = df[df['tract'].notna() & (df['tract'] != '')]
    hemi_raw = df.get('hemi', pd.Series(['unknown'] * len(df)))
    df['hemi_norm'] = [unify_hemi(h) for h in hemi_raw]
    rows = []
    for _, r in df.iterrows():
        name = canonicalize_tract(strip_hemi_tokens(r['tract']))
        hemi = r['hemi_norm']
        # Infer hemisphere from variable_name/regions if unknown/bilateral
        if hemi in ('unknown', 'bilateral'):
            texts = [str(r.get('variable_name', '')).lower(),
                     str(r.get('ggseg_region', '')).lower(),
                     str(r.get('brain_region', '')).lower(),
                     str(r.get('tract', '')).lower()]
            if any(re.search(r'\bleft\b', t) or re.search(r'\blh\b', t) or re.search(r'left\.hemisphere', t) for t in texts):
                hemi = 'left'
            elif any(re.search(r'\bright\b', t) or re.search(r'\brh\b', t) or re.search(r'right\.hemisphere', t) for t in texts):
                hemi = 'right'
        z = float(r.get('z_change', 0.0))
        base_row = dict(tract=name, z_value=z)
        if name in MIDLINE_TRACTS and hemi in ('bilateral', 'unknown'):
            rows.append({**base_row, 'hemi': 'left'})
            rows.append({**base_row, 'hemi': 'right'})
        else:
            rows.append({**base_row, 'hemi': hemi if hemi in ('left', 'right') else 'bilateral'})
    tidy = pd.DataFrame(rows)
    tidy = tidy.groupby(['tract', 'hemi'], as_index=False)['z_value'].max()
    return tidy

def _render_jhu_group(df_group: pd.DataFrame,
                      out_dir_group: str,
                      atlas_dir: Optional[str],
                      atlas_nii: Optional[str],
                      labels_txt: Optional[str],
                      bg_img: Optional[nib.Nifti1Image]) -> None:
    os.makedirs(out_dir_group, exist_ok=True)
    tidy = _prepare_wm_tidy(df_group)
    try:
        base_img, labels, name_to_code, src = fetch_jhu(atlas_dir=atlas_dir, atlas_nii=atlas_nii, labels_txt=labels_txt)
        cov_rows = []
        match_rows = []
        image_paths_overlay = []
        glass_paths = []
        for _, r in tidy.iterrows():
            hemi_req = r['hemi']
            tract = r['tract']
            hemi_list = [hemi_req] if hemi_req in ('left','right') else ['left','right']
            for hemi in hemi_list:
                lab = find_label(labels, tract, hemi)
                matched = bool(lab)
                cov_rows.append({'tract': tract, 'hemi': hemi, 'matched': int(matched)})
                match_rows.append({'atlas': 'JHU-ICBM-tracts', 'tract': tract, 'hemi_request': hemi_req, 'hemi_assigned': hemi, 'label': lab, 'matched': int(matched)})
                if not matched:
                    continue
                stat = build_stat_map(base_img, labels, lab, float(r['z_value']), name_to_code=name_to_code)
                fname = f"WM_JHU_Overlay_{tract.replace(' ', '_')}_{hemi}.png"
                outp = os.path.join(out_dir_group, fname)
                ok = save_overlay(stat, outp,
                                   title=f"{tract} ({hemi}) | Z={r['z_value']:.2f}",
                                   bg_img=bg_img,
                                   threshold=0.01,
                                   cmap='RdBu_r',
                                   alpha=0.9,
                                   black_bg=True,
                                   symmetric_cbar=True,
                                   dpi=300)
                if ok:
                    image_paths_overlay.append(outp)
                # Glass brain
                glass_png = os.path.join(out_dir_group, f"WM_JHU_Glass_{tract.replace(' ', '_')}_{hemi}.png")
                ok_g = save_glass_brain(stat, glass_png,
                                        title=f"Glass: {tract} ({hemi}) | Z={r['z_value']:.2f}",
                                        threshold=0.01,
                                        cmap='RdBu_r',
                                        plot_abs=False,
                                        black_bg=True,
                                        dpi=300)
                if ok_g:
                    glass_paths.append(glass_png)
        cov_df = pd.DataFrame(cov_rows)
        cov_df['group'] = cov_df['tract'].apply(lambda t: (
            'Commissural' if re.search(r'forceps minor|forceps major|corpus callosum|splenium|genu', t, re.I) else
            'Projection' if re.search(r'corticospinal|corona radiata|thalamic radiation', t, re.I) else
            'Association' if re.search(r'uncinate|fronto[- ]?occipital|inferior longitudinal|superior longitudinal|middle longitudinal|cingulum', t, re.I) else
            'Cerebellar' if re.search(r'cerebellar peduncle|mcp|scp|icp', t, re.I) else
            'Limbic' if re.search(r'fornix', t, re.I) else
            'Brainstem' if re.search(r'lemniscus|brainstem|pons|medulla', t, re.I) else
            'Other'
        ))
        gsum = cov_df.groupby('group', as_index=False).agg(total=('matched', 'size'), matched=('matched', 'sum'))
        gsum['unmatched'] = gsum['total'] - gsum['matched']
        gsum['coverage_rate'] = np.where(gsum['total'] > 0, gsum['matched'] / gsum['total'], np.nan)
        gsum.to_csv(os.path.join(out_dir_group, 'WM_JHU_Voxel_Coverage_Summary.csv'), index=False)
        # Emit label matching diagnostics
        try:
            match_df = pd.DataFrame(match_rows)
            match_df.to_csv(os.path.join(out_dir_group, 'WM_JHU_Label_Matches.csv'), index=False)
            umap_df = match_df[match_df['matched'] == 0].copy()
            umap_df['suggest_label'] = ''
            if not umap_df.empty:
                def suggest(row):
                    b = canonicalize_tract(strip_hemi_tokens(row['tract']))
                    def norm(x: str) -> str:
                        x = re.sub(r'[_\-]', ' ', x.lower())
                        return re.sub(r'\s{2,}', ' ', x).strip()
                    options = [lab for lab in labels if norm(b) in norm(lab)]
                    options_hemi = [lab for lab in options if row['hemi_assigned'] in lab.lower()]
                    return options_hemi[0] if options_hemi else (options[0] if options else '')
                umap_df['suggest_label'] = [suggest(r) for _, r in umap_df.iterrows()]
            umap_df.to_csv(os.path.join(out_dir_group, 'WM_JHU_Unmapped.csv'), index=False)
        except Exception:
            pass
        montage_path = os.path.join(out_dir_group, 'WM_JHU_Overlay_Montage.png')
        compose_montage(image_paths_overlay, montage_path, cols=4)
        print(f"Saved {len(image_paths_overlay)} overlays; montage: {montage_path} | atlas_source={src}")
        if glass_paths:
            glass_montage = os.path.join(out_dir_group, 'WM_JHU_Glass_Montage.png')
            compose_montage(glass_paths, glass_montage, cols=4)
            print(f"Saved {len(glass_paths)} glass brains; montage: {glass_montage} | atlas_source={src}")
    except Exception as e:
        print(f"JHU atlas not available; skipped JHU voxel overlays for this group. Reason: {e}")

def main(csv_path: str = CSV_DEFAULT_PARSED,
         atlas_dir: Optional[str] = None,
         atlas_nii: Optional[str] = None,
         labels_txt: Optional[str] = None,
         uofm_dir: Optional[str] = None,
         bg_nii: Optional[str] = None,
         prob_threshold: float = 0.1,
         sig_idp_csv: Optional[str] = None,
         uofm_sig_only: bool = False,
         selected_only: bool = False,
         nitrc_cookie_jar: Optional[str] = None,
         download_uofm_v2: bool = False,
         groups: Optional[List[str]] = None):
    # Load CSV once and derive tract/hemisphere base columns
    os.makedirs(OUT_DIR, exist_ok=True)
    if download_uofm_v2:
        if not uofm_dir:
            uofm_dir = os.path.join(ROOT, 'ischemic', 'atlases', 'uofm_jhu_atlas_v2')
        if not nitrc_cookie_jar:
            raise ValueError("nitrc_cookie_jar is required when download_uofm_v2 is enabled.")
        if not os.path.isdir(uofm_dir) or not discover_uofm_prob_maps(uofm_dir):
            install_uofm_atlas_from_nitrc(NITRC_UOFM_V2_LINK_ID, nitrc_cookie_jar, uofm_dir)
    path_to_use = csv_path
    if not os.path.isfile(path_to_use):
        if os.path.isfile(CSV_DEFAULT_PARSED):
            path_to_use = CSV_DEFAULT_PARSED
        else:
            raise FileNotFoundError(f"CSV not found: {csv_path} (and default {CSV_DEFAULT_PARSED} missing)")
    df_all = pd.read_csv(path_to_use)
    if 'hemi' not in df_all.columns and 'hemisphere' in df_all.columns:
        df_all['hemi'] = df_all['hemisphere']
    if selected_only:
        df_all = filter_selected_idps(df_all)
        if df_all.empty:
            print(f"No rows selected for mapping after filtering: {path_to_use}")
            return

    st = df_all['structure_type'].fillna('').astype(str)
    is_wm = st.str.lower().eq('white matter')

    def _col(c: str) -> pd.Series:
        if c in df_all.columns:
            return df_all[c].fillna('').astype(str)
        return pd.Series([''] * len(df_all))

    wm_tract = _col('wm_tract')
    wm_struct = _col('anatomical_structure')
    wm_var = _col('variable_name')
    wm_br = _col('brain_region')
    wm_gg = _col('ggseg_region')

    tract_wm = wm_tract
    tract_wm = tract_wm.mask(tract_wm.str.strip().eq(''), wm_struct)
    tract_wm = tract_wm.mask(tract_wm.str.strip().eq(''), wm_var)
    tract_wm = tract_wm.mask(tract_wm.str.strip().eq(''), wm_br)
    tract_wm = tract_wm.mask(tract_wm.str.strip().eq(''), wm_gg)

    non_wm = _col('ggseg_region')
    non_wm = non_wm.mask(non_wm.str.strip().eq(''), _col('brain_region'))
    non_wm = non_wm.mask(non_wm.str.strip().eq(''), _col('anatomical_structure'))
    non_wm = non_wm.mask(non_wm.str.strip().eq(''), _col('variable_name'))

    df_all['tract'] = np.where(is_wm.values, tract_wm.values, non_wm.values)
    df_all['tract'] = pd.Series(df_all['tract']).fillna('').astype(str).str.strip()

    # Background image (used both for UofM and WM overlays)
    bg_img = load_bg_img(bg_nii, uofm_dir)
    sig_effects_by_group = load_sig_wm_tract_effects(sig_idp_csv) if sig_idp_csv else {}
    uofm_sig_only = bool(uofm_sig_only) or bool(sig_idp_csv)
    if uofm_sig_only and not sig_idp_csv:
        raise ValueError("uofm_sig_only requires sig_idp_csv.")

    # Determine groups to process
    if groups is None or len(groups) == 0:
        wanted = {'ischemic_vs_control', 'mi_vs_chronic', 'mi_vs_control', 'chronic_vs_control'}
        present = set(df_all['model_comparison'].unique())
        groups = [g for g in ['ischemic_vs_control', 'mi_vs_control', 'chronic_vs_control', 'mi_vs_chronic'] if g in present] or list(present)

    for g in groups:
        print(f"\n=== Processing group: {g} ===")
        df_g = df_all[df_all['model_comparison'] == g].copy()
        if df_g.empty:
            print(f"Skipped group {g}: no rows after filtering.")
            continue

        # UofM probability maps (per-group output folder for organizational clarity)
        if uofm_dir and os.path.isdir(uofm_dir):
            try:
                effects_g = sig_effects_by_group.get(g, {})
                if uofm_sig_only and not effects_g:
                    print(f"Skipped UofM overlays for group {g}: no significant WM tracts found in {sig_idp_csv}")
                    continue
                out_dir_net_group = os.path.join(OUT_DIR_NETWORK, g, 'UofM')
                _, montage_path = generate_uofm_overlays(
                    uofm_dir,
                    out_dir_net_group,
                    bg_img,
                    threshold=prob_threshold,
                    make_interactive=True,
                    effects_by_tract=effects_g,
                    sig_only=uofm_sig_only
                )
                if montage_path:
                    print(f"Saved UofM network overlays; montage: {montage_path}")
                if GLASS_BRAIN_ENABLED:
                    glass_dir = os.path.join(OUT_DIR_NETWORK, g, 'UofM_Glass')
                    _, glass_montage = generate_uofm_glass_overlays(
                        uofm_dir,
                        glass_dir,
                        threshold=prob_threshold,
                        make_interactive=True,
                        effects_by_tract=effects_g,
                        sig_only=uofm_sig_only
                    )
                    if glass_montage:
                        print(f"Saved UofM glass-brain overlays; montage: {glass_montage}")
                    if sig_idp_csv:
                        sum_path = write_uofm_glass_idp_summary(sig_idp_csv, g, glass_dir)
                        if sum_path:
                            print(f"Saved UofM glass IDP summary: {sum_path}")
            except Exception as e:
                print(f"UofM probability map rendering failed for group {g}: {e}")

        # JHU voxel overlays for white matter (group-specific)
        out_dir_group = os.path.join(OUT_DIR, g)
        _render_jhu_group(df_g, out_dir_group, atlas_dir, atlas_nii, labels_txt, bg_img)
        # Ensure PDF companions exist for all PNG outputs in group directory
        ensure_pdf_companions(out_dir_group)

        # Hemisphere diagnostics
        _write_hemi_diagnostics(df_g, out_dir_group)

        # Harvard-Oxford overlays (Cortical & Subcortical)
        out_dir_cort_group = os.path.join(OUT_DIR_CORTICAL, g, 'HarvardOxford')
        _render_ho_group(df_g, out_dir_cort_group, atlas_kind='cortical', bg_img=bg_img)
        # Cortical surface renders (fsaverage)
        _render_ho_surface_cortical(df_g, out_dir_cort_group, bg_img=bg_img)
        ensure_pdf_companions(out_dir_cort_group)
        out_dir_sub_group = os.path.join(OUT_DIR_SUBCORTICAL, g, 'HarvardOxford')
        _render_ho_group(df_g, out_dir_sub_group, atlas_kind='subcortical', bg_img=bg_img)
        ensure_pdf_companions(out_dir_sub_group)

    if GLASS_BRAIN_ENABLED and sig_idp_csv and os.path.exists(sig_idp_csv):
        try:
            df_sig = pd.read_csv(sig_idp_csv)
            if not df_sig.empty and 'model_comparison' in df_sig.columns:
                present = list(dict.fromkeys(df_sig['model_comparison'].astype(str).tolist()))
                preferred = [g for g in ['ischemic_vs_control', 'mi_vs_control', 'chronic_vs_control', 'mi_vs_chronic'] if g in present]
                all_groups = preferred + [g for g in present if g not in preferred]
                paths: List[str] = []
                for g in all_groups:
                    out_dir = os.path.join(OUT_DIR_NETWORK, g, 'UofM_Glass')
                    p = write_uofm_glass_idp_summary(sig_idp_csv, g, out_dir)
                    if p:
                        paths.append(p)
                if paths:
                    dfs = [pd.read_csv(p) for p in paths if p and os.path.exists(p)]
                    if dfs:
                        os.makedirs(OUT_DIR_NETWORK, exist_ok=True)
                        out_all = pd.concat(dfs, ignore_index=True)
                        out_all_path = os.path.join(OUT_DIR_NETWORK, 'UofM_Glass_IDPs_AllGroups.csv')
                        out_all.to_csv(out_all_path, index=False)
                        print(f"Saved UofM glass IDP summary (all groups): {out_all_path}")
        except Exception:
            pass

if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser(description='Generate WM JHU voxel overlays and UofM probability maps.')
    parser.add_argument('csv', nargs='?', default=CSV_DEFAULT_PARSED,
                        help='Path to Brain_IDPs_Parsed_for_Mapping.csv')
    parser.add_argument('--atlas-dir', default=None,
                        help='Directory containing JHU-ICBM-tracts-maxprob-1mm.nii.gz and JHU-ICBM-tracts-labels.txt')
    parser.add_argument('--atlas-nii', default=None,
                        help='Path to JHU-ICBM-tracts-maxprob-1mm.nii.gz')
    parser.add_argument('--labels-txt', default=None,
                        help='Path to JHU-ICBM-tracts-labels.txt')
    parser.add_argument('--uofm-dir', default=os.path.join(ROOT, 'ischemic', 'atlases', 'uofm_jhu_atlas_v1'),
                        help='Root directory of UofM probability maps')
    parser.add_argument('--bg-nii', default=None,
                        help='Background T1 NIfTI for overlays; if omitted, use UofM template or MNI152')
    parser.add_argument('--prob-threshold', type=float, default=0.1,
                        help='Threshold for UofM probability maps (e.g., 0.1)')
    parser.add_argument('--sig-idp-csv', default=None,
                        help='Path to significant IDP CSV; filters UofM maps and encodes z in overlays')
    parser.add_argument('--uofm-sig-only', action='store_true',
                        help='Only render UofM maps that match significant WM tract(s) from --sig-idp-csv')
    parser.add_argument('--selected-only', action='store_true',
                        help='Only process rows selected/significant in the input CSV (selected_for_mapping or p<0.05 logic)')
    parser.add_argument('--nitrc-cookie-jar', default=None,
                        help='Path to a Netscape/Mozilla cookies.txt for nitrc.org (used for download options)')
    parser.add_argument('--download-uofm-v2', action='store_true',
                        help='Download and install UofM_JHU_Atlas_v2 into --uofm-dir using --nitrc-cookie-jar')
    parser.add_argument('--glass-brain', action='store_true',
                        help='Also render glass brain view of UofM probability maps')
    parser.add_argument('--groups', nargs='+', default=None,
                        help='Model comparison groups to process (e.g., ischemic_vs_control mi_vs_chronic mi_vs_control). Default: present in CSV')
    args = parser.parse_args()
    GLASS_BRAIN_ENABLED = bool(getattr(args, 'glass_brain', False))
    main(args.csv,
         atlas_dir=args.atlas_dir,
         atlas_nii=args.atlas_nii,
         labels_txt=args.labels_txt,
         uofm_dir=args.uofm_dir,
         bg_nii=args.bg_nii,
         prob_threshold=args.prob_threshold,
         sig_idp_csv=args.sig_idp_csv,
         uofm_sig_only=bool(getattr(args, 'uofm_sig_only', False)),
         selected_only=bool(getattr(args, 'selected_only', False)),
         nitrc_cookie_jar=getattr(args, 'nitrc_cookie_jar', None),
         download_uofm_v2=bool(getattr(args, 'download_uofm_v2', False)),
         groups=args.groups)
