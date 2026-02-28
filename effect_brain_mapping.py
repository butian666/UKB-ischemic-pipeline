# ============================================
# IDP脑图映射系统 - 基于nilearn的精确解剖定位
# idp_brain_mapping_precise.py
# ============================================

import os
import re
import numpy as np
import pandas as pd
import nibabel as nib
import matplotlib.pyplot as plt
import seaborn as sns
from typing import Dict, List, Tuple, Optional, Union
from dataclasses import dataclass
import warnings
warnings.filterwarnings('ignore')

from nilearn import datasets, plotting, image
from nilearn.image import resample_to_img
from nilearn.plotting import cm
from nilearn.maskers import NiftiLabelsMasker, NiftiSpheresMasker

# ============================================
# 模块 1: 脑区字典（基于您提供的R字典）
# ============================================

@dataclass
class BrainRegion:
    """脑区信息数据类"""
    pattern: str
    brain_region: str
    structure_type: str
    atlas_type: str
    ggseg_region: str
    priority: int
    hemisphere: str = 'bilateral'
    atlas_label: Optional[int] = None
    mni_coords: Optional[Tuple[float, float, float]] = None

class BrainAtlasDictionary:
    """脑区映射字典 - Python版本"""
    
    def __init__(self):
        self.create_brain_regions()
        self.create_measure_types()
        self.load_atlases()
    
    def create_brain_regions(self):
        """创建脑区映射规则（基于您的R字典）"""
        self.brain_regions = [
            # === 额叶区域 ===
            BrainRegion("superiorfrontal|superior\\.frontal|front\\.sup", "Superior Frontal", "Cortical", "cortical", "superiorfrontal", 1),
            BrainRegion("rostralmiddlefrontal|rostral\\.middle\\.frontal", "Rostral Middle Frontal", "Cortical", "cortical", "rostralmiddlefrontal", 1),
            BrainRegion("caudalmiddlefrontal|caudal\\.middle\\.frontal", "Caudal Middle Frontal", "Cortical", "cortical", "caudalmiddlefrontal", 1),
            BrainRegion("middlefrontal|middle\\.frontal|front\\.middle", "Middle Frontal", "Cortical", "cortical", "rostralmiddlefrontal", 2),
            BrainRegion("lateralorbitofrontal|lateral\\.orbitofrontal", "Lateral Orbitofrontal", "Cortical", "cortical", "lateralorbitofrontal", 1),
            BrainRegion("medialorbitofrontal|medial\\.orbitofrontal", "Medial Orbitofrontal", "Cortical", "cortical", "medialorbitofrontal", 1),
            BrainRegion("parstriangularis|pars\\.triangularis", "Pars Triangularis", "Cortical", "cortical", "parstriangularis", 1),
            BrainRegion("parsopercularis|pars\\.opercularis", "Pars Opercularis", "Cortical", "cortical", "parsopercularis", 1),
            BrainRegion("parsorbitalis|pars\\.orbitalis", "Pars Orbitalis", "Cortical", "cortical", "parsorbitalis", 1),
            BrainRegion("precentral|g\\.precentral", "Precentral", "Cortical", "cortical", "precentral", 1),
            BrainRegion("paracentral|g\\.s\\.paracentral", "Paracentral", "Cortical", "cortical", "paracentral", 1),
            BrainRegion("frontalpole|frontal\\.pole", "Frontal Pole", "Cortical", "cortical", "frontalpole", 1),
            
            # === 顶叶区域 ===
            BrainRegion("superiorparietal|superior\\.parietal", "Superior Parietal", "Cortical", "cortical", "superiorparietal", 1),
            BrainRegion("inferiorparietal|inferior\\.parietal", "Inferior Parietal", "Cortical", "cortical", "inferiorparietal", 1),
            BrainRegion("postcentral|g\\.postcentral", "Postcentral", "Cortical", "cortical", "postcentral", 1),
            BrainRegion("precuneus|g\\.precuneus", "Precuneus", "Cortical", "cortical", "precuneus", 1),
            BrainRegion("supramarginal", "Supramarginal", "Cortical", "cortical", "supramarginal", 1),
            BrainRegion("angular", "Angular", "Cortical", "cortical", "inferiorparietal", 1),
            
            # === 颞叶区域 ===
            BrainRegion("superiortemporal|superior\\.temporal", "Superior Temporal", "Cortical", "cortical", "superiortemporal", 1),
            BrainRegion("middletemporal|middle\\.temporal", "Middle Temporal", "Cortical", "cortical", "middletemporal", 1),
            BrainRegion("inferiortemporal|inferior\\.temporal", "Inferior Temporal", "Cortical", "cortical", "inferiortemporal", 1),
            BrainRegion("transversetemporal|transverse\\.temporal", "Transverse Temporal", "Cortical", "cortical", "transversetemporal", 1),
            BrainRegion("temporalpole|temporal\\.pole|pole\\.temporal", "Temporal Pole", "Cortical", "cortical", "temporalpole", 1),
            BrainRegion("parahippocampal", "Parahippocampal", "Cortical", "cortical", "parahippocampal", 1),
            BrainRegion("entorhinal", "Entorhinal", "Cortical", "cortical", "entorhinal", 1),
            BrainRegion("fusiform", "Fusiform", "Cortical", "cortical", "fusiform", 1),
            BrainRegion("bankssts|banks\\.sts", "Banks STS", "Cortical", "cortical", "bankssts", 1),
            
            # === 枕叶区域 ===
            BrainRegion("lateraloccipital|lateral\\.occipital", "Lateral Occipital", "Cortical", "cortical", "lateraloccipital", 1),
            BrainRegion("pericalcarine", "Pericalcarine", "Cortical", "cortical", "pericalcarine", 1),
            BrainRegion("cuneus|g\\.cuneus", "Cuneus", "Cortical", "cortical", "cuneus", 1),
            BrainRegion("lingual", "Lingual", "Cortical", "cortical", "lingual", 1),
            
            # === 扣带回 ===
            BrainRegion("rostralanteriorcingulate", "Rostral Anterior Cingulate", "Cortical", "cortical", "rostralanteriorcingulate", 1),
            BrainRegion("caudalanteriorcingulate", "Caudal Anterior Cingulate", "Cortical", "cortical", "caudalanteriorcingulate", 1),
            BrainRegion("posteriorcingulate", "Posterior Cingulate", "Cortical", "cortical", "posteriorcingulate", 1),
            BrainRegion("isthmuscingulate", "Isthmus Cingulate", "Cortical", "cortical", "isthmuscingulate", 1),
            
            # === 脑岛 ===
            BrainRegion("insula", "Insula", "Cortical", "cortical", "insula", 1),
            
            # === 皮下核团 ===
            BrainRegion("caudate", "Caudate", "Subcortical", "subcortical", "caudate", 1),
            BrainRegion("putamen", "Putamen", "Subcortical", "subcortical", "putamen", 1),
            BrainRegion("pallidum", "Pallidum", "Subcortical", "subcortical", "pallidum", 1),
            BrainRegion("accumbens", "Nucleus Accumbens", "Subcortical", "subcortical", "accumbens", 1),
            
            # === 海马体和杏仁核 ===
            BrainRegion("whole\\.amygdala|amygdala$", "Amygdala", "Subcortical", "subcortical", "amygdala", 1),
            BrainRegion("whole\\.hippocampus|hippocampus$", "Hippocampus", "Subcortical", "subcortical", "hippocampus", 1),
            
            # === 丘脑 ===
            BrainRegion("whole\\.thalamus|thalamus$", "Thalamus", "Subcortical", "subcortical", "thalamus", 1),

            # === 丘脑亚核与相关缩写（增强覆盖） ===
            BrainRegion("VLa|VL\\.a|Ventral\\.?Lateral\\.?anterior", "Ventral Lateral Anterior (Thalamus)", "Subcortical", "subcortical", "thalamus", 1),
            BrainRegion("VLp|VL\\.p|Ventral\\.?Lateral\\.?posterior", "Ventral Lateral Posterior (Thalamus)", "Subcortical", "subcortical", "thalamus", 1),
            BrainRegion("VA\\b|VAmc|VA\\.mc|Ventral\\.?Anterior(.*magnocellular)?", "Ventral Anterior (Thalamus)", "Subcortical", "subcortical", "thalamus", 1),
            BrainRegion("MDl|MD\\.l|Mediodorsal.*lateral", "Mediodorsal Lateral (Thalamus)", "Subcortical", "subcortical", "thalamus", 1),
            BrainRegion("MDm|MD\\.m|Mediodorsal.*medial", "Mediodorsal Medial (Thalamus)", "Subcortical", "subcortical", "thalamus", 2),
            BrainRegion("CL\\b|central\\s+lateral", "Central Lateral (Thalamus)", "Subcortical", "subcortical", "thalamus", 1),
            BrainRegion("CM\\b|centromedian", "Centromedian (Thalamus)", "Subcortical", "subcortical", "thalamus", 1),
            BrainRegion("LP\\b|lateral\\s+posterior", "Lateral Posterior (Thalamus)", "Subcortical", "subcortical", "thalamus", 1),
            BrainRegion("LGN|lateral\\s+geniculate", "Lateral Geniculate Nucleus", "Subcortical", "subcortical", "thalamus", 2),
            BrainRegion("MGN|medial\\s+geniculate", "Medial Geniculate Nucleus", "Subcortical", "subcortical", "thalamus", 2),
            BrainRegion("PuA|pulvinar\\.?anterior", "Pulvinar Anterior", "Subcortical", "subcortical", "thalamus", 1),
            BrainRegion("PuM|pulvinar\\.?medial", "Pulvinar Medial", "Subcortical", "subcortical", "thalamus", 1),
            BrainRegion("PuL|pulvinar\\.?lateral", "Pulvinar Lateral", "Subcortical", "subcortical", "thalamus", 1),
            BrainRegion("PuI|pulvinar\\.?inferior", "Pulvinar Inferior", "Subcortical", "subcortical", "thalamus", 1),

            # === 海马/杏仁体相关细分（增强覆盖） ===
            BrainRegion("HATA", "Hippocampal Amygdaloid Transition Area", "Subcortical", "subcortical", "amygdala", 1),
            BrainRegion(r"Corticoamygdaloid.*transition|Cortico-amygdaloid.*transition|Corticoamygdaloid\.transitio", "Corticoamygdaloid Transition", "Subcortical", "subcortical", "amygdala", 1),
            BrainRegion(r"GC[-_]?ML[-_]?DG", "Dentate Gyrus (GC-ML-DG)", "Subcortical", "subcortical", "hippocampus", 1),
            BrainRegion("CA1|CA2|CA3|CA4|subiculum|presubiculum|hippocampal\\.head|hippocampal\\.body", "Hippocampal Subfields", "Subcortical", "subcortical", "hippocampus", 2),
            
            # === 小脑 ===
            BrainRegion("cerebellum|cerebell", "Cerebellum", "Cerebellum", "cerebellum", "", 2),
            
            # === 脑干 ===
            BrainRegion("whole\\.brainstem|brain\\.stem|brainstem", "Brainstem", "Brainstem", "brainstem", "", 1),
            BrainRegion("midbrain", "Midbrain", "Brainstem", "brainstem", "", 1),
            BrainRegion("pons", "Pons", "Brainstem", "brainstem", "", 1),
            BrainRegion("medulla", "Medulla", "Brainstem", "brainstem", "", 1),
            
            # === 白质纤维束 ===
            BrainRegion("corpus\\.callosum", "Corpus Callosum", "White Matter", "white_matter", "", 2),
            BrainRegion("forceps\\.?minor|genu", "Genu of corpus callosum", "White Matter", "white_matter", "", 1),
            BrainRegion("forceps\\.?major|splenium", "Splenium of corpus callosum", "White Matter", "white_matter", "", 1),
            BrainRegion("body\\.?of\\.?corpus\\.?callosum", "Body of corpus callosum", "White Matter", "white_matter", "", 1),
            BrainRegion("anterior\\.?thalamic\\.?radiation|\\batr\\b", "Anterior thalamic radiation", "White Matter", "white_matter", "", 1),
            BrainRegion("corticospinal", "Corticospinal tract", "White Matter", "white_matter", "", 1),
            BrainRegion("superior\\.longitudinal\\.fasciculus", "Superior Longitudinal Fasciculus", "White Matter", "white_matter", "", 1),
            BrainRegion("inferior\\.longitudinal\\.fasciculus", "Inferior Longitudinal Fasciculus", "White Matter", "white_matter", "", 1),
            BrainRegion("uncinate\\.fasciculus", "Uncinate Fasciculus", "White Matter", "white_matter", "", 1),
            BrainRegion("cingulum.*hippocampus", "Cingulum (hippocampus)", "White Matter", "white_matter", "", 1),
            BrainRegion("cingulum.*cingulate|cingulate\\.?gyrus", "Cingulum (cingulate gyrus)", "White Matter", "white_matter", "", 1),
            BrainRegion("cingulum", "Cingulum", "White Matter", "white_matter", "", 2),
            BrainRegion("inferior\\.?fronto\\.?occipital|\\bifo\\b", "Sagittal stratum", "White Matter", "white_matter", "", 1),
            BrainRegion("superior\\.?fronto\\.?occipital|\\bsfof\\b", "Superior fronto-occipital fasciculus", "White Matter", "white_matter", "", 1),
            BrainRegion("internal\\.capsule", "Internal Capsule", "White Matter", "white_matter", "", 2),
            BrainRegion("corona\\.radiata", "Corona Radiata", "White Matter", "white_matter", "", 2),
            
            # === 通用匹配模式 ===
            BrainRegion("frontal|front", "Frontal", "Cortical", "cortical", "superiorfrontal", 3),
            BrainRegion("parietal|pariet", "Parietal", "Cortical", "cortical", "superiorparietal", 3),
            BrainRegion("temporal|temp", "Temporal", "Cortical", "cortical", "middletemporal", 3),
            BrainRegion("occipital|occip", "Occipital", "Cortical", "cortical", "lateraloccipital", 3),
            BrainRegion("cingulat|cingul", "Cingulate", "Cortical", "cortical", "rostralanteriorcingulate", 3),
        ]
    
    def create_measure_types(self):
        """创建测量类型映射"""
        self.measure_patterns = [
            (r"^area\.", "Surface Area"),
            (r"^volume\.", "Volume"),
            (r"mean\.thickness", "Cortical Thickness"),
            (r"grey\.white\.contrast", "Grey-White Contrast"),
            (r"weighted\.mean\.(fa|md|l1|l2|l3|mo|isovf|icvf|od)", "DTI Tract Weighted"),
            (r"mean\.(fa|md|l1|l2|l3|mo|isovf|icvf|od).*on\.fa\.skeleton", "DTI Skeleton"),
            (r"mean\.(fa|md|l1|l2|l3|mo|isovf|icvf|od)", "DTI Measure"),
            (r"median\.bold|90th\.percentile.*z\.statistic|median\.z\.statistic", "fMRI Activation"),
            (r"head\.motion|tfmri\.head\.motion", "Head Motion"),
        ]
    
    def load_atlases(self):
        """加载标准脑图谱"""
        print("Loading brain atlases...")
        
        # 加载多个图谱以获得更精确的定位
        self.atlases = {}
        
        # AAL图谱
        try:
            self.atlases['AAL'] = datasets.fetch_atlas_aal()
            print("✓ Loaded AAL atlas")
        except:
            print("! Could not load AAL atlas")
        
        # Harvard-Oxford图谱
        try:
            self.atlases['Harvard-Oxford'] = datasets.fetch_atlas_harvard_oxford('cort-maxprob-thr25-2mm')
            self.atlases['Harvard-Oxford-sub'] = datasets.fetch_atlas_harvard_oxford('sub-maxprob-thr25-2mm')
            print("✓ Loaded Harvard-Oxford atlas")
        except:
            print("! Could not load Harvard-Oxford atlas")
        
        # Destrieux图谱
        try:
            self.atlases['Destrieux'] = datasets.fetch_atlas_destrieux_2009()
            print("✓ Loaded Destrieux atlas")
        except:
            print("! Could not load Destrieux atlas")

# ============================================
# 模块 2: IDP解析器
# ============================================

class IDPParser:
    """IDP变量名解析器"""
    
    def __init__(self, dictionary: BrainAtlasDictionary):
        self.dictionary = dictionary
        self.atlas_coordinates = self.load_atlas_coordinates()
    
    def load_atlas_coordinates(self):
        """加载标准脑区MNI坐标"""
        # 这里使用预定义的MNI坐标，可以从图谱中提取或使用标准值
        coords = {
            # 额叶
            "Superior Frontal": (-20, 35, 45),
            "Middle Frontal": (-35, 30, 35),
            "Inferior Frontal": (-45, 30, 15),
            "Precentral": (-40, -5, 50),
            "Frontal Pole": (-10, 65, 20),
            "Orbitofrontal": (-25, 35, -15),
            "Frontal": (-25, 35, 35),
            
            # 顶叶
            "Superior Parietal": (-25, -60, 60),
            "Inferior Parietal": (-45, -45, 45),
            "Postcentral": (-40, -25, 55),
            "Precuneus": (-5, -60, 40),
            "Supramarginal": (-55, -35, 30),
            "Angular": (-45, -60, 35),
            "Parietal": (-30, -55, 50),
            
            # 颞叶
            "Superior Temporal": (-55, -20, 5),
            "Middle Temporal": (-60, -35, -5),
            "Inferior Temporal": (-50, -30, -20),
            "Temporal Pole": (-35, 10, -35),
            "Fusiform": (-35, -40, -20),
            "Parahippocampal": (-25, -20, -20),
            "Entorhinal": (-20, -10, -30),
            "Temporal": (-55, -25, -5),
            
            # 枕叶
            "Lateral Occipital": (-35, -85, 5),
            "Cuneus": (-5, -80, 25),
            "Lingual": (-15, -70, -5),
            "Pericalcarine": (-10, -85, 5),
            "Occipital": (-20, -85, 10),
            
            # 边缘系统
            "Anterior Cingulate": (-5, 35, 15),
            "Posterior Cingulate": (-5, -45, 30),
            "Isthmus Cingulate": (-5, -45, 20),
            "Insula": (-35, 0, 5),
            "Cingulate": (-5, 20, 25),
            
            # 皮下结构
            "Hippocampus": (-25, -20, -15),
            "Amygdala": (-20, -5, -20),
            "Thalamus": (-10, -15, 10),
            "Caudate": (-12, 10, 10),
            "Putamen": (-25, 5, 0),
            "Pallidum": (-20, -2, -2),
            "Nucleus Accumbens": (-10, 10, -8),
            
            # 小脑和脑干
            "Cerebellum": (-20, -60, -35),
            "Brainstem": (0, -30, -40),
            
            # 白质
            "Corpus Callosum": (0, 0, 25),
            "Genu of corpus callosum": (0, 25, 20),
            "Body of corpus callosum": (0, 0, 25),
            "Splenium of corpus callosum": (0, -45, 20),
            "Internal Capsule": (-20, 0, 10),
            "Corona Radiata": (-25, -10, 30),
            "Cingulum": (-10, -10, 30),
            "Cingulum (cingulate gyrus)": (-8, -5, 30),
            "Cingulum (hippocampus)": (-20, -25, 10),
            "Superior Longitudinal Fasciculus": (-35, -20, 35),
            "Inferior Longitudinal Fasciculus": (-40, -35, -10),
            "Uncinate Fasciculus": (-30, 10, -15),
            "Anterior thalamic radiation": (-12, 2, 12),
            "Corticospinal tract": (-18, -18, 40),
            "Sagittal stratum": (-30, -35, 5),
            "Superior fronto-occipital fasciculus": (-28, -20, 30),
        }
        return coords

    
    def parse_variable_name(self, variable_name: str) -> Dict:
        """解析单个IDP变量名"""
        result = {
            "original_name": variable_name,
            "brain_region": "Unknown",
            "hemisphere": "bilateral",
            "measure_type": "Unknown",
            "structure_type": "Unknown",
            "mni_coords": None,
            "confidence": 0.0,
            "atlas_label": None
        }
        
        if not variable_name or pd.isna(variable_name):
            return result
        
        var_lower = variable_name.lower()
        
        # 1. 解析半球
        if "left" in var_lower or "_l_" in var_lower or ".l." in var_lower:
            result["hemisphere"] = "left"
        elif "right" in var_lower or "_r_" in var_lower or ".r." in var_lower:
            result["hemisphere"] = "right"
        
        # 2. 解析测量类型
        for pattern, measure_type in self.dictionary.measure_patterns:
            if re.search(pattern, var_lower):
                result["measure_type"] = measure_type
                break
        
        # 3. 解析脑区（按优先级排序）
        best_match = None
        best_priority = float('inf')
        
        for region in sorted(self.dictionary.brain_regions, key=lambda x: x.priority):
            if re.search(region.pattern, var_lower):
                if region.priority < best_priority:
                    best_priority = region.priority
                    best_match = region
        
        if best_match:
            result["brain_region"] = best_match.brain_region
            result["structure_type"] = best_match.structure_type
            result["confidence"] = 1.0 / best_match.priority
            
            # 获取MNI坐标
            coord_key = best_match.brain_region
            if coord_key not in self.atlas_coordinates:
                coord_alias = {
                    "Rostral Middle Frontal": "Middle Frontal",
                    "Caudal Middle Frontal": "Middle Frontal",
                    "Lateral Orbitofrontal": "Orbitofrontal",
                    "Medial Orbitofrontal": "Orbitofrontal",
                    "Pars Triangularis": "Inferior Frontal",
                    "Pars Opercularis": "Inferior Frontal",
                    "Pars Orbitalis": "Inferior Frontal",
                    "Transverse Temporal": "Superior Temporal",
                    "Banks STS": "Superior Temporal",
                    "Rostral Anterior Cingulate": "Anterior Cingulate",
                    "Caudal Anterior Cingulate": "Anterior Cingulate",
                    "Ventral Lateral Anterior (Thalamus)": "Thalamus",
                    "Ventral Lateral Posterior (Thalamus)": "Thalamus",
                    "Ventral Anterior (Thalamus)": "Thalamus",
                    "Mediodorsal Lateral (Thalamus)": "Thalamus",
                    "Mediodorsal Medial (Thalamus)": "Thalamus",
                    "Central Lateral (Thalamus)": "Thalamus",
                    "Centromedian (Thalamus)": "Thalamus",
                    "Lateral Posterior (Thalamus)": "Thalamus",
                    "Pulvinar Anterior": "Thalamus",
                    "Pulvinar Medial": "Thalamus",
                    "Pulvinar Lateral": "Thalamus",
                    "Pulvinar Inferior": "Thalamus",
                }.get(coord_key)
                if coord_alias is not None:
                    coord_key = coord_alias

            if coord_key in self.atlas_coordinates:
                coords = self.atlas_coordinates[coord_key]
                # 根据半球调整坐标
                if result["hemisphere"] == "left":
                    coords = (-abs(coords[0]), coords[1], coords[2])
                elif result["hemisphere"] == "right":
                    coords = (abs(coords[0]), coords[1], coords[2])
                result["mni_coords"] = coords
        
        return result
    
    def parse_dataframe(self, df: pd.DataFrame, variable_column: str = "variable_name") -> pd.DataFrame:
        """解析整个数据框"""
        print(f"\nParsing {len(df)} IDP variables...")
        
        parsed_results = []
        for idx, row in df.iterrows():
            if variable_column in df.columns:
                var_name = row[variable_column]
                parsed = self.parse_variable_name(var_name)
                
                # 合并原始数据和解析结果
                result_row = row.to_dict()
                result_row.update(parsed)
                parsed_results.append(result_row)
            
            if (idx + 1) % 100 == 0:
                print(f"  Processed {idx + 1}/{len(df)} variables...")
        
        result_df = pd.DataFrame(parsed_results)
        
        # 打印统计
        self._print_stats(result_df)
        
        return result_df
    
    def _print_stats(self, df: pd.DataFrame):
        """打印解析统计"""
        total = len(df)
        identified = len(df[df['brain_region'] != 'Unknown'])
        
        print(f"\n{'='*50}")
        print(f"Total variables: {total}")
        print(f"Successfully identified: {identified} ({identified/total*100:.1f}%)")
        print(f"Unknown regions: {total - identified}")
        
        print("\nTop 10 Brain Regions:")
        top_regions = df[df['brain_region'] != 'Unknown']['brain_region'].value_counts().head(10)
        for region, count in top_regions.items():
            print(f"  {region}: {count}")

# ============================================
# 模块 3: 精确脑图映射器
# ============================================

class PreciseBrainMapper:
    """使用多个图谱创建精确的脑图映射"""
    
    def __init__(self):
        # 加载MNI模板
        self.template = datasets.load_mni152_template(resolution=2)
        self.brain_mask = datasets.load_mni152_brain_mask(resolution=2)
        
        # 加载图谱
        self.load_precise_atlases()
        
    def load_precise_atlases(self):
        """加载精确的脑图谱"""
        print("\nLoading precise brain atlases...")
        
        # AAL图谱
        try:
            self.aal = datasets.fetch_atlas_aal()
            print("✓ AAL atlas loaded")
        except:
            self.aal = None
        
        # Harvard-Oxford图谱
        try:
            self.ho_cort = datasets.fetch_atlas_harvard_oxford('cort-maxprob-thr25-2mm')
            self.ho_sub = datasets.fetch_atlas_harvard_oxford('sub-maxprob-thr25-2mm')
            print("✓ Harvard-Oxford atlas loaded")
        except:
            self.ho_cort = None
            self.ho_sub = None
    
    def create_roi_from_atlas(self, region_name: str, hemisphere: str = "bilateral") -> Optional[np.ndarray]:
        """从图谱中创建ROI"""
        roi = None
        
        # 尝试从不同图谱中获取ROI
        if self.aal is not None:
            roi = self._get_roi_from_aal(region_name, hemisphere)
        
        if roi is None and self.ho_cort is not None:
            roi = self._get_roi_from_harvard_oxford(region_name, hemisphere)
        
        return roi
    
    def _get_roi_from_aal(self, region_name: str, hemisphere: str) -> Optional[np.ndarray]:
        """从AAL图谱获取ROI"""
        # 这里需要根据AAL的标签映射
        # 简化示例，实际需要完整的映射表
        return None
    
    def _get_roi_from_harvard_oxford(self, region_name: str, hemisphere: str) -> Optional[np.ndarray]:
        """从Harvard-Oxford图谱获取ROI"""
        # 简化示例
        return None
    
    def create_brain_map(self, parsed_df: pd.DataFrame, 
                        value_column: Optional[str] = None,
                        use_effect_weighting: bool = True) -> nib.Nifti1Image:
        """创建精确的脑图"""
        print("\nCreating precise brain map with Gaussian weighting...")
        
        # 初始化空白脑图
        brain_data = np.zeros(self.template.shape[:3])
        weight_data = np.zeros(self.template.shape[:3])
        
        # 确定数值列
        if value_column is None:
            parsed_df['demo_value'] = np.random.randn(len(parsed_df))
            value_column = 'demo_value'
        
        # 映射每个IDP到脑图
        mapped_count = 0
        for idx, row in parsed_df.iterrows():
            if row['mni_coords'] is not None:
                # 获取效应值（如果有）
                effect_size = abs(row.get('effect_magnitude', row.get(value_column, 1.0)))
                
                # 根据效应值动态调整半径
                base_radius = 4  # 基础半径 8mm
                if use_effect_weighting:
                    # 效应值越大，影响范围越大（但有上限）
                    adjusted_radius = base_radius * (1 + min(effect_size, 2.0) * 0.3)
                else:
                    adjusted_radius = base_radius
                
                # 创建高斯加权ROI
                roi = self._create_sphere_roi(  # 现在使用高斯加权版本
                    row['mni_coords'],
                    radius=adjusted_radius,
                    structure_type=row['structure_type']
                )
                
                # 获取值和权重
                value = row[value_column] if pd.notna(row[value_column]) else 1.0
                weight = row.get('confidence', 1.0)
                
                # 使用最大值合并策略（避免重叠区域的过度累加）
                weighted_roi = roi * value * weight
                brain_data = np.maximum(brain_data, weighted_roi)
                
                # 记录权重（用于后续归一化）
                weight_data = np.maximum(weight_data, roi * weight)
                
                mapped_count += 1
                
                # 显示进度
                if (idx + 1) % 50 == 0:
                    print(f"  Processed {idx + 1}/{len(parsed_df)} IDPs...")
        
        # 应用额外的平滑（使整体效果更自然）
        brain_img = nib.Nifti1Image(brain_data, self.template.affine)
        brain_img = image.smooth_img(brain_img, fwhm=4)  # 4mm FWHM平滑
        
        # 应用脑mask
        brain_img = image.math_img("img * mask", img=brain_img, mask=self.brain_mask)
        
        print(f"✓ Mapped {mapped_count}/{len(parsed_df)} IDPs to brain using Gaussian weighting")
        
        return brain_img

    
    def _create_sphere_roi(self, mni_coords: Tuple, radius: float, 
                        structure_type: str) -> np.ndarray:
        """创建高斯加权的ROI（替代原来的球形ROI）"""
        # 转换MNI坐标到体素坐标
        affine = self.template.affine
        shape = self.template.shape[:3]
        voxel_coords = nib.affines.apply_affine(np.linalg.inv(affine), mni_coords).astype(float)
        if not np.all(np.isfinite(voxel_coords)):
            return np.zeros(shape, dtype=float)
        
        # 计算体素大小
        voxel_size = np.abs(np.diag(affine)[:3]).mean()
        radius_voxels = radius / voxel_size
        
        # 根据结构类型调整参数
        structure_params = {
            "Cortical": {
                "radius_scale": 1.2,      # 皮层区域稍大
                "sigma_factor": 2.0,       # 较平滑的衰减
                "cutoff_multiplier": 2.5   # 截断距离
            },
            "Subcortical": {
                "radius_scale": 0.8,       # 皮下结构较小
                "sigma_factor": 1.5,       # 较陡的衰减
                "cutoff_multiplier": 2.0
            },
            "White Matter": {
                "radius_scale": 1.5,       # 白质束较大
                "sigma_factor": 2.5,       # 平滑衰减
                "cutoff_multiplier": 3.0
            },
            "Cerebellum": {
                "radius_scale": 1.0,
                "sigma_factor": 1.8,
                "cutoff_multiplier": 2.2
            },
            "Brainstem": {
                "radius_scale": 0.9,
                "sigma_factor": 1.6,
                "cutoff_multiplier": 2.0
            }
        }
        
        # 获取结构特定参数，如果没有则使用默认值
        params = structure_params.get(structure_type, {
            "radius_scale": 1.0,
            "sigma_factor": 2.0,
            "cutoff_multiplier": 2.5
        })
        
        # 应用缩放
        scaled_radius = radius_voxels * params["radius_scale"]
        
        # 高斯函数参数
        sigma = scaled_radius / params["sigma_factor"]  # 标准差
        if sigma <= 0 or not np.isfinite(sigma):
            return np.zeros(shape, dtype=float)
        
        # 截断：超过一定距离的设为0（节省计算和避免过度扩散）
        cutoff_distance = scaled_radius * params["cutoff_multiplier"]
        if cutoff_distance <= 0 or not np.isfinite(cutoff_distance):
            return np.zeros(shape, dtype=float)

        x0, y0, z0 = voxel_coords
        x_min = max(int(np.floor(x0 - cutoff_distance)), 0)
        x_max = min(int(np.ceil(x0 + cutoff_distance)) + 1, shape[0])
        y_min = max(int(np.floor(y0 - cutoff_distance)), 0)
        y_max = min(int(np.ceil(y0 + cutoff_distance)) + 1, shape[1])
        z_min = max(int(np.floor(z0 - cutoff_distance)), 0)
        z_max = min(int(np.ceil(z0 + cutoff_distance)) + 1, shape[2])

        if x_min >= x_max or y_min >= y_max or z_min >= z_max:
            return np.zeros(shape, dtype=float)

        x, y, z = np.ogrid[x_min:x_max, y_min:y_max, z_min:z_max]
        distances = np.sqrt((x - x0)**2 + (y - y0)**2 + (z - z0)**2)
        roi_local = np.exp(-(distances**2) / (2 * sigma**2))
        roi_local[distances > cutoff_distance] = 0

        roi = np.zeros(shape, dtype=float)
        roi[x_min:x_max, y_min:y_max, z_min:z_max] = roi_local
        
        # 可选：归一化到[0,1]范围（中心点为1）
        # 这已经自动满足，因为高斯函数在中心点(距离=0)时值为1
        
        return roi


# ============================================
# 模块 4: 多视角可视化器
# ============================================

class MultiViewVisualizer:
    """多视角脑图可视化"""
    
    def __init__(self, output_dir: str = "brain_visualizations"):
        self.output_dir = output_dir
        os.makedirs(output_dir, exist_ok=True)
        
    def create_comprehensive_views(self, brain_img: nib.Nifti1Image, 
                                 parsed_df: pd.DataFrame,
                                 title: str = "IDP Brain Mapping"):
        """创建综合视图"""
        print("\nCreating comprehensive brain visualizations...")
        
        # 1. 三个解剖平面视图
        self.plot_anatomical_planes(brain_img, title)
        
        # 2. 3D玻璃脑视图
        self.plot_glass_brain_views(brain_img, title)
        
        # 3. 表面投影视图
        self.plot_surface_projections(brain_img, title)
        
        # 4. 交互式HTML视图
        self.create_interactive_html(brain_img, title)
        
        # 5. 统计报告
        self.create_statistical_report(parsed_df)
        
        print(f"✓ All visualizations saved to: {self.output_dir}")
    
    def plot_anatomical_planes(self, brain_img: nib.Nifti1Image, title: str):
        """绘制解剖学三个平面的脑图"""
        from matplotlib.gridspec import GridSpec
        
        fig = plt.figure(figsize=(18, 16))
        gs = GridSpec(3, 3, figure=fig, hspace=0.3, wspace=0.2)
        
        # 计算显示阈值
        data = brain_img.get_fdata()
        non_zero = data[data != 0]
        
        if len(non_zero) > 0:
            vmax = np.percentile(np.abs(non_zero), 95)
            threshold = np.percentile(np.abs(non_zero), 70)
        else:
            vmax = 1
            threshold = 0
        
        # 轴位视图 (Axial) - 占据第一行
        ax1 = fig.add_subplot(gs[0, :])  # 第一行所有列
        display_axial = plotting.plot_stat_map(
            brain_img,
            display_mode='z',
            cut_coords=[-20, -10, 0, 10, 20, 30],  # 多个切片
            threshold=threshold,
            vmax=vmax,
            title="Axial View (Transverse)",
            axes=ax1,
            cmap='RdBu_r',
            colorbar=True
        )
        
        # 矢状位视图 (Sagittal) - 占据第二行
        ax2 = fig.add_subplot(gs[1, :])  # 第二行所有列
        display_sagittal = plotting.plot_stat_map(
            brain_img,
            display_mode='x',
            cut_coords=[-30, -15, 0, 15, 30, 45],  # 多个切片
            threshold=threshold,
            vmax=vmax,
            title="Sagittal View",
            axes=ax2,
            cmap='RdBu_r',
            colorbar=True
        )
        
        # 冠状位视图 (Coronal) - 占据第三行
        ax3 = fig.add_subplot(gs[2, :])  # 第三行所有列
        display_coronal = plotting.plot_stat_map(
            brain_img,
            display_mode='y',
            cut_coords=[-40, -25, -10, 5, 20, 35],  # 多个切片
            threshold=threshold,
            vmax=vmax,
            title="Coronal View",
            axes=ax3,
            cmap='RdBu_r',
            colorbar=True
        )
        
        plt.suptitle(f"{title} - Anatomical Planes", fontsize=16, fontweight='bold', y=0.98)
        
        output_file = os.path.join(self.output_dir, "anatomical_planes.png")
        plt.savefig(output_file, dpi=300, bbox_inches='tight')
        plt.show()
        print(f"  ✓ Saved anatomical planes view: {output_file}")
  
    
    def plot_glass_brain_views(self, brain_img: nib.Nifti1Image, title: str):
        """绘制3D玻璃脑视图"""
        fig, axes = plt.subplots(2, 2, figsize=(16, 16))
        
        # 计算阈值
        data = brain_img.get_fdata()
        non_zero = data[data != 0]
        if len(non_zero) > 0:
            threshold = np.percentile(np.abs(non_zero), 30)
            vmax = np.percentile(np.abs(non_zero), 95)
        else:
            threshold = 0
            vmax = 1
        
        # 正交视图
        plotting.plot_glass_brain(
            brain_img,
            display_mode='ortho',
            threshold=threshold,
            vmax=vmax,
            title="Orthogonal Glass Brain",
            axes=axes[0, 0],
            cmap='RdBu_r',
            colorbar=True
        )
        
        # 左右视图
        plotting.plot_glass_brain(
            brain_img,
            display_mode='lr',
            threshold=threshold,
            vmax=vmax,
            title="Left-Right View",
            axes=axes[0, 1],
            cmap='RdBu_r'
        )
        
        # 前后视图
        plotting.plot_glass_brain(
            brain_img,
            display_mode='yz',
            threshold=threshold,
            vmax=vmax,
            title="Anterior-Posterior View",
            axes=axes[1, 0],
            cmap='RdBu_r'
        )
        
        # 多平面视图
        plotting.plot_glass_brain(
            brain_img,
            display_mode='lyrz',
            threshold=threshold,
            vmax=vmax,
            title="Multi-plane Glass Brain",
            axes=axes[1, 1],
            cmap='RdBu_r'
        )
        
        plt.suptitle(f"{title} - 3D Glass Brain Views", fontsize=16, fontweight='bold')
        plt.tight_layout()
        
        output_file = os.path.join(self.output_dir, "glass_brain_3d.png")
        plt.savefig(output_file, dpi=300, bbox_inches='tight')
        plt.show()
        print(f"  ✓ Saved 3D glass brain: {output_file}")
    
    def plot_surface_projections(self, brain_img: nib.Nifti1Image, title: str):
        """绘制脑表面投影"""
        try:
            from nilearn import surface
            
            # 获取fsaverage表面
            fsaverage = datasets.fetch_surf_fsaverage()
            
            # 投影到表面
            texture_left = surface.vol_to_surf(brain_img, fsaverage.pial_left)
            texture_right = surface.vol_to_surf(brain_img, fsaverage.pial_right)
            
            # 计算显示范围
            vmax = max(
                np.percentile(np.abs(texture_left[~np.isnan(texture_left)]), 95) if len(texture_left[~np.isnan(texture_left)]) > 0 else 1,
                np.percentile(np.abs(texture_right[~np.isnan(texture_right)]), 95) if len(texture_right[~np.isnan(texture_right)]) > 0 else 1
            )
            
            fig = plt.figure(figsize=(20, 15))
            
            # 左半球 - 外侧面
            ax1 = fig.add_subplot(2, 3, 1, projection='3d')
            plotting.plot_surf_stat_map(
                fsaverage.pial_left, texture_left,
                hemi='left', view='lateral',
                threshold=vmax*0.1,
                vmax=vmax,
                title="Left Hemisphere - Lateral",
                axes=ax1,
                cmap='RdBu_r',
                bg_map=fsaverage.sulc_left
            )
            
            # 左半球 - 内侧面
            ax2 = fig.add_subplot(2, 3, 2, projection='3d')
            plotting.plot_surf_stat_map(
                fsaverage.pial_left, texture_left,
                hemi='left', view='medial',
                threshold=vmax*0.1,
                vmax=vmax,
                title="Left Hemisphere - Medial",
                axes=ax2,
                cmap='RdBu_r',
                bg_map=fsaverage.sulc_left
            )
            
            # 左半球 - 背侧面
            ax3 = fig.add_subplot(2, 3, 3, projection='3d')
            plotting.plot_surf_stat_map(
                fsaverage.pial_left, texture_left,
                hemi='left', view='dorsal',
                threshold=vmax*0.1,
                vmax=vmax,
                title="Left Hemisphere - Dorsal",
                axes=ax3,
                cmap='RdBu_r',
                bg_map=fsaverage.sulc_left
            )
            
            # 右半球 - 外侧面
            ax4 = fig.add_subplot(2, 3, 4, projection='3d')
            plotting.plot_surf_stat_map(
                fsaverage.pial_right, texture_right,
                hemi='right', view='lateral',
                threshold=vmax*0.1,
                vmax=vmax,
                title="Right Hemisphere - Lateral",
                axes=ax4,
                cmap='RdBu_r',
                bg_map=fsaverage.sulc_right
            )
            
            # 右半球 - 内侧面
            ax5 = fig.add_subplot(2, 3, 5, projection='3d')
            plotting.plot_surf_stat_map(
                fsaverage.pial_right, texture_right,
                hemi='right', view='medial',
                threshold=vmax*0.1,
                vmax=vmax,
                title="Right Hemisphere - Medial",
                axes=ax5,
                cmap='RdBu_r',
                bg_map=fsaverage.sulc_right
            )
            
            # 右半球 - 腹侧面
            ax6 = fig.add_subplot(2, 3, 6, projection='3d')
            plotting.plot_surf_stat_map(
                fsaverage.pial_right, texture_right,
                hemi='right', view='ventral',
                threshold=vmax*0.1,
                vmax=vmax,
                title="Right Hemisphere - Ventral",
                axes=ax6,
                cmap='RdBu_r',
                bg_map=fsaverage.sulc_right
            )
            
            plt.suptitle(f"{title} - Surface Projections", fontsize=16, fontweight='bold')
            plt.tight_layout()
            
            output_file = os.path.join(self.output_dir, "surface_projections.png")
            plt.savefig(output_file, dpi=300, bbox_inches='tight')
            plt.show()
            print(f"  ✓ Saved surface projections: {output_file}")
            
        except Exception as e:
            print(f"  ! Could not create surface projections: {e}")
    
    def create_interactive_html(self, brain_img: nib.Nifti1Image, title: str):
        """创建交互式HTML视图"""
        try:
            data = brain_img.get_fdata()
            non_zero = data[data != 0]
            if non_zero.size > 0:
                threshold = float(np.percentile(np.abs(non_zero), 70))
            else:
                threshold = 0.0

            # 创建交互式视图
            html_view = plotting.view_img(
                brain_img,
                threshold=threshold,
                cmap='RdBu_r',
                symmetric_cmap=True,
                title=title
            )
            
            # 保存HTML
            html_file = os.path.join(self.output_dir, "interactive_brain_map.html")
            html_view.save_as_html(html_file)
            print(f"  ✓ Saved interactive HTML: {html_file}")
            
            return html_view
            
        except Exception as e:
            print(f"  ! Could not create interactive HTML: {e}")
    
    def create_statistical_report(self, parsed_df: pd.DataFrame):
        """创建统计报告"""
        fig, axes = plt.subplots(3, 3, figsize=(18, 15))
        
        # 1. 脑区分布
        ax = axes[0, 0]
        region_counts = parsed_df[parsed_df['brain_region'] != 'Unknown']['brain_region'].value_counts().head(15)
        ax.barh(range(len(region_counts)), region_counts.values)
        ax.set_yticks(range(len(region_counts)))
        ax.set_yticklabels(region_counts.index, fontsize=8)
        ax.set_xlabel("Count")
        ax.set_title("Top 15 Brain Regions")
        ax.invert_yaxis()
        
        # 2. 结构类型分布
        ax = axes[0, 1]
        struct_counts = parsed_df['structure_type'].value_counts()
        colors = plt.cm.Set3(range(len(struct_counts)))
        ax.pie(struct_counts.values, labels=struct_counts.index, autopct='%1.1f%%', colors=colors)
        ax.set_title("Structure Types Distribution")
        
        # 3. 半球分布
        ax = axes[0, 2]
        hem_counts = parsed_df['hemisphere'].value_counts()
        hem_colors = {'left': '#ff9999', 'right': '#66b3ff', 'bilateral': '#99ff99'}
        colors = [hem_colors.get(h, '#dddddd') for h in hem_counts.index]
        ax.pie(hem_counts.values, labels=hem_counts.index, autopct='%1.1f%%', colors=colors)
        ax.set_title("Hemisphere Distribution")
        
        # 4. 测量类型分布
        ax = axes[1, 0]
        measure_counts = parsed_df[parsed_df['measure_type'] != 'Unknown']['measure_type'].value_counts().head(10)
        ax.bar(range(len(measure_counts)), measure_counts.values)
        ax.set_xticks(range(len(measure_counts)))
        ax.set_xticklabels(measure_counts.index, rotation=45, ha='right', fontsize=8)
        ax.set_ylabel("Count")
        ax.set_title("Top 10 Measurement Types")
        
        # 5. 置信度分布
        ax = axes[1, 1]
        if 'confidence' in parsed_df.columns:
            ax.hist(parsed_df['confidence'], bins=20, edgecolor='black', alpha=0.7)
            ax.set_xlabel("Confidence Score")
            ax.set_ylabel("Count")
            ax.set_title("Parsing Confidence Distribution")
            ax.axvline(parsed_df['confidence'].mean(), color='red', 
                      linestyle='--', label=f'Mean: {parsed_df["confidence"].mean():.2f}')
            ax.legend()
        
        # 6. 脑叶分布
        ax = axes[1, 2]
        lobe_mapping = {
            'Frontal': ['frontal', 'precentral', 'orbitofrontal'],
            'Parietal': ['parietal', 'postcentral', 'precuneus'],
            'Temporal': ['temporal', 'fusiform', 'parahippocampal'],
            'Occipital': ['occipital', 'cuneus', 'lingual', 'calcarine'],
            'Subcortical': ['caudate', 'putamen', 'pallidum', 'thalamus', 'hippocampus', 'amygdala'],
            'Other': []
        }
        
        lobe_counts = {}
        for lobe, keywords in lobe_mapping.items():
            count = 0
            for keyword in keywords:
                count += parsed_df['brain_region'].str.lower().str.contains(keyword, na=False).sum()
            lobe_counts[lobe] = count
        
        lobes = list(lobe_counts.keys())
        counts = list(lobe_counts.values())
        ax.bar(lobes, counts, color=plt.cm.Set2(range(len(lobes))))
        ax.set_xlabel("Brain Lobe")
        ax.set_ylabel("Count")
        ax.set_title("Distribution by Brain Lobe")
        ax.tick_params(axis='x', rotation=45)
        
        # 7. 数据完整性
        ax = axes[2, 0]
        completeness = pd.DataFrame({
            'Identified': [len(parsed_df[parsed_df['brain_region'] != 'Unknown'])],
            'Unknown': [len(parsed_df[parsed_df['brain_region'] == 'Unknown'])]
        })
        completeness.plot(kind='bar', stacked=True, ax=ax, color=['#90EE90', '#FFB6C1'])
        ax.set_title("Data Parsing Completeness")
        ax.set_xticklabels(['IDPs'], rotation=0)
        ax.set_ylabel("Count")
        ax.legend(title="Status")
        
        # 8. 结构类型 x 半球 热图
        ax = axes[2, 1]
        pivot_table = pd.crosstab(parsed_df['structure_type'], parsed_df['hemisphere'])
        sns.heatmap(pivot_table, annot=True, fmt='d', cmap='YlOrRd', ax=ax)
        ax.set_title("Structure Type x Hemisphere")
        
        # 9. 统计摘要
        ax = axes[2, 2]
        ax.axis('off')
        
        total = len(parsed_df)
        identified = len(parsed_df[parsed_df['brain_region'] != 'Unknown'])
        
        stats_text = f"""
PARSING SUMMARY
{'='*40}
Total IDPs: {total}
Successfully mapped: {identified} ({identified/total*100:.1f}%)
Unknown regions: {total - identified} ({(total-identified)/total*100:.1f}%)

Unique brain regions: {parsed_df['brain_region'].nunique()}
Unique structure types: {parsed_df['structure_type'].nunique()}
Unique measure types: {parsed_df['measure_type'].nunique()}

Most common region:
{region_counts.index[0] if len(region_counts) > 0 else 'N/A'}
({region_counts.values[0] if len(region_counts) > 0 else 0} occurrences)

Hemisphere distribution:
Left: {hem_counts.get('left', 0)}
Right: {hem_counts.get('right', 0)}
Bilateral: {hem_counts.get('bilateral', 0)}
        """
        
        ax.text(0.1, 0.5, stats_text, transform=ax.transAxes,
                fontsize=10, verticalalignment='center', family='monospace',
                bbox=dict(boxstyle='round', facecolor='lightgray', alpha=0.5))
        
        plt.suptitle("IDP Brain Mapping Statistical Report", fontsize=16, fontweight='bold')
        plt.tight_layout()
        
        output_file = os.path.join(self.output_dir, "statistical_report.png")
        plt.savefig(output_file, dpi=300, bbox_inches='tight')
        plt.show()
        print(f"  ✓ Saved statistical report: {output_file}")

# ============================================
# 模块 5: 主处理流程
# ============================================

class IDPBrainMappingPipeline:
    """完整的IDP脑图映射流程"""
    
    def __init__(self):
        # 初始化所有组件
        self.dictionary = BrainAtlasDictionary()
        self.parser = IDPParser(self.dictionary)
        self.mapper = PreciseBrainMapper()
        self.visualizer = None

    def _safe_dirname(self, name: str) -> str:
        safe = re.sub(r'[^A-Za-z0-9._-]+', '_', str(name)).strip('_')
        return safe if safe else "unknown_model"

    def process_dataframe(self, df: pd.DataFrame, *, tag: str, output_dir: str):
        os.makedirs(output_dir, exist_ok=True)
        print("=" * 80)
        print(f"MODEL: {tag}")
        print(f"Output dir: {output_dir}")
        print("=" * 80)

        print("\n1. Parsing IDP variables...")
        parsed_df = self.parser.parse_dataframe(df, variable_column='variable_name')

        parsed_file = os.path.join(output_dir, "parsed_idp_results.csv")
        parsed_df.to_csv(parsed_file, index=False)
        print(f"\n2. Saved parsed results to: {parsed_file}")

        print("\n3. Creating brain map...")
        value_col = 'effect_magnitude' if 'effect_magnitude' in parsed_df.columns else None
        brain_img = self.mapper.create_brain_map(parsed_df, value_column=value_col, use_effect_weighting=True)

        brain_file = os.path.join(output_dir, "idp_brain_map.nii.gz")
        nib.save(brain_img, brain_file)
        print(f"   ✓ Saved brain map to: {brain_file}")

        print("\n4. Generating visualizations...")
        visualizer = MultiViewVisualizer(output_dir=output_dir)
        visualizer.create_comprehensive_views(
            brain_img,
            parsed_df,
            title=f"MI-COG2 IDP Brain Mapping - {tag}"
        )

        self._print_final_summary(parsed_df, output_dir=output_dir)
        return parsed_df, brain_img
    
    def process_micog_data(self, file_path: str):
        """处理MI-COG2数据的主函数"""
        print("="*80)
        print("IDP BRAIN MAPPING PIPELINE - MI-COG2 DATA")
        print("="*80)
        print(f"Data file: {file_path}")
        
        # 1. 加载数据
        print("\n1. Loading data...")
        df = pd.read_csv(file_path)
        print(f"   ✓ Loaded {len(df)} records")

        if 'model_comparison' in df.columns:
            model_values = [m for m in df['model_comparison'].dropna().unique().tolist() if str(m).strip()]
            model_values_clean = [m for m in model_values if "no_sex_ethnicity" not in str(m)]
            if len(model_values_clean) > 0:
                model_values = model_values_clean
            if len(model_values) > 0:
                outputs = []
                for model_name in model_values:
                    model_df = df[df['model_comparison'] == model_name].copy()
                    safe_model = self._safe_dirname(model_name)
                    out_dir = os.path.join("brain_visualizations", safe_model)
                    outputs.append(self.process_dataframe(model_df, tag=str(model_name), output_dir=out_dir))

                print("\n" + "=" * 80)
                print("PIPELINE COMPLETED SUCCESSFULLY!")
                print("=" * 80)
                return outputs

        out_dir = "brain_visualizations"
        parsed_df, brain_img = self.process_dataframe(df, tag="all_models", output_dir=out_dir)

        print("\n" + "=" * 80)
        print("PIPELINE COMPLETED SUCCESSFULLY!")
        print("=" * 80)

        return parsed_df, brain_img
    
    def _print_final_summary(self, parsed_df: pd.DataFrame, *, output_dir: str):
        """打印最终摘要"""
        print("\nFINAL SUMMARY:")
        print("-"*40)
        
        total = len(parsed_df)
        identified = len(parsed_df[parsed_df['brain_region'] != 'Unknown'])
        
        print(f"Total IDPs processed: {total}")
        print(f"Successfully mapped: {identified} ({identified/total*100:.1f}%)")
        
        # 按结构类型统计
        print("\nBy Structure Type:")
        for struct_type, count in parsed_df['structure_type'].value_counts().items():
            print(f"  {struct_type}: {count} ({count/total*100:.1f}%)")
        
        # 按测量类型统计
        print("\nTop Measurement Types:")
        measure_counts = parsed_df[parsed_df['measure_type'] != 'Unknown']['measure_type'].value_counts().head(5)
        for measure, count in measure_counts.items():
            print(f"  {measure}: {count}")
        
        print("\nOutput files created:")
        print(f"  - {os.path.join(output_dir, 'parsed_idp_results.csv')}")
        print(f"  - {os.path.join(output_dir, 'idp_brain_map.nii.gz')}")
        print(f"  - {os.path.join(output_dir, 'anatomical_planes.png')}")
        print(f"  - {os.path.join(output_dir, 'glass_brain_3d.png')}")
        print(f"  - {os.path.join(output_dir, 'surface_projections.png')}")
        print(f"  - {os.path.join(output_dir, 'statistical_report.png')}")
        print(f"  - {os.path.join(output_dir, 'interactive_brain_map.html')}")

# ============================================
# 主执行函数
# ============================================

def main():
    """主执行函数"""
    # MI-COG1数据文件路径
    root_dir = os.getcwd()
    candidates = [
        os.path.join(root_dir, "Output_Tables", "Four_Model", "All_Independent_Effect_IDPs_for_Visualization.csv"),
        os.path.join(root_dir, "All_Independent_Effect_IDPs_for_Visualization.csv"),
    ]
    data_file = next((p for p in candidates if os.path.exists(p)), candidates[0])

    # 创建并运行流程
    pipeline = IDPBrainMappingPipeline()
    results = pipeline.process_micog_data(data_file)
    return results

if __name__ == "__main__":
    # 运行主程序
    results = main()
