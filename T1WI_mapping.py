import pandas as pd
import numpy as np
import re
from typing import Dict, List, Tuple, Optional
import warnings
warnings.filterwarnings('ignore')

class UltraEnhancedBrainRegionParser:
    """超级增强版脑区解析器 - 包含更全面的脑区映射"""
    
    def __init__(self):
        # 完整的脑区映射字典（按优先级排序）
        self.brain_mapping_dict = [
            
            # ============= Brodmann区域 (BA) - 最高优先级 =============
            {'pattern': r'BA1\b|BA\.1\b', 'brain_region': 'BA1 (Primary Somatosensory)', 
             'structure_type': 'Cortical', 'atlas_type': 'brodmann', 'ggseg_region': 'postcentral', 'priority': 0.5},
            {'pattern': r'BA2\b|BA\.2\b', 'brain_region': 'BA2 (Primary Somatosensory)', 
             'structure_type': 'Cortical', 'atlas_type': 'brodmann', 'ggseg_region': 'postcentral', 'priority': 0.5},
            {'pattern': r'BA3\b|BA\.3\b', 'brain_region': 'BA3 (Primary Somatosensory)', 
             'structure_type': 'Cortical', 'atlas_type': 'brodmann', 'ggseg_region': 'postcentral', 'priority': 0.5},
            {'pattern': r'BA4a\b|BA\.4a\b', 'brain_region': 'BA4a (Primary Motor - Anterior)', 
             'structure_type': 'Cortical', 'atlas_type': 'brodmann', 'ggseg_region': 'precentral', 'priority': 0.5},
            {'pattern': r'BA4p\b|BA\.4p\b', 'brain_region': 'BA4p (Primary Motor - Posterior)', 
             'structure_type': 'Cortical', 'atlas_type': 'brodmann', 'ggseg_region': 'precentral', 'priority': 0.5},
            {'pattern': r'BA4\b|BA\.4\b', 'brain_region': 'BA4 (Primary Motor)', 
             'structure_type': 'Cortical', 'atlas_type': 'brodmann', 'ggseg_region': 'precentral', 'priority': 0.5},
            {'pattern': r'BA5\b|BA\.5\b', 'brain_region': 'BA5 (Somatosensory Association)', 
             'structure_type': 'Cortical', 'atlas_type': 'brodmann', 'ggseg_region': 'superiorparietal', 'priority': 0.5},
            {'pattern': r'BA6\b|BA\.6\b', 'brain_region': 'BA6 (Premotor/SMA)', 
             'structure_type': 'Cortical', 'atlas_type': 'brodmann', 'ggseg_region': 'superiorfrontal', 'priority': 0.5},
            {'pattern': r'BA7\b|BA\.7\b', 'brain_region': 'BA7 (Somatosensory Association)', 
             'structure_type': 'Cortical', 'atlas_type': 'brodmann', 'ggseg_region': 'superiorparietal', 'priority': 0.5},
            {'pattern': r'BA8\b|BA\.8\b', 'brain_region': 'BA8 (Frontal Eye Fields)', 
             'structure_type': 'Cortical', 'atlas_type': 'brodmann', 'ggseg_region': 'superiorfrontal', 'priority': 0.5},
            {'pattern': r'BA9\b|BA\.9\b', 'brain_region': 'BA9 (Dorsolateral Prefrontal)', 
             'structure_type': 'Cortical', 'atlas_type': 'brodmann', 'ggseg_region': 'superiorfrontal', 'priority': 0.5},
            {'pattern': r'BA10\b|BA\.10\b', 'brain_region': 'BA10 (Anterior Prefrontal)', 
             'structure_type': 'Cortical', 'atlas_type': 'brodmann', 'ggseg_region': 'frontalpole', 'priority': 0.5},
            {'pattern': r'BA11\b|BA\.11\b', 'brain_region': 'BA11 (Orbitofrontal)', 
             'structure_type': 'Cortical', 'atlas_type': 'brodmann', 'ggseg_region': 'lateralorbitofrontal', 'priority': 0.5},
            {'pattern': r'BA17\b|BA\.17\b', 'brain_region': 'BA17 (Primary Visual)', 
             'structure_type': 'Cortical', 'atlas_type': 'brodmann', 'ggseg_region': 'pericalcarine', 'priority': 0.5},
            {'pattern': r'BA18\b|BA\.18\b', 'brain_region': 'BA18 (Secondary Visual)', 
             'structure_type': 'Cortical', 'atlas_type': 'brodmann', 'ggseg_region': 'lateraloccipital', 'priority': 0.5},
            {'pattern': r'BA19\b|BA\.19\b', 'brain_region': 'BA19 (Associative Visual)', 
             'structure_type': 'Cortical', 'atlas_type': 'brodmann', 'ggseg_region': 'lateraloccipital', 'priority': 0.5},
            {'pattern': r'BA20\b|BA\.20\b', 'brain_region': 'BA20 (Inferior Temporal)', 
             'structure_type': 'Cortical', 'atlas_type': 'brodmann', 'ggseg_region': 'inferiortemporal', 'priority': 0.5},
            {'pattern': r'BA21\b|BA\.21\b', 'brain_region': 'BA21 (Middle Temporal)', 
             'structure_type': 'Cortical', 'atlas_type': 'brodmann', 'ggseg_region': 'middletemporal', 'priority': 0.5},
            {'pattern': r'BA22\b|BA\.22\b', 'brain_region': 'BA22 (Superior Temporal)', 
             'structure_type': 'Cortical', 'atlas_type': 'brodmann', 'ggseg_region': 'superiortemporal', 'priority': 0.5},
            {'pattern': r'BA37\b|BA\.37\b', 'brain_region': 'BA37 (Fusiform)', 
             'structure_type': 'Cortical', 'atlas_type': 'brodmann', 'ggseg_region': 'fusiform', 'priority': 0.5},
            {'pattern': r'BA38\b|BA\.38\b', 'brain_region': 'BA38 (Temporal Pole)', 
             'structure_type': 'Cortical', 'atlas_type': 'brodmann', 'ggseg_region': 'temporalpole', 'priority': 0.5},
            {'pattern': r'BA39\b|BA\.39\b', 'brain_region': 'BA39 (Angular)', 
             'structure_type': 'Cortical', 'atlas_type': 'brodmann', 'ggseg_region': 'inferiorparietal', 'priority': 0.5},
            {'pattern': r'BA40\b|BA\.40\b', 'brain_region': 'BA40 (Supramarginal)', 
             'structure_type': 'Cortical', 'atlas_type': 'brodmann', 'ggseg_region': 'supramarginal', 'priority': 0.5},
            {'pattern': r'BA41\b|BA\.41\b', 'brain_region': 'BA41 (Primary Auditory)', 
             'structure_type': 'Cortical', 'atlas_type': 'brodmann', 'ggseg_region': 'transversetemporal', 'priority': 0.5},
            {'pattern': r'BA42\b|BA\.42\b', 'brain_region': 'BA42 (Primary/Secondary Auditory)', 
             'structure_type': 'Cortical', 'atlas_type': 'brodmann', 'ggseg_region': 'superiortemporal', 'priority': 0.5},
            {'pattern': r'BA43\b|BA\.43\b', 'brain_region': 'BA43 (Gustatory)', 
             'structure_type': 'Cortical', 'atlas_type': 'brodmann', 'ggseg_region': 'postcentral', 'priority': 0.5},
            {'pattern': r'BA44\b|BA\.44\b', 'brain_region': 'BA44 (Broca - Pars Opercularis)', 
             'structure_type': 'Cortical', 'atlas_type': 'brodmann', 'ggseg_region': 'parsopercularis', 'priority': 0.5},
            {'pattern': r'BA45\b|BA\.45\b', 'brain_region': 'BA45 (Broca - Pars Triangularis)', 
             'structure_type': 'Cortical', 'atlas_type': 'brodmann', 'ggseg_region': 'parstriangularis', 'priority': 0.5},
            {'pattern': r'BA46\b|BA\.46\b', 'brain_region': 'BA46 (Dorsolateral Prefrontal)', 
             'structure_type': 'Cortical', 'atlas_type': 'brodmann', 'ggseg_region': 'rostralmiddlefrontal', 'priority': 0.5},
            {'pattern': r'BA47\b|BA\.47\b', 'brain_region': 'BA47 (Inferior Frontal)', 
             'structure_type': 'Cortical', 'atlas_type': 'brodmann', 'ggseg_region': 'parsorbitalis', 'priority': 0.5},
            
            # ============= 海马体子区域 =============
            {'pattern': r'CA1\.head|CA1\s+head', 'brain_region': 'CA1 Head', 
             'structure_type': 'Hippocampus', 'atlas_type': 'hippocampus', 'ggseg_region': 'hippocampus', 'priority': 0.5},
            {'pattern': r'CA1\.body|CA1\s+body', 'brain_region': 'CA1 Body', 
             'structure_type': 'Hippocampus', 'atlas_type': 'hippocampus', 'ggseg_region': 'hippocampus', 'priority': 0.5},
            {'pattern': r'CA1\.tail|CA1\s+tail', 'brain_region': 'CA1 Tail', 
             'structure_type': 'Hippocampus', 'atlas_type': 'hippocampus', 'ggseg_region': 'hippocampus', 'priority': 0.5},
            {'pattern': r'CA2\.head|CA2\s+head', 'brain_region': 'CA2 Head', 
             'structure_type': 'Hippocampus', 'atlas_type': 'hippocampus', 'ggseg_region': 'hippocampus', 'priority': 0.5},
            {'pattern': r'CA2\.body|CA2\s+body', 'brain_region': 'CA2 Body', 
             'structure_type': 'Hippocampus', 'atlas_type': 'hippocampus', 'ggseg_region': 'hippocampus', 'priority': 0.5},
            {'pattern': r'CA3\.head|CA3\s+head', 'brain_region': 'CA3 Head', 
             'structure_type': 'Hippocampus', 'atlas_type': 'hippocampus', 'ggseg_region': 'hippocampus', 'priority': 0.5},
            {'pattern': r'CA3\.body|CA3\s+body', 'brain_region': 'CA3 Body', 
             'structure_type': 'Hippocampus', 'atlas_type': 'hippocampus', 'ggseg_region': 'hippocampus', 'priority': 0.5},
            {'pattern': r'CA4\.head|CA4\s+head', 'brain_region': 'CA4 Head', 
             'structure_type': 'Hippocampus', 'atlas_type': 'hippocampus', 'ggseg_region': 'hippocampus', 'priority': 0.5},
            {'pattern': r'CA4\.body|CA4\s+body', 'brain_region': 'CA4 Body', 
             'structure_type': 'Hippocampus', 'atlas_type': 'hippocampus', 'ggseg_region': 'hippocampus', 'priority': 0.5},
            {'pattern': r'subiculum\.head|subiculum\s+head', 'brain_region': 'Subiculum Head', 
             'structure_type': 'Hippocampus', 'atlas_type': 'hippocampus', 'ggseg_region': 'hippocampus', 'priority': 0.5},
            {'pattern': r'subiculum\.body|subiculum\s+body', 'brain_region': 'Subiculum Body', 
             'structure_type': 'Hippocampus', 'atlas_type': 'hippocampus', 'ggseg_region': 'hippocampus', 'priority': 0.5},
            {'pattern': r'presubiculum\.head|presubiculum\s+head', 'brain_region': 'Presubiculum Head', 
             'structure_type': 'Hippocampus', 'atlas_type': 'hippocampus', 'ggseg_region': 'hippocampus', 'priority': 0.5},
            {'pattern': r'presubiculum\.body|presubiculum\s+body', 'brain_region': 'Presubiculum Body', 
             'structure_type': 'Hippocampus', 'atlas_type': 'hippocampus', 'ggseg_region': 'hippocampus', 'priority': 0.5},
            {'pattern': r'whole\.hippocampal\.head|hippocampal\.head', 'brain_region': 'Hippocampal Head', 
             'structure_type': 'Hippocampus', 'atlas_type': 'hippocampus', 'ggseg_region': 'hippocampus', 'priority': 0.5},
            {'pattern': r'whole\.hippocampus|hippocampus$', 'brain_region': 'Hippocampus', 
             'structure_type': 'Subcortical', 'atlas_type': 'subcortical', 'ggseg_region': 'hippocampus', 'priority': 1},
            {'pattern': r'CA1\b', 'brain_region': 'CA1', 
             'structure_type': 'Hippocampus', 'atlas_type': 'hippocampus', 'ggseg_region': 'hippocampus', 'priority': 1},
            {'pattern': r'CA2\b', 'brain_region': 'CA2', 
             'structure_type': 'Hippocampus', 'atlas_type': 'hippocampus', 'ggseg_region': 'hippocampus', 'priority': 1},
            {'pattern': r'CA3\b', 'brain_region': 'CA3', 
             'structure_type': 'Hippocampus', 'atlas_type': 'hippocampus', 'ggseg_region': 'hippocampus', 'priority': 1},
            {'pattern': r'CA4\b|DG\b|dentate', 'brain_region': 'CA4/Dentate Gyrus', 
             'structure_type': 'Hippocampus', 'atlas_type': 'hippocampus', 'ggseg_region': 'hippocampus', 'priority': 1},
            {'pattern': r'subiculum', 'brain_region': 'Subiculum', 
             'structure_type': 'Hippocampus', 'atlas_type': 'hippocampus', 'ggseg_region': 'hippocampus', 'priority': 1},
            {'pattern': r'presubiculum', 'brain_region': 'Presubiculum', 
             'structure_type': 'Hippocampus', 'atlas_type': 'hippocampus', 'ggseg_region': 'hippocampus', 'priority': 1},
            {'pattern': r'perirhinal', 'brain_region': 'Perirhinal Cortex', 
             'structure_type': 'Cortical', 'atlas_type': 'cortical', 'ggseg_region': 'entorhinal', 'priority': 1},
            
            # ============= 杏仁核子区域 =============
            {'pattern': r'lateral\.nucleus|lateral\s+nucleus', 'brain_region': 'Lateral Nucleus (Amygdala)', 
             'structure_type': 'Amygdala', 'atlas_type': 'subcortical', 'ggseg_region': 'amygdala', 'priority': 0.5},
            {'pattern': r'basal\.nucleus|basal\s+nucleus', 'brain_region': 'Basal Nucleus (Amygdala)', 
             'structure_type': 'Amygdala', 'atlas_type': 'subcortical', 'ggseg_region': 'amygdala', 'priority': 0.5},
            {'pattern': r'accessory\.basal\.nucleus', 'brain_region': 'Accessory Basal Nucleus (Amygdala)', 
             'structure_type': 'Amygdala', 'atlas_type': 'subcortical', 'ggseg_region': 'amygdala', 'priority': 0.5},
            {'pattern': r'anterior\.amygdaloid\.area|AAA', 'brain_region': 'Anterior Amygdaloid Area', 
              'structure_type': 'Amygdala', 'atlas_type': 'subcortical', 'ggseg_region': 'amygdala', 'priority': 0.5},
            {'pattern': r'central\.nucleus|\bCeM\b', 'brain_region': 'Central Nucleus (Amygdala)', 
             'structure_type': 'Amygdala', 'atlas_type': 'subcortical', 'ggseg_region': 'amygdala', 'priority': 0.5},
            {'pattern': r'medial\.nucleus|\bMeA\b', 'brain_region': 'Medial Nucleus (Amygdala)', 
             'structure_type': 'Amygdala', 'atlas_type': 'subcortical', 'ggseg_region': 'amygdala', 'priority': 0.5},
            {'pattern': r'cortical\.nucleus|\bCoA\b', 'brain_region': 'Cortical Nucleus (Amygdala)', 
             'structure_type': 'Amygdala', 'atlas_type': 'subcortical', 'ggseg_region': 'amygdala', 'priority': 0.5},
            {'pattern': r'paralaminar\.nucleus', 'brain_region': 'Paralaminar Nucleus (Amygdala)', 
             'structure_type': 'Amygdala', 'atlas_type': 'subcortical', 'ggseg_region': 'amygdala', 'priority': 0.5},
            
            # ============= 丘脑核团 =============
            {'pattern': r'LGN|lateral\.geniculate', 'brain_region': 'Lateral Geniculate Nucleus', 
             'structure_type': 'Thalamus', 'atlas_type': 'subcortical', 'ggseg_region': 'thalamus', 'priority': 0.5},
            {'pattern': r'MGN|medial\.geniculate', 'brain_region': 'Medial Geniculate Nucleus', 
             'structure_type': 'Thalamus', 'atlas_type': 'subcortical', 'ggseg_region': 'thalamus', 'priority': 0.5},
            {'pattern': r'VPL|ventral\.posterior\.lateral', 'brain_region': 'Ventral Posterior Lateral', 
             'structure_type': 'Thalamus', 'atlas_type': 'subcortical', 'ggseg_region': 'thalamus', 'priority': 0.5},
            {'pattern': r'VPM|ventral\.posterior\.medial', 'brain_region': 'Ventral Posterior Medial', 
             'structure_type': 'Thalamus', 'atlas_type': 'subcortical', 'ggseg_region': 'thalamus', 'priority': 0.5},
            {'pattern': r'VA\b|ventral\.anterior', 'brain_region': 'Ventral Anterior', 
             'structure_type': 'Thalamus', 'atlas_type': 'subcortical', 'ggseg_region': 'thalamus', 'priority': 0.5},
            {'pattern': r'VL\b|ventral\.lateral', 'brain_region': 'Ventral Lateral', 
             'structure_type': 'Thalamus', 'atlas_type': 'subcortical', 'ggseg_region': 'thalamus', 'priority': 0.5},
            {'pattern': r'VLa|VL\.a|ventral\.lateral\.anterior', 'brain_region': 'Ventral Lateral Anterior', 
             'structure_type': 'Thalamus', 'atlas_type': 'subcortical', 'ggseg_region': 'thalamus', 'priority': 0.5},
            {'pattern': r'VLp|VL\.p|ventral\.lateral\.posterior', 'brain_region': 'Ventral Lateral Posterior', 
             'structure_type': 'Thalamus', 'atlas_type': 'subcortical', 'ggseg_region': 'thalamus', 'priority': 0.5},
            {'pattern': r'LD\b|lateral\.dorsal', 'brain_region': 'Lateral Dorsal', 
             'structure_type': 'Thalamus', 'atlas_type': 'subcortical', 'ggseg_region': 'thalamus', 'priority': 0.5},
            {'pattern': r'LP\b|lateral\.posterior', 'brain_region': 'Lateral Posterior', 
             'structure_type': 'Thalamus', 'atlas_type': 'subcortical', 'ggseg_region': 'thalamus', 'priority': 0.5},
            {'pattern': r'MD\b|mediodorsal|medial\.dorsal', 'brain_region': 'Mediodorsal', 
             'structure_type': 'Thalamus', 'atlas_type': 'subcortical', 'ggseg_region': 'thalamus', 'priority': 0.5},
            {'pattern': r'MDl|MD\.l|mediodorsal.*lateral', 'brain_region': 'Mediodorsal Lateral', 
             'structure_type': 'Thalamus', 'atlas_type': 'subcortical', 'ggseg_region': 'thalamus', 'priority': 0.5},
            {'pattern': r'VA\.mc|VAmc|ventral\.anterior.*magnocellular', 'brain_region': 'Ventral Anterior Magnocellular', 
             'structure_type': 'Thalamus', 'atlas_type': 'subcortical', 'ggseg_region': 'thalamus', 'priority': 0.5},
            {'pattern': r'CL\b|central\s+lateral', 'brain_region': 'Central Lateral', 
             'structure_type': 'Thalamus', 'atlas_type': 'subcortical', 'ggseg_region': 'thalamus', 'priority': 0.5},
            {'pattern': r'CM\b|centromedian', 'brain_region': 'Centromedian', 
             'structure_type': 'Thalamus', 'atlas_type': 'subcortical', 'ggseg_region': 'thalamus', 'priority': 0.5},
            {'pattern': r'PuA|pulvinar\.anterior', 'brain_region': 'Pulvinar Anterior', 
             'structure_type': 'Thalamus', 'atlas_type': 'subcortical', 'ggseg_region': 'thalamus', 'priority': 0.5},
            {'pattern': r'PuM|pulvinar\.medial', 'brain_region': 'Pulvinar Medial', 
             'structure_type': 'Thalamus', 'atlas_type': 'subcortical', 'ggseg_region': 'thalamus', 'priority': 0.5},
            {'pattern': r'PuL|pulvinar\.lateral', 'brain_region': 'Pulvinar Lateral', 
             'structure_type': 'Thalamus', 'atlas_type': 'subcortical', 'ggseg_region': 'thalamus', 'priority': 0.5},
            {'pattern': r'PuI|pulvinar\.inferior', 'brain_region': 'Pulvinar Inferior', 
             'structure_type': 'Thalamus', 'atlas_type': 'subcortical', 'ggseg_region': 'thalamus', 'priority': 0.5},
            {'pattern': r'AV\b|anteroventral', 'brain_region': 'Anteroventral', 
             'structure_type': 'Thalamus', 'atlas_type': 'subcortical', 'ggseg_region': 'thalamus', 'priority': 0.5},
            {'pattern': r'MV\.Re|MV\b', 'brain_region': 'Medial Ventral/Reuniens', 
             'structure_type': 'Thalamus', 'atlas_type': 'subcortical', 'ggseg_region': 'thalamus', 'priority': 0.5},
            {'pattern': r'L\.Sg|L\-Sg', 'brain_region': 'Limitans/Suprageniculate', 
             'structure_type': 'Thalamus', 'atlas_type': 'subcortical', 'ggseg_region': 'thalamus', 'priority': 0.5},
            
            # ============= 皮质沟回 (Sulci and Gyri) =============
            {'pattern': r'S\.central|s\.central', 'brain_region': 'Central Sulcus', 
             'structure_type': 'Cortical', 'atlas_type': 'cortical', 'ggseg_region': 'precentral', 'priority': 0.8},
            {'pattern': r'G\.rectus|g\.rectus', 'brain_region': 'Gyrus Rectus', 
             'structure_type': 'Cortical', 'atlas_type': 'cortical', 'ggseg_region': 'medialorbitofrontal', 'priority': 0.8},
            {'pattern': r'S\.collat\.transv\.post', 'brain_region': 'Collateral Sulcus Transverse Posterior', 
             'structure_type': 'Cortical', 'atlas_type': 'cortical', 'ggseg_region': 'fusiform', 'priority': 0.8},
            {'pattern': r'S\.collat\.transv\.ant', 'brain_region': 'Collateral Sulcus Transverse Anterior', 
             'structure_type': 'Cortical', 'atlas_type': 'cortical', 'ggseg_region': 'parahippocampal', 'priority': 0.8},
            {'pattern': r'S\.oc\.sup\.transversal', 'brain_region': 'Superior Occipital Transverse Sulcus', 
             'structure_type': 'Cortical', 'atlas_type': 'cortical', 'ggseg_region': 'lateraloccipital', 'priority': 0.8},
            {'pattern': r'S\.orbital\.H\.Shaped|S\.orbital\.H\-Shaped', 'brain_region': 'H-Shaped Orbital Sulcus', 
             'structure_type': 'Cortical', 'atlas_type': 'cortical', 'ggseg_region': 'lateralorbitofrontal', 'priority': 0.8},
            {'pattern': r'Lat\.Fis\.ant\.Vertical|lat\.fis\.ant\.vertical', 'brain_region': 'Anterior Vertical Lateral Fissure', 
             'structure_type': 'Cortical', 'atlas_type': 'cortical', 'ggseg_region': 'parstriangularis', 'priority': 0.8},
            {'pattern': r'Lat\.Fis\.ant\.Horizont|lat\.fis\.ant\.horizont', 'brain_region': 'Anterior Horizontal Lateral Fissure', 
             'structure_type': 'Cortical', 'atlas_type': 'cortical', 'ggseg_region': 'parstriangularis', 'priority': 0.8},
            {'pattern': r'S\.precentral|s\.precentral', 'brain_region': 'Precentral Sulcus', 
             'structure_type': 'Cortical', 'atlas_type': 'cortical', 'ggseg_region': 'precentral', 'priority': 0.8},
            {'pattern': r'S\.postcentral|s\.postcentral', 'brain_region': 'Postcentral Sulcus', 
             'structure_type': 'Cortical', 'atlas_type': 'cortical', 'ggseg_region': 'postcentral', 'priority': 0.8},
            {'pattern': r'S\.frontal\.sup|s\.frontal\.sup', 'brain_region': 'Superior Frontal Sulcus', 
             'structure_type': 'Cortical', 'atlas_type': 'cortical', 'ggseg_region': 'superiorfrontal', 'priority': 0.8},
            {'pattern': r'S\.frontal\.inf|s\.frontal\.inf', 'brain_region': 'Inferior Frontal Sulcus', 
             'structure_type': 'Cortical', 'atlas_type': 'cortical', 'ggseg_region': 'rostralmiddlefrontal', 'priority': 0.8},
            {'pattern': r'S\.temporal\.sup|s\.temporal\.sup', 'brain_region': 'Superior Temporal Sulcus', 
             'structure_type': 'Cortical', 'atlas_type': 'cortical', 'ggseg_region': 'superiortemporal', 'priority': 0.8},
            {'pattern': r'S\.temporal\.inf|s\.temporal\.inf', 'brain_region': 'Inferior Temporal Sulcus', 
             'structure_type': 'Cortical', 'atlas_type': 'cortical', 'ggseg_region': 'inferiortemporal', 'priority': 0.8},
            {'pattern': r'S\.intraparietal|s\.intraparietal', 'brain_region': 'Intraparietal Sulcus', 
             'structure_type': 'Cortical', 'atlas_type': 'cortical', 'ggseg_region': 'superiorparietal', 'priority': 0.8},
            {'pattern': r'S\.calcarine|s\.calcarine', 'brain_region': 'Calcarine Sulcus', 
             'structure_type': 'Cortical', 'atlas_type': 'cortical', 'ggseg_region': 'pericalcarine', 'priority': 0.8},
            {'pattern': r'S\.parieto\.occipital|s\.parieto\.occipital', 'brain_region': 'Parieto-occipital Sulcus', 
             'structure_type': 'Cortical', 'atlas_type': 'cortical', 'ggseg_region': 'cuneus', 'priority': 0.8},
            {'pattern': r'S\.cingulate|s\.cingulate', 'brain_region': 'Cingulate Sulcus', 
             'structure_type': 'Cortical', 'atlas_type': 'cortical', 'ggseg_region': 'rostralanteriorcingulate', 'priority': 0.8},
            
            # ============= 白质纤维束特定位置 =============
            {'pattern': r'posterior\.thalamic\.radiation', 'brain_region': 'Posterior Thalamic Radiation', 
             'structure_type': 'White Matter', 'atlas_type': 'white_matter', 'ggseg_region': '', 'priority': 0.7},
            {'pattern': r'anterior\.thalamic\.radiation', 'brain_region': 'Anterior Thalamic Radiation', 
             'structure_type': 'White Matter', 'atlas_type': 'white_matter', 'ggseg_region': '', 'priority': 0.7},
            {'pattern': r'superior\.thalamic\.radiation', 'brain_region': 'Superior Thalamic Radiation', 
             'structure_type': 'White Matter', 'atlas_type': 'white_matter', 'ggseg_region': '', 'priority': 0.7},
            {'pattern': r'medial\.lemniscus', 'brain_region': 'Medial Lemniscus', 
             'structure_type': 'White Matter', 'atlas_type': 'white_matter', 'ggseg_region': '', 'priority': 0.7},
            {'pattern': r'sagittal\.stratum', 'brain_region': 'Sagittal Stratum', 
             'structure_type': 'White Matter', 'atlas_type': 'white_matter', 'ggseg_region': '', 'priority': 0.7},
            {'pattern': r'tapetum', 'brain_region': 'Tapetum', 
             'structure_type': 'White Matter', 'atlas_type': 'white_matter', 'ggseg_region': '', 'priority': 0.7},
            {'pattern': r'fornix\.cres|fornix\s+cres', 'brain_region': 'Fornix Cres', 
             'structure_type': 'White Matter', 'atlas_type': 'white_matter', 'ggseg_region': '', 'priority': 0.7},
            {'pattern': r'fornix\.stria\.terminalis|stria\.terminalis', 'brain_region': 'Stria Terminalis', 
             'structure_type': 'White Matter', 'atlas_type': 'white_matter', 'ggseg_region': '', 'priority': 0.7},
            {'pattern': r'fornix', 'brain_region': 'Fornix', 
             'structure_type': 'White Matter', 'atlas_type': 'white_matter', 'ggseg_region': '', 'priority': 1},
            {'pattern': r'acoustic\.radiation', 'brain_region': 'Acoustic Radiation', 
             'structure_type': 'White Matter', 'atlas_type': 'white_matter', 'ggseg_region': '', 'priority': 0.7},
            {'pattern': r'optic\.radiation', 'brain_region': 'Optic Radiation', 
             'structure_type': 'White Matter', 'atlas_type': 'white_matter', 'ggseg_region': '', 'priority': 0.7},
            {'pattern': r'cerebellar\.peduncle\.inferior|inferior\.cerebellar\.peduncle|ICP', 'brain_region': 'Inferior Cerebellar Peduncle', 
             'structure_type': 'White Matter', 'atlas_type': 'white_matter', 'ggseg_region': '', 'priority': 0.7},
            {'pattern': r'cerebellar\.peduncle\.middle|middle\.cerebellar\.peduncle|MCP', 'brain_region': 'Middle Cerebellar Peduncle', 
             'structure_type': 'White Matter', 'atlas_type': 'white_matter', 'ggseg_region': '', 'priority': 0.7},
            {'pattern': r'cerebellar\.peduncle\.superior|superior\.cerebellar\.peduncle|SCP', 'brain_region': 'Superior Cerebellar Peduncle', 
             'structure_type': 'White Matter', 'atlas_type': 'white_matter', 'ggseg_region': '', 'priority': 0.7},
            
            # ============= 其他特定区域 =============
            {'pattern': r'supracalcarine\.cortex|supracalcarine', 'brain_region': 'Supracalcarine Cortex', 
             'structure_type': 'Cortical', 'atlas_type': 'cortical', 'ggseg_region': 'cuneus', 'priority': 0.8},
            {'pattern': r'intracalcarine\.cortex|intracalcarine', 'brain_region': 'Intracalcarine Cortex', 
             'structure_type': 'Cortical', 'atlas_type': 'cortical', 'ggseg_region': 'pericalcarine', 'priority': 0.8},
            {'pattern': r'planum\.polare', 'brain_region': 'Planum Polare', 
             'structure_type': 'Cortical', 'atlas_type': 'cortical', 'ggseg_region': 'superiortemporal', 'priority': 0.8},
            {'pattern': r'planum\.temporale', 'brain_region': 'Planum Temporale', 
             'structure_type': 'Cortical', 'atlas_type': 'cortical', 'ggseg_region': 'superiortemporal', 'priority': 0.8},
            {'pattern': r'heschl.*gyrus|transverse\.temporal', 'brain_region': 'Heschls Gyrus', 
             'structure_type': 'Cortical', 'atlas_type': 'cortical', 'ggseg_region': 'transversetemporal', 'priority': 0.8},
            
            # ============= 原有的额叶区域 ===
            {'pattern': r'superiorfrontal|superior\.frontal|front\.sup', 'brain_region': 'Superior Frontal', 
             'structure_type': 'Cortical', 'atlas_type': 'cortical', 'ggseg_region': 'superiorfrontal', 'priority': 1},
            {'pattern': r'rostralmiddlefrontal|rostral\.middle\.frontal', 'brain_region': 'Rostral Middle Frontal', 
             'structure_type': 'Cortical', 'atlas_type': 'cortical', 'ggseg_region': 'rostralmiddlefrontal', 'priority': 1},
            {'pattern': r'caudalmiddlefrontal|caudal\.middle\.frontal', 'brain_region': 'Caudal Middle Frontal', 
             'structure_type': 'Cortical', 'atlas_type': 'cortical', 'ggseg_region': 'caudalmiddlefrontal', 'priority': 1},
            {'pattern': r'middlefrontal|middle\.frontal|front\.middle', 'brain_region': 'Middle Frontal', 
             'structure_type': 'Cortical', 'atlas_type': 'cortical', 'ggseg_region': 'rostralmiddlefrontal', 'priority': 2},
            {'pattern': r'lateralorbitofrontal|lateral\.orbitofrontal', 'brain_region': 'Lateral Orbitofrontal', 
             'structure_type': 'Cortical', 'atlas_type': 'cortical', 'ggseg_region': 'lateralorbitofrontal', 'priority': 1},
            {'pattern': r'medialorbitofrontal|medial\.orbitofrontal', 'brain_region': 'Medial Orbitofrontal', 
             'structure_type': 'Cortical', 'atlas_type': 'cortical', 'ggseg_region': 'medialorbitofrontal', 'priority': 1},
            {'pattern': r'parstriangularis|pars\.triangularis', 'brain_region': 'Pars Triangularis', 
             'structure_type': 'Cortical', 'atlas_type': 'cortical', 'ggseg_region': 'parstriangularis', 'priority': 1},
            {'pattern': r'parsopercularis|pars\.opercularis', 'brain_region': 'Pars Opercularis', 
             'structure_type': 'Cortical', 'atlas_type': 'cortical', 'ggseg_region': 'parsopercularis', 'priority': 1},
            {'pattern': r'parsorbitalis|pars\.orbitalis', 'brain_region': 'Pars Orbitalis', 
             'structure_type': 'Cortical', 'atlas_type': 'cortical', 'ggseg_region': 'parsorbitalis', 'priority': 1},
            {'pattern': r'precentral|g\.precentral', 'brain_region': 'Precentral', 
             'structure_type': 'Cortical', 'atlas_type': 'cortical', 'ggseg_region': 'precentral', 'priority': 1},
            {'pattern': r'paracentral|g\.s\.paracentral', 'brain_region': 'Paracentral', 
             'structure_type': 'Cortical', 'atlas_type': 'cortical', 'ggseg_region': 'paracentral', 'priority': 1},
            {'pattern': r'frontalpole|frontal\.pole', 'brain_region': 'Frontal Pole', 
             'structure_type': 'Cortical', 'atlas_type': 'cortical', 'ggseg_region': 'frontalpole', 'priority': 1},
            
            # === 顶叶区域 ===
            {'pattern': r'superiorparietal|superior\.parietal', 'brain_region': 'Superior Parietal', 
             'structure_type': 'Cortical', 'atlas_type': 'cortical', 'ggseg_region': 'superiorparietal', 'priority': 1},
            {'pattern': r'inferiorparietal|inferior\.parietal', 'brain_region': 'Inferior Parietal', 
             'structure_type': 'Cortical', 'atlas_type': 'cortical', 'ggseg_region': 'inferiorparietal', 'priority': 1},
            {'pattern': r'postcentral|g\.postcentral', 'brain_region': 'Postcentral', 
             'structure_type': 'Cortical', 'atlas_type': 'cortical', 'ggseg_region': 'postcentral', 'priority': 1},
            {'pattern': r'precuneus|g\.precuneus', 'brain_region': 'Precuneus', 
             'structure_type': 'Cortical', 'atlas_type': 'cortical', 'ggseg_region': 'precuneus', 'priority': 1},
            {'pattern': r'supramarginal', 'brain_region': 'Supramarginal', 
             'structure_type': 'Cortical', 'atlas_type': 'cortical', 'ggseg_region': 'supramarginal', 'priority': 1},
            {'pattern': r'angular', 'brain_region': 'Angular', 
             'structure_type': 'Cortical', 'atlas_type': 'cortical', 'ggseg_region': 'inferiorparietal', 'priority': 1},
            
            # === 颞叶区域 ===
            {'pattern': r'superiortemporal|superior\.temporal', 'brain_region': 'Superior Temporal', 
             'structure_type': 'Cortical', 'atlas_type': 'cortical', 'ggseg_region': 'superiortemporal', 'priority': 1},
            {'pattern': r'middletemporal|middle\.temporal', 'brain_region': 'Middle Temporal', 
             'structure_type': 'Cortical', 'atlas_type': 'cortical', 'ggseg_region': 'middletemporal', 'priority': 1},
            {'pattern': r'inferiortemporal|inferior\.temporal', 'brain_region': 'Inferior Temporal', 
             'structure_type': 'Cortical', 'atlas_type': 'cortical', 'ggseg_region': 'inferiortemporal', 'priority': 1},
            {'pattern': r'transversetemporal|transverse\.temporal', 'brain_region': 'Transverse Temporal', 
             'structure_type': 'Cortical', 'atlas_type': 'cortical', 'ggseg_region': 'transversetemporal', 'priority': 1},
            {'pattern': r'temporalpole|temporal\.pole|pole\.temporal', 'brain_region': 'Temporal Pole', 
             'structure_type': 'Cortical', 'atlas_type': 'cortical', 'ggseg_region': 'temporalpole', 'priority': 1},
            {'pattern': r'parahippocampal', 'brain_region': 'Parahippocampal', 
             'structure_type': 'Cortical', 'atlas_type': 'cortical', 'ggseg_region': 'parahippocampal', 'priority': 1},
            {'pattern': r'entorhinal', 'brain_region': 'Entorhinal', 
             'structure_type': 'Cortical', 'atlas_type': 'cortical', 'ggseg_region': 'entorhinal', 'priority': 1},
            {'pattern': r'fusiform', 'brain_region': 'Fusiform', 
             'structure_type': 'Cortical', 'atlas_type': 'cortical', 'ggseg_region': 'fusiform', 'priority': 1},
            {'pattern': r'bankssts|banks\.sts', 'brain_region': 'Banks STS', 
             'structure_type': 'Cortical', 'atlas_type': 'cortical', 'ggseg_region': 'bankssts', 'priority': 1},
            
            # === 枕叶区域 ===
            {'pattern': r'lateraloccipital|lateral\.occipital', 'brain_region': 'Lateral Occipital', 
             'structure_type': 'Cortical', 'atlas_type': 'cortical', 'ggseg_region': 'lateraloccipital', 'priority': 1},
            {'pattern': r'pericalcarine', 'brain_region': 'Pericalcarine', 
             'structure_type': 'Cortical', 'atlas_type': 'cortical', 'ggseg_region': 'pericalcarine', 'priority': 1},
            {'pattern': r'cuneus|g\.cuneus', 'brain_region': 'Cuneus', 
             'structure_type': 'Cortical', 'atlas_type': 'cortical', 'ggseg_region': 'cuneus', 'priority': 1},
            {'pattern': r'lingual', 'brain_region': 'Lingual', 
             'structure_type': 'Cortical', 'atlas_type': 'cortical', 'ggseg_region': 'lingual', 'priority': 1},
            
            # === 扣带回 ===
            {'pattern': r'rostralanteriorcingulate', 'brain_region': 'Rostral Anterior Cingulate', 
             'structure_type': 'Cortical', 'atlas_type': 'cortical', 'ggseg_region': 'rostralanteriorcingulate', 'priority': 1},
            {'pattern': r'caudalanteriorcingulate', 'brain_region': 'Caudal Anterior Cingulate', 
             'structure_type': 'Cortical', 'atlas_type': 'cortical', 'ggseg_region': 'caudalanteriorcingulate', 'priority': 1},
            {'pattern': r'posteriorcingulate', 'brain_region': 'Posterior Cingulate', 
             'structure_type': 'Cortical', 'atlas_type': 'cortical', 'ggseg_region': 'posteriorcingulate', 'priority': 1},
            {'pattern': r'isthmuscingulate', 'brain_region': 'Isthmus Cingulate', 
             'structure_type': 'Cortical', 'atlas_type': 'cortical', 'ggseg_region': 'isthmuscingulate', 'priority': 1},
            
            # === 脑岛 ===
            {'pattern': r'insula', 'brain_region': 'Insula', 
             'structure_type': 'Cortical', 'atlas_type': 'cortical', 'ggseg_region': 'insula', 'priority': 1},
            
            # === 皮下核团 ===
            {'pattern': r'caudate', 'brain_region': 'Caudate', 
             'structure_type': 'Subcortical', 'atlas_type': 'subcortical', 'ggseg_region': 'caudate', 'priority': 1},
            {'pattern': r'putamen', 'brain_region': 'Putamen', 
             'structure_type': 'Subcortical', 'atlas_type': 'subcortical', 'ggseg_region': 'putamen', 'priority': 1},
            {'pattern': r'pallidum', 'brain_region': 'Pallidum', 
             'structure_type': 'Subcortical', 'atlas_type': 'subcortical', 'ggseg_region': 'pallidum', 'priority': 1},
            {'pattern': r'accumbens', 'brain_region': 'Nucleus Accumbens', 
             'structure_type': 'Subcortical', 'atlas_type': 'subcortical', 'ggseg_region': 'accumbens', 'priority': 1},
            
            # === 海马体和杏仁核 ===
            {'pattern': r'whole\.amygdala|amygdala$', 'brain_region': 'Amygdala', 
             'structure_type': 'Subcortical', 'atlas_type': 'subcortical', 'ggseg_region': 'amygdala', 'priority': 1},
            {'pattern': r'HATA', 'brain_region': 'Hippocampal Amygdaloid Transition Area', 
             'structure_type': 'Subcortical', 'atlas_type': 'subcortical', 'ggseg_region': 'amygdala', 'priority': 1},
            {'pattern': r'Corticoamygdaloid.*transition|Cortico\-amygdaloid.*transition|Corticoamygdaloid\.transitio', 'brain_region': 'Corticoamygdaloid Transition', 
             'structure_type': 'Subcortical', 'atlas_type': 'subcortical', 'ggseg_region': 'amygdala', 'priority': 1},
            {'pattern': r'GC[\-_]?ML[\-_]?DG', 'brain_region': 'Dentate Gyrus (GC-ML-DG)', 
             'structure_type': 'Subcortical', 'atlas_type': 'subcortical', 'ggseg_region': 'hippocampus', 'priority': 1},
            {'pattern': r'CA1|CA2|CA3|CA4|subiculum|presubiculum|hippocampal\.head|hippocampal\.body', 'brain_region': 'Hippocampal Subfields', 
             'structure_type': 'Subcortical', 'atlas_type': 'subcortical', 'ggseg_region': 'hippocampus', 'priority': 2},
            
            # === 丘脑 ===
            {'pattern': r'whole\.thalamus|thalamus$', 'brain_region': 'Thalamus', 
             'structure_type': 'Subcortical', 'atlas_type': 'subcortical', 'ggseg_region': 'thalamus', 'priority': 1},
            
            # === 小脑 ===
            {'pattern': r'cerebellum|cerebell', 'brain_region': 'Cerebellum', 
             'structure_type': 'Cerebellum', 'atlas_type': 'cerebellum', 'ggseg_region': '', 'priority': 2},
            
            # === 脑干 ===
            {'pattern': r'whole\.brainstem|brain\.stem|brainstem', 'brain_region': 'Brainstem', 
             'structure_type': 'Brainstem', 'atlas_type': 'brainstem', 'ggseg_region': '', 'priority': 1},
            {'pattern': r'midbrain', 'brain_region': 'Midbrain', 
             'structure_type': 'Brainstem', 'atlas_type': 'brainstem', 'ggseg_region': '', 'priority': 1},
            {'pattern': r'pons', 'brain_region': 'Pons', 
             'structure_type': 'Brainstem', 'atlas_type': 'brainstem', 'ggseg_region': '', 'priority': 1},
            {'pattern': r'medulla', 'brain_region': 'Medulla', 
             'structure_type': 'Brainstem', 'atlas_type': 'brainstem', 'ggseg_region': '', 'priority': 1},
            
            # === 白质纤维束 ===
            {'pattern': r'corpus\.callosum', 'brain_region': 'Corpus Callosum', 
             'structure_type': 'White Matter', 'atlas_type': 'white_matter', 'ggseg_region': '', 'priority': 2},
            {'pattern': r'superior\.longitudinal\.fasciculus', 'brain_region': 'Superior Longitudinal Fasciculus', 
             'structure_type': 'White Matter', 'atlas_type': 'white_matter', 'ggseg_region': '', 'priority': 1},
            {'pattern': r'inferior\.longitudinal\.fasciculus', 'brain_region': 'Inferior Longitudinal Fasciculus', 
             'structure_type': 'White Matter', 'atlas_type': 'white_matter', 'ggseg_region': '', 'priority': 1},
            {'pattern': r'uncinate\.fasciculus', 'brain_region': 'Uncinate Fasciculus', 
             'structure_type': 'White Matter', 'atlas_type': 'white_matter', 'ggseg_region': '', 'priority': 1},
            {'pattern': r'cingulum', 'brain_region': 'Cingulum', 
             'structure_type': 'White Matter', 'atlas_type': 'white_matter', 'ggseg_region': '', 'priority': 2},
            {'pattern': r'internal\.capsule', 'brain_region': 'Internal Capsule', 
             'structure_type': 'White Matter', 'atlas_type': 'white_matter', 'ggseg_region': '', 'priority': 2},
            {'pattern': r'corona\.radiata', 'brain_region': 'Corona Radiata', 
             'structure_type': 'White Matter', 'atlas_type': 'white_matter', 'ggseg_region': '', 'priority': 2},
            {'pattern': r'tract\.', 'brain_region': 'White Matter Tract', 
             'structure_type': 'White Matter', 'atlas_type': 'white_matter', 'ggseg_region': '', 'priority': 5},
            
            # === 通用匹配模式（低优先级） ===
            {'pattern': r'frontal|front', 'brain_region': 'Frontal', 
             'structure_type': 'Cortical', 'atlas_type': 'cortical', 'ggseg_region': 'superiorfrontal', 'priority': 3},
            {'pattern': r'parietal|pariet', 'brain_region': 'Parietal', 
             'structure_type': 'Cortical', 'atlas_type': 'cortical', 'ggseg_region': 'superiorparietal', 'priority': 3},
            {'pattern': r'temporal|temp', 'brain_region': 'Temporal', 
             'structure_type': 'Cortical', 'atlas_type': 'cortical', 'ggseg_region': 'middletemporal', 'priority': 3},
            {'pattern': r'occipital|occip', 'brain_region': 'Occipital', 
             'structure_type': 'Cortical', 'atlas_type': 'cortical', 'ggseg_region': 'lateraloccipital', 'priority': 3},
            {'pattern': r'cingulat|cingul', 'brain_region': 'Cingulate', 
             'structure_type': 'Cortical', 'atlas_type': 'cortical', 'ggseg_region': 'rostralanteriorcingulate', 'priority': 3},
            
            # === 额外的通用模式 ===
            {'pattern': r'totalsurface', 'brain_region': 'Total Surface', 
             'structure_type': 'Measurement', 'atlas_type': 'measurement', 'ggseg_region': '', 'priority': 4},
            {'pattern': r'whole\.brain|brain\.grey\.white', 'brain_region': 'Whole Brain', 
             'structure_type': 'Global', 'atlas_type': 'global', 'ggseg_region': '', 'priority': 4},
            {'pattern': r'white\.matter\b|white\s+matter', 'brain_region': 'White Matter', 
             'structure_type': 'White Matter', 'atlas_type': 'white_matter', 'ggseg_region': '', 'priority': 3},
            {'pattern': r'grey\.matter|gray\.matter', 'brain_region': 'Grey Matter', 
             'structure_type': 'Grey Matter', 'atlas_type': 'grey_matter', 'ggseg_region': '', 'priority': 3},
            {'pattern': r'ventricular|ventricle', 'brain_region': 'Ventricles', 
             'structure_type': 'CSF', 'atlas_type': 'csf', 'ggseg_region': '', 'priority': 2},
            {'pattern': r'hyperintensit', 'brain_region': 'White Matter Hyperintensities', 
             'structure_type': 'Pathology', 'atlas_type': 'pathology', 'ggseg_region': '', 'priority': 2},
        ]
        
        # 测量类型字典（扩展）
        self.measure_type_dict = [
            {'pattern': r'^area\.|area\.of', 'measure_type': 'Surface Area'},
            {'pattern': r'^volume\.|volume\.of', 'measure_type': 'Volume'},
            {'pattern': r'mean\.thickness|thickness\.of', 'measure_type': 'Cortical Thickness'},
            {'pattern': r'grey\.white\.contrast', 'measure_type': 'Grey-White Contrast'},
            {'pattern': r'weighted\.mean\.(fa|md|l1|l2|l3|isovf|icvf|od|mo)', 'measure_type': 'DTI Tract Weighted'},
            {'pattern': r'mean\.(fa|md|l1|l2|l3|isovf|icvf|od|mo).*on\.fa\.skeleton', 'measure_type': 'DTI Skeleton'},
            {'pattern': r'mean\.fa\b', 'measure_type': 'FA (Fractional Anisotropy)'},
            {'pattern': r'mean\.md\b', 'measure_type': 'MD (Mean Diffusivity)'},
            {'pattern': r'mean\.l1\b', 'measure_type': 'L1 (Axial Diffusivity)'},
            {'pattern': r'mean\.l2\b', 'measure_type': 'L2 (Radial Diffusivity 1)'},
            {'pattern': r'mean\.l3\b', 'measure_type': 'L3 (Radial Diffusivity 2)'},
            {'pattern': r'mean\.icvf\b', 'measure_type': 'ICVF (Intracellular Volume Fraction)'},
            {'pattern': r'mean\.isovf\b', 'measure_type': 'ISOVF (Isotropic Volume Fraction)'},
            {'pattern': r'mean\.od\b', 'measure_type': 'OD (Orientation Dispersion)'},
            {'pattern': r'mean\.mo\b', 'measure_type': 'MO (Mode of Anisotropy)'},
            {'pattern': r'median\.bold|90th\.percentile.*z\.statistic|median\.z\.statistic', 'measure_type': 'fMRI Activation'},
            {'pattern': r'head\.motion|tfmri\.head\.motion', 'measure_type': 'Head Motion'},
            {'pattern': r'percentile|bold|motion', 'measure_type': 'fMRI/Motion'},
            {'pattern': r'number\.of.*matches', 'measure_type': 'Behavioral'},
                        {'pattern': r'normalised\.for\.head\.size', 'measure_type': 'Normalized Volume'},
            {'pattern': r'fast\.(?:grey|gray)\.matter', 'measure_type': 'Grey Matter Segmentation'},
            {'pattern': r'fast\.white\.matter', 'measure_type': 'White Matter Segmentation'},
            {'pattern': r'fast\.csf', 'measure_type': 'CSF Segmentation'},
            {'pattern': r'first\.all\.fast\.volumes', 'measure_type': 'FIRST Subcortical Volumes'},
            {'pattern': r'aseg\.stats', 'measure_type': 'FreeSurfer Segmentation'},
            {'pattern': r'aparc\.stats', 'measure_type': 'FreeSurfer Parcellation'},
            {'pattern': r'connectivity', 'measure_type': 'Connectivity'},
            {'pattern': r'network', 'measure_type': 'Network Measure'},
            {'pattern': r'centrality', 'measure_type': 'Centrality'},
            {'pattern': r'modularity', 'measure_type': 'Modularity'},
            {'pattern': r'path\.length', 'measure_type': 'Path Length'},
            {'pattern': r'clustering', 'measure_type': 'Clustering Coefficient'},
            {'pattern': r'betweenness', 'measure_type': 'Betweenness Centrality'},
            {'pattern': r'degree', 'measure_type': 'Degree'},
            {'pattern': r'efficiency', 'measure_type': 'Efficiency'},
            {'pattern': r'surface\.area', 'measure_type': 'Surface Area'},
            {'pattern': r'curvature', 'measure_type': 'Curvature'},
            {'pattern': r'folding\.index', 'measure_type': 'Folding Index'},
            {'pattern': r'gyrification', 'measure_type': 'Gyrification Index'},
            {'pattern': r'sulcal\.depth', 'measure_type': 'Sulcal Depth'},
            {'pattern': r'cortical\.complexity', 'measure_type': 'Cortical Complexity'},
            {'pattern': r'jacobian', 'measure_type': 'Jacobian Determinant'},
            {'pattern': r'tensor', 'measure_type': 'Tensor'},
            {'pattern': r'eigenvalue', 'measure_type': 'Eigenvalue'},
            {'pattern': r'eigenvector', 'measure_type': 'Eigenvector'},
            {'pattern': r'trace', 'measure_type': 'Trace'},
            {'pattern': r'mode', 'measure_type': 'Mode'},
            {'pattern': r'kurtosis', 'measure_type': 'Kurtosis'},
            {'pattern': r'skewness', 'measure_type': 'Skewness'},
        ]
        
        # 编译所有正则表达式
        for item in self.brain_mapping_dict:
            item['compiled_pattern'] = re.compile(item['pattern'], re.IGNORECASE)
        
        for item in self.measure_type_dict:
            item['compiled_pattern'] = re.compile(item['pattern'], re.IGNORECASE)
    
    def get_hemisphere(self, variable_name):
        """提取半球信息，增强模式识别"""
        # 清理变量名
        clean_name = variable_name.lower()
        
        # 定义半球模式
        hemisphere_patterns = {
            'left': [
                r'\.l\.|\.l$|^l\.',  # .l. 或 .l 结尾或 l. 开头
                r'\.lh\.|\.lh$|^lh\.',  # .lh. 或 .lh 结尾或 lh. 开头
                r'_l_|_l$|^l_',  # _l_ 或 _l 结尾或 l_ 开头
                r'_lh_|_lh$|^lh_',  # _lh_ 或 _lh 结尾或 lh_ 开头
                r'left',  # 明确的 left
                r'\bl\b',  # 独立的 l
                r'\blh\b'  # 独立的 lh
            ],
            'right': [
                r'\.r\.|\.r$|^r\.',  # .r. 或 .r 结尾或 r. 开头
                r'\.rh\.|\.rh$|^rh\.',  # .rh. 或 .rh 结尾或 rh. 开头
                r'_r_|_r$|^r_',  # _r_ 或 _r 结尾或 r_ 开头
                r'_rh_|_rh$|^rh_',  # _rh_ 或 _rh 结尾或 rh_ 开头
                r'right',  # 明确的 right
                r'\br\b',  # 独立的 r
                r'\brh\b'  # 独立的 rh
            ]
        }
        
        # 检查左半球
        for pattern in hemisphere_patterns['left']:
            if re.search(pattern, clean_name):
                return 'left'
        
        # 检查右半球
        for pattern in hemisphere_patterns['right']:
            if re.search(pattern, clean_name):
                return 'right'
        
        # 检查是否为双侧结构
        bilateral_keywords = ['bilateral', 'both', 'total', 'whole', 'mean', 'average']
        for keyword in bilateral_keywords:
            if keyword in clean_name:
                return 'bilateral'
        
        # 某些结构默认为双侧
        midline_structures = ['brainstem', 'corpus.callosum', 'ventricle', 'thalamus', 
                            'hypothalamus', 'pituitary', 'pineal', 'fornix', 
                            'septum', 'vermis', 'aqueduct']
        for structure in midline_structures:
            if structure in clean_name:
                return 'bilateral'
        
        return None
    
    def parse_brain_region(self, variable_name):
        """
        解析脑区信息，返回匹配的区域及其属性
        """
        clean_name = variable_name.lower()
        clean_name = re.sub(r'^(?:weighted\.)?mean\.', '', clean_name)
        clean_name = re.sub(r'\.mean\.', '.', clean_name)
        matches = []
        
        # 遍历所有脑区模式
        for region_info in self.brain_mapping_dict:
            if region_info['compiled_pattern'].search(clean_name):
                match = {
                    'brain_region': region_info['brain_region'],
                    'structure_type': region_info['structure_type'],
                    'atlas_type': region_info['atlas_type'],
                    'ggseg_region': region_info['ggseg_region'],
                    'priority': region_info['priority']
                }
                matches.append(match)
        
        # 如果有多个匹配，选择优先级最高的（数值最小）
        if matches:
            matches.sort(key=lambda x: x['priority'])
            return matches[0]
        
        return None
    
    def parse_measure_type(self, variable_name):
        """解析测量类型"""
        clean_name = variable_name.lower()
        
        for measure_info in self.measure_type_dict:
            if measure_info['compiled_pattern'].search(clean_name):
                return measure_info['measure_type']
        
        # 默认测量类型推断
        if 'volume' in clean_name:
            return 'Volume'
        elif 'area' in clean_name:
            return 'Surface Area'
        elif 'thickness' in clean_name:
            return 'Cortical Thickness'
        elif 'fa' in clean_name:
            return 'FA (Fractional Anisotropy)'
        elif 'md' in clean_name:
            return 'MD (Mean Diffusivity)'
        elif 'connectivity' in clean_name:
            return 'Connectivity'
        
        return None
    
    def parse_modality(self, variable_name):
        """推断成像模态"""
        clean_name = variable_name.lower()
        
        modality_patterns = {
            'T1-weighted MRI': ['t1', 'mprage', 'spgr', 'structural', 'anatomical', 
                              'volume', 'thickness', 'area', 'morphometry'],
            'DTI/DWI': ['dti', 'dwi', 'diffusion', 'fa', 'md', 'ad', 'rd', 'tract', 
                       'fiber', 'fibre', 'tensor', 'kurtosis', 'noddi'],
            'fMRI': ['fmri', 'bold', 'tfmri', 'rfmri', 'resting', 'task', 
                    'activation', 'connectivity', 'network'],
            'PET': ['pet', 'fdg', 'amyloid', 'tau', 'dopamine', 'suvr', 'binding'],
            'CT': ['ct', 'hounsfield', 'attenuation'],
            'MEG': ['meg', 'magnetoencephalography', 'magnetic'],
            'EEG': ['eeg', 'electroencephalography', 'electric', 'erp'],
            'SPECT': ['spect', 'perfusion', 'cbf'],
            'MRS': ['mrs', 'spectroscopy', 'metabolite', 'naa', 'cho', 'cr'],
            'ASL': ['asl', 'arterial', 'perfusion', 'cbf'],
            'SWI': ['swi', 'susceptibility', 'microbleed'],
            'QSM': ['qsm', 'quantitative.susceptibility'],
            'T2-weighted MRI': ['t2', 'flair', 'hyperintensit'],
            'Multi-modal': ['multimodal', 'combined', 'fusion']
        }
        
        # 检查每种模态
        for modality, keywords in modality_patterns.items():
            for keyword in keywords:
                if keyword in clean_name:
                    return modality
        
        # 基于测量类型推断
        measure_type = self.parse_measure_type(variable_name)
        if measure_type:
            if any(keyword in measure_type.lower() for keyword in ['volume', 'area', 'thickness']):
                return 'T1-weighted MRI'
            elif any(keyword in measure_type.lower() for keyword in ['fa', 'md', 'diffusivity', 'tract']):
                return 'DTI/DWI'
            elif any(keyword in measure_type.lower() for keyword in ['activation', 'fmri', 'bold']):
                return 'fMRI'
        
        return 'Unknown'
    
    def parse_network_measure(self, variable_name):
        """解析网络测量相关信息"""
        clean_name = variable_name.lower()
        
        network_measures = {
            'degree': ['degree', 'deg'],
            'betweenness': ['betweenness', 'between'],
            'closeness': ['closeness', 'close'],
            'eigenvector': ['eigenvector', 'eigen'],
            'pagerank': ['pagerank', 'page.rank'],
            'clustering': ['clustering', 'cluster.coefficient'],
            'modularity': ['modularity', 'module'],
            'efficiency': ['efficiency', 'eff'],
            'path_length': ['path.length', 'shortest.path', 'distance'],
            'small_world': ['small.world', 'smallworld'],
            'assortativity': ['assortativity', 'assort'],
            'rich_club': ['rich.club', 'richclub'],
            'centrality': ['centrality', 'central'],
            'hub': ['hub', 'hubs'],
            'community': ['community', 'communities'],
            'participation': ['participation', 'participation.coefficient'],
            'diversity': ['diversity', 'shannon', 'entropy']
        }
        
        for measure, keywords in network_measures.items():
            for keyword in keywords:
                if keyword in clean_name:
                    return measure
        
        return None
    
    def parse_behavioral_domain(self, variable_name):
        """解析行为/认知域信息"""
        clean_name = variable_name.lower()
        
        behavioral_domains = {
            'memory': ['memory', 'recall', 'recognition', 'episodic', 'semantic', 
                      'working.memory', 'wm', 'encoding', 'retrieval'],
            'attention': ['attention', 'attentional', 'vigilance', 'sustained', 
                         'selective', 'divided', 'focus'],
            'executive': ['executive', 'inhibition', 'switching', 'flexibility', 
                         'planning', 'decision', 'control', 'stroop', 'flanker'],
            'language': ['language', 'verbal', 'semantic', 'phonological', 
                        'naming', 'fluency', 'comprehension', 'reading'],
            'visuospatial': ['visuospatial', 'visual', 'spatial', 'navigation', 
                           'rotation', 'imagery', 'perception'],
            'motor': ['motor', 'movement', 'coordination', 'dexterity', 
                     'tapping', 'grasping', 'reaching'],
            'emotion': ['emotion', 'affect', 'mood', 'fear', 'anxiety', 
                       'depression', 'happiness', 'sadness', 'anger'],
            'social': ['social', 'theory.of.mind', 'tom', 'empathy', 
                      'mentalizing', 'face', 'interpersonal'],
            'reward': ['reward', 'punishment', 'reinforcement', 'value', 
                      'gambling', 'risk', 'delay.discounting'],
            'sensory': ['sensory', 'somatosensory', 'tactile', 'pain', 
                       'temperature', 'proprioception', 'auditory', 'olfactory'],
            'processing_speed': ['speed', 'reaction.time', 'rt', 'latency', 
                               'response.time', 'processing'],
            'intelligence': ['intelligence', 'iq', 'reasoning', 'problem.solving', 
                           'fluid', 'crystallized', 'g.factor']
        }
        
        for domain, keywords in behavioral_domains.items():
            for keyword in keywords:
                if keyword in clean_name:
                    return domain
        
        return None
    
    def parse_statistical_measure(self, variable_name):
        """解析统计测量类型"""
        clean_name = variable_name.lower()
        
        stat_patterns = {
            'mean': r'\bmean\b|average',
            'median': r'\bmedian\b',
            'std': r'\bstd\b|standard.deviation',
            'variance': r'\bvar\b|variance',
            'min': r'\bmin\b|minimum',
            'max': r'\bmax\b|maximum',
            'sum': r'\bsum\b|total',
            'count': r'\bcount\b|number',
            'percentile': r'percentile|pct|quantile',
            'zscore': r'z.score|z\.statistic|zscore',
            'tscore': r't.score|t\.statistic|tscore',
            'pvalue': r'p.value|pval|p\.val',
            'effect_size': r'effect.size|cohens.d|eta',
            'correlation': r'corr|correlation|r\.value',
            'beta': r'\bbeta\b|coefficient',
            'se': r'\bse\b|standard.error',
            'ci': r'\bci\b|confidence.interval'
        }
        
        for stat_type, pattern in stat_patterns.items():
            if re.search(pattern, clean_name):
                return stat_type
        
        return None
    
    def parse_time_point(self, variable_name):
        """解析时间点信息"""
        clean_name = variable_name.lower()
        
        # 时间点模式
        time_patterns = [
            (r'baseline|bl\b|t0|time0|wave0|w0', 'baseline'),
            (r't1\b|time1|wave1|w1|visit1|v1', 'time1'),
            (r't2\b|time2|wave2|w2|visit2|v2', 'time2'),
            (r't3\b|time3|wave3|w3|visit3|v3', 'time3'),
            (r'followup|follow.up|fu', 'followup'),
            (r'pre\b|before', 'pre'),
            (r'post\b|after', 'post'),
            (r'year(\d+)|y(\d+)', 'year'),
            (r'month(\d+)|m(\d+)', 'month'),
            (r'week(\d+)|wk(\d+)', 'week'),
            (r'day(\d+)|d(\d+)', 'day')
        ]
        
        for pattern, time_type in time_patterns:
            match = re.search(pattern, clean_name)
            if match:
                if time_type in ['year', 'month', 'week', 'day'] and match.groups():
                    return f"{time_type}_{match.group(1)}"
                return time_type
        
        return None
    
    def parse_age_range(self, variable_name):
        """解析年龄范围信息"""
        clean_name = variable_name.lower()
        
        age_patterns = {
            'pediatric': ['pediatric', 'child', 'infant', 'toddler', 'adolescent', 'teen'],
            'young_adult': ['young.adult', 'ya', '18.25', '20.30'],
            'middle_aged': ['middle.age', 'ma', '40.60', '45.65'],
            'elderly': ['elderly', 'older', 'senior', 'geriatric', '65plus', '70plus'],
            'lifespan': ['lifespan', 'all.ages', 'across.ages']
        }
        
        for age_range, keywords in age_patterns.items():
            for keyword in keywords:
                if keyword in clean_name:
                    return age_range
        
        # 查找具体年龄数字
        age_match = re.search(r'age(\d+)', clean_name)
        if age_match:
            age = int(age_match.group(1))
            if age < 18:
                return 'pediatric'
            elif age < 30:
                return 'young_adult'
            elif age < 60:
                return 'middle_aged'
            else:
                return 'elderly'
        
        return None
    
    def parse_full(self, variable_name):
        """
        完整解析变量名，返回所有解析出的信息
        
        Parameters:
        -----------
        variable_name : str
            要解析的变量名
        
        Returns:
        --------
        dict : 包含所有解析信息的字典
        """
        result = {
            'original_name': variable_name,
            'hemisphere': self.get_hemisphere(variable_name),
            'brain_region': None,
            'structure_type': None,
            'atlas_type': None,
            'ggseg_region': None,
            'measure_type': self.parse_measure_type(variable_name),
            'modality': self.parse_modality(variable_name),
            'network_measure': self.parse_network_measure(variable_name),
            'behavioral_domain': self.parse_behavioral_domain(variable_name),
            'statistical_measure': self.parse_statistical_measure(variable_name),
            'time_point': self.parse_time_point(variable_name),
            'age_range': self.parse_age_range(variable_name),
            'description': None
        }
        
        # 解析脑区信息
        region_info = self.parse_brain_region(variable_name)
        if region_info:
            result['brain_region'] = region_info['brain_region']
            result['structure_type'] = region_info['structure_type']
            result['atlas_type'] = region_info['atlas_type']
            result['ggseg_region'] = region_info['ggseg_region']
        
        # 生成描述
        result['description'] = self._generate_description(result)
        
        return result

    def parse_variable_name_enhanced(self, variable_name):
        result = {
            'brain_region': 'Unknown',
            'structure_type': 'Unknown',
            'atlas_type': 'Unknown',
            'ggseg_region': 'Unknown',
            'lobe': 'Unknown',
            'hemisphere': 'Bilateral',
            'measure_type': 'Unknown',
            'tissue_type': 'Unknown',
            'match_confidence': 0.0
        }

        if pd.isna(variable_name):
            return result

        var_lower = str(variable_name).lower().strip()

        hemisphere = self.get_hemisphere(variable_name)
        if hemisphere == 'left':
            result['hemisphere'] = 'Left'
        elif hemisphere == 'right':
            result['hemisphere'] = 'Right'
        elif hemisphere == 'bilateral':
            result['hemisphere'] = 'Bilateral'

        measure_type = self.parse_measure_type(variable_name)
        if measure_type:
            result['measure_type'] = measure_type
        else:
            if 'volume' in var_lower:
                result['measure_type'] = 'Volume'
            elif 'area' in var_lower:
                result['measure_type'] = 'Surface Area'
            elif 'thickness' in var_lower:
                result['measure_type'] = 'Cortical Thickness'

        if any(k in var_lower for k in ['grey matter', 'gray matter', ' fast.grey', ' fast.gray', 'gm ']):
            result['tissue_type'] = 'Grey Matter'
        elif any(k in var_lower for k in ['white matter', ' fast.white', 'wm ']):
            result['tissue_type'] = 'White Matter'
        elif any(k in var_lower for k in ['csf', 'cerebrospinal fluid', ' fast.csf']):
            result['tissue_type'] = 'CSF'

        region_info = self.parse_brain_region(variable_name)
        if region_info:
            result['brain_region'] = region_info.get('brain_region', 'Unknown')
            result['structure_type'] = region_info.get('structure_type', 'Unknown')
            result['atlas_type'] = region_info.get('atlas_type', 'Unknown')
            ggseg_region = region_info.get('ggseg_region', None)
            if ggseg_region:
                result['ggseg_region'] = ggseg_region

            priority = region_info.get('priority', None)
            if isinstance(priority, (int, float)) and priority > 0:
                result['match_confidence'] = min(1.0, 1.0 / float(priority))
            else:
                result['match_confidence'] = 0.5

            if result['structure_type'] == 'Cortical':
                region_lower = str(result['brain_region']).lower()
                if 'frontal' in region_lower:
                    result['lobe'] = 'Frontal'
                elif 'parietal' in region_lower:
                    result['lobe'] = 'Parietal'
                elif 'temporal' in region_lower:
                    result['lobe'] = 'Temporal'
                elif 'occipital' in region_lower:
                    result['lobe'] = 'Occipital'
                elif 'cingulate' in region_lower:
                    result['lobe'] = 'Cingulate'
                elif 'insula' in region_lower:
                    result['lobe'] = 'Insula'
            else:
                result['lobe'] = result['structure_type']

        if 'fusiform' in var_lower:
            if 'temporal' in var_lower:
                result['brain_region'] = 'Temporal Fusiform'
                result['lobe'] = 'Temporal'
            elif 'occipital' in var_lower:
                result['brain_region'] = 'Occipital Fusiform'
                result['lobe'] = 'Occipital'

        if 'operculum' in var_lower or 'opercular' in var_lower:
            if 'frontal' in var_lower:
                result['brain_region'] = 'Frontal Operculum'
                result['lobe'] = 'Frontal'
            elif 'parietal' in var_lower:
                result['brain_region'] = 'Parietal Operculum'
                result['lobe'] = 'Parietal'
            elif 'central' in var_lower:
                result['brain_region'] = 'Central Operculum'
                result['lobe'] = 'Central'

        abbreviation_mapping = {
            'sma': 'Supplementary Motor Area',
            'dlpfc': 'Dorsolateral Prefrontal Cortex',
            'vmpfc': 'Ventromedial Prefrontal Cortex',
            'acc': 'Anterior Cingulate Cortex',
            'pcc': 'Posterior Cingulate Cortex',
            'sts': 'Superior Temporal Sulcus',
            'ips': 'Intraparietal Sulcus',
            'tpj': 'Temporoparietal Junction'
        }

        for abbrev, full_name in abbreviation_mapping.items():
            if abbrev in var_lower.split():
                result['brain_region'] = full_name
                full_lower = full_name.lower()
                if 'frontal' in full_lower:
                    result['lobe'] = 'Frontal'
                elif 'temporal' in full_lower:
                    result['lobe'] = 'Temporal'
                elif 'parietal' in full_lower:
                    result['lobe'] = 'Parietal'
                elif 'cingulate' in full_lower:
                    result['lobe'] = 'Cingulate'
                break

        return result
    
    def _generate_description(self, parsed_info):
        """
        根据解析的信息生成人类可读的描述
        """
        parts = []
        
        # 添加测量类型
        if parsed_info['measure_type']:
            parts.append(parsed_info['measure_type'])
        
        # 添加半球信息
        if parsed_info['hemisphere']:
            if parsed_info['hemisphere'] == 'left':
                parts.append('Left')
            elif parsed_info['hemisphere'] == 'right':
                parts.append('Right')
            elif parsed_info['hemisphere'] == 'bilateral':
                parts.append('Bilateral')
        
        # 添加脑区
        if parsed_info['brain_region']:
            parts.append(parsed_info['brain_region'])
        
        # 添加模态
        if parsed_info['modality'] and parsed_info['modality'] != 'Unknown':
            parts.append(f"({parsed_info['modality']})")
        
        # 添加网络测量
        if parsed_info['network_measure']:
            parts.append(f"Network: {parsed_info['network_measure']}")
        
        # 添加行为域
        if parsed_info['behavioral_domain']:
            parts.append(f"Domain: {parsed_info['behavioral_domain']}")
        
        # 添加时间点
        if parsed_info['time_point']:
            parts.append(f"Time: {parsed_info['time_point']}")
        
        # 添加年龄范围
        if parsed_info['age_range']:
            parts.append(f"Age: {parsed_info['age_range']}")
        
        # 添加统计测量
        if parsed_info['statistical_measure']:
            parts.append(f"Stat: {parsed_info['statistical_measure']}")
        
        return ' '.join(parts) if parts else 'Unknown variable'
    
    def batch_parse(self, variable_names):
        """
        批量解析变量名列表
        
        Parameters:
        -----------
        variable_names : list
            变量名列表
        
        Returns:
        --------
        list : 解析结果列表
        """
        results = []
        for var_name in variable_names:
            results.append(self.parse_full(var_name))
        return results
    
    def to_dataframe(self, parsed_results):
        """
        将解析结果转换为DataFrame
        
        Parameters:
        -----------
        parsed_results : list
            parse_full 或 batch_parse 的结果
        
        Returns:
        --------
        pd.DataFrame : 解析结果的DataFrame
        """
        import pandas as pd
        
        if isinstance(parsed_results, dict):
            parsed_results = [parsed_results]
        
        return pd.DataFrame(parsed_results)
    
    def export_mapping(self, parsed_results, output_file):
        """
        导出解析结果到文件
        
        Parameters:
        -----------
        parsed_results : list
            解析结果
        output_file : str
            输出文件路径（支持.csv, .xlsx, .json）
        """
        import pandas as pd
        import json
        
        if output_file.endswith('.csv'):
            df = self.to_dataframe(parsed_results)
            df.to_csv(output_file, index=False)
        elif output_file.endswith('.xlsx'):
            df = self.to_dataframe(parsed_results)
            df.to_excel(output_file, index=False)
        elif output_file.endswith('.json'):
            with open(output_file, 'w') as f:
                json.dump(parsed_results, f, indent=2)
        else:
            raise ValueError(f"Unsupported file format: {output_file}")
    
    def generate_summary_stats(self, parsed_results):
        """
        生成解析结果的统计摘要
        
        Parameters:
        -----------
        parsed_results : list
            解析结果
        
        Returns:
        --------
        dict : 统计摘要
        """
        import pandas as pd
        from collections import Counter
        
        df = self.to_dataframe(parsed_results)
        
        summary = {
            'total_variables': len(df),
            'hemisphere_distribution': dict(df['hemisphere'].value_counts()),
            'structure_type_distribution': dict(df['structure_type'].value_counts()),
            'modality_distribution': dict(df['modality'].value_counts()),
            'measure_type_distribution': dict(df['measure_type'].value_counts()),
            'unique_brain_regions': df['brain_region'].nunique(),
            'top_brain_regions': dict(df['brain_region'].value_counts().head(10)),
            'has_time_info': (df['time_point'].notna()).sum(),
            'has_behavioral_info': (df['behavioral_domain'].notna()).sum(),
            'has_network_info': (df['network_measure'].notna()).sum()
        }
        
        return summary
        
        if pd.isna(variable_name):
            return result
        
        # 预处理：转换为小写，保留原始用于其他匹配
        var_lower = str(variable_name).lower().strip()
        var_original = str(variable_name).strip()
        
        # 1. 提取半球信息
        hemisphere_patterns = {
            'Left': [r'$left$', r'left hemisphere', r'_lh_', r'\.lh\.', r'\bleft\b'],
            'Right': [r'$right$', r'right hemisphere', r'_rh_', r'\.rh\.', r'\bright\b'],
            'Bilateral': [r'bilateral', r'\bboth\b']
        }
        
        for hemi, patterns in hemisphere_patterns.items():
            for pattern in patterns:
                if re.search(pattern, var_lower):
                    result['hemisphere'] = hemi
                    # 清理半球信息
                    var_lower = re.sub(pattern, '', var_lower).strip()
                    break
            if result['hemisphere'] != 'Bilateral':
                break
        
        # 2. 提取测量类型
        for measure_dict in self.measure_type_dict:
            if re.search(measure_dict['pattern'], var_lower):
                result['measure_type'] = measure_dict['measure_type']
                break
        
        # 额外的测量类型检测（基于常见词汇）
        if result['measure_type'] == 'Unknown':
            measure_keywords = {
                'volume': 'Volume',
                'area': 'Surface Area',
                'thickness': 'Cortical Thickness',
                'intensity': 'Mean Intensity',
                'fa': 'FA (Fractional Anisotropy)',
                'md': 'MD (Mean Diffusivity)',
                'icvf': 'ICVF',
                'isovf': 'ISOVF',
                'od': 'OD (Orientation Dispersion)'
            }
            for keyword, measure_type in measure_keywords.items():
                if keyword in var_lower:
                    result['measure_type'] = measure_type
                    break
        
        # 3. 提取组织类型
        tissue_patterns = {
            'Grey Matter': [r'grey matter', r'gray matter', r'\bgm\b', r'grey', r'gray'],
            'White Matter': [r'white matter', r'\bwm\b', r'white'],
            'CSF': [r'csf', r'cerebrospinal fluid']
        }
        
        for tissue, patterns in tissue_patterns.items():
            for pattern in patterns:
                if re.search(pattern, var_lower):
                    result['tissue_type'] = tissue
                    break
            if result['tissue_type'] != 'Unknown':
                break
        
        # 4. 脑区匹配（使用优先级系统）
        best_match = None
        best_priority = float('inf')
        
        for mapping in self.brain_mapping_dict:
            if re.search(mapping['pattern'], var_lower):
                if mapping['priority'] < best_priority:
                    best_match = mapping
                    best_priority = mapping['priority']
                    # 如果找到优先级为1的匹配，直接使用
                    if best_priority == 1:
                        break
        
        if best_match:
            result['brain_region'] = best_match['brain_region']
            result['structure_type'] = best_match['structure_type']
            result['atlas_type'] = best_match['atlas_type']
            result['ggseg_region'] = best_match['ggseg_region']
            result['match_confidence'] = 1.0 / best_priority  # 优先级越低，置信度越高
            
            # 从映射推断脑叶
            if best_match['brain_region'] in self.lobe_mapping:
                result['lobe'] = self.lobe_mapping[best_match['brain_region']]
            elif best_match['structure_type'] == 'Cortical':
                # 尝试从脑区名称推断
                if 'frontal' in best_match['brain_region'].lower():
                    result['lobe'] = 'Frontal'
                elif 'parietal' in best_match['brain_region'].lower():
                    result['lobe'] = 'Parietal'
                elif 'temporal' in best_match['brain_region'].lower():
                    result['lobe'] = 'Temporal'
                elif 'occipital' in best_match['brain_region'].lower():
                    result['lobe'] = 'Occipital'
                elif 'cingulate' in best_match['brain_region'].lower():
                    result['lobe'] = 'Cingulate'
                elif 'insula' in best_match['brain_region'].lower():
                    result['lobe'] = 'Insula'
            else:
                result['lobe'] = best_match['structure_type']
        
        # 5. 特殊情况处理
        # 处理包含多个脑区关键词的情况
        if 'fusiform' in var_lower:
            if 'temporal' in var_lower:
                result['brain_region'] = 'Temporal Fusiform'
                result['lobe'] = 'Temporal'
            elif 'occipital' in var_lower:
                result['brain_region'] = 'Occipital Fusiform'
                result['lobe'] = 'Occipital'
        
        # 处理operculum相关
        if 'operculum' in var_lower or 'opercular' in var_lower:
            if 'frontal' in var_lower:
                result['brain_region'] = 'Frontal Operculum'
                result['lobe'] = 'Frontal'
            elif 'parietal' in var_lower:
                result['brain_region'] = 'Parietal Operculum'
                result['lobe'] = 'Parietal'
            elif 'central' in var_lower:
                result['brain_region'] = 'Central Operculum'
                result['lobe'] = 'Central'
        
        # 处理特殊的缩写
        abbreviation_mapping = {
            'sma': 'Supplementary Motor Area',
            'dlpfc': 'Dorsolateral Prefrontal Cortex',
            'vmpfc': 'Ventromedial Prefrontal Cortex',
            'acc': 'Anterior Cingulate Cortex',
            'pcc': 'Posterior Cingulate Cortex',
            'sts': 'Superior Temporal Sulcus',
            'ips': 'Intraparietal Sulcus',
            'tpj': 'Temporoparietal Junction'
        }
        
        for abbrev, full_name in abbreviation_mapping.items():
            if abbrev in var_lower.split():
                result['brain_region'] = full_name
                if 'frontal' in full_name.lower():
                    result['lobe'] = 'Frontal'
                elif 'temporal' in full_name.lower():
                    result['lobe'] = 'Temporal'
                elif 'parietal' in full_name.lower():
                    result['lobe'] = 'Parietal'
                elif 'cingulate' in full_name.lower():
                    result['lobe'] = 'Cingulate'
                break
        
        return result

EnhancedBrainRegionParser = UltraEnhancedBrainRegionParser

def process_brain_data_with_enhanced_parser(file_path: str, save_outputs: bool = True) -> pd.DataFrame:
    """
    使用增强解析器处理脑区数据
    """
    
    print("\n" + "="*80)
    print(" ENHANCED BRAIN REGION PARSING PIPELINE ")
    print("="*80)
    
    # 读取数据
    print("\n[Step 1] Loading data...")
    df = pd.read_csv(file_path)
    print(f"✓ Loaded: {df.shape[0]} rows, {df.shape[1]} columns")
    
    # 初始化增强解析器
    print("\n[Step 2] Initializing enhanced parser...")
    parser = UltraEnhancedBrainRegionParser()
    print(f"✓ Loaded {len(parser.brain_mapping_dict)} brain region patterns")
    print(f"✓ Loaded {len(parser.measure_type_dict)} measurement type patterns")
    
    # 解析所有记录
    print("\n[Step 3] Parsing variable names...")
    parsed_records = []
    
    for idx, row in df.iterrows():
        # 解析
        parsed = parser.parse_variable_name_enhanced(row['variable_name'])
        
        # 创建完整记录
        record = {
            'index': idx,
            'variable_name': row['variable_name'],
            'brain_region': parsed['brain_region'],
            'structure_type': parsed['structure_type'],
            'atlas_type': parsed['atlas_type'],
            'ggseg_region': parsed['ggseg_region'],
            'lobe': parsed['lobe'],
            'hemisphere': parsed['hemisphere'],
            'measure_type': parsed['measure_type'],
            'tissue_type': parsed['tissue_type'],
            'match_confidence': parsed['match_confidence']
        }
        
        # 添加其他列
        for col in df.columns:
            if col != 'variable_name' and col not in record:
                record[col] = row[col]
        
        parsed_records.append(record)
        
        # 进度提示
        if (idx + 1) % 100 == 0:
            progress = (idx + 1) / len(df) * 100
            print(f"  Progress: {idx + 1}/{len(df)} ({progress:.1f}%)")
    
    # 创建结果DataFrame
    result_df = pd.DataFrame(parsed_records)
    
    # 统计分析
    print("\n[Step 4] Analysis Results:")
    print("-" * 60)
    
    total = len(result_df)
    identified = len(result_df[result_df['brain_region'] != 'Unknown'])
    unknown = total - identified
    
    print(f"Total records: {total}")
    print(f"Successfully identified: {identified} ({identified/total*100:.1f}%)")
    print(f"Unknown regions: {unknown} ({unknown/total*100:.1f}%)")
    
    # 按结构类型统计
    print("\nBy Structure Type:")
    structure_counts = result_df['structure_type'].value_counts()
    for struct, count in structure_counts.items():
        print(f"  {struct:20s}: {count:5d} ({count/total*100:5.1f}%)")
    
    # 按脑叶统计
    print("\nBy Brain Lobe:")
    lobe_counts = result_df[result_df['lobe'] != 'Unknown']['lobe'].value_counts()
    for lobe, count in lobe_counts.items():
        print(f"  {lobe:20s}: {count:5d} ({count/total*100:5.1f}%)")
    
    # 按置信度统计
    print("\nBy Match Confidence:")
    high_conf = len(result_df[result_df['match_confidence'] >= 0.5])
    med_conf = len(result_df[(result_df['match_confidence'] > 0) & (result_df['match_confidence'] < 0.5)])
    no_match = len(result_df[result_df['match_confidence'] == 0])
    print(f"  High confidence (≥0.5): {high_conf} ({high_conf/total*100:.1f}%)")
    print(f"  Medium confidence (<0.5): {med_conf} ({med_conf/total*100:.1f}%)")
    print(f"  No match: {no_match} ({no_match/total*100:.1f}%)")
    
    # Top识别的脑区
    print("\nTop 15 Identified Brain Regions:")
    top_regions = result_df[result_df['brain_region'] != 'Unknown']['brain_region'].value_counts().head(15)
    for region, count in top_regions.items():
        print(f"  {region:35s}: {count:4d}")
    
    # 保存结果
    if save_outputs:
        # 保存完整解析结果
        output_file = file_path.replace('.csv', '_enhanced_parsed.csv')
        result_df.to_csv(output_file, index=False)
        print(f"\n✓ Saved parsed results: {output_file}")
        
        # 保存未识别记录
        unknown_df = result_df[result_df['brain_region'] == 'Unknown']
        if len(unknown_df) > 0:
            unknown_file = file_path.replace('.csv', '_unknown_regions.csv')
            unknown_df[['variable_name', 'brain_region', 'measure_type']].to_csv(unknown_file, index=False)
            print(f"✓ Saved unknown regions: {unknown_file}")
        
        # 保存汇总统计
        summary_file = file_path.replace('.csv', '_summary_stats.txt')
        with open(summary_file, 'w') as f:
            f.write("BRAIN REGION PARSING SUMMARY\n")
            f.write("="*60 + "\n\n")
            f.write(f"Total records: {total}\n")
            f.write(f"Identified: {identified} ({identified/total*100:.1f}%)\n")
            f.write(f"Unknown: {unknown} ({unknown/total*100:.1f}%)\n\n")
            
            f.write("Structure Type Distribution:\n")
            for struct, count in structure_counts.items():
                f.write(f"  {struct}: {count} ({count/total*100:.1f}%)\n")
            
            f.write("\nTop Brain Regions:\n")
            for region, count in top_regions.items():
                f.write(f"  {region}: {count}\n")
        
        print(f"✓ Saved summary statistics: {summary_file}")
    
    return result_df

# 使用示例
if __name__ == "__main__":
    # 处理您的数据文件
    import os
    root_dir = os.getcwd()
    candidates = [
        os.path.join(root_dir, "Output_Tables", "Four_Model", "All_Independent_Effect_IDPs_for_Visualization.csv"),
        os.path.join(root_dir, "All_Independent_Effect_IDPs_for_Visualization.csv"),
    ]
    csv_file = next((p for p in candidates if os.path.exists(p)), candidates[0])

    # 使用增强解析器处理
    enhanced_df = process_brain_data_with_enhanced_parser(csv_file, save_outputs=True)

    # 显示一些解析示例
    print("\n" + "="*80)
    print("Sample Parsed Results (First 20):")
    print("="*80)

    display_cols = ['variable_name', 'brain_region', 'lobe',
                    'hemisphere', 'measure_type', 'match_confidence']

    for idx, row in enhanced_df[display_cols].head(20).iterrows():
        print(f"\n[{idx+1}]")
        print(f"  Original: {row['variable_name'][:60]}...")
        print(f"  Region: {row['brain_region']}")
        print(f"  Lobe: {row['lobe']}")
        print(f"  Hemisphere: {row['hemisphere']}")
        print(f"  Measure: {row['measure_type']}")
        print(f"  Confidence: {row['match_confidence']:.2f}")
        print(f"  Measure: {row['measure_type']}")
        print(f"  Confidence: {row['match_confidence']:.2f}")
