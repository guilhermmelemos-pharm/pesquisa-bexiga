# Lemos Lambda v2.1

MODELOS = {
    "flash": "gemini-2.5-flash",
    "pro":   "gemini-2.5-pro"
}

# Padrões de lixo genômico — bloqueados ANTES de qualquer análise
LIXO_REGEX = [
    r"^LINC\d+", r"^LNC-", r"^MIR\d+", r"^MIRNA", r"^hsa-mir",
    r"^SNHG\d+", r"^NEAT\d*$", r"\w+-AS\d*$", r"^RPL\d+P\d*",
    r"^RPS\d+P\d*", r"^GAPDHP\d*", r"^C\dorf\d+", r"^C\d+orf\d+",
    r"^AATBC$", r"^CSTP\d*", r"^LOC\d+", r"^ENSG\d+",
    r"^AC\d{6}\.\d+", r"^AL\d{6}\.\d+", r"^LINC", r"^lnc",
]

BLACKLIST = set([
    "DEMONSTRATE","ELUCIDATE","FACILITATE","ATTENUATE","INVESTIGATE","INDICATE",
    "DETERMINE","EVALUATE","ASSOCIATED","RELEASE","SPECULATE","SUGGEST","OBSERVED",
    "RESULTED","COMPARED","REVEALED","ANALYZED","REPORTED","EXPRESSION","ACTIVITY",
    "FUNCTION","LEVEL","MECHANISM","PATHWAY","POTENTIAL","THERAPEUTIC","TARGET",
    "EFFECT","ROLE","ACTION","STUDY","APPROPRIATE","ALONGSIDE","IMMEDIATE","GUIDELINE",
    "WORLDWIDE","NATIONWIDE","OVERALL","SIGNIFICANT","IMPORTANT","NOVEL","CURRENT",
    "SPECIFIC","PREVIOUS","SIMILAR","ADDITIONAL","VARIOUS","MULTIPLE","MAJOR",
    "RECENT","INITIAL","FINAL","GENERAL","PRIMARY","SECONDARY","TERTIARY",
    "SUBSTRATE","CANDIDATE","PROTOCOL","INTERMEDIATE","CONSISTENT","ESSENTIAL",
    "INTESTINE","LIVER","KIDNEY","BLADDER","HEART","BRAIN","TISSUE","ORGAN",
    "MUSCLE","BONE","SKIN","LUNG","COLON","STOMACH","PANCREAS","SPLEEN",
    "THYROID","ADRENAL","OVARY","UTERUS","TESTIS","PROSTATE","VESSEL","ARTERY",
    "PATIENT","CLINICAL","HOSPITAL","SURGERY","THERAPY","TREATMENT","DISEASE",
    "SYNDROME","DISORDER","SYMPTOM","DIAGNOSIS","PROGNOSIS","OUTCOME","SURVIVAL",
    "MORTALITY","MORBIDITY","PREVALENCE","INCIDENCE","COHORT","TRIAL","PLACEBO",
    "PROTEIN","GENE","RECEPTOR","CHANNEL","INHIBITOR","AGONIST","ANTAGONIST",
    "DRUG","COMPOUND","CELLS","CULTURE","MEDIA","BUFFER","ASSAY","WESTERN",
    "ELISA","FACS","DMSO","DMEM","RPMI","PBS","SALINE","ETHANOL",
    "RESPONSE","REGULATION","SIGNALING","MEDIATED","INDUCED","HYPOTHESIS",
    "CHROMATIN","GENOME","TRANSCRIPTION","TRANSLATION","REPLICATION",
    "MICE","RAT","HUMAN","ANIMAL","MODEL","MOUSE","RABBIT","MONKEY","ZEBRAFISH",
    "GLUCOSE","INSULIN","TESTOSTERONE","ESTRADIOL","CORTISOL","ALDOSTERONE",
    "CHOLESTEROL","ATP","ADP","AMP","ADENOSINE","OXYGEN","ROS","NITRIC",
    "SODIUM","POTASSIUM","CALCIUM","CHLORIDE","WATER","UREA","CREATININE",
    "ALBUMIN","LACTATE","PYRUVATE","ACETATE","CITRATE","SUCCINATE",
    "HCT116","HEK293","CHO","HELA","MCF7","PC12","RAW264","JURKAT","U937",
    "A549","MDCK","VERO","COS7","NIH3T3","C2C12","L929","THP1",
    "AAV","AAV9","LENTIVIRUS","ADENOVIRUS","RETROVIRUS","PLASMID",
    "REVIEW","DATA","RESULT","CONCLUSION","INTRODUCTION","METHOD","MATERIAL",
    "FIGURE","TABLE","SUPPLEMENT","APPENDIX","ABSTRACT","KEYWORD",
])

REGEX_PATTERNS = [
    r"\b[A-Z]{2,6}\d{1,4}[A-Z]?\b",
    r"\b[A-Za-z]{3,}\d{1,3}\b",
    r"\b[A-Za-z]{6,}(?:mab|ib|ol|one|ide|ate|ase|pril|afil|tin|ine|arin|oside)\b"
]

TIPOS_LABEL = {
    "Receptor":             "🔵 Receptor",
    "Drug":                 "💊 Droga",
    "Modulator":            "🟡 Modulador",
    "Gene":                 "🧬 Gene/TF",
    "Ion Channel":          "⚡ Canal Iônico",
    "Enzyme":               "⚙️ Enzima",
    "Pathway":              "🔀 Via",
    "Biomarker":            "📊 Biomarcador",
    "Target (OpenTargets)": "🎯 OpenTargets",
    "Molecule":             "🔬 Molécula",
}

PRESETS = {
    "🧪 Chaperonas & Estresse":    ["Trehalose","TMAO","4-PBA","Taurine","Betaine","HSP70","HSP90"],
    "🔥 Gasotransmissores":        ["HO-1","CORM-2","CORM-401","GYY4137","AP39","CBS","CSE"],
    "⚡ Mecanossensores & Canais": ["Piezo1","Piezo2","TRPA1","TRPM8","TRPV4","TREK-1","TMEM16A","P2X3","P2X7"],
    "🧬 Epigenética & RNAs":       ["TET2","HDAC6","DNMT1","METTL3","YTHDF2","FTO"],
    "🥑 Resolução & Lipídios":     ["Resolvin D1","Resolvin E1","Maresin-1","Lipoxin A4","GPR55","GPR18","GPX4"],
    "👻 Receptores Órfãos":        ["GPR35","GPR84","GPR68","GPR183","TAS2R14","TAS2R38","CRTH2"],
    "🎯 Alvos Emergentes":         ["NLRP3","ACSL4","HMGB1","DDR1","HAVCR2","STING","cGAS"],
    "🔄 Autofagia & Metabolismo":  ["AMPK","ULK1","Beclin1","Rapamycin","Metformin","Spermidine"],
    "🩸 Vasculatura & Fibrose":    ["VEGFR2","Nintedanib","Pirfenidone","LOX","CTGF","Relaxin"],
    "🦠 Microbioma & Metabólitos": ["Butyrate","Propionate","Urolithin A","Indole","SCFA","Trimethylamine"],
}

SISTEMAS = {
    "Selecionar Sistema...":  "",
    "Urinary System":         "Bladder or Kidney or Ureter or Urothelium or Detrusor",
    "Nervous System":         "Brain or Spinal Cord or Neurons or Glia or Astrocytes",
    "Cardiovascular System":  "Heart or Endothelium or Artery or Vein or Myocardium",
    "Digestive System":       "Stomach or Intestine or Liver or Pancreas or Colon",
    "Respiratory System":     "Lung or Bronchi or Alveoli or Trachea or Airway",
    "Endocrine System":       "Thyroid or Adrenal or Pituitary or Islets",
    "Immune System":          "Lymphocytes or Macrophages or Spleen or Inflammation",
    "Musculoskeletal System": "Skeletal Muscle or Bone or Cartilage or Synovium",
    "Male Reproductive":      "Testis or Epididymis or Prostate or Seminal",
    "Female Reproductive":    "Ovary or Uterus or Endometrium or Fallopian",
    "Skin":                   "Dermis or Epidermis or Keratinocytes or Fibroblasts",
}

TAGS = {
    "blue_ocean": {"label": "💎 Blue Ocean",  "color": "#00d4ff"},
    "embryonic":  {"label": "🌱 Embrionário", "color": "#7fff00"},
    "gold":       {"label": "🥇 Ouro",        "color": "#ffd700"},
    "trending":   {"label": "🚀 Tendência",   "color": "#ff8c00"},
    "ghost":      {"label": "👻 Fantasma",    "color": "#888888"},
    "saturated":  {"label": "🔴 Saturado",    "color": "#ff4444"},
}

# Mapeamento órgão → tecido HPA
HPA_TISSUE_MAP = {
    "bladder": "urinary bladder", "urothelium": "urinary bladder",
    "detrusor": "urinary bladder", "overactive bladder": "urinary bladder",
    "heart": "heart muscle", "myocardium": "heart muscle", "cardiac": "heart muscle",
    "lung": "lung", "bronchi": "lung", "kidney": "kidney", "renal": "kidney",
    "liver": "liver", "colon": "colon", "intestine": "colon",
    "brain": "cerebral cortex", "neuron": "cerebral cortex",
    "stomach": "stomach", "pancreas": "pancreas",
    "prostate": "prostate", "testis": "testis",
    "ovary": "ovary", "uterus": "uterus", "endometrium": "endometrium",
    "skin": "skin", "muscle": "skeletal muscle",
    "thyroid": "thyroid gland", "adrenal": "adrenal gland",
}
