Search.setIndex({"docnames": ["1-Introduction/Julia", "1-Introduction/QuantumManyBody", "1-Introduction/TensorNetworks", "2-Tensors/TensorKit", "2-Tensors/TensorOperations", "3-MatrixProductStates/Algorithms", "3-MatrixProductStates/Applications", "3-MatrixProductStates/InfiniteMPS", "3-MatrixProductStates/MatrixProductOperators", "3-MatrixProductStates/MatrixProductStates", "References", "intro"], "filenames": ["1-Introduction/Julia.md", "1-Introduction/QuantumManyBody.md", "1-Introduction/TensorNetworks.md", "2-Tensors/TensorKit.ipynb", "2-Tensors/TensorOperations.ipynb", "3-MatrixProductStates/Algorithms.md", "3-MatrixProductStates/Applications.md", "3-MatrixProductStates/InfiniteMPS.md", "3-MatrixProductStates/MatrixProductOperators.md", "3-MatrixProductStates/MatrixProductStates.md", "References.md", "intro.md"], "titles": ["<span class=\"section-number\">3. </span>Getting started with Julia", "<span class=\"section-number\">1. </span>Quantum Many-Body Theory", "<span class=\"section-number\">2. </span>Tensor Network Theory", "<span class=\"section-number\">4. </span>TensorKit.jl", "<span class=\"section-number\">5. </span>TensorOperations.jl", "<span class=\"section-number\">9. </span>Algorithms", "<span class=\"section-number\">10. </span>Applications", "<span class=\"section-number\">8. </span>Infinite Matrix Product States", "<span class=\"section-number\">7. </span>Matrix Product Operators", "<span class=\"section-number\">6. </span>Matrix Product States", "<span class=\"section-number\">11. </span>References", "Tensor Network Methods with Julia"], "terms": {"introduct": [1, 2], "context": [1, 2], "histori": 1, "purpos": 1, "relev": [1, 2], "second": 1, "quantiz": 1, "tensor": 1, "product": [1, 10], "classic": [1, 6], "In": 2, "thi": [2, 11], "lectur": [2, 10], "we": 2, "introduc": 2, "basic": [2, 5], "concept": 2, "start": [2, 11], "brief": 2, "modern": 2, "physic": [2, 10, 11], "mathemat": 2, "includ": 2, "graphic": 2, "notat": 2, "final": 2, "discuss": [2, 7], "comput": 2, "complex": 2, "quantum": [2, 11], "mani": [2, 11], "bodi": [2, 11], "The": 2, "fascin": 2, "journei": 2, "through": 2, "evolut": 2, "profound": 2, "theoret": 2, "idea": 2, "well": [2, 11], "develop": [2, 11], "method": [2, 10], "tool": 2, "These": 2, "have": [2, 11], "been": [2, 11], "varieti": 2, "especi": 2, "studi": 2, "machin": 2, "learn": 2, "earli": 2, "foundat": 2, "root": 2, "can": 2, "trace": 2, "back": 2, "matrix": [2, 10], "19th": 2, "centuri": 2, "pioneer": 2, "mathematician": 2, "like": 2, "arthur": 2, "caylei": 2, "jame": 2, "sylvest": 2, "dimension": 2, "arrai": 2, "number": 2, "began": 2, "emerg": 2, "late": 2, "20th": 2, "state": [2, 10], "dmrg": 2, "birth": 2, "attribut": 2, "mp": [2, 9], "1960": 2, "One": 2, "earliest": 2, "still": 2, "most": 2, "wide": 2, "us": 2, "algorithm": [2, 11], "It": [2, 11], "wa": 2, "steven": 2, "white": 2, "1992": 2, "provid": [2, 11], "one": 2, "effici": 2, "simul": 2, "system": 2, "inform": 2, "1980": 2, "1990": 2, "field": [2, 6], "driven": 2, "add": 2, "name": 2, "here": 2, "entangl": [2, 9], "becam": 2, "central": 2, "higher": 2, "As": 2, "progress": 2, "were": 2, "extend": 2, "lead": 2, "more": 2, "gener": 2, "tn": 2, "two": 2, "project": 2, "pair": 2, "pep": 2, "scale": 2, "renorm": 2, "ansatz": [2, 9], "mera": 2, "2000": 2, "other": 2, "disciplin": 2, "appli": 2, "promin": 2, "being": 2, "unsuprisingli": 2, "thei": 2, "also": 2, "plai": 2, "role": 2, "where": 2, "natur": 2, "languag": 2, "explor": 2, "circuit": 2, "ongo": 2, "research": 2, "applic": [2, 11], "continu": 2, "vibrant": 2, "evolv": 2, "variou": 2, "direct": 2, "contract": 2, "understand": 2, "phase": 2, "matter": 2, "befor": 2, "necessari": 2, "what": 2, "ar": 2, "furthermor": 2, "realli": 2, "instruct": 2, "reiter": 2, "some": 2, "case": 2, "which": 2, "noth": 2, "specif": 2, "fact": 2, "defin": 2, "term": 2, "think": 2, "follow": 2, "thought": 2, "from": 2, "viewpoint": 2, "repres": 2, "regular": 2, "either": 2, "real": 2, "nevertheless": 2, "much": 2, "trivial": 2, "arbitrari": 2, "space": [2, 10], "an": [2, 8], "object": 2, "describ": 2, "list": 2, "correspond": 2, "compon": 2, "basi": 2, "For": 2, "exampl": [2, 11], "its": 2, "form": [2, 8], "vec": 2, "v": 2, "left": 2, "v_1": 2, "v_2": 2, "right": 2, "t": 2, "remind": 2, "properti": 2, "make": 2, "sure": 2, "oper": [2, 11], "ad": 2, "togeth": 2, "i": 2, "e": 2, "w": 2, "multipli": 2, "scalar": 2, "alpha": 2, "behav": 2, "expect": [2, 7, 8, 9], "notion": 2, "associ": 2, "commut": 2, "distribut": 2, "given": 2, "necessarili": 2, "distinct": 2, "possibl": 2, "between": 2, "them": 2, "just": 2, "function": 2, "preserv": 2, "structur": [2, 7], "word": 2, "A": 2, "colon": 2, "anoth": 2, "becaus": 2, "requir": 2, "complet": 2, "determin": 2, "action": 2, "veri": 2, "wai": 2, "consid": 2, "construct": 2, "v_i": 2, "w_i": 2, "leftarrow": 2, "sum_j": 2, "a_": 2, "ij": 2, "v_j": 2, "w_j": 2, "abstract": 2, "concret": 2, "usual": 2, "particular": 2, "column": 2, "label": 2, "input": 2, "while": 2, "row": 2, "output": 2, "same": 2, "logic": 2, "abov": 2, "combin": 2, "new": 2, "otim": 2, "origin": 2, "equal": 2, "hold": 2, "all": 2, "lambda": 2, "equip": 2, "canon": 2, "take": 2, "respect": 2, "combinatino": 2, "when": 2, "how": 2, "written": 2, "sum_": 2, "i_1": 2, "i_2": 2, "t_": 2, "i_1i_2": 2, "v_": 2, "w_": 2, "tent": 2, "denot": 2, "induc": 2, "common": 2, "express": 2, "reshap": 2, "th": 2, "built": 2, "similarli": 2, "than": 2, "definit": 2, "element": 2, "up": 2, "addition": 2, "laid": 2, "out": 2, "slight": 2, "misus": 2, "terminolog": 2, "call": 2, "indic": 2, "cartesian": 2, "cdot": 2, "i_n": 2, "trick": 2, "allow": 2, "reinterpret": 2, "vice": 2, "versa": 2, "due": 2, "itself": 2, "again": 2, "keep": 2, "mind": 2, "2": 2, "1": 2, "now": 2, "howev": 2, "themselv": 2, "compris": 2, "If": 2, "order": 2, "establish": 2, "begin": 2, "lcr": 2, "w_1": 2, "w_2": 2, "w_m": 2, "v_n": 2, "j_1": 2, "j_2": 2, "j_n": 2, "i_m": 2, "j": 2, "n": 2, "attent": 2, "reader": 2, "might": 2, "alreadi": 2, "note": [2, 10], "strongli": 2, "resembl": 2, "coincid": 2, "easili": 2, "identifi": 2, "identif": 2, "isomorph": 2, "cong": 2, "finit": 2, "without": 2, "addit": 2, "choic": 2, "uniqu": 2, "differ": 2, "major": 2, "set": 2, "equival": 2, "mean": 2, "relat": 2, "entir": 2, "summar": 2, "along": 2, "constitu": 2, "lift": 2, "import": 2, "facet": 2, "test": [3, 4], "simpl": 5, "updat": 5, "trotter": 5, "tebd": 5, "transvers": 6, "Ising": 6, "model": 6, "larg": 7, "base": 7, "vhv19": [7, 10], "gaug": [7, 9], "revisit": 7, "valu": [7, 8, 9], "factor": 7, "statist": [8, 11], "mechan": [8, 11], "mpo": 8, "jordan": 8, "block": 8, "polynomi": 9, "area": 9, "law": 9, "freedom": 9, "correl": 9, "lauren": 10, "vanderstraeten": 10, "jutho": 10, "haegeman": 10, "frank": 10, "verstraet": 10, "tangent": 10, "uniform": 10, "scipost": 10, "page": 10, "007": 10, "januari": 10, "2019": 10, "url": 10, "http": 10, "org": 10, "scipostphyslectnot": 10, "7": 10, "doi": 10, "10": 10, "21468": 10, "seri": 11, "tutori": 11, "illustr": 11, "theori": 11, "aim": 11, "give": 11, "hand": 11, "practis": 11, "code": 11, "showcas": 11, "softwar": 11, "librari": 11, "get": 11, "tensorkit": 11, "jl": 11, "tensoroper": 11, "infinit": 11, "refer": 11}, "objects": {}, "objtypes": {}, "objnames": {}, "titleterms": {"get": 0, "start": 0, "julia": [0, 11], "quantum": 1, "mani": 1, "bodi": 1, "theori": [1, 2], "tensor": [2, 11], "network": [2, 11], "content": 2, "overview": 2, "histori": 2, "linear": 2, "algebra": 2, "vector": 2, "matric": 2, "multi": 2, "product": [2, 7, 8, 9, 11], "map": 2, "conclus": 2, "tensorkit": 3, "jl": [3, 4], "tensoroper": 4, "algorithm": 5, "applic": 6, "infinit": 7, "matrix": [7, 8, 9, 11], "state": [7, 9, 11], "oper": 8, "refer": 10, "method": 11, "introduct": 11, "other": 11}, "envversion": {"sphinx.domains.c": 2, "sphinx.domains.changeset": 1, "sphinx.domains.citation": 1, "sphinx.domains.cpp": 6, "sphinx.domains.index": 1, "sphinx.domains.javascript": 2, "sphinx.domains.math": 2, "sphinx.domains.python": 3, "sphinx.domains.rst": 2, "sphinx.domains.std": 2, "sphinx.ext.intersphinx": 1, "sphinxcontrib.bibtex": 9, "sphinx": 56}})