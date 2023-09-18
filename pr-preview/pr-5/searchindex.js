Search.setIndex({"docnames": ["1-Introduction/LinearAlgebra", "1-Introduction/QuantumManyBody", "1-Introduction/Software", "1-Introduction/Symmetries", "1-Introduction/TensorNetworks", "2-Tensors/TensorKit", "2-Tensors/TensorOperations", "3-MatrixProductStates/Algorithms", "3-MatrixProductStates/Applications", "3-MatrixProductStates/InfiniteMPS", "3-MatrixProductStates/MatrixProductOperators", "3-MatrixProductStates/MatrixProductStates", "4-Algorithms/FixedpointAlgorithms", "References", "intro"], "filenames": ["1-Introduction/LinearAlgebra.md", "1-Introduction/QuantumManyBody.md", "1-Introduction/Software.md", "1-Introduction/Symmetries.md", "1-Introduction/TensorNetworks.md", "2-Tensors/TensorKit.ipynb", "2-Tensors/TensorOperations.ipynb", "3-MatrixProductStates/Algorithms.md", "3-MatrixProductStates/Applications.md", "3-MatrixProductStates/InfiniteMPS.md", "3-MatrixProductStates/MatrixProductOperators.md", "3-MatrixProductStates/MatrixProductStates.md", "4-Algorithms/FixedpointAlgorithms.md", "References.md", "intro.md"], "titles": ["<span class=\"section-number\">3. </span>(Multi-) Linear Algebra", "<span class=\"section-number\">1. </span>Quantum Many-Body Theory", "<span class=\"section-number\">2. </span>Getting started with Julia", "<span class=\"section-number\">5. </span>Symmetries in quantum many-body physics", "<span class=\"section-number\">4. </span>Tensor Network Theory", "<span class=\"section-number\">6. </span>TensorKit.jl", "<span class=\"section-number\">7. </span>TensorOperations.jl", "<span class=\"section-number\">11. </span>Algorithms", "<span class=\"section-number\">12. </span>Applications", "<span class=\"section-number\">10. </span>Infinite Matrix Product States", "<span class=\"section-number\">9. </span>Matrix Product Operators", "<span class=\"section-number\">8. </span>Matrix Product States", "<span class=\"section-number\">13. </span>Fixed-Point algorithms", "<span class=\"section-number\">14. </span>References", "Tensor Network Methods with Julia"], "terms": {"thi": [0, 3, 4, 9, 12, 14], "lectur": [0, 4, 12, 13], "cover": 0, "some": [0, 3, 4, 9], "basic": [0, 4, 7, 9], "concept": [0, 3, 4, 9], "oper": [0, 3, 9, 14], "serv": 0, "foundat": [0, 4], "most": [0, 3, 4, 9, 12], "what": [0, 3, 4], "follow": [0, 3, 4, 9, 12], "The": [0, 3, 4, 9, 12], "goal": [0, 3], "provid": [0, 4, 9, 14], "intuit": [0, 3, 4, 9], "understand": [0, 4], "without": [0, 3, 9, 12], "insist": 0, "too": 0, "much": [0, 3, 9], "mathemat": [0, 3, 13], "rigour": 0, "import": [0, 3, 4], "introduc": [0, 3, 4, 9, 12], "defin": [0, 3, 4, 9, 12], "resort": 0, "usual": [0, 12], "definit": [0, 4, 9], "which": [0, 3, 4, 9, 12], "veri": [0, 3, 4, 9, 12], "simultan": [0, 3], "also": [0, 3, 4, 9], "showcas": [0, 14], "featur": [0, 3], "tensorkit": [0, 3, 9, 12, 14], "jl": [0, 9, 14], "julia": [0, 9], "packag": [0, 9], "extrem": [0, 3], "well": [0, 3, 4, 9, 12, 14], "suit": [0, 3], "demonstr": [0, 12], "ar": [0, 3, 4, 9, 12], "discuss": [0, 3, 4, 9, 12], "us": [0, 3, 4, 9, 12], "befor": [0, 3, 9], "network": [0, 3, 9, 13], "necessari": [0, 4], "furthermor": [0, 4], "realli": [0, 3, 4, 12], "instruct": 0, "reiter": 0, "case": [0, 3, 9], "noth": [0, 9], "specif": [0, 3, 4, 9], "In": [0, 3, 4, 9, 12], "fact": [0, 3, 4], "mani": [0, 4, 9, 14], "idea": [0, 3, 4], "term": [0, 3, 9], "think": [0, 3], "thought": [0, 3, 4], "from": [0, 3, 4, 9, 12], "viewpoint": 0, "comput": [0, 3, 4, 9, 12], "where": [0, 3, 4, 9, 12], "thei": [0, 3, 4, 9], "repres": [0, 3, 4, 9, 12], "regular": [0, 3], "one": [0, 3, 4, 9, 12], "two": [0, 3, 4, 9, 12], "dimension": [0, 3, 4, 9, 12], "arrai": [0, 4], "either": 0, "real": [0, 3, 12], "complex": [0, 9], "number": [0, 3, 4, 9, 12], "nevertheless": [0, 3, 4, 12], "can": [0, 3, 4, 9, 12], "readili": [0, 9], "gener": [0, 4, 9, 12], "arbitrari": [0, 9, 12], "space": [0, 3, 4, 9, 13], "an": [0, 3, 4, 9, 10, 12, 13], "object": [0, 3, 9, 12], "describ": [0, 3], "list": 0, "correspond": [0, 3, 4, 9, 12], "compon": [0, 3, 4], "basi": [0, 3, 4], "For": [0, 3, 4, 9, 12], "exampl": [0, 4, 14], "its": [0, 3, 4, 9], "form": [0, 3, 4, 9, 10, 12], "vec": [0, 3], "v": [0, 4, 9], "left": [0, 3, 9, 12], "v_1": 0, "v_2": 0, "right": [0, 3, 4, 9, 12], "t": [0, 3, 12, 13], "As": [0, 3, 4, 9, 12], "remind": 0, "properti": [0, 3, 4, 9], "make": [0, 4, 12], "sure": 0, "ad": [0, 3], "togeth": [0, 3], "i": [0, 3, 4, 9, 12], "e": [0, 3, 4, 9], "w": [0, 4], "multipli": [0, 3], "scalar": [0, 4], "alpha": [0, 3, 4, 9], "These": [0, 3, 4, 9], "behav": [0, 3], "expect": [0, 3, 10, 11], "notion": [0, 3, 4], "associ": [0, 3, 9], "commut": [0, 3], "distribut": 0, "given": [0, 3, 4, 9, 12], "necessarili": [0, 3, 4], "distinct": [0, 3, 4], "possibl": [0, 4], "between": [0, 4, 9], "them": [0, 3, 4], "just": [0, 4, 9], "function": 0, "preserv": [0, 4], "structur": [0, 3, 4], "other": [0, 3, 4], "word": [0, 3, 4], "A": [0, 3, 4, 9, 12, 13], "colon": 0, "anoth": [0, 3], "becaus": [0, 4, 9, 12], "requir": [0, 4, 9], "complet": [0, 9], "determin": [0, 3, 4, 9], "action": [0, 3], "lead": [0, 3, 4, 9], "natur": [0, 4, 9], "wai": [0, 3, 4, 12], "matrix": [0, 3, 4, 12, 13], "consid": [0, 3, 4, 9, 12], "construct": [0, 3, 4, 9], "v_i": 0, "w_i": 0, "begin": [0, 3], "rcl": 0, "leftarrow": [0, 9, 12], "sum_j": 0, "a_": [0, 12], "ij": 0, "v_j": 0, "end": [0, 3, 9], "w_j": 0, "abstract": 0, "concret": [0, 3, 4], "particular": [0, 3, 4, 9], "column": [0, 4, 9], "label": [0, 3, 4], "input": 0, "while": [0, 3, 4, 9], "row": 0, "output": 0, "context": [0, 1, 3, 4], "we": [0, 3, 4, 9, 12], "creat": 0, "through": [0, 4, 9, 12], "syntax": 0, "close": [0, 3, 12], "\u2102": [0, 9, 12], "2": [0, 3, 4, 9, 12], "type": 0, "bbc": 0, "tab": 0, "complexspac": [0, 9], "3": [0, 3, 4, 9, 12], "equival": [0, 3, 4], "tensormap": [0, 3, 9], "rand": 0, "float64": [0, 3], "1": [0, 3, 4, 9, 12], "true": [0, 9], "same": [0, 3, 4, 9, 12], "logic": 0, "abov": [0, 3, 4, 9, 12], "combin": [0, 4], "new": [0, 3, 9], "otim": [0, 3, 9], "origin": [0, 9], "equal": [0, 3, 4, 9], "hold": 0, "all": [0, 3, 4, 9, 12], "lambda": [0, 3, 4], "\u03bb": 0, "equip": 0, "canon": [0, 9, 12], "take": [0, 3, 4, 9, 12], "respect": [0, 3, 12], "combinatino": 0, "when": [0, 3, 4, 9, 12], "how": [0, 3, 4, 12], "written": [0, 4], "sum_": [0, 3, 9, 12], "i_1": [0, 4], "i_2": [0, 4], "t_": 0, "i_1i_2": 0, "v_": [0, 12], "w_": 0, "shorthand": 0, "extract": 0, "stridedview": 0, "typeof": 0, "ident": [0, 3, 9], "0": [0, 3, 4, 9, 12], "734924": 0, "754518": 0, "800525": 0, "937594": 0, "535333": 0, "248799": 0, "here": [0, 3, 4, 9, 12], "tent": 0, "name": [0, 3, 4, 12], "wa": [0, 3, 4, 12], "denot": [0, 3], "induc": 0, "more": [0, 3, 4, 9, 12], "common": [0, 4], "express": [0, 3, 4, 9], "reshap": 0, "th": [0, 9], "built": 0, "similarli": [0, 3, 4, 9, 12], "than": [0, 3, 9], "final": [0, 3, 4, 9, 12], "element": [0, 3, 4, 9, 12], "up": [0, 3, 4, 9, 12], "addition": [0, 4], "laid": 0, "out": [0, 3, 4, 12], "slight": 0, "misus": 0, "terminolog": 0, "call": [0, 3, 4, 9, 12], "indic": [0, 9], "cartesian": 0, "cdot": [0, 3, 4], "i_n": 0, "trick": 0, "allow": [0, 3, 4, 9], "reinterpret": [0, 4], "vice": [0, 3], "versa": [0, 3], "linearindic": 0, "tupl": 0, "unitrang": 0, "int64": 0, "5": [0, 4, 9, 12], "4": [0, 3, 4, 9, 12], "6": [0, 3, 4, 12], "collect": 0, "cartesianindic": 0, "forc": 0, "print": 0, "cartesianindex": 0, "due": [0, 3, 4], "itself": [0, 3, 12], "again": [0, 3, 4, 9, 12], "keep": [0, 3], "mind": 0, "now": [0, 3, 9, 12], "howev": [0, 3, 4, 9, 12], "themselv": [0, 3], "compris": 0, "If": [0, 3, 4, 9], "order": 0, "establish": [0, 4], "w_1": 0, "w_2": 0, "w_m": 0, "v_n": 0, "mapsto": 0, "j_1": 0, "j_2": 0, "j_n": 0, "i_m": 0, "n": [0, 3, 4, 9], "j": [0, 3, 4, 12, 13], "v1": 0, "v2": 0, "w1": 0, "w2": 0, "attent": 0, "reader": [0, 3], "might": [0, 3, 12], "have": [0, 3, 4, 9, 12, 14], "alreadi": [0, 3, 9], "note": [0, 3, 4, 9, 12, 13], "strongli": 0, "resembl": 0, "coincid": 0, "easili": [0, 4, 9, 12], "identifi": [0, 3], "identif": 0, "isomorph": [0, 4], "cong": 0, "b": [0, 3, 4], "finit": [0, 3, 9, 12], "addit": [0, 3, 4], "trivial": [0, 3, 9], "choic": [0, 3, 4, 9], "uniqu": [0, 3, 4, 9], "differ": [0, 3, 4, 9, 12], "major": [0, 4], "set": [0, 3, 4, 12], "mean": [0, 3, 9], "still": [0, 3, 4, 9], "relat": [0, 3, 9], "entir": [0, 9], "summar": [0, 3], "along": [0, 4], "constitu": [0, 4], "like": [0, 3, 4], "lift": 0, "facet": 0, "introduct": [1, 3, 4], "histori": 1, "purpos": 1, "relev": [1, 4], "second": [1, 3, 4, 9, 12], "quantiz": 1, "tensor": [1, 9, 12, 13], "product": [1, 12, 13], "classic": [1, 3, 8], "section": [3, 4, 9, 12], "give": [3, 4, 9, 14], "gentl": 3, "framework": [3, 9], "least": [3, 4], "restrict": [3, 9, 12], "our": [3, 9, 12], "illustr": [3, 4, 12, 14], "rather": [3, 12], "first": [3, 4, 9, 12], "coupl": 3, "model": [3, 8, 12], "gradual": 3, "build": 3, "finish": 3, "present": 3, "It": [3, 4, 12, 14], "goe": 3, "sai": 3, "onli": [3, 4, 9, 12], "scratch": 3, "surfac": [3, 4], "vast": 3, "topic": 3, "interest": 3, "refer": [3, 9, 14], "immens": 3, "literatur": 3, "special": 3, "cours": [3, 4, 13], "recal": [3, 9], "transvers": [3, 8, 12], "field": [3, 4, 8, 12], "Ising": [3, 8, 12], "Its": 3, "degre": 3, "freedom": [3, 4, 9, 11], "qubit": 3, "lattic": [3, 9, 12, 13], "hamiltonian": [3, 12], "read": 3, "h": [3, 12], "sigma": [3, 4], "z_i": [3, 12], "z_": 3, "h_x": [3, 12], "sum_i": [3, 9, 12], "x_i": [3, 12], "let": [3, 9, 12], "simpli": [3, 4, 12], "period": 3, "boundari": [3, 9, 12], "condit": [3, 9, 12], "besid": 3, "obviou": 3, "translat": [3, 9, 12], "below": [3, 12], "invari": [3, 9, 12], "under": [3, 9, 12], "flip": [3, 4], "spin": [3, 12], "z": [3, 9, 12], "pauli": [3, 12], "ket": [3, 9, 12], "uparrow": 3, "leftrightarrow": 3, "downarrow": 3, "That": 3, "constitut": 3, "clear": [3, 4, 12], "energi": [3, 9, 12], "depend": [3, 9, 13], "neighbour": [3, 4], "being": [3, 4, 12], "anti": 3, "align": 3, "clearli": [3, 9], "extern": 3, "magnet": 3, "orthogon": [3, 4, 9, 12], "implement": [3, 4, 9, 12], "correctli": 3, "unitari": [3, 4, 9, 12], "p": [3, 4, 12], "bigotimes_i": 3, "notic": [3, 12], "accord": [3, 12], "twice": 3, "leav": [3, 4, 9], "untouch": 3, "dagger": [3, 4, 9, 12], "hp": 3, "everi": [3, 9, 12], "thu": [3, 4, 12], "even": [3, 4, 12], "though": [3, 12], "ha": [3, 4, 9, 12], "regardless": 3, "valu": [3, 10, 11, 12], "you": [3, 12], "know": 3, "previou": [3, 4, 9, 12], "ground": [3, 9, 12], "state": [3, 4, 12, 13], "subspac": 3, "phenomenon": 3, "known": [3, 4, 12], "spontanu": 3, "ssb": 3, "short": 3, "investig": 3, "vanish": 3, "infinit": [3, 14], "rightarrow": 3, "infti": [3, 9], "effectli": 3, "reduc": [3, 4, 9, 12], "paramagnet": 3, "psi_": 3, "frac": [3, 9, 12], "sqrt": [3, 9], "eigenvalu": [3, 9, 12], "eigenvector": [3, 4, 9, 12], "x": [3, 4, 12], "reason": [3, 4, 9], "mention": [3, 12], "disord": 3, "minim": [3, 4, 12], "ferromagnet": 3, "obvious": [3, 12], "contrari": [3, 12], "span": 3, "get": [3, 14], "map": [3, 4, 9], "onto": [3, 9], "broken": 3, "sinc": [3, 9, 12], "degeneraci": 3, "integ": [3, 4], "chang": [3, 4], "smoothli": 3, "slowli": 3, "turn": [3, 12], "therefor": [3, 9], "small": [3, 12], "larg": [3, 4, 12], "said": 3, "belong": 3, "transit": 3, "abruptli": 3, "place": [3, 4, 9], "happen": [3, 4], "becom": [3, 4, 12], "critic": [3, 9, 12], "inspir": 3, "credo": 3, "local": [3, 9, 12], "probe": 3, "wit": 3, "magnetis": 3, "site": [3, 9, 12], "o": [3, 9], "anticommut": 3, "op": 3, "expecti": 3, "braket": [3, 12], "latter": 3, "could": 3, "chosen": [3, 4], "so": [3, 4, 9, 12], "seem": 3, "ill": 3, "remedi": 3, "sign": 3, "select": 3, "after": [3, 4], "limit": [3, 12], "taken": [3, 12], "synopsi": 3, "singl": [3, 4, 9], "particl": 3, "mechan": [3, 10, 14], "multipl": [3, 4], "part": [3, 9, 12], "pattern": 3, "hallmark": 3, "doe": [3, 12], "paradigm": 3, "classifi": 3, "base": [3, 4, 9], "principl": [3, 13], "put": 3, "forward": 3, "landau": 3, "bear": 3, "hi": 3, "rememb": 3, "s": [3, 9, 12], "theorem": 3, "continu": [3, 4], "system": [3, 4, 9, 12], "often": [3, 4], "via": [3, 4], "lagrangian": 3, "rise": [3, 9], "current": [3, 12], "almost": 3, "impli": [3, 4, 9], "d": [3, 4, 9, 12], "dt": 3, "psi": [3, 9, 12], "proof": [3, 4], "simpl": [3, 4, 7, 9, 12], "exercis": [3, 4], "simplest": [3, 4], "hamiltonain": 3, "consequ": 3, "total": [3, 12], "act": [3, 9, 12], "o_i": 3, "o_it": 3, "o_": 3, "exp": 3, "pi": [3, 13], "ip": 3, "momentum": 3, "By": [3, 9, 12], "virtu": 3, "understood": 3, "good": [3, 4], "eigenst": 3, "translation": 3, "non": [3, 9], "implic": 3, "xxz": 3, "heisenberg": 3, "whose": 3, "x_": 3, "y_i": 3, "y_": 3, "delta": [3, 4], "2s": 3, "satisfi": [3, 9], "mathfrak": [3, 12], "su": [3, 12], "a_i": [3, 12], "b_j": 3, "delta_": 3, "sum_c": 3, "varepsilon_": 3, "abc": 3, "c_i": [3, 9], "xxx": [3, 12], "y": 3, "neq": 3, "full": [3, 4, 12], "half": [3, 9], "see": [3, 9, 12], "wherea": [3, 12], "simeq": 3, "u": [3, 4, 9, 12], "retain": 3, "automat": [3, 9], "theta": 3, "interpret": [3, 4, 9, 13], "rotat": [3, 4], "around": [3, 4], "axi": 3, "angl": 3, "quantiti": 3, "m_z": 3, "motiv": 3, "gentli": 3, "backbon": 3, "roughli": 3, "speak": 3, "g": 3, "rule": 3, "compos": 3, "step": [3, 12], "time": [3, 4, 9, 12, 13], "discret": 3, "consist": [3, 4, 9], "carri": [3, 4], "transform": [3, 4, 9], "should": [3, 4, 9, 12], "about": [3, 12], "result": [3, 4, 9, 12], "ani": [3, 4, 9], "next": [3, 4, 9], "over": [3, 4, 12], "theta_2": 3, "theta_1": 3, "g_1": 3, "g_2": 3, "endow": 3, "There": 3, "exist": [3, 4], "1g": 3, "g1": 3, "foral": 3, "abelian": 3, "composit": 3, "k": [3, 4], "hk": 3, "gh": 3, "would": [3, 9], "formal": [3, 9], "undon": 3, "opposit": 3, "henc": [3, 12], "invers": [3, 4], "gg": 3, "subgroup": 3, "suggest": [3, 12], "subset": 3, "li": 3, "heart": 3, "mathbb": [3, 4, 9], "_2": [3, 12], "explan": [3, 12], "notat": [3, 9], "undergo": 3, "symbol": 3, "law": [3, 11], "_n": 3, "modulo": 3, "c": [3, 4, 9, 12, 13], "encount": 3, "unimodular": 3, "matric": [3, 4, 9, 12], "det": 3, "uu": 3, "geq": 3, "none": 3, "3d": 3, "unit": [3, 9], "m": [3, 9], "r": [3, 4, 9, 12], "mm": 3, "tm": 3, "were": [3, 4], "deal": [3, 9], "question": 3, "absenc": 3, "invert": 3, "linear": [3, 4, 12, 14], "singular": [3, 9, 12], "do": 3, "wonder": 3, "come": 3, "exactli": [3, 9, 12], "underli": 3, "linearli": 3, "vector": [3, 4, 9, 12], "immedi": 3, "rais": 3, "plethora": 3, "kind": [3, 4], "ones": 3, "answer": 3, "sake": 3, "index": [3, 9], "x_g": 3, "x_gx_h": 3, "alwai": [3, 9], "dimens": [3, 4, 9, 12], "probabl": 3, "x_0": 3, "x_1": 3, "inde": [3, 4, 9], "x_1x_1": 3, "squar": [3, 4], "bar": [3, 9], "x_h": 3, "overlin": 3, "equiv": 3, "y_g": 3, "kroneck": [3, 4], "check": [3, 9], "oplu": 3, "observ": [3, 9], "unitarili": 3, "ux_gu": 3, "independ": [3, 4], "pmatrix": 3, "show": [3, 9], "crux": 3, "appropri": 3, "brought": 3, "block": [3, 10], "diagon": [3, 4, 9], "fbox": 3, "1_g": 3, "2_g": 3, "vdot": 3, "ddot": 3, "amongst": 3, "irrep": 3, "shown": [3, 4, 9], "d_": 3, "One": [3, 4, 9], "kei": 3, "decompos": [3, 4], "sometim": 3, "fusion": 3, "clebsch": 3, "gordan": 3, "coeffici": 3, "explicitli": [3, 4, 9], "schur": 3, "lemma": 3, "x_gy": 3, "yx_g": 3, "proport": 3, "pose": 3, "s_1": 3, "s_2": 3, "bigoplus_": 3, "been": [3, 4, 14], "analyt": [3, 12], "low": [3, 4, 9], "tabul": 3, "z_g": 3, "strong": 3, "didn": 3, "assum": [3, 9, 12], "argu": 3, "bring": [3, 9, 12], "appear": [3, 12], "write": [3, 9], "bigoplus_c": 3, "b_c": 3, "_c": 3, "decomposit": [3, 9, 12], "store": [3, 9], "effici": [3, 4, 9, 12], "track": 3, "particularli": [3, 4, 9], "paragraph": 3, "herebi": 3, "drastic": 3, "amount": [3, 4], "memori": 3, "abl": 3, "manipul": 3, "exploit": 3, "maximum": 3, "rank": [3, 4], "su\u2082spac": 3, "l": [3, 9, 12], "rep": 3, "su\u2082": 3, "essenti": 3, "copi": 3, "summand": 3, "want": 3, "ss": 3, "productspac": [3, 9], "data": 3, "fusiontre": 3, "fals": 3, "7": [3, 12, 13], "0025861102e": 3, "313": 3, "inspect": 3, "domain": [3, 4], "codomain": [3, 4], "assert": [3, 9], "dim": 3, "sortedvectordict": 3, "su2irrep": 3, "entri": 3, "00259e": 3, "fill": 3, "compat": 3, "cannot": 3, "fuse": 3, "third": 3, "92162e": 3, "310": 3, "four": [3, 4], "conclud": 3, "global": [3, 9, 12], "familiar": 3, "gaug": [3, 11, 12], "ubiquit": 3, "phenomena": [3, 9], "actual": 3, "each": [3, 4, 9, 12], "redund": 3, "descript": 3, "brief": [3, 4, 12], "overview": [3, 9], "mostli": [3, 4, 9], "neglect": 3, "spatial": 3, "reflect": [3, 4], "don": [3, 12], "anymor": 3, "classif": 3, "notori": 3, "rich": 3, "beauti": [3, 12], "especi": [3, 4], "higher": [3, 4], "algorithm": [3, 4, 9], "tremend": 3, "speedup": 3, "stabil": 3, "benefit": 3, "uniform": [3, 9, 12, 13], "repeat": [3, 4, 9, 12], "indefinit": 3, "discoveri": 3, "topolog": 3, "matter": [3, 4], "anyon": 3, "excit": 3, "grow": 3, "fascin": [3, 4], "explor": [3, 4], "categor": 3, "beyond": [3, 12], "scope": [3, 12], "intric": 3, "algebra": [3, 4, 12, 14], "categori": 3, "chain": [3, 9, 12], "storag": 3, "start": [4, 9, 12, 14], "modern": 4, "physic": [4, 9, 12, 13, 14], "simplifi": 4, "quantum": [4, 9, 13, 14], "bodi": [4, 14], "bc17": [4, 12, 13], "journei": 4, "evolut": 4, "profound": 4, "theoret": [4, 13], "develop": [4, 14], "method": [4, 9, 13], "tool": [4, 9], "varieti": [4, 9], "studi": 4, "machin": 4, "learn": [4, 9], "earli": 4, "root": 4, "back": [4, 12], "19th": 4, "centuri": 4, "pioneer": 4, "mathematician": 4, "arthur": 4, "caylei": 4, "jame": 4, "sylvest": 4, "multi": [4, 9, 14], "began": 4, "emerg": 4, "late": 4, "20th": 4, "dmrg": 4, "birth": 4, "attribut": 4, "mp": [4, 9, 11, 12], "1960": 4, "earliest": 4, "wide": 4, "steven": 4, "white": 4, "1992": 4, "simul": 4, "inform": 4, "1980": 4, "1990": 4, "driven": 4, "add": 4, "entangl": [4, 11], "becam": 4, "central": [4, 9], "progress": 4, "extend": [4, 9], "tn": 4, "project": 4, "pair": 4, "pep": 4, "scale": [4, 9], "renorm": [4, 12], "ansatz": [4, 9, 11, 12], "mera": 4, "2000": 4, "disciplin": 4, "appli": [4, 9], "promin": 4, "unsuprisingli": 4, "plai": [4, 9], "role": [4, 9], "languag": 4, "circuit": 4, "ongo": 4, "research": 4, "applic": [4, 14], "vibrant": 4, "evolv": 4, "variou": 4, "direct": [4, 12], "phase": 4, "main": 4, "advantag": [4, 12], "admit": 4, "greatli": 4, "involv": 4, "numer": [4, 9], "node": 4, "graph": 4, "depict": 4, "leg": [4, 9], "stick": [4, 12], "individu": [4, 12], "recoverd": 4, "fix": [4, 9, 14], "open": [4, 12], "diagram": [4, 9], "r_": [4, 12], "i_3": 4, "i_4": 4, "eq": [4, 9], "tensor_isomorph": 4, "freeli": [4, 9], "move": [4, 9, 12], "long": 4, "shape": 4, "certain": [4, 9], "explicit": 4, "represent": [4, 12], "ow": 4, "precis": [4, 12], "detail": [4, 9, 12], "convent": 4, "d_1": 4, "d_2": 4, "d_3": 4, "d_i": 4, "immateri": 4, "Of": 4, "eachoth": 4, "three": [4, 9], "complic": 4, "join": 4, "similar": [4, 9], "einstein": 4, "summat": 4, "signifi": 4, "partial": 4, "cyclic": 4, "slide": 4, "loop": 4, "placement": 4, "tr": [4, 9], "ab": 4, "ba": 4, "found": [4, 9], "drawn": 4, "famililiar": 4, "inner": 4, "arbitrarili": 4, "evalu": [4, 9, 12], "sequenc": 4, "wise": 4, "specifi": [4, 9, 12], "text": 4, "indix": 4, "assign": 4, "implicitli": 4, "sum": [4, 9], "beta": [4, 9], "gamma": 4, "epsilon": [4, 12], "zeta": 4, "eta": [4, 12], "f": 4, "quickli": 4, "unwieldi": 4, "wish": [4, 9], "pairwis": 4, "spirit": 4, "minor": 4, "modif": 4, "ncon": 4, "increas": 4, "neg": 4, "posit": [4, 9], "vari": 4, "wildli": 4, "problem": [4, 12], "both": 4, "size": [4, 9, 12], "length": [4, 9, 12], "2n": 4, "float": 4, "point": [4, 9, 14], "flop": 4, "substanti": 4, "prefer": 4, "typic": [4, 9], "area": [4, 11], "boil": [4, 12], "down": [4, 9, 12], "heurist": 4, "cut": 4, "bubbl": 4, "ineffici": 4, "infeas": 4, "optim": [4, 9, 12], "ladder": 4, "highlight": 4, "few": 4, "np": 4, "hard": 4, "larger": 4, "find": [4, 9, 12], "30": [4, 12], "40": 4, "instrument": 4, "approxim": [4, 9, 12], "partit": 4, "everyth": [4, 9], "piec": 4, "need": 4, "eigen": 4, "ax": 4, "imag": 4, "normal": [4, 12], "made": 4, "spectral": [4, 9], "equat": [4, 12], "av": 4, "diagrammat": [4, 12], "svd": 4, "seen": [4, 9], "eigendecomposit": 4, "rectangular": 4, "isometr": [4, 12], "best": [4, 12], "truncat": 4, "largest": [4, 9], "tensori": 4, "version": [4, 9], "semi": 4, "hermitian": 4, "decopmosit": 4, "q": 4, "upper": 4, "triangular": 4, "solv": [4, 12], "solut": [4, 12], "easi": 4, "gaussian": 4, "elimin": 4, "overdetermin": 4, "variant": 4, "qrpo": 4, "transpos": 4, "rq": 4, "ql": 4, "lq": 4, "reveal": 4, "trade": 4, "off": 4, "expens": 4, "practic": 4, "zero": 4, "commonli": [4, 9], "perform": [4, 9], "test": [5, 6], "updat": [7, 9, 12], "trotter": 7, "tebd": 7, "excel": [9, 12], "review": [9, 13], "vhv19": [9, 12, 13], "thorough": 9, "technic": 9, "tangent": [9, 13], "exposit": [9, 12], "supplement": 9, "work": [9, 12], "routin": [9, 12], "tutori": [9, 14], "hilbert": 9, "rangl": 9, "boldsymbol": 9, "_l": 9, "prod_": 9, "s_m": 9, "_r": 9, "altern": 9, "view": 9, "bond": [9, 12], "control": 9, "10": [9, 12, 13], "accuraci": 9, "class": [9, 12], "gap": [9, 12], "accur": 9, "smaller": 9, "never": 9, "safe": 9, "ignor": 9, "bulk": 9, "faithfulli": 9, "captur": 9, "impos": 9, "transat": 9, "diagramat": 9, "instead": 9, "cell": 9, "techniqu": [9, 12], "unform": 9, "calcul": 9, "transfer": 9, "character": [9, 12], "norm": [9, 12], "contract": 9, "power": 9, "magnitud": 9, "lambda_0": 9, "lambda_i": 9, "remain": 9, "mangitud": 9, "projector": 9, "To": 9, "ensur": 9, "properli": 9, "rescal": 9, "trace": 9, "With": 9, "overlap": 9, "effect": [9, 12], "choos": [9, 12], "langl": 9, "middl": 9, "suppos": 9, "extens": 9, "o_n": 9, "dictat": 9, "look": 9, "bra": 9, "beta_m": 9, "alpha_n": 9, "abritrari": 9, "locat": 9, "insert": 9, "disconnect": 9, "rest": 9, "exponenti": 9, "decai": 9, "connect": 9, "why": 9, "harder": 9, "xi": 9, "lambda_1": 9, "log": 9, "lambda_": 9, "mathrm": 9, "max": [9, 12], "sublead": 9, "focuss": 9, "symmetri": [9, 14], "sector": 9, "target": 9, "exit": 9, "crucial": 9, "despit": 9, "inher": 9, "convers": 9, "mai": [9, 13], "parametr": 9, "orthonorm": [9, 12], "a_l": [9, 12], "iter": [9, 12], "procedur": 9, "qr": 9, "docomposit": 9, "initi": [9, 12], "guess": [9, 12], "repeatedli": [9, 12], "bound": 9, "converg": [9, 12], "room": 9, "a_r": [9, 12], "mix": [9, 12], "center": [9, 12], "a_c": [9, 12], "obtain": [9, 12], "contrast": 9, "lr": 9, "therebi": 9, "link": 9, "usv": 9, "absorb": 9, "residu": 9, "rm": 9, "arriv": 9, "straightforwardli": 9, "schmidt": 9, "across": 9, "i_l": 9, "i_r": 9, "bipartit": 9, "enabl": [9, 12], "sens": 9, "maxim": 9, "isometri": [9, 12], "correspondingli": 9, "guarante": 9, "lower": 9, "variat": [9, 12, 13], "cost": 9, "tild": [9, 12], "go": 9, "aspect": 9, "virtual": 9, "standard": 9, "cr": 9, "al": 9, "linearalgebra": 9, "ac": 9, "0000000000000002": 9, "verifi": 9, "tensoroper": [9, 14], "macro": 9, "al_id": 9, "conj": 9, "ar_id": 9, "id": 9, "lh": 9, "rh": 9, "randn": 9, "expectation_valu": [9, 12], "complexf64": [9, 12], "12594290002782813": 9, "04590018022880266im": 9, "encod": [9, 12], "correlation_length": 9, "5165490768548735": 9, "export": 9, "statist": [10, 14], "mpo": [10, 12], "jordan": 10, "polynomi": 11, "correl": 11, "manner": 12, "min_": 12, "simplic": 12, "consider": 12, "interact": 12, "rang": 12, "densiti": 12, "group": 12, "thermodynam": 12, "break": 12, "random": 12, "tri": 12, "sequenti": 12, "sweep": 12, "until": 12, "reach": 12, "bit": 12, "a_1": 12, "a_2": 12, "seemingli": 12, "daunt": 12, "those": 12, "denomin": 12, "mathcal": 12, "h_i": 12, "_i": 12, "smallest": 12, "forth": 12, "At": 12, "reus": 12, "manifestli": 12, "suffer": 12, "artefact": 12, "surprisingli": 12, "proven": 12, "success": 12, "variation": 12, "mpskit": 12, "mpskitmodel": 12, "z_j": 12, "h_z": 12, "free": 12, "paramet": 12, "factor": 12, "16": 12, "12": 12, "default": 12, "hand": [12, 13, 14], "transverse_field_is": 12, "summon": 12, "\u03c8": 12, "finitemp": 12, "\u03c8\u2080": 12, "_": 12, "find_groundst": 12, "36m": 12, "1m": 12, "22m": 12, "39m": 12, "1minfo": 12, "39miteraton": 12, "error": 12, "011595463511624361": 12, "582283370855894e": 12, "9": 12, "460475630158445e": 12, "468399699008025e": 12, "277891267746376e": 12, "71255942169506e": 12, "8": 12, "7399667830750483e": 12, "431380902833998e": 12, "3795766268948445e": 12, "812476228906224e": 12, "2665001950294204e": 12, "11": [12, 13], "211856986221164e": 12, "499906405675467e": 12, "13": 12, "6724368716204494e": 12, "14": 12, "260336721279763e": 12, "15": 12, "3174101534118947e": 12, "672167805144742e": 12, "directli": 12, "unbound": 12, "f_l": 12, "f_r": 12, "obei": 12, "offer": 12, "design": 12, "reli": 12, "intial": 12, "h_": 12, "h_c": 12, "a_lc": 12, "ca_r": 12, "deriv": 12, "minimum": 12, "manifold": 12, "last": 12, "intertwin": 12, "explain": 12, "chose": 12, "toler": 12, "arnoldi": 12, "yield": 12, "epsilon_l": 12, "epsilon_r": 12, "comment": 12, "soltuion": 12, "u_lv_l": 12, "u_l": 12, "v_l": 12, "aris": 12, "sigma_lv_l": 12, "u_rv_r": 12, "u_r": 12, "sigma_rv_r": 12, "approach": 12, "exact": 12, "s_c": 12, "s_lc": 12, "ca": 12, "s_r": 12, "sigma_": 12, "arithmet": 12, "u_": 12, "poor": 12, "robust": 12, "l_": 12, "l_c": 12, "qquad": 12, "r_c": 12, "polar": 12, "l_cp": 12, "r_cu": 12, "estim": 12, "cdot1": 12, "heisenberg_xyz": 12, "infinitemp": 12, "env": 12, "39mvump": 12, "galerkin": 12, "4220089820747899": 12, "2020580471723": 12, "2400733002372661": 12, "26210352292048267": 12, "02955160285293315": 12, "015641648430033624": 12, "006166599355735042": 12, "002299243149539928": 12, "0007714101040879418": 12, "0002828370948092264": 12, "937353857428135e": 12, "7604713611165196e": 12, "355150803673915e": 12, "226940877724936e": 12, "916784999497858e": 12, "466905367258202e": 12, "17": 12, "7757989783168754e": 12, "18": 12, "0867105636325642e": 12, "19": 12, "0852361661610255e": 12, "20": 12, "60329468214558e": 12, "21": 12, "083354427897799e": 12, "22": [12, 13], "390405673351967e": 12, "23": 12, "139456873879459e": 12, "24": 12, "593549878924305e": 12, "25": 12, "382509596649082e": 12, "26": 12, "437551009210185e": 12, "27": 12, "1029682785138377e": 12, "28": 12, "270699926054731e": 12, "29": 12, "231874137792546e": 12, "2832422910470052e": 12, "31": 12, "029856999474449e": 12, "periodicarrai": 12, "4013806435182499": 12, "627715291671279e": 12, "17im": 12, "compar": 12, "quasi": 12, "401": 12, "484": 12, "038": 12, "971": 12, "hco": [12, 13], "decim": 12, "jacob": 13, "bridgeman": 13, "christoph": 13, "chubb": 13, "wave": 13, "danc": 13, "introductori": 13, "journal": 13, "50": 13, "223001": 13, "2017": 13, "url": 13, "http": 13, "dx": 13, "doi": 13, "org": 13, "1088": 13, "1751": 13, "8121": 13, "aa6dc3": 13, "jutho": 13, "haegeman": 13, "ignacio": 13, "cirac": 13, "tobia": 13, "osborn": 13, "iztok": 13, "\u017e": 13, "orn": 13, "henri": 13, "verscheld": 13, "frank": 13, "verstraet": 13, "letter": 13, "107": 13, "070601": 13, "2011": 13, "lauren": 13, "vanderstraeten": 13, "scipost": 13, "page": 13, "007": 13, "januari": 13, "2019": 13, "scipostphyslectnot": 13, "arxiv": 13, "1810": 13, "07006": 13, "21468": 13, "seri": 14, "theori": 14, "aim": 14, "practis": 14, "code": 14, "softwar": 14, "librari": 14}, "objects": {}, "objtypes": {}, "objnames": {}, "titleterms": {"multi": 0, "linear": 0, "algebra": 0, "content": [0, 4, 9], "overview": [0, 4], "vector": 0, "matric": 0, "tensor": [0, 3, 4, 14], "product": [0, 3, 4, 9, 10, 11, 14], "map": 0, "conclus": [0, 4], "quantum": [1, 3], "mani": [1, 3], "bodi": [1, 3], "theori": [1, 3, 4], "get": 2, "start": 2, "julia": [2, 14], "symmetri": 3, "physic": 3, "exampl": [3, 9, 12], "applic": [3, 8], "break": 3, "order": [3, 4], "paramet": 3, "phase": 3, "noether": 3, "conserv": 3, "quantitit": 3, "group": [3, 4], "represent": [3, 9], "definit": 3, "complex": [3, 4], "conjug": 3, "direct": 3, "sum": 3, "irreduc": 3, "symmetr": 3, "outlook": 3, "gener": 3, "network": [4, 14], "histori": 4, "graphic": 4, "notat": 4, "oper": [4, 10], "index": 4, "split": 4, "indic": 4, "outer": 4, "trace": 4, "contract": 4, "factor": 4, "eigenvalu": 4, "decomposit": 4, "singular": 4, "valu": [4, 9], "polar": 4, "qr": 4, "nullspac": 4, "tensorkit": 5, "jl": [5, 6], "tensoroper": 6, "algorithm": [7, 12, 14], "infinit": 9, "matrix": [9, 10, 11, 14], "state": [9, 11, 14], "thermodynam": 9, "limit": 9, "normal": 9, "expect": 9, "correl": 9, "function": 9, "gaug": 9, "revisit": 9, "entangl": 9, "entropi": 9, "truncat": 9, "code": 9, "mpskit": 9, "infinitemp": 9, "fix": 12, "point": 12, "dmrg": 12, "vump": 12, "refer": 13, "method": 14, "introduct": 14, "other": 14}, "envversion": {"sphinx.domains.c": 2, "sphinx.domains.changeset": 1, "sphinx.domains.citation": 1, "sphinx.domains.cpp": 6, "sphinx.domains.index": 1, "sphinx.domains.javascript": 2, "sphinx.domains.math": 2, "sphinx.domains.python": 3, "sphinx.domains.rst": 2, "sphinx.domains.std": 2, "sphinx.ext.intersphinx": 1, "sphinxcontrib.bibtex": 9, "sphinx": 56}})