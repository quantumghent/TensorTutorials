Search.setIndex({"docnames": ["1-Introduction/Julia", "1-Introduction/QuantumManyBody", "1-Introduction/Symmetries", "1-Introduction/TensorNetworks", "2-Tensors/TensorKit", "2-Tensors/TensorOperations", "3-MatrixProductStates/Algorithms", "3-MatrixProductStates/Applications", "3-MatrixProductStates/InfiniteMPS", "3-MatrixProductStates/MatrixProductOperators", "3-MatrixProductStates/MatrixProductStates", "4-Algorithms/FixedpointAlgorithms", "References", "intro"], "filenames": ["1-Introduction/Julia.md", "1-Introduction/QuantumManyBody.md", "1-Introduction/Symmetries.md", "1-Introduction/TensorNetworks.md", "2-Tensors/TensorKit.ipynb", "2-Tensors/TensorOperations.ipynb", "3-MatrixProductStates/Algorithms.md", "3-MatrixProductStates/Applications.md", "3-MatrixProductStates/InfiniteMPS.md", "3-MatrixProductStates/MatrixProductOperators.md", "3-MatrixProductStates/MatrixProductStates.md", "4-Algorithms/FixedpointAlgorithms.md", "References.md", "intro.md"], "titles": ["<span class=\"section-number\">3. </span>Getting started with Julia", "<span class=\"section-number\">1. </span>Quantum Many-Body Theory", "<span class=\"section-number\">4. </span>Symmetries in quantum many-body physics", "<span class=\"section-number\">2. </span>Tensor Network Theory", "<span class=\"section-number\">5. </span>TensorKit.jl", "<span class=\"section-number\">6. </span>TensorOperations.jl", "<span class=\"section-number\">10. </span>Algorithms", "<span class=\"section-number\">11. </span>Applications", "<span class=\"section-number\">9. </span>Infinite Matrix Product States", "<span class=\"section-number\">8. </span>Matrix Product Operators", "<span class=\"section-number\">7. </span>Matrix Product States", "<span class=\"section-number\">12. </span>Fixed-Point algorithms", "<span class=\"section-number\">13. </span>References", "Tensor Network Methods with Julia"], "terms": {"introduct": [1, 2], "context": [1, 2, 3], "histori": [1, 3], "purpos": [1, 3], "relev": [1, 3], "second": [1, 2, 8, 11], "quantiz": 1, "tensor": [1, 8, 11, 12], "product": [1, 11, 12], "classic": [1, 2, 7], "The": [2, 8, 11], "goal": 2, "thi": [2, 8, 11, 13], "section": [2, 8, 11], "give": [2, 8, 13], "veri": [2, 8, 11], "gentl": 2, "concept": [2, 8], "notion": 2, "mathemat": [2, 12], "framework": [2, 8], "least": 2, "we": [2, 8, 11], "restrict": [2, 8, 11], "our": [2, 8, 11], "take": [2, 8, 11], "given": [2, 8, 11], "illustr": [2, 11, 13], "rather": [2, 11], "first": [2, 8, 11], "discuss": [2, 8, 11], "coupl": 2, "import": 2, "some": [2, 8], "concret": 2, "model": [2, 7, 11, 12], "gradual": 2, "build": 2, "up": [2, 8, 11], "more": [2, 8, 11], "finish": 2, "an": [2, 8, 9, 11, 12], "present": 2, "here": [2, 8, 11], "It": [2, 11, 13], "goe": 2, "without": [2, 8, 11], "sai": 2, "onli": [2, 8, 11], "scratch": 2, "surfac": 2, "vast": 2, "topic": 2, "interest": 2, "reader": 2, "refer": [2, 8, 13], "immens": 2, "literatur": 2, "special": 2, "cours": [2, 12], "recal": [2, 8], "one": [2, 8, 11], "dimension": [2, 8, 11], "transvers": [2, 7, 11], "field": [2, 7, 11], "Ising": [2, 7, 11], "defin": [2, 8, 11], "abov": [2, 8, 11], "Its": 2, "degre": 2, "freedom": [2, 8, 10], "ar": [2, 8, 11], "qubit": 2, "lattic": [2, 8, 11, 12], "its": [2, 8], "hamiltonian": [2, 11], "read": 2, "h": [2, 11, 12], "sum_": [2, 8, 11], "i": [2, 8, 11], "sigma": 2, "z_i": [2, 11], "z_": 2, "1": [2, 8, 11], "h_x": [2, 11], "sum_i": [2, 8, 11], "x_i": [2, 11], "let": [2, 8, 11], "us": [2, 8, 11], "simpli": [2, 11], "consid": [2, 8, 11], "period": 2, "boundari": [2, 8, 11], "condit": [2, 8, 11], "besid": 2, "obviou": 2, "translat": [2, 8, 11], "which": [2, 8, 11], "below": [2, 11], "also": [2, 8], "invari": [2, 8, 11], "under": [2, 8, 11], "flip": 2, "all": [2, 8, 11], "spin": [2, 11], "simultan": 2, "z": [2, 8, 11], "e": [2, 8], "pauli": [2, 11], "basi": 2, "ket": [2, 8, 11], "uparrow": 2, "leftrightarrow": 2, "downarrow": 2, "That": 2, "oper": [2, 8, 13], "constitut": 2, "clear": [2, 11], "from": [2, 8, 11], "energi": [2, 8, 11], "term": [2, 8], "depend": [2, 8, 12], "neighbour": 2, "being": [2, 11], "anti": 2, "align": 2, "clearli": [2, 8], "trivial": [2, 8], "extern": 2, "magnet": 2, "orthogon": [2, 8, 11], "implement": [2, 8, 11], "correctli": 2, "repres": [2, 8, 11], "unitari": [2, 8, 11], "p": [2, 11], "bigotimes_i": 2, "notic": [2, 11], "2": [2, 8, 11], "accord": [2, 11], "intuit": [2, 8], "twice": 2, "equival": 2, "leav": [2, 8], "untouch": 2, "fact": 2, "0": [2, 8, 11], "dagger": [2, 8, 11], "hp": 2, "ident": [2, 8], "everi": [2, 8, 11], "thu": [2, 11], "set": [2, 11], "close": [2, 11], "even": [2, 11], "though": [2, 11], "ha": [2, 8, 11], "regardless": 2, "valu": [2, 9, 10, 11], "you": [2, 11], "might": [2, 11], "know": 2, "previou": [2, 8, 11], "ground": [2, 8, 11], "state": [2, 11, 12], "subspac": 2, "necessarili": 2, "phenomenon": 2, "known": [2, 11], "spontanu": 2, "ssb": 2, "short": 2, "investig": 2, "extrem": 2, "case": [2, 8], "vanish": 2, "infinit": [2, 13], "rightarrow": 2, "infti": [2, 8], "In": [2, 8, 11], "effectli": 2, "reduc": [2, 8, 11], "paramagnet": 2, "uniqu": [2, 8], "psi_": 2, "otim": [2, 8], "n": [2, 8], "where": [2, 8, 11], "frac": [2, 8, 11], "sqrt": [2, 8], "eigenvalu": [2, 8, 11], "eigenvector": [2, 8, 11], "x": [2, 11, 12], "other": 2, "word": 2, "For": [2, 8, 11], "reason": [2, 8], "mention": [2, 11], "disord": 2, "minim": [2, 11], "behav": 2, "ferromagnet": 2, "obvious": [2, 11], "two": [2, 8, 11], "distinct": 2, "contrari": [2, 11], "thei": [2, 8], "span": 2, "action": 2, "get": [2, 13], "map": [2, 8], "onto": [2, 8], "vice": 2, "versa": 2, "broken": 2, "sinc": [2, 8, 11], "degeneraci": 2, "integ": 2, "can": [2, 8, 11], "chang": 2, "smoothli": 2, "when": [2, 8, 11], "slowli": 2, "turn": [2, 11], "therefor": [2, 8], "small": [2, 11], "larg": [2, 11], "said": 2, "belong": 2, "differ": [2, 8, 11], "finit": [2, 8, 11], "transit": 2, "abruptli": 2, "expect": [2, 9, 10], "place": [2, 8], "As": [2, 8, 11], "out": [2, 11], "happen": 2, "becom": [2, 11], "critic": [2, 8, 11], "inspir": 2, "credo": 2, "introduc": [2, 8, 11], "local": [2, 8, 11], "probe": 2, "wit": 2, "magnetis": 2, "site": [2, 8, 11], "o": [2, 8], "anticommut": 2, "op": 2, "follow": [2, 8, 11], "expecti": 2, "braket": [2, 11], "while": [2, 8], "howev": [2, 8, 11], "latter": 2, "could": 2, "have": [2, 8, 11, 13], "chosen": 2, "so": [2, 8, 11], "seem": 2, "ill": 2, "remedi": 2, "ad": 2, "lambda": 2, "sign": 2, "select": 2, "after": 2, "limit": [2, 11], "taken": [2, 11], "synopsi": 2, "singl": [2, 8], "particl": 2, "mechan": [2, 9, 13], "multipl": 2, "part": [2, 8, 11], "pattern": 2, "hallmark": 2, "featur": 2, "doe": [2, 11], "commut": 2, "paradigm": 2, "classifi": 2, "base": [2, 8], "principl": [2, 12], "wa": [2, 11], "put": 2, "forward": 2, "landau": 2, "bear": 2, "hi": 2, "name": [2, 11], "rememb": 2, "s": [2, 8, 11], "theorem": 2, "continu": 2, "system": [2, 8, 11], "most": [2, 8, 11], "often": 2, "via": 2, "lagrangian": 2, "rise": [2, 8], "current": [2, 11], "almost": 2, "impli": [2, 8], "d": [2, 8, 11], "dt": 2, "psi": [2, 8, 11], "t": [2, 11, 12], "proof": 2, "left": [2, 8, 11], "simpl": [2, 6, 8, 11], "exercis": 2, "simplest": 2, "hamiltonain": 2, "itself": [2, 11], "consequ": 2, "total": [2, 11], "anoth": 2, "act": [2, 8, 11], "o_i": 2, "o_it": 2, "o_": 2, "exp": 2, "pi": [2, 12], "ip": 2, "number": [2, 8, 11], "momentum": 2, "By": [2, 8, 11], "virtu": 2, "understood": 2, "good": 2, "eigenst": 2, "translation": 2, "non": [2, 8], "implic": 2, "xxz": [2, 12], "heisenberg": 2, "whose": 2, "j": [2, 11, 12], "x_": 2, "y_i": 2, "y_": 2, "delta": 2, "2s": 2, "satisfi": [2, 8], "mathfrak": [2, 11], "su": [2, 11], "relat": [2, 8], "a_i": [2, 11], "b_j": 2, "delta_": 2, "sum_c": 2, "varepsilon_": 2, "abc": 2, "c_i": [2, 8], "comput": [2, 3, 8, 11], "xxx": [2, 11], "y": 2, "neq": 2, "compon": 2, "same": [2, 8, 11], "mean": [2, 8], "full": [2, 11], "half": [2, 8], "3": [2, 8, 11], "see": [2, 8, 11], "wherea": [2, 11], "simeq": 2, "u": [2, 8, 11], "retain": 2, "If": [2, 8], "automat": [2, 8], "theta": 2, "interpret": [2, 8, 12], "rotat": 2, "around": 2, "axi": 2, "angl": 2, "quantiti": 2, "associ": [2, 8], "particular": [2, 8], "m_z": 2, "vec": 2, "cdot": 2, "label": 2, "motiv": 2, "gentli": 2, "form": [2, 8, 9, 11], "backbon": 2, "roughli": 2, "speak": 2, "g": 2, "rule": 2, "how": [2, 11], "compos": 2, "them": 2, "step": [2, 11], "time": [2, 8, 11, 12], "discret": 2, "consist": [2, 8], "carri": 2, "transform": [2, 8], "should": [2, 8, 11], "realli": [2, 11], "think": 2, "about": [2, 11], "These": [2, 8], "multipli": 2, "new": [2, 8], "result": [2, 8, 11], "ani": [2, 8], "next": [2, 8], "over": [2, 11], "theta_2": 2, "theta_1": 2, "lead": [2, 8], "what": 2, "A": [2, 8, 11, 12], "g_1": 2, "g_2": 2, "endow": 2, "There": 2, "exist": 2, "1g": 2, "g1": 2, "foral": 2, "note": [2, 8, 11, 12], "abelian": 2, "composit": 2, "still": [2, 8], "element": [2, 8, 11], "k": 2, "hk": 2, "gh": 2, "properti": [2, 8], "would": [2, 8], "like": 2, "formal": [2, 8], "undon": 2, "opposit": 2, "henc": [2, 11], "invers": 2, "gg": 2, "togeth": 2, "befor": [2, 8], "subgroup": 2, "suggest": [2, 11], "subset": 2, "li": 2, "heart": 2, "mathbb": [2, 8], "_2": [2, 11], "explan": [2, 11], "notat": [2, 3, 8], "undergo": 2, "symbol": 2, "express": [2, 8], "keep": 2, "than": [2, 8], "law": [2, 10], "denot": 2, "_n": 2, "addit": 2, "modulo": 2, "correspond": [2, 8, 11], "c": [2, 8, 11, 12], "right": [2, 8, 11], "encount": 2, "unimodular": 2, "matric": [2, 8, 11], "det": 2, "uu": 2, "similarli": [2, 8, 11], "geq": 2, "none": 2, "3d": 2, "real": [2, 11], "unit": [2, 8], "determin": [2, 8], "m": [2, 8, 12], "r": [2, 8, 11], "mm": 2, "tm": 2, "were": 2, "deal": [2, 8], "question": 2, "absenc": 2, "invert": 2, "linear": [2, 3, 11], "singular": [2, 8, 11], "structur": 2, "identifi": 2, "now": [2, 8, 11], "do": 2, "wonder": 2, "come": 2, "exactli": [2, 8, 11], "underli": 2, "idea": 2, "linearli": 2, "vector": [2, 8, 11], "space": [2, 8, 12], "immedi": 2, "rais": 2, "plethora": 2, "kind": 2, "construct": [2, 8], "ones": 2, "answer": 2, "sake": 2, "index": [2, 8], "x_g": 2, "x_gx_h": 2, "alwai": [2, 8], "matrix": [2, 11, 12], "call": [2, 8, 11], "dimens": [2, 8, 11], "probabl": 2, "x_0": 2, "x_1": 2, "inde": [2, 8], "x_1x_1": 2, "squar": 2, "regular": 2, "bar": [2, 8], "x_h": 2, "overlin": 2, "equiv": 2, "y_g": 2, "wai": [2, 11], "kroneck": 2, "check": [2, 8], "oplu": 2, "observ": [2, 8], "choic": [2, 8], "unitarili": 2, "ux_gu": 2, "independ": 2, "again": [2, 8, 11], "begin": 2, "pmatrix": 2, "end": [2, 8], "show": [2, 8], "crux": 2, "appropri": 2, "brought": 2, "block": [2, 9], "diagon": [2, 8], "fbox": 2, "1_g": 2, "2_g": 2, "vdot": 2, "ddot": 2, "amongst": 2, "themselv": 2, "irrep": 2, "shown": [2, 8], "equal": [2, 8], "alpha": [2, 8], "d_": 2, "respect": [2, 11], "One": [2, 8], "kei": 2, "decompos": 2, "sometim": 2, "fusion": 2, "clebsch": 2, "gordan": 2, "coeffici": 2, "explicitli": [2, 8], "due": 2, "schur": 2, "lemma": 2, "x_gy": 2, "yx_g": 2, "proport": 2, "pose": 2, "well": [2, 8, 11, 13], "summar": 2, "s_1": 2, "s_2": 2, "bigoplus_": 2, "been": [2, 13], "analyt": [2, 11], "low": [2, 8], "tabul": 2, "z_g": 2, "strong": 2, "didn": 2, "assum": [2, 8, 11], "argu": 2, "bring": [2, 8, 11], "appear": [2, 11], "write": [2, 8], "bigoplus_c": 2, "b_c": 2, "_c": 2, "decomposit": [2, 8, 11], "store": [2, 8], "much": [2, 8], "effici": [2, 8, 11], "track": 2, "tensorkit": [2, 8, 11, 13], "particularli": [2, 8], "suit": 2, "describ": 2, "paragraph": 2, "herebi": 2, "drastic": 2, "amount": 2, "memori": 2, "object": [2, 8, 11], "abl": 2, "manipul": 2, "exploit": 2, "maximum": 2, "rank": 2, "su\u2082spac": 2, "l": [2, 8, 11], "36m": [2, 8, 11], "1m": [2, 8, 11], "22m": [2, 8, 11], "39m": [2, 8, 11], "1minfo": [2, 8, 11], "39mprecompil": [2, 8, 11], "07d1fe3e": 2, "3e46": 2, "537d": 2, "9eac": 2, "e9e13d0d4cec": 2, "rep": 2, "su\u2082": 2, "essenti": 2, "copi": 2, "summand": 2, "want": 2, "ss": 2, "tensormap": [2, 8], "productspac": [2, 8], "data": 2, "fusiontre": 2, "fals": 2, "inspect": 2, "domain": 2, "codomain": 2, "assert": [2, 8], "dim": 2, "4": [2, 8, 11, 12], "sortedvectordict": 2, "su2irrep": 2, "float64": 2, "entri": 2, "fill": 2, "b": 2, "compat": 2, "cannot": 2, "fuse": 2, "third": 2, "final": [2, 8, 11], "5": [2, 8, 11], "0e": 2, "324": 2, "6": [2, 11], "91119e": 2, "310": 2, "323": 2, "four": 2, "conclud": 2, "global": [2, 8, 11], "familiar": 2, "gaug": [2, 10, 11], "ubiquit": 2, "phenomena": [2, 8], "thought": 2, "actual": 2, "each": [2, 8, 11], "redund": 2, "descript": 2, "nevertheless": [2, 11], "brief": [2, 11], "overview": [2, 8], "mostli": [2, 8], "neglect": 2, "spatial": 2, "reflect": 2, "don": [2, 11], "anymor": 2, "classif": 2, "notori": 2, "rich": 2, "beauti": [2, 11], "especi": 2, "higher": 2, "algorithm": [2, 8], "tremend": 2, "speedup": 2, "stabil": 2, "alreadi": [2, 8], "benefit": 2, "network": [2, 8, 12], "uniform": [2, 8, 11, 12], "repeat": [2, 8, 11], "indefinit": 2, "discoveri": 2, "topolog": [2, 12], "matter": 2, "anyon": [2, 12], "excit": 2, "grow": 2, "fascin": 2, "explor": 2, "categor": 2, "beyond": [2, 11], "scope": [2, 11], "intric": 2, "algebra": [2, 3, 11], "categori": 2, "specif": [2, 8], "chain": [2, 8, 11, 12], "ftl": [2, 12], "07": [2, 12], "allow": [2, 8], "storag": 2, "multi": [3, 8], "graphic": 3, "complex": [3, 8], "test": [4, 5], "basic": [6, 8], "updat": [6, 8, 11], "trotter": 6, "tebd": 6, "mp": [8, 10, 11], "excel": [8, 11], "review": [8, 12], "vhv19": [8, 11, 12], "provid": [8, 13], "thorough": 8, "technic": 8, "tangent": [8, 12], "method": [8, 12], "exposit": [8, 11], "supplement": 8, "work": [8, 11], "jl": [8, 13], "detail": [8, 11], "numer": 8, "routin": [8, 11], "julia": 8, "version": 8, "tutori": [8, 13], "readili": 8, "extend": 8, "physic": [8, 11, 12, 13], "hilbert": 8, "quantum": [8, 12, 13], "rangl": 8, "boldsymbol": 8, "v": 8, "_l": 8, "prod_": 8, "s_m": 8, "_r": 8, "altern": 8, "view": 8, "three": 8, "indic": 8, "bond": [8, 11], "control": 8, "9": [8, 11], "arbitrari": [8, 11], "accuraci": 8, "certain": 8, "class": [8, 11], "gap": [8, 11], "accur": 8, "approxim": [8, 11], "smaller": 8, "eq": 8, "never": 8, "safe": 8, "ignor": 8, "bulk": 8, "faithfulli": 8, "captur": 8, "natur": 8, "impos": 8, "transat": 8, "diagramat": 8, "instead": 8, "cell": 8, "size": [8, 11], "techniqu": [8, 11], "appli": 8, "just": 8, "central": 8, "unform": 8, "calcul": 8, "transfer": 8, "leg": 8, "leftarrow": [8, 11], "complet": 8, "posit": 8, "character": [8, 11], "length": [8, 11], "evalu": [8, 11], "norm": [8, 11], "contract": 8, "noth": 8, "spectral": 8, "th": 8, "power": 8, "fix": [8, 13], "point": [8, 13], "largest": 8, "magnitud": 8, "lambda_0": 8, "lambda_i": 8, "remain": 8, "mangitud": 8, "projector": 8, "To": 8, "ensur": 8, "properli": 8, "rescal": 8, "requir": 8, "trace": 8, "With": 8, "overlap": 8, "between": 8, "effect": [8, 11], "choos": [8, 11], "langl": 8, "middl": 8, "suppos": 8, "wish": 8, "extens": 8, "o_n": 8, "dictat": 8, "everyth": 8, "look": 8, "beta": 8, "bra": 8, "beta_m": 8, "alpha_n": 8, "abritrari": 8, "locat": 8, "becaus": [8, 11], "insert": 8, "learn": 8, "disconnect": 8, "rest": 8, "exponenti": 8, "decai": 8, "connect": 8, "why": 8, "gener": [8, 11], "harder": 8, "xi": 8, "lambda_1": 8, "log": 8, "lambda_": 8, "mathrm": 8, "max": [8, 11], "sublead": 8, "typic": 8, "focuss": 8, "symmetri": [8, 13], "sector": 8, "target": 8, "exit": 8, "plai": 8, "crucial": 8, "role": 8, "scale": 8, "rcc18": [8, 12], "despit": 8, "ansatz": [8, 10, 11], "inher": 8, "convers": 8, "true": 8, "mai": [8, 12], "easili": [8, 11], "seen": 8, "parametr": 8, "canon": [8, 11], "start": [8, 11, 13], "orthonorm": [8, 11], "a_l": [8, 11], "find": [8, 11], "iter": [8, 11], "procedur": 8, "qr": 8, "docomposit": 8, "initi": [8, 11], "guess": [8, 11], "repeatedli": [8, 11], "perform": 8, "bound": 8, "converg": [8, 11], "room": 8, "a_r": [8, 11], "found": 8, "similar": 8, "mix": [8, 11], "center": [8, 11], "a_c": [8, 11], "obtain": [8, 11], "contrast": 8, "origin": 8, "commonli": 8, "lr": 8, "therebi": 8, "freeli": 8, "move": [8, 11], "through": [8, 11], "link": 8, "usv": 8, "absorb": 8, "definit": 8, "residu": 8, "entir": 8, "rm": 8, "tr": 8, "arriv": 8, "straightforwardli": 8, "down": [8, 11], "schmidt": 8, "across": 8, "i_l": 8, "i_r": 8, "bipartit": 8, "enabl": [8, 11], "sum": 8, "optim": [8, 11], "sens": 8, "maxim": 8, "column": 8, "isometri": [8, 11], "correspondingli": 8, "guarante": 8, "lower": 8, "variat": [8, 11, 12], "cost": 8, "tild": [8, 11], "packag": 8, "mani": [8, 13], "tool": 8, "go": 8, "aspect": 8, "specifi": [8, 11], "virtual": 8, "standard": 8, "complexspac": 8, "\u2102": [8, 11], "bb1c41ca": 8, "d63c": 8, "52ed": 8, "829e": 8, "0820dda26502": 8, "cr": 8, "al": 8, "linearalgebra": 8, "ac": 8, "9999999999999999": 8, "verifi": 8, "diagram": 8, "tensoroper": [8, 13], "macro": 8, "al_id": 8, "conj": 8, "ar_id": 8, "id": 8, "lh": 8, "rh": 8, "randn": 8, "expectation_valu": [8, 11], "complexf64": [8, 11], "9253889853121735": 8, "05692568523573335im": 8, "encod": [8, 11], "correlation_length": 8, "5050041892918671": 8, "export": 8, "varieti": 8, "statist": [9, 13], "mpo": [9, 11], "jordan": 9, "polynomi": 10, "entangl": 10, "area": 10, "correl": [10, 12], "manner": 11, "boil": 11, "min_": 11, "simplic": 11, "consider": 11, "represent": 11, "interact": [11, 12], "rang": 11, "densiti": 11, "renorm": 11, "group": 11, "direct": 11, "thermodynam": 11, "break": 11, "bc17": [11, 12], "lectur": [11, 12], "random": 11, "tri": 11, "sequenti": 11, "sweep": 11, "until": 11, "reach": 11, "bit": 11, "a_1": 11, "a_2": 11, "seemingli": 11, "daunt": 11, "problem": 11, "make": 11, "those": 11, "denomin": 11, "mathcal": 11, "h_i": 11, "solv": 11, "_i": 11, "smallest": 11, "back": 11, "forth": 11, "At": 11, "reus": 11, "manifestli": 11, "suffer": 11, "artefact": 11, "surprisingli": 11, "proven": 11, "success": 11, "variation": 11, "mpskit": 11, "mpskitmodel": 11, "z_j": 11, "h_z": 11, "free": 11, "paramet": 11, "usual": 11, "factor": 11, "16": [11, 12], "12": 11, "open": 11, "stick": 11, "default": 11, "hand": [11, 12, 13], "transverse_field_is": 11, "summon": 11, "\u03c8": 11, "finitemp": 11, "\u03c8\u2080": 11, "_": 11, "find_groundst": 11, "ca635005": 11, "6f8c": 11, "4cd1": 11, "b51d": 11, "8491250ef2ab": 11, "39miteraton": 11, "error": 11, "012416077282990786": 11, "096580055657124e": 11, "99774835027203e": 11, "7": [11, 12], "452126130298275e": 11, "21619849666518e": 11, "8": [11, 12], "2253741288248538e": 11, "6585061708822792e": 11, "382970261218168e": 11, "217298314883982e": 11, "10": [11, 12], "965743561640723e": 11, "401238662513454e": 11, "11": [11, 12], "7860459770692195e": 11, "0489342641394526e": 11, "13": 11, "950846732206643e": 11, "14": 11, "487316766748785e": 11, "15": 11, "0199505177377904e": 11, "311759751328411e": 11, "directli": 11, "unbound": 11, "diagrammat": 11, "f_l": 11, "f_r": 11, "obei": 11, "normal": 11, "offer": 11, "advantag": 11, "design": 11, "reli": 11, "individu": 11, "intial": 11, "solut": 11, "equat": 11, "h_": 11, "h_c": 11, "a_lc": 11, "ca_r": 11, "deriv": 11, "minimum": 11, "manifold": 11, "last": 11, "intertwin": 11, "explain": 11, "chose": 11, "toler": 11, "eta": 11, "arnoldi": 11, "yield": 11, "epsilon_l": 11, "epsilon_r": 11, "isometr": 11, "comment": 11, "soltuion": 11, "epsilon": 11, "a_": 11, "u_lv_l": 11, "u_l": 11, "v_l": 11, "aris": 11, "sigma_lv_l": 11, "u_rv_r": 11, "u_r": 11, "sigma_rv_r": 11, "approach": 11, "best": 11, "exact": 11, "s_c": 11, "s_lc": 11, "ca": 11, "s_r": 11, "sigma_": 11, "precis": [11, 12], "arithmet": 11, "u_": 11, "v_": 11, "poor": 11, "robust": 11, "l_": 11, "l_c": 11, "qquad": 11, "r_c": 11, "r_": 11, "polar": 11, "l_cp": 11, "r_cu": 11, "demonstr": 11, "estim": 11, "cdot1": 11, "heisenberg_xyz": 11, "infinitemp": 11, "env": 11, "39mvump": 11, "galerkin": 11, "0836871681973854": 11, "43919484489605565": 11, "07043631724717614": 11, "02970868195222154": 11, "021688260172916175": 11, "009884968971296924": 11, "0038863737616032967": 11, "0011919008029591459": 11, "0004028626930456042": 11, "00013994254153142723": 11, "9484740738776824e": 11, "7949596092087722e": 11, "51561664938521e": 11, "4134476269623677e": 11, "889508439248023e": 11, "3308476656957004e": 11, "17": 11, "238811952663072e": 11, "18": 11, "676241715984628e": 11, "19": 11, "752457723175491e": 11, "20": 11, "651857635716631e": 11, "21": 11, "509222783987944e": 11, "22": [11, 12], "568143817948816e": 11, "23": 11, "6307171028007104e": 11, "24": 11, "389956720484358e": 11, "25": 11, "302789585332627e": 11, "26": 11, "037506472937299e": 11, "27": 11, "811911733229695e": 11, "28": 11, "018643441914258e": 11, "29": 11, "1819841791771657e": 11, "30": 11, "842991318561591e": 11, "periodicarrai": 11, "40138064351825": 11, "3877787807814457e": 11, "17im": 11, "compar": 11, "quasi": 11, "401": 11, "484": 11, "038": 11, "971": 11, "hco": [11, 12], "decim": 11, "jacob": 12, "bridgeman": 12, "christoph": 12, "chubb": 12, "wave": 12, "danc": 12, "introductori": 12, "journal": 12, "theoret": 12, "50": 12, "223001": 12, "2017": 12, "url": 12, "http": 12, "dx": 12, "doi": 12, "org": 12, "1088": 12, "1751": 12, "8121": 12, "aa6dc3": 12, "adrian": 12, "feiguin": 12, "simon": 12, "trebst": 12, "andrea": 12, "ww": 12, "ludwig": 12, "matthia": 12, "troyer": 12, "alexei": 12, "kitaev": 12, "zhenghan": 12, "wang": 12, "michael": 12, "freedman": 12, "liquid": 12, "golden": 12, "letter": 12, "98": 12, "160409": 12, "2007": 12, "jutho": 12, "haegeman": 12, "ignacio": 12, "cirac": 12, "tobia": 12, "osborn": 12, "iztok": 12, "\u017e": 12, "orn": 12, "henri": 12, "verscheld": 12, "frank": 12, "verstraet": 12, "107": 12, "070601": 12, "2011": 12, "marek": 12, "ram": 12, "piotr": 12, "czarnik": 12, "lukasz": 12, "cincio": 12, "extrapol": 12, "function": 12, "asymptot": 12, "applic": [12, 13], "bose": 12, "hubbard": 12, "041033": 12, "novemb": 12, "2018": 12, "arxiv": 12, "ab": 12, "1801": 12, "08554": 12, "1103": 12, "physrevx": 12, "lauren": 12, "vanderstraeten": 12, "scipost": 12, "page": 12, "007": 12, "januari": 12, "2019": 12, "scipostphyslectnot": 12, "1810": 12, "07006": 12, "21468": 12, "seri": 13, "theori": 13, "aim": 13, "practis": 13, "code": 13, "exampl": 13, "showcas": 13, "softwar": 13, "librari": 13, "develop": 13, "bodi": 13}, "objects": {}, "objtypes": {}, "objnames": {}, "titleterms": {"get": 0, "start": 0, "julia": [0, 13], "quantum": [1, 2], "mani": [1, 2], "bodi": [1, 2], "theori": [1, 2, 3], "symmetri": 2, "physic": 2, "exampl": [2, 8, 11], "applic": [2, 7], "break": 2, "order": 2, "paramet": 2, "phase": 2, "noether": 2, "conserv": 2, "quantitit": 2, "group": 2, "represent": [2, 8], "definit": 2, "complex": 2, "conjug": 2, "tensor": [2, 3, 13], "product": [2, 8, 9, 10, 13], "direct": 2, "sum": 2, "irreduc": 2, "symmetr": 2, "outlook": 2, "gener": 2, "network": [3, 13], "tensorkit": 4, "jl": [4, 5], "tensoroper": 5, "algorithm": [6, 11, 13], "infinit": 8, "matrix": [8, 9, 10, 13], "state": [8, 10, 13], "content": 8, "thermodynam": 8, "limit": 8, "normal": 8, "expect": 8, "valu": 8, "correl": 8, "function": 8, "gaug": 8, "revisit": 8, "entangl": 8, "entropi": 8, "truncat": 8, "code": 8, "mpskit": 8, "infinitemp": 8, "oper": 9, "fix": 11, "point": 11, "dmrg": 11, "vump": 11, "refer": 12, "method": 13, "introduct": 13, "other": 13}, "envversion": {"sphinx.domains.c": 2, "sphinx.domains.changeset": 1, "sphinx.domains.citation": 1, "sphinx.domains.cpp": 6, "sphinx.domains.index": 1, "sphinx.domains.javascript": 2, "sphinx.domains.math": 2, "sphinx.domains.python": 3, "sphinx.domains.rst": 2, "sphinx.domains.std": 2, "sphinx.ext.intersphinx": 1, "sphinxcontrib.bibtex": 9, "sphinx": 56}})