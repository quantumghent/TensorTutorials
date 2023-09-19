Search.setIndex({"docnames": ["1-Introduction/LinearAlgebra", "1-Introduction/QuantumManyBody", "1-Introduction/Software", "1-Introduction/Symmetries", "1-Introduction/TensorNetworks", "2-MatrixProductStates/Algorithms", "2-MatrixProductStates/Algorithms-backup", "2-MatrixProductStates/Applications", "2-MatrixProductStates/InfiniteMPS", "2-MatrixProductStates/MatrixProductOperators", "2-MatrixProductStates/MatrixProductStates", "3-Algorithms/FixedpointAlgorithms", "References", "intro"], "filenames": ["1-Introduction/LinearAlgebra.md", "1-Introduction/QuantumManyBody.md", "1-Introduction/Software.md", "1-Introduction/Symmetries.md", "1-Introduction/TensorNetworks.md", "2-MatrixProductStates/Algorithms.md", "2-MatrixProductStates/Algorithms-backup.md", "2-MatrixProductStates/Applications.md", "2-MatrixProductStates/InfiniteMPS.md", "2-MatrixProductStates/MatrixProductOperators.md", "2-MatrixProductStates/MatrixProductStates.md", "3-Algorithms/FixedpointAlgorithms.md", "References.md", "intro.md"], "titles": ["<span class=\"section-number\">3. </span>(Multi-) Linear Algebra", "<span class=\"section-number\">1. </span>Quantum Many-Body Theory", "<span class=\"section-number\">2. </span>Getting Started with Numerics", "<span class=\"section-number\">5. </span>Symmetries in Quantum Many-Body Physics", "<span class=\"section-number\">4. </span>Tensor Network Theory", "<span class=\"section-number\">9. </span>Tensor Network Algorithms", "Tensor Network Algorithms", "<span class=\"section-number\">10. </span>Applications", "<span class=\"section-number\">8. </span>Infinite Matrix Product States", "<span class=\"section-number\">7. </span>Matrix Product Operators", "<span class=\"section-number\">6. </span>Matrix Product States", "<span class=\"section-number\">11. </span>Fixed-Point algorithms", "<span class=\"section-number\">12. </span>References", "Tensor Network Methods with Julia"], "terms": {"thi": [0, 2, 3, 4, 5, 6, 8, 11, 13], "lectur": [0, 4, 5, 11, 12, 13], "cover": 0, "some": [0, 2, 3, 4, 5, 6, 8], "basic": [0, 4, 6, 8], "concept": [0, 3, 4, 8], "oper": [0, 3, 5, 8, 13], "serv": [0, 4, 6], "foundat": [0, 4], "most": [0, 2, 3, 4, 8, 11], "what": [0, 3, 4, 6], "follow": [0, 3, 4, 5, 6, 8, 11], "The": [0, 2, 3, 4, 5, 8, 11, 13], "goal": [0, 3], "provid": [0, 4, 8, 13], "intuit": [0, 3, 4, 8], "understand": [0, 4, 5], "without": [0, 3, 8, 11], "insist": 0, "too": 0, "much": [0, 2, 3, 5, 6, 8], "mathemat": [0, 3, 12], "rigour": 0, "import": [0, 3, 4, 5, 6], "introduc": [0, 3, 4, 5, 6, 8, 11], "defin": [0, 3, 4, 6, 8, 11], "resort": [0, 5, 6], "usual": [0, 6, 11], "definit": [0, 4, 8], "which": [0, 2, 3, 4, 5, 6, 8, 11], "veri": [0, 2, 3, 4, 5, 6, 8, 11], "simultan": [0, 3], "also": [0, 2, 3, 4, 5, 6, 8], "showcas": [0, 4, 13], "featur": [0, 3, 4], "tensorkit": [0, 2, 3, 4, 8, 11], "jl": [0, 2, 4, 8], "julia": [0, 4, 8], "packag": [0, 4, 8], "extrem": [0, 3], "well": [0, 2, 3, 4, 6, 8, 11, 13], "suit": [0, 3], "demonstr": [0, 11], "ar": [0, 2, 3, 4, 5, 6, 8, 11], "discuss": [0, 3, 4, 5, 6, 8, 11], "us": [0, 2, 3, 4, 5, 6, 8, 11], "36m": [0, 8, 11], "1m": [0, 8, 11], "22m": [0, 4, 8, 11], "39m": [0, 4, 8, 11], "1minfo": [0, 8, 11], "39mprecompil": [0, 8, 11], "07d1fe3e": 0, "3e46": 0, "537d": 0, "9eac": 0, "e9e13d0d4cec": 0, "befor": [0, 3, 5, 8], "network": [0, 3, 8, 12], "necessari": [0, 4], "furthermor": [0, 4], "realli": [0, 3, 4, 11], "instruct": 0, "reiter": 0, "case": [0, 3, 4, 5, 8], "noth": [0, 4, 8], "specif": [0, 3, 4, 5, 6, 8], "In": [0, 2, 3, 4, 5, 6, 8, 11], "fact": [0, 3, 4, 5], "mani": [0, 2, 4, 5, 6, 8, 13], "idea": [0, 3, 4, 5], "term": [0, 3, 5, 6, 8], "think": [0, 3], "thought": [0, 3, 4], "from": [0, 2, 3, 4, 6, 8, 11], "viewpoint": 0, "comput": [0, 2, 3, 4, 5, 6, 8, 11, 12], "where": [0, 3, 4, 5, 6, 8, 11], "thei": [0, 3, 4, 5, 6, 8], "repres": [0, 3, 4, 8, 11], "regular": [0, 3], "one": [0, 3, 4, 5, 6, 8, 11], "two": [0, 3, 4, 5, 6, 8, 11, 12], "dimension": [0, 3, 4, 8, 11], "arrai": [0, 4], "either": [0, 4, 5, 6], "real": [0, 3, 5, 6, 11], "complex": [0, 8], "number": [0, 3, 4, 5, 6, 8, 11], "nevertheless": [0, 3, 4, 11], "can": [0, 2, 3, 4, 5, 6, 8, 11], "readili": [0, 8], "gener": [0, 2, 4, 5, 6, 8, 11], "arbitrari": [0, 8, 11], "space": [0, 3, 4, 6, 8, 12], "an": [0, 2, 3, 4, 5, 6, 8, 9, 11, 12], "object": [0, 3, 6, 8, 11], "describ": [0, 3, 6], "list": [0, 2], "correspond": [0, 3, 4, 6, 8, 11], "compon": [0, 3, 4], "basi": [0, 3, 4], "For": [0, 3, 4, 5, 6, 8, 11], "exampl": [0, 2, 4, 6, 13], "its": [0, 3, 4, 5, 6, 8], "form": [0, 3, 4, 5, 6, 8, 9, 11], "vec": [0, 3], "v": [0, 4, 8, 12], "left": [0, 3, 4, 5, 6, 8, 11], "v_1": 0, "v_2": 0, "right": [0, 3, 4, 5, 6, 8, 11], "t": [0, 3, 5, 6, 11, 12], "As": [0, 3, 4, 5, 6, 8, 11], "remind": 0, "properti": [0, 3, 4, 5, 6, 8], "make": [0, 4, 5, 6, 11], "sure": 0, "ad": [0, 3], "togeth": [0, 3], "i": [0, 3, 4, 5, 6, 8, 11, 12], "e": [0, 3, 4, 5, 6, 8, 12], "w": [0, 4], "multipli": [0, 3], "scalar": [0, 4], "alpha": [0, 3, 8], "These": [0, 2, 3, 4, 8], "behav": [0, 3], "expect": [0, 3, 5, 6, 9, 10], "notion": [0, 3, 4], "associ": [0, 3, 8], "commut": [0, 3, 5, 6], "distribut": 0, "given": [0, 2, 3, 4, 5, 6, 8, 11], "necessarili": [0, 3, 4, 6], "distinct": [0, 3, 4], "possibl": [0, 4, 5, 6], "between": [0, 4, 8], "them": [0, 3, 4, 5], "just": [0, 4, 6, 8], "function": [0, 4, 5, 6, 12], "preserv": [0, 4], "structur": [0, 3, 4, 5, 6], "other": [0, 2, 3, 4, 5, 6], "word": [0, 3, 4, 5], "A": [0, 3, 4, 5, 6, 8, 11, 12], "colon": 0, "anoth": [0, 3, 5], "becaus": [0, 4, 8, 11], "requir": [0, 4, 5, 6, 8], "complet": [0, 8], "determin": [0, 3, 4, 5, 8, 12], "action": [0, 3], "lead": [0, 3, 4, 8], "natur": [0, 4, 6, 8], "wai": [0, 3, 4, 5, 6, 11], "matrix": [0, 3, 4, 5, 6, 11, 12], "consid": [0, 3, 4, 5, 6, 8, 11], "construct": [0, 3, 4, 6, 8], "v_i": 0, "w_i": 0, "begin": [0, 3, 4, 5], "rcl": 0, "leftarrow": [0, 8, 11], "equiv": [0, 3], "sum_j": 0, "a_": [0, 11], "ij": [0, 5], "v_j": 0, "end": [0, 3, 4, 5, 6, 8], "base": [0, 3, 4, 6, 8], "abstract": 0, "concret": [0, 3, 4], "particular": [0, 3, 4, 6, 8], "column": [0, 4, 8], "label": [0, 3, 4], "input": 0, "call": [0, 3, 4, 5, 6, 8, 11], "domain": [0, 3, 4], "while": [0, 3, 4, 5, 6, 8], "row": [0, 6], "output": 0, "codomain": [0, 3, 4], "context": [0, 1, 3, 4], "we": [0, 2, 3, 4, 5, 6, 8, 11], "creat": [0, 4], "through": [0, 4, 8, 11], "syntax": 0, "close": [0, 3, 5, 6, 11], "\u2102": [0, 4, 8, 11], "2": [0, 3, 4, 5, 6, 8, 11], "type": 0, "bbc": 0, "tab": 0, "complexspac": [0, 4, 8], "3": [0, 3, 4, 8, 11], "equival": [0, 3, 4], "tensormap": [0, 3, 4, 8], "rand": [0, 4], "float64": [0, 3, 4], "1": [0, 3, 4, 5, 6, 8, 11], "true": [0, 4, 8], "same": [0, 3, 4, 5, 6, 8, 11], "logic": 0, "abov": [0, 3, 4, 5, 6, 8, 11], "combin": [0, 4], "new": [0, 3, 4, 8], "otim": [0, 3, 6, 8], "origin": [0, 6, 8], "equal": [0, 3, 4, 8], "hold": 0, "all": [0, 3, 4, 5, 6, 8, 11], "lambda": [0, 3, 4, 6], "\u03bb": 0, "equip": 0, "canon": [0, 8, 11], "take": [0, 3, 4, 5, 6, 8, 11], "w_j": 0, "respect": [0, 3, 4, 6, 11], "when": [0, 3, 4, 8, 11], "how": [0, 3, 4, 5, 6, 11], "written": [0, 4, 13], "sum_": [0, 3, 5, 6, 8, 11], "i_1": [0, 4], "i_2": [0, 4], "t_": 0, "i_1i_2": 0, "v_": [0, 11], "w_": 0, "shorthand": 0, "extract": 0, "stridedview": 0, "typeof": 0, "ident": [0, 3, 8], "0": [0, 3, 4, 5, 6, 8, 11], "449467": 0, "684866": 0, "773791": 0, "442901": 0, "86258": 0, "143055": 0, "here": [0, 2, 3, 4, 5, 6, 8, 11], "tent": 0, "name": [0, 3, 4, 5, 11], "wa": [0, 3, 4, 11], "denot": [0, 3, 5], "induc": 0, "more": [0, 2, 3, 4, 5, 6, 8, 11], "common": [0, 4, 6], "express": [0, 3, 4, 8], "reshap": [0, 4], "th": [0, 8], "built": [0, 4], "similarli": [0, 3, 4, 8, 11], "than": [0, 3, 5, 8], "final": [0, 3, 4, 8, 11], "element": [0, 3, 4, 8, 11], "up": [0, 3, 4, 5, 6, 8, 11], "addition": [0, 2, 4], "laid": 0, "out": [0, 2, 3, 4, 5, 6, 11], "slight": 0, "misus": 0, "terminolog": 0, "indic": [0, 8], "cartesian": 0, "cdot": [0, 3, 4], "i_n": 0, "trick": [0, 5], "allow": [0, 2, 3, 4, 8], "reinterpret": [0, 4], "vice": [0, 3], "versa": [0, 3], "linearindic": 0, "tupl": [0, 4], "unitrang": [0, 4], "int64": [0, 4], "5": [0, 3, 4, 8, 11], "4": [0, 3, 4, 8, 11, 12], "6": [0, 4, 11], "collect": 0, "cartesianindic": 0, "forc": 0, "print": 0, "cartesianindex": 0, "due": [0, 3, 4, 5, 6], "itself": [0, 3, 6, 11], "again": [0, 2, 3, 4, 5, 8, 11], "keep": [0, 2, 3], "mind": 0, "now": [0, 3, 5, 6, 8, 11], "howev": [0, 3, 4, 5, 6, 8, 11], "themselv": [0, 3], "compris": 0, "If": [0, 3, 4, 5, 6, 8], "order": [0, 2, 5, 6], "establish": [0, 4], "w_1": 0, "w_2": 0, "w_m": 0, "v_n": 0, "mapsto": 0, "j_1": 0, "j_2": 0, "j_n": 0, "i_m": 0, "n": [0, 3, 4, 5, 6, 8, 12], "j": [0, 3, 4, 5, 11, 12], "v1": 0, "v2": 0, "w1": 0, "w2": 0, "attent": 0, "reader": [0, 3], "might": [0, 3, 6, 11], "have": [0, 2, 3, 4, 5, 6, 8, 11, 13], "alreadi": [0, 2, 3, 4, 5, 8], "note": [0, 3, 4, 8, 11, 12], "strongli": 0, "resembl": [0, 6], "coincid": 0, "easili": [0, 4, 8, 11], "identifi": [0, 3], "identif": [0, 12], "isomorph": [0, 4], "cong": 0, "b": [0, 3, 4, 5, 6, 12], "finit": [0, 3, 8, 11], "addit": [0, 2, 3, 4, 5], "trivial": [0, 3, 4, 8], "choic": [0, 3, 4, 8], "uniqu": [0, 3, 4, 8], "differ": [0, 2, 3, 4, 8, 11], "major": [0, 4], "set": [0, 2, 3, 4, 6, 11], "mean": [0, 3, 5, 6, 8], "still": [0, 3, 4, 5, 6, 8], "relat": [0, 3, 8], "entir": [0, 5, 6, 8], "summar": [0, 3], "along": [0, 4], "constitu": [0, 4, 5], "thu": [0, 3, 4, 11], "like": [0, 3, 4, 6], "lift": 0, "facet": 0, "introduct": [1, 3, 4, 6], "histori": [1, 2], "purpos": [1, 2, 5, 6], "relev": [1, 2, 4, 5, 6], "second": [1, 3, 4, 6, 8, 11], "quantiz": 1, "tensor": [1, 8, 11, 12], "product": [1, 5, 6, 11, 12], "classic": [1, 3, 7, 12], "On": [2, 6], "page": [2, 12], "link": [2, 6, 8, 12], "inform": [2, 4], "point": [2, 4, 6, 8, 13], "refer": [2, 3, 6, 8, 13], "program": 2, "languag": [2, 4], "resourc": 2, "learn": [2, 4, 8], "about": [2, 3, 11], "tool": [2, 4, 8], "develop": [2, 4, 5, 6, 13], "sometim": [2, 3], "field": [2, 3, 4, 6, 7, 11], "manag": [2, 5, 6], "track": [2, 3], "chang": [2, 3, 4, 6], "made": [2, 4, 5], "project": [2, 4, 5], "s": [2, 3, 4, 6, 8, 11, 12], "sourc": 2, "code": [2, 13], "document": 2, "ani": [2, 3, 4, 5, 6, 8], "file": 2, "It": [2, 3, 4, 5, 6, 11, 13], "multipl": [2, 3, 4], "contributor": 2, "work": [2, 5, 6, 8, 11], "collabor": 2, "facilit": [2, 4], "organ": 2, "synchron": 2, "popular": 2, "system": [2, 3, 4, 8, 11, 12], "git": 2, "free": [2, 11], "linu": 2, "torvald": 2, "2005": 2, "ha": [2, 3, 4, 5, 6, 8, 11], "becom": [2, 3, 4, 5, 11], "de": [2, 13], "facto": 2, "standard": [2, 8], "industri": 2, "avail": 2, "offici": 2, "websit": 2, "book": 2, "pro": 2, "good": [2, 3, 4], "place": [2, 3, 4, 5, 6, 8], "full": [2, 3, 4, 5, 6, 11], "exposit": [2, 8, 11], "found": [2, 4, 5, 8], "There": [2, 3, 5, 6], "tutori": [2, 8, 13], "topic": [2, 3], "activ": 2, "forum": 2, "ask": 2, "question": [2, 3, 5], "slack": 2, "channel": 2, "stack": 2, "overflow": 2, "open": [2, 4, 11], "commun": 2, "typic": [2, 4, 8], "own": [2, 6], "host": [2, 6], "github": 2, "incomplet": 2, "cours": [2, 3, 4, 6, 12], "below": [2, 3, 11], "tensoroper": [2, 4, 8], "krylovkit": 2, "optimkit": 2, "mpskit": [2, 11], "pepskit": 2, "check": [2, 3, 8], "our": [2, 3, 5, 6, 8, 11], "repositori": 2, "librari": [2, 13], "quantum": [2, 4, 8, 12, 13], "physic": [2, 4, 5, 6, 8, 11, 12, 13], "research": [2, 4], "you": [2, 3, 11], "find": [2, 4, 5, 6, 8, 11], "itensor": 2, "c": [2, 3, 4, 6, 8, 11, 12], "calcul": [2, 8], "tenpi": 2, "python": 2, "quimb": 2, "bodi": [2, 4, 6, 13], "section": [3, 4, 6, 8, 11], "give": [3, 4, 6, 8, 13], "gentl": 3, "framework": [3, 6, 8], "least": [3, 4], "restrict": [3, 8, 11], "illustr": [3, 4, 5, 6, 11, 13], "rather": [3, 5, 11], "first": [3, 4, 5, 6, 8, 11], "coupl": [3, 6], "model": [3, 7, 11, 12], "gradual": 3, "build": [3, 5], "finish": 3, "present": 3, "goe": 3, "sai": [3, 6], "onli": [3, 4, 5, 6, 8, 11], "scratch": 3, "surfac": [3, 4], "vast": 3, "interest": [3, 6], "immens": 3, "literatur": 3, "special": 3, "recal": [3, 8], "transvers": [3, 7, 11], "Ising": [3, 7, 11], "Its": 3, "degre": 3, "freedom": [3, 4, 8, 10], "qubit": 3, "lattic": [3, 5, 6, 8, 11, 12], "hamiltonian": [3, 6, 11], "read": 3, "h": [3, 5, 6, 11, 12], "sigma": [3, 4], "z_i": [3, 11], "z_": 3, "h_x": [3, 11], "sum_i": [3, 5, 8, 11], "x_i": [3, 11], "let": [3, 5, 8, 11], "simpli": [3, 4, 5, 6, 11], "period": [3, 5, 6], "boundari": [3, 5, 6, 8, 11], "condit": [3, 5, 6, 8, 11], "besid": 3, "obviou": 3, "translat": [3, 6, 8, 11], "invari": [3, 8, 11], "under": [3, 6, 8, 11], "flip": [3, 4], "spin": [3, 6, 11], "z": [3, 6, 8, 11, 12], "pauli": [3, 11], "ket": [3, 5, 6, 8, 11], "uparrow": 3, "leftrightarrow": 3, "downarrow": 3, "That": [3, 6], "constitut": 3, "clear": [3, 4, 6, 11], "energi": [3, 5, 6, 8, 11], "depend": [3, 6, 8, 12], "neighbour": [3, 4, 5], "being": [3, 4, 11], "anti": 3, "align": [3, 5], "clearli": [3, 8], "extern": [3, 6], "magnet": 3, "orthogon": [3, 4, 8, 11], "implement": [3, 4, 5, 8, 11], "correctli": 3, "unitari": [3, 4, 8, 11], "p": [3, 4, 11], "bigotimes_i": 3, "notic": [3, 11], "accord": [3, 6, 11], "twice": 3, "leav": [3, 4, 8], "untouch": 3, "dagger": [3, 4, 8, 11], "hp": 3, "everi": [3, 6, 8, 11], "even": [3, 4, 5, 6, 11], "though": [3, 11], "regardless": 3, "valu": [3, 5, 6, 9, 10, 11], "know": 3, "previou": [3, 8, 11], "ground": [3, 6, 8, 11], "state": [3, 4, 5, 6, 11, 12], "subspac": 3, "phenomenon": 3, "known": [3, 4, 5, 6, 11], "spontanu": 3, "ssb": 3, "short": 3, "investig": 3, "vanish": 3, "infinit": [3, 12, 13], "rightarrow": 3, "infti": [3, 5, 8], "effectli": 3, "reduc": [3, 4, 8, 11], "paramagnet": 3, "psi_": 3, "frac": [3, 5, 6, 8, 11], "sqrt": [3, 8], "eigenvalu": [3, 5, 8, 11], "eigenvector": [3, 4, 5, 8, 11], "x": [3, 4, 11, 12], "reason": [3, 4, 6, 8], "mention": [3, 5, 11], "disord": 3, "minim": [3, 4, 5, 11], "ferromagnet": 3, "obvious": [3, 11], "contrari": [3, 11], "span": 3, "get": [3, 6, 13], "map": [3, 4, 5, 6, 8], "onto": [3, 8], "broken": 3, "sinc": [3, 5, 6, 8, 11], "degeneraci": 3, "integ": [3, 4], "smoothli": 3, "slowli": 3, "turn": [3, 5, 6, 11], "therefor": [3, 5, 6, 8], "small": [3, 5, 11], "larg": [3, 4, 5, 6, 11], "said": 3, "belong": 3, "transit": 3, "abruptli": 3, "happen": [3, 4], "critic": [3, 8, 11], "inspir": 3, "credo": 3, "local": [3, 5, 6, 8, 11], "probe": 3, "wit": 3, "magnetis": 3, "site": [3, 5, 6, 8, 11], "o": [3, 5, 6, 8], "anticommut": 3, "op": 3, "expecti": 3, "braket": [3, 11], "latter": [3, 4], "could": [3, 5], "chosen": [3, 4, 6], "so": [3, 4, 5, 6, 8, 11], "seem": 3, "ill": [3, 5, 6], "remedi": 3, "sign": 3, "select": 3, "after": [3, 4, 5, 6], "limit": [3, 5, 11], "taken": [3, 5, 6, 11], "synopsi": 3, "singl": [3, 4, 6, 8], "particl": [3, 5], "mechan": [3, 9, 13], "part": [3, 5, 6, 8, 11], "pattern": 3, "hallmark": 3, "doe": [3, 4, 5, 6, 11], "paradigm": 3, "classifi": 3, "principl": [3, 12], "put": [3, 5, 6], "forward": 3, "landau": 3, "bear": 3, "hi": 3, "rememb": 3, "theorem": 3, "continu": [3, 4, 5, 6], "often": [3, 4, 6], "via": [3, 4], "lagrangian": 3, "rise": [3, 6, 8], "current": [3, 11], "almost": 3, "impli": [3, 4, 8], "d": [3, 4, 6, 8, 11], "dt": 3, "psi": [3, 5, 6, 8, 11], "proof": [3, 4], "simpl": [3, 4, 5, 6, 8, 11], "exercis": [3, 4, 5], "simplest": [3, 4, 5, 6], "hamiltonain": 3, "consequ": 3, "total": [3, 6, 11], "act": [3, 5, 6, 8, 11], "o_i": 3, "o_it": 3, "o_": [3, 6], "exp": 3, "pi": [3, 12], "ip": 3, "momentum": 3, "By": [3, 8, 11], "virtu": 3, "understood": 3, "eigenst": [3, 5, 6], "translation": 3, "non": [3, 5, 6, 8], "implic": 3, "xxz": [3, 12], "heisenberg": 3, "whose": [3, 6], "x_": 3, "y_i": 3, "y_": 3, "delta": [3, 5, 6], "2s": 3, "satisfi": [3, 8], "mathfrak": [3, 11], "su": [3, 11], "a_i": [3, 11], "b_j": 3, "delta_": [3, 6], "sum_c": 3, "varepsilon_": 3, "abc": 3, "c_i": [3, 8], "xxx": [3, 11], "y": [3, 12], "neq": 3, "half": [3, 8], "see": [3, 5, 6, 8, 11], "wherea": [3, 11], "simeq": 3, "u": [3, 4, 8, 11], "retain": [3, 5, 6], "automat": [3, 6, 8], "theta": 3, "interpret": [3, 4, 8, 12], "rotat": [3, 4], "around": [3, 4], "axi": 3, "angl": 3, "quantiti": 3, "m_z": 3, "motiv": 3, "gentli": 3, "backbon": 3, "roughli": 3, "speak": [3, 6], "g": [3, 6, 12], "rule": 3, "compos": 3, "step": [3, 5, 6, 11], "time": [3, 4, 8, 11, 12], "discret": [3, 6], "consist": [3, 4, 5, 8], "carri": [3, 4], "transform": [3, 4, 8], "should": [3, 4, 8, 11], "result": [3, 4, 5, 6, 8, 11], "next": [3, 4, 8], "over": [3, 4, 6, 11], "theta_2": 3, "theta_1": 3, "g_1": 3, "g_2": 3, "endow": 3, "exist": [3, 4, 5], "1g": 3, "g1": 3, "foral": 3, "abelian": 3, "composit": 3, "k": [3, 4, 6, 12], "hk": 3, "gh": 3, "would": [3, 5, 8], "formal": [3, 5, 8], "undon": 3, "opposit": 3, "henc": [3, 11], "invers": [3, 4], "gg": 3, "subgroup": 3, "suggest": [3, 11], "subset": 3, "li": 3, "heart": 3, "mathbb": [3, 4, 6, 8], "_2": [3, 11], "explan": [3, 11], "notat": [3, 8], "undergo": 3, "symbol": 3, "law": [3, 10], "_n": 3, "modulo": 3, "encount": 3, "unimodular": 3, "matric": [3, 4, 8, 11], "det": 3, "uu": 3, "geq": 3, "none": 3, "3d": 3, "unit": [3, 8], "m": [3, 4, 5, 6, 8, 12], "r": [3, 4, 8, 11, 12], "mm": 3, "tm": 3, "were": [3, 4], "deal": [3, 5, 6, 8], "absenc": 3, "invert": 3, "linear": [3, 4, 11, 13], "singular": [3, 5, 8, 11], "do": [3, 5, 6], "wonder": 3, "come": [3, 4, 5, 6], "exactli": [3, 4, 5, 6, 8, 11], "underli": [3, 4], "linearli": 3, "vector": [3, 4, 8, 11], "immedi": 3, "rais": 3, "plethora": 3, "kind": [3, 4, 5, 6], "ones": [3, 5], "answer": 3, "sake": 3, "index": [3, 6, 8], "x_g": 3, "x_gx_h": 3, "alwai": [3, 6, 8], "dimens": [3, 4, 5, 6, 8, 11, 12], "probabl": 3, "x_0": 3, "x_1": 3, "inde": [3, 4, 5, 6, 8], "x_1x_1": 3, "squar": [3, 4], "bar": [3, 8], "x_h": 3, "overlin": 3, "y_g": 3, "kroneck": [3, 4], "oplu": 3, "observ": [3, 6, 8], "unitarili": 3, "ux_gu": 3, "independ": [3, 4, 6], "pmatrix": 3, "show": [3, 5, 8], "crux": 3, "appropri": [3, 5, 6], "brought": [3, 5, 6], "block": [3, 6, 9], "diagon": [3, 4, 8], "fbox": 3, "1_g": 3, "2_g": 3, "vdot": 3, "ddot": 3, "amongst": 3, "irrep": 3, "shown": [3, 4, 8], "d_": 3, "One": [3, 4, 8], "kei": 3, "decompos": [3, 4], "fusion": 3, "clebsch": 3, "gordan": 3, "coeffici": 3, "explicitli": [3, 4, 8], "schur": 3, "lemma": 3, "x_gy": 3, "yx_g": 3, "proport": 3, "pose": 3, "s_1": 3, "s_2": 3, "bigoplus_": 3, "been": [3, 4, 5, 6, 13], "analyt": [3, 11], "low": [3, 4, 6, 8], "tabul": 3, "z_g": 3, "strong": 3, "didn": 3, "assum": [3, 5, 6, 8, 11], "argu": 3, "bring": [3, 5, 6, 8, 11], "appear": [3, 11], "write": [3, 8], "bigoplus_c": 3, "b_c": 3, "_c": 3, "decomposit": [3, 5, 6, 8, 11], "store": [3, 4, 8], "effici": [3, 4, 5, 6, 8, 11, 12], "particularli": [3, 4, 5, 6, 8], "paragraph": 3, "herebi": 3, "drastic": 3, "amount": [3, 4], "memori": [3, 4], "abl": [3, 6], "manipul": 3, "exploit": [3, 5, 6], "maximum": 3, "rank": [3, 4], "su\u2082spac": 3, "l": [3, 4, 8, 11, 12], "rep": 3, "su\u2082": 3, "essenti": 3, "copi": [3, 4], "summand": 3, "want": [3, 5, 6], "ss": 3, "productspac": [3, 8], "data": [3, 4], "fusiontre": 3, "fals": 3, "5e": 3, "323": 3, "inspect": 3, "assert": [3, 8], "dim": 3, "sortedvectordict": 3, "su2irrep": 3, "entri": 3, "fill": 3, "compat": 3, "cannot": [3, 5, 6], "fuse": 3, "third": 3, "0e": 3, "324": 3, "four": [3, 4], "conclud": [3, 5], "global": [3, 8, 11], "familiar": 3, "gaug": [3, 5, 6, 10, 11], "ubiquit": 3, "phenomena": [3, 8], "actual": [3, 5, 6], "each": [3, 4, 5, 6, 8, 11], "redund": 3, "descript": [3, 5, 6], "brief": [3, 4, 11], "overview": [3, 6, 8], "mostli": [3, 4, 8], "neglect": 3, "spatial": [3, 5, 6, 12], "reflect": [3, 4, 6], "don": [3, 6, 11], "anymor": 3, "classif": 3, "notori": 3, "rich": 3, "beauti": [3, 11], "especi": [3, 4], "higher": [3, 4, 5, 6], "algorithm": [3, 4, 8, 12], "tremend": 3, "speedup": 3, "stabil": [3, 5, 6], "benefit": 3, "uniform": [3, 8, 11, 12], "repeat": [3, 4, 8, 11], "indefinit": 3, "discoveri": 3, "topolog": [3, 6, 12], "matter": [3, 4, 6], "anyon": [3, 12], "excit": 3, "grow": [3, 5], "fascin": [3, 4], "explor": [3, 4], "categor": 3, "beyond": [3, 11], "scope": [3, 4, 11], "intric": [3, 6], "algebra": [3, 4, 11, 13], "categori": 3, "chain": [3, 8, 11, 12], "ftl": [3, 12], "07": [3, 12], "storag": 3, "start": [4, 6, 8, 11, 13], "modern": 4, "simplifi": [4, 5], "bc17": [4, 11, 12], "re": 4, "export": [4, 8], "macro": [4, 8], "separ": 4, "load": 4, "test": 4, "journei": 4, "evolut": 4, "profound": 4, "theoret": [4, 6, 12], "method": [4, 6, 8, 12], "varieti": [4, 5, 8], "studi": 4, "machin": 4, "earli": 4, "root": 4, "back": [4, 6, 11], "19th": 4, "centuri": 4, "pioneer": 4, "mathematician": 4, "arthur": 4, "caylei": 4, "jame": 4, "sylvest": 4, "multi": [4, 8, 13], "began": 4, "emerg": 4, "late": 4, "20th": 4, "dmrg": [4, 6], "birth": 4, "attribut": 4, "mp": [4, 5, 6, 8, 10, 11], "1960": 4, "earliest": 4, "wide": 4, "steven": 4, "white": 4, "1992": 4, "simul": [4, 12], "1980": 4, "1990": 4, "driven": 4, "add": [4, 6], "entangl": [4, 10, 12], "becam": 4, "central": [4, 5, 6, 8], "progress": 4, "extend": [4, 8], "tn": 4, "pair": [4, 6], "pep": 4, "scale": [4, 5, 6, 8], "renorm": [4, 6, 11, 12], "ansatz": [4, 6, 8, 10, 11], "mera": 4, "2000": 4, "disciplin": 4, "appli": [4, 5, 6, 8], "promin": 4, "unsuprisingli": 4, "plai": [4, 8], "role": [4, 8], "circuit": 4, "ongo": 4, "applic": [4, 6, 12, 13], "vibrant": 4, "evolv": [4, 6], "variou": 4, "direct": [4, 11], "phase": 4, "main": 4, "advantag": [4, 11], "admit": [4, 6], "greatli": 4, "involv": [4, 5], "numer": [4, 6, 8, 13], "node": 4, "graph": 4, "depict": 4, "leg": [4, 8], "stick": [4, 11], "individu": [4, 6, 11], "recoverd": 4, "fix": [4, 6, 8, 13], "diagram": [4, 8], "r_": [4, 11], "i_3": 4, "i_4": 4, "freeli": [4, 8], "move": [4, 8, 11], "long": [4, 5], "shape": 4, "certain": [4, 5, 8], "explicit": 4, "represent": [4, 5, 6, 11], "ow": 4, "precis": [4, 6, 11, 12], "detail": [4, 8, 11], "convent": 4, "d_1": 4, "d_2": 4, "d_3": 4, "d_i": 4, "immateri": 4, "conveni": 4, "perform": [4, 6, 8], "size": [4, 6, 8, 11, 12], "need": 4, "adress": 4, "implicit": 4, "eltyp": 4, "7": [4, 11, 12], "8": [4, 8, 11, 12], "9": [4, 11, 12], "11": [4, 11, 12], "13": [4, 11], "15": [4, 11], "10": [4, 11, 12], "12": [4, 11], "14": [4, 11, 12], "16": [4, 11, 12], "Of": 4, "three": [4, 6, 8], "complic": 4, "join": 4, "similar": [4, 5, 6, 8], "einstein": 4, "summat": 4, "signifi": 4, "partial": 4, "cyclic": 4, "slide": 4, "loop": [4, 5, 6], "placement": 4, "text": [4, 6], "tr": [4, 6, 8], "ab": [4, 12], "ba": 4, "drawn": 4, "famililiar": 4, "inner": 4, "arbitrarili": 4, "evalu": [4, 5, 6, 8, 11], "sequenc": [4, 5, 6, 12], "wise": 4, "specifi": [4, 8, 11], "assign": [4, 6], "implicitli": 4, "sum": [4, 5, 6, 8], "obtain": [4, 5, 6, 8, 11], "must": 4, "less": 4, "alloc": 4, "quickli": [4, 5], "unwieldi": 4, "wish": [4, 6, 8], "pairwis": 4, "spirit": 4, "minor": 4, "modif": 4, "ncon": 4, "increas": [4, 5, 6], "neg": 4, "posit": [4, 8], "f": [4, 12], "32095": 4, "80467": 4, "80061": 4, "12631": 4, "vari": 4, "wildli": 4, "problem": [4, 5, 6, 11], "both": 4, "length": [4, 8, 11], "2n": [4, 5, 6], "float": 4, "flop": 4, "substanti": [4, 5, 6], "prefer": 4, "area": [4, 6, 10], "boil": [4, 11], "down": [4, 8, 11], "heurist": 4, "cut": 4, "bubbl": 4, "ineffici": 4, "infeas": 4, "optim": [4, 5, 8, 11, 12], "ladder": 4, "highlight": 4, "few": [4, 5], "np": 4, "hard": [4, 5], "larger": [4, 6], "30": [4, 11], "40": 4, "phv14": [4, 12], "process": 4, "opt": 4, "keyword": 4, "enabl": [4, 8, 11], "done": [4, 6], "compil": 4, "onc": [4, 6], "\u03b1": 4, "\u03b2": 4, "\u03b3": 4, "\u03f5": 4, "\u03b6": 4, "\u03b7": 4, "\u03b4": 4, "instrument": 4, "approxim": [4, 5, 8, 11], "partit": [4, 6], "everyth": [4, 8], "s1": 4, "s2": 4, "eigen": 4, "ax": 4, "imag": 4, "normal": [4, 5, 11], "spectral": [4, 8], "equat": [4, 6, 11], "av": 4, "diagrammat": [4, 11], "randn": [4, 8], "complexf64": [4, 8, 11], "eigendecomposit": 4, "eig": 4, "32m": 4, "1mtest": 4, "pass": [4, 6], "svd": [4, 5, 6], "seen": [4, 5, 8], "rectangular": 4, "isometr": [4, 11], "best": [4, 11], "truncat": [4, 5, 6], "largest": [4, 5, 8], "tensori": 4, "version": [4, 8], "tsvd": 4, "permut": 4, "id": [4, 8], "methoderror": 4, "match": 4, "closest": 4, "candid": 4, "abstracttensormap": 4, "trunc": 4, "alg": 4, "mikvc": 4, "src": 4, "253": 4, "vararg": 4, "kwarg": 4, "50": [4, 12], "stacktrac": 4, "top": 4, "level": 4, "semi": 4, "hermitian": [4, 5, 6], "q": 4, "leftorth": 4, "upper": 4, "triangular": 4, "solv": [4, 5, 11], "solut": [4, 5, 6, 11], "easi": [4, 5, 6], "gaussian": 4, "elimin": 4, "overdetermin": 4, "variant": 4, "qrpo": 4, "transpos": 4, "rq": 4, "ql": 4, "lq": 4, "reveal": 4, "trade": 4, "off": 4, "expens": [4, 5], "practic": [4, 5, 6], "zero": 4, "leftnul": 4, "norm": [4, 8, 11], "commonli": [4, 8], "focu": [5, 6], "advanc": 5, "importantli": 5, "effect": [5, 6, 8, 11], "correl": [5, 6, 10, 12], "deriv": [5, 11], "desir": [5, 6], "initi": [5, 6, 8, 11], "psi_0": [5, 6], "later": 5, "naiv": [5, 6], "try": 5, "exponenti": [5, 6, 8], "prohibit": 5, "random": [5, 11], "exhibit": 5, "power": [5, 8], "interact": [5, 6, 11, 12], "h_": [5, 6, 11], "although": [5, 6], "unfeas": 5, "smaller": [5, 6, 8], "subsystem": [5, 6], "ih_": 5, "instead": [5, 6, 8], "suzuki": [5, 6], "trotter": [5, 6], "a4": 5, "mathcal": [5, 6, 11], "split": [5, 6], "interv": [5, 6], "h_e": [5, 6], "h_o": [5, 6], "error": [5, 6, 11], "choos": [5, 6, 8, 11], "suffici": [5, 6], "famili": [5, 6], "hs05": [5, 6, 12], "procedur": [5, 6, 8], "dynam": [5, 6], "aforement": 5, "eq": [5, 6, 8], "group": [5, 6, 11, 12], "within": [5, 6], "odd": [5, 6], "sum_n": [5, 6], "qquad": [5, 6, 11], "overlap": [5, 6, 8], "bond": [5, 6, 8, 11], "layer": [5, 6], "chi": [5, 6], "bra": [5, 8], "face": 5, "difficulti": [5, 6], "strategi": 5, "smallest": [5, 11], "replac": 5, "tau": [5, 6], "lowest": [5, 6], "damp": 5, "lim_": 5, "e_i": 5, "psi_i": 5, "e_0": [5, 6], "approx": 5, "regard": 5, "tackl": [5, 6], "until": [5, 11], "converg": [5, 8, 11], "reach": [5, 6, 11], "tip": 5, "iceberg": 5, "To": [5, 6, 8], "briefli": [5, 6], "comment": [5, 6, 11], "contain": [5, 6], "parallel": [5, 6], "approach": [5, 6, 11], "updat": [5, 6, 8, 11], "tebd_trunc": [5, 6], "account": [5, 6], "rest": [5, 6, 8], "surround": [5, 6], "captur": [5, 6, 8], "techniqu": [5, 6, 8, 11], "environ": [5, 6], "qualiti": [5, 6], "directli": [5, 6, 11], "affect": [5, 6], "jwx08": [5, 6, 12], "ignor": [5, 6, 8], "accur": [5, 6, 8, 12], "jorusv": [5, 6, 12], "08": [5, 6, 12], "gain": [5, 6], "accuraci": [5, 6, 8], "cost": [5, 6, 8], "spent": 6, "worth": 6, "character": [6, 8, 11], "symmetri": [6, 8, 13], "protect": 6, "aspect": [6, 8], "intens": 6, "hilbert": [6, 8], "locat": [6, 8], "bigotimes_": 6, "fig": 6, "ref": 6, "Such": 6, "geometri": 6, "systemat": 6, "contini": 6, "theori": [6, 13], "cleverli": 6, "amen": 6, "adapt": 6, "cmp": 6, "ctn": 6, "measur": 6, "langl": [6, 8], "rangl": [6, 8], "rm": [6, 8], "rho": 6, "connect": [6, 8], "mu": 6, "hand": [6, 11, 12, 13], "dictat": [6, 8], "govern": 6, "thermal": 6, "equilibrium": 6, "geometr": 6, "subregion": 6, "leq": 6, "bound": [6, 8], "constant": 6, "neighbor": 6, "nearest": 6, "asid": 6, "microscop": 6, "relativist": 6, "guarante": [6, 8], "condens": 6, "highli": 6, "class": [6, 8, 11], "qantum": 6, "footnot": 6, "protocol": 6, "prepar": 6, "estim": [6, 11], "varepsilon": 6, "setup": 6, "temperatur": 6, "mix": [6, 8, 11], "densiti": [6, 11], "beta": [6, 8], "strength": 6, "drawback": 6, "tack": 6, "rel": 6, "succunctli": 6, "vump": 6, "gradient": 6, "proven": [6, 11], "encod": [6, 8, 11], "ensembl": 6, "wich": 6, "analyz": 6, "sens": [6, 8], "those": [6, 11], "fewer": 6, "s_i": 6, "run": 6, "configur": 6, "rewritten": 6, "contract": [6, 8, 12], "look": [6, 8], "someth": 6, "weight": 6, "zaunerstaubervf": [6, 12], "18": [6, 11, 12], "transfer": [6, 8, 12], "pull": 6, "throught": 6, "parititon": 6, "altern": [6, 8], "patch": 6, "corner": [6, 12], "ctmrg": 6, "no96": [6, 12], "anyth": 6, "decim": [6, 11], "vid03": [6, 12], "design": [6, 11, 13], "reliabl": 6, "achiev": 6, "schr\u00f6dinger": 6, "intract": 6, "suitabl": 6, "difficult": 6, "ne": 6, "prod_": [6, 8], "gate": 6, "reli": [6, 11], "exi": 6, "tentir": 6, "propto": 6, "unless": 6, "disconnect": [6, 8], "factor": [6, 11], "exponenenti": 6, "efficien": 6, "manual": 6, "sensibl": 6, "imaginari": 6, "excel": [8, 11], "review": [8, 12], "vhv19": [8, 11, 12], "thorough": 8, "technic": 8, "tangent": [8, 12], "supplement": 8, "routin": [8, 11], "boldsymbol": 8, "_l": 8, "s_m": 8, "_r": 8, "view": 8, "control": 8, "gap": [8, 11], "never": 8, "safe": 8, "bulk": 8, "faithfulli": 8, "impos": 8, "transat": 8, "diagramat": 8, "cell": 8, "unform": 8, "magnitud": 8, "lambda_0": 8, "lambda_i": 8, "remain": 8, "mangitud": 8, "projector": 8, "ensur": 8, "properli": 8, "rescal": 8, "trace": 8, "With": 8, "middl": 8, "suppos": 8, "extens": 8, "o_n": 8, "beta_m": 8, "alpha_n": 8, "abritrari": 8, "insert": 8, "decai": 8, "why": 8, "harder": 8, "xi": 8, "lambda_1": 8, "log": 8, "lambda_": 8, "mathrm": 8, "max": [8, 11], "sublead": 8, "focuss": 8, "sector": 8, "target": 8, "exit": 8, "crucial": 8, "rcc18": [8, 12], "despit": 8, "inher": 8, "convers": 8, "mai": [8, 12], "parametr": 8, "orthonorm": [8, 11], "a_l": [8, 11], "iter": [8, 11], "qr": 8, "docomposit": 8, "guess": [8, 11], "repeatedli": [8, 11], "room": 8, "a_r": [8, 11], "center": [8, 11], "a_c": [8, 11], "contrast": 8, "lr": 8, "therebi": 8, "usv": 8, "absorb": 8, "residu": 8, "arriv": 8, "straightforwardli": 8, "schmidt": 8, "across": 8, "i_l": 8, "i_r": 8, "bipartit": 8, "maxim": 8, "isometri": [8, 11], "correspondingli": 8, "lower": 8, "variat": [8, 11, 12], "tild": [8, 11], "go": 8, "virtual": 8, "bb1c41ca": 8, "d63c": 8, "52ed": 8, "829e": 8, "0820dda26502": 8, "cr": 8, "al": 8, "linearalgebra": 8, "ac": 8, "0000000000000002": 8, "verifi": 8, "al_id": 8, "conj": 8, "ar_id": 8, "lh": 8, "rh": 8, "expectation_valu": [8, 11], "2965360780687573": 8, "024684087674045976im": 8, "correlation_length": 8, "5210339009216689": 8, "statist": [9, 13], "mpo": [9, 11], "jordan": [9, 12], "polynomi": 10, "manner": 11, "min_": 11, "simplic": 11, "consider": 11, "rang": 11, "thermodynam": 11, "break": 11, "tri": 11, "sequenti": 11, "sweep": 11, "bit": 11, "a_1": 11, "a_2": 11, "seemingli": 11, "daunt": 11, "denomin": 11, "h_i": 11, "_i": 11, "forth": 11, "At": 11, "reus": 11, "manifestli": 11, "suffer": 11, "artefact": 11, "surprisingli": 11, "success": 11, "variation": 11, "mpskitmodel": 11, "z_j": 11, "h_z": 11, "paramet": 11, "default": 11, "transverse_field_is": 11, "summon": 11, "\u03c8": 11, "finitemp": 11, "\u03c8\u2080": 11, "_": 11, "find_groundst": 11, "ca635005": 11, "6f8c": 11, "4cd1": 11, "b51d": 11, "8491250ef2ab": 11, "39miteraton": 11, "004918527829644078": 11, "6018489268067254e": 11, "172740472500253e": 11, "7112592674678824e": 11, "8748929872927394e": 11, "447243030289668e": 11, "9588728321554764e": 11, "1763451120881838e": 11, "679898216523513e": 11, "863024577370902e": 11, "421004073641341e": 11, "957687044386289e": 11, "1794092881999024e": 11, "705208493519156e": 11, "877880521695336e": 11, "497530684181754e": 11, "1354652128300114e": 11, "17": 11, "191727427841558e": 11, "878483208839837e": 11, "unbound": 11, "f_l": 11, "f_r": 11, "obei": 11, "offer": 11, "intial": 11, "h_c": 11, "a_lc": 11, "ca_r": 11, "minimum": 11, "manifold": 11, "last": 11, "intertwin": 11, "explain": 11, "chose": 11, "toler": 11, "eta": 11, "arnoldi": 11, "yield": 11, "epsilon_l": 11, "epsilon_r": 11, "soltuion": 11, "epsilon": 11, "u_lv_l": 11, "u_l": 11, "v_l": 11, "aris": 11, "sigma_lv_l": 11, "u_rv_r": 11, "u_r": 11, "sigma_rv_r": 11, "exact": 11, "s_c": 11, "s_lc": 11, "ca": 11, "s_r": 11, "sigma_": 11, "arithmet": 11, "u_": 11, "poor": 11, "robust": 11, "l_": 11, "l_c": 11, "r_c": 11, "polar": 11, "l_cp": 11, "r_cu": 11, "cdot1": 11, "heisenberg_xyz": 11, "infinitemp": 11, "env": 11, "39mvump": 11, "galerkin": 11, "07246103469777": 11, "2666467311860786": 11, "03674145848146899": 11, "01630290541053547": 11, "006755801125032585": 11, "00213024746120037": 11, "0007360704395154456": 11, "0002556331026616265": 11, "97742113007201e": 11, "213758322737539e": 11, "1525599364714219e": 11, "1966267301429994e": 11, "5313217925607893e": 11, "647102155657192e": 11, "0890138861170584e": 11, "782151029816918e": 11, "9100058380751707e": 11, "0930912593699508e": 11, "19": 11, "122387382308404e": 11, "20": 11, "559261655394076e": 11, "21": 11, "920596218088465e": 11, "22": [11, 12], "2525120558260656e": 11, "23": 11, "599293928545564e": 11, "24": 11, "288214116430126e": 11, "25": [11, 12], "2613085514875523e": 11, "26": 11, "849412800016547e": 11, "27": 11, "8785021698245177e": 11, "28": 11, "564909882632748e": 11, "periodicarrai": 11, "4013806435182494": 11, "3201926491189795e": 11, "17im": 11, "compar": 11, "quasi": 11, "401": 11, "484": 11, "038": 11, "971": 11, "hco": [11, 12], "jacob": [12, 13], "bridgeman": [12, 13], "christoph": 12, "chubb": 12, "wave": 12, "danc": 12, "introductori": 12, "journal": 12, "223001": 12, "2017": 12, "url": 12, "http": 12, "dx": 12, "doi": 12, "org": 12, "1088": 12, "1751": 12, "8121": 12, "aa6dc3": 12, "adrian": 12, "feiguin": 12, "simon": 12, "trebst": 12, "andrea": 12, "ww": 12, "ludwig": 12, "matthia": 12, "troyer": 12, "alexei": 12, "kitaev": 12, "zhenghan": 12, "wang": 12, "michael": 12, "freedman": 12, "liquid": 12, "golden": 12, "letter": 12, "98": 12, "160409": 12, "2007": 12, "jutho": [12, 13], "haegeman": [12, 13], "ignacio": 12, "cirac": 12, "tobia": 12, "osborn": 12, "iztok": 12, "\u017e": 12, "orn": 12, "henri": 12, "verscheld": 12, "frank": 12, "verstraet": 12, "107": 12, "070601": 12, "2011": 12, "miss": 12, "booktitl": 12, "hatano2005find": 12, "jiang": 12, "weng": 12, "xiang": 12, "101": 12, "090603": 12, "august": 12, "2008": 12, "arxiv": 12, "0806": 12, "3719": 12, "1103": 12, "physrevlett": 12, "Or": 12, "\u00fa": 12, "vidal": 12, "250602": 12, "decemb": 12, "ap": 12, "nishino": 12, "okunishi": 12, "societi": 12, "japan": 12, "65": 12, "891": 12, "894": 12, "april": 12, "1996": 12, "cond": 12, "mat": 12, "9507087": 12, "1143": 12, "jpsj": 12, "robert": 12, "pfeifer": 12, "faster": 12, "phy": 12, "rev": 12, "90": 12, "033315": 12, "sep": 12, "2014": 12, "physrev": 12, "marek": 12, "ram": 12, "piotr": 12, "czarnik": 12, "lukasz": 12, "cincio": 12, "extrapol": 12, "asymptot": 12, "bose": 12, "hubbard": 12, "041033": 12, "novemb": 12, "2018": 12, "1801": 12, "08554": 12, "physrevx": 12, "lauren": 12, "vanderstraeten": 12, "scipost": 12, "007": 12, "januari": 12, "2019": 12, "scipostphyslectnot": 12, "1810": 12, "07006": 12, "21468": 12, "guifr": 12, "slightli": 12, "91": 12, "147902": 12, "octob": 12, "2003": 12, "quant": 12, "ph": 12, "0301063": 12, "zauner": 12, "stauber": 12, "fishman": 12, "97": 12, "045145": 12, "1701": 12, "07035": 12, "physrevb": 12, "seri": 13, "aim": 13, "practis": 13, "softwar": 13, "lander": 13, "burgelman": 13, "luka": 13, "devo": 13, "daan": 13, "maerten": 13, "bram": 13, "vancraeynest": 13, "cuiper": 13, "kevin": 13, "vervoort": 13}, "objects": {}, "objtypes": {}, "objnames": {}, "titleterms": {"multi": 0, "linear": 0, "algebra": 0, "content": [0, 4, 8], "overview": [0, 4], "vector": 0, "matric": 0, "tensor": [0, 2, 3, 4, 5, 6, 13], "product": [0, 3, 4, 8, 9, 10, 13], "map": 0, "conclus": [0, 4, 5], "quantum": [1, 3, 5, 6], "mani": [1, 3], "bodi": [1, 3], "theori": [1, 3, 4], "get": 2, "start": 2, "numer": 2, "version": 2, "control": 2, "julia": [2, 13], "packag": 2, "noteworthi": 2, "network": [2, 4, 5, 6, 13], "softwar": 2, "symmetri": 3, "physic": 3, "exampl": [3, 5, 8, 11], "applic": [3, 7], "break": 3, "order": [3, 4], "paramet": 3, "phase": 3, "noether": 3, "conserv": 3, "quantitit": 3, "group": [3, 4], "represent": [3, 8], "definit": 3, "complex": [3, 4], "conjug": 3, "direct": 3, "sum": 3, "irreduc": 3, "symmetr": 3, "outlook": [3, 5, 6], "gener": 3, "histori": 4, "graphic": 4, "notat": 4, "oper": [4, 6, 9], "index": 4, "split": 4, "indic": 4, "outer": 4, "trace": 4, "contract": 4, "factor": 4, "eigenvalu": 4, "decomposit": 4, "singular": 4, "valu": [4, 8], "polar": 4, "qr": 4, "nullspac": 4, "algorithm": [5, 6, 11, 13], "simul": [5, 6], "system": [5, 6], "time": [5, 6], "evolv": 5, "block": 5, "decim": 5, "tebd": [5, 6], "One": [5, 6], "dimension": [5, 6], "nearest": 5, "neighbor": 5, "hamiltonian": 5, "groundstat": 5, "search": 5, "imaginari": 5, "evolut": [5, 6], "todo": 6, "classic": 6, "statist": 6, "mechan": 6, "case": 6, "studi": 6, "approxim": 6, "The": 6, "infinit": 8, "matrix": [8, 9, 10, 13], "state": [8, 10, 13], "thermodynam": 8, "limit": 8, "normal": 8, "expect": 8, "correl": 8, "function": 8, "gaug": 8, "revisit": 8, "entangl": 8, "entropi": 8, "truncat": 8, "code": 8, "mpskit": 8, "infinitemp": 8, "fix": 11, "point": 11, "dmrg": 11, "vump": 11, "refer": 12, "method": 13, "introduct": 13, "other": 13}, "envversion": {"sphinx.domains.c": 2, "sphinx.domains.changeset": 1, "sphinx.domains.citation": 1, "sphinx.domains.cpp": 6, "sphinx.domains.index": 1, "sphinx.domains.javascript": 2, "sphinx.domains.math": 2, "sphinx.domains.python": 3, "sphinx.domains.rst": 2, "sphinx.domains.std": 2, "sphinx.ext.intersphinx": 1, "sphinxcontrib.bibtex": 9, "sphinx": 56}})