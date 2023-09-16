Search.setIndex({"docnames": ["1-Introduction/Julia", "1-Introduction/QuantumManyBody", "1-Introduction/TensorNetworks", "2-Tensors/TensorKit", "2-Tensors/TensorOperations", "3-MatrixProductStates/Algorithms", "3-MatrixProductStates/Applications", "3-MatrixProductStates/InfiniteMPS", "3-MatrixProductStates/MatrixProductOperators", "3-MatrixProductStates/MatrixProductStates", "References", "intro"], "filenames": ["1-Introduction/Julia.md", "1-Introduction/QuantumManyBody.md", "1-Introduction/TensorNetworks.md", "2-Tensors/TensorKit.ipynb", "2-Tensors/TensorOperations.ipynb", "3-MatrixProductStates/Algorithms.md", "3-MatrixProductStates/Applications.md", "3-MatrixProductStates/InfiniteMPS.md", "3-MatrixProductStates/MatrixProductOperators.md", "3-MatrixProductStates/MatrixProductStates.md", "References.md", "intro.md"], "titles": ["<span class=\"section-number\">3. </span>Getting started with Julia", "<span class=\"section-number\">1. </span>Quantum Many-Body Theory", "<span class=\"section-number\">2. </span>Tensor Network Theory", "<span class=\"section-number\">4. </span>TensorKit.jl", "<span class=\"section-number\">5. </span>TensorOperations.jl", "<span class=\"section-number\">9. </span>Algorithms", "<span class=\"section-number\">10. </span>Applications", "<span class=\"section-number\">8. </span>Infinite Matrix Product States", "<span class=\"section-number\">7. </span>Matrix Product Operators", "<span class=\"section-number\">6. </span>Matrix Product States", "<span class=\"section-number\">11. </span>References", "Tensor Network Methods with Julia"], "terms": {"introduct": 1, "context": [1, 2], "histori": [1, 2], "purpos": [1, 2], "relev": [1, 2], "second": [1, 7], "quantiz": 1, "tensor": [1, 7, 10], "product": [1, 10], "classic": [1, 6], "multi": 2, "linear": 2, "algebra": 2, "graphic": 2, "notat": [2, 7], "comput": [2, 7], "complex": [2, 7], "test": [3, 4], "basic": [5, 7], "simpl": [5, 7], "updat": [5, 7], "trotter": 5, "tebd": 5, "transvers": 6, "field": 6, "Ising": 6, "model": [6, 10], "thi": [7, 11], "section": 7, "discuss": 7, "mp": [7, 9], "properti": 7, "our": 7, "mostli": 7, "base": 7, "excel": 7, "review": [7, 10], "vhv19": [7, 10], "which": 7, "provid": [7, 11], "thorough": 7, "technic": 7, "overview": 7, "tangent": [7, 10], "space": [7, 10], "method": [7, 10], "uniform": [7, 10], "The": 7, "formal": 7, "exposit": 7, "supplement": 7, "some": 7, "veri": 7, "work": 7, "us": 7, "jl": [7, 11], "end": 7, "For": 7, "more": 7, "detail": 7, "numer": 7, "implement": 7, "routin": 7, "we": 7, "refer": [7, 11], "julia": 7, "version": 7, "tutori": [7, 11], "again": 7, "finit": 7, "introduc": 7, "previou": 7, "can": 7, "readili": 7, "extend": 7, "consid": 7, "an": [7, 8], "one": 7, "dimension": 7, "chain": 7, "local": 7, "physic": [7, 10, 11], "hilbert": 7, "mathbb": 7, "c": 7, "d": 7, "dimens": 7, "everi": 7, "site": 7, "repres": 7, "quantum": [7, 11], "system": 7, "ha": 7, "form": [7, 8], "left": 7, "psi": 7, "A": 7, "right": 7, "rangl": 7, "sum_": 7, "s": 7, "boldsymbol": 7, "v": 7, "_l": 7, "dagger": 7, "prod_": 7, "m": [7, 10], "z": 7, "s_m": 7, "_r": 7, "here": 7, "each": 7, "time": 7, "index": 7, "As": 7, "befor": 7, "altern": 7, "view": 7, "three": 7, "indic": 7, "where": 7, "so": 7, "call": 7, "bond": 7, "assum": 7, "same": 7, "case": 7, "control": 7, "correspond": 7, "infti": 7, "8": [7, 10], "1": 7, "ani": 7, "up": 7, "arbitrari": 7, "accuraci": 7, "certain": 7, "class": 7, "low": 7, "energi": 7, "gap": 7, "howev": 7, "accur": 7, "approxim": 7, "much": 7, "smaller": 7, "note": [7, 10], "while": 7, "eq": 7, "have": [7, 11], "also": 7, "two": 7, "boundari": 7, "vector": 7, "condit": 7, "never": 7, "mean": 7, "These": 7, "therefor": 7, "safe": 7, "ignor": 7, "follow": 7, "all": 7, "bulk": 7, "ar": 7, "faithfulli": 7, "captur": 7, "invari": 7, "under": 7, "translat": 7, "natur": 7, "impos": 7, "transat": 7, "lead": 7, "In": 7, "diagramat": 7, "first": 7, "instead": 7, "non": 7, "trivial": 7, "repeat": 7, "unit": 7, "cell": 7, "size": 7, "would": 7, "One": 7, "central": 7, "object": 7, "unform": 7, "calcul": 7, "transfer": 7, "oper": [7, 11], "defin": 7, "act": 7, "matric": 7, "interpret": 7, "4": [7, 10], "leg": 7, "otim": 7, "leftarrow": 7, "shown": 7, "complet": 7, "posit": 7, "map": 7, "its": 7, "eigenvalu": 7, "number": 7, "eigenvaluesof": 7, "character": 7, "length": 7, "eigenvector": 7, "evalu": 7, "observ": 7, "norm": 7, "contract": 7, "clearli": 7, "noth": 7, "than": 7, "abov": 7, "spectral": 7, "decomposit": 7, "n": 7, "th": 7, "power": 7, "e": 7, "l": 7, "r": 7, "fix": 7, "point": 7, "largest": 7, "magnitud": 7, "lambda_0": 7, "lambda_i": 7, "remain": 7, "mangitud": 7, "take": 7, "express": 7, "reduc": 7, "projector": 7, "onto": 7, "To": 7, "ensur": 7, "properli": 7, "should": 7, "rescal": 7, "sqrt": 7, "well": [7, 11], "requir": 7, "trace": 7, "equal": 7, "With": 7, "place": 7, "overlap": 7, "between": 7, "sinc": 7, "effect": 7, "alwai": 7, "choos": 7, "langl": 7, "bar": 7, "middl": 7, "suppos": 7, "wish": 7, "extens": 7, "o": 7, "frac": 7, "o_n": 7, "If": 7, "singl": 7, "dictat": 7, "given": 7, "everyth": 7, "similarli": 7, "let": 7, "look": 7, "alpha": 7, "beta": 7, "bra": 7, "beta_m": 7, "alpha_n": 7, "ket": 7, "abritrari": 7, "locat": 7, "becaus": 7, "onli": 7, "depend": 7, "differ": 7, "insert": 7, "from": 7, "learn": 7, "determin": 7, "ground": 7, "inde": 7, "recal": 7, "now": 7, "0": 7, "see": 7, "part": 7, "just": 7, "disconnect": 7, "rest": 7, "exponenti": 7, "decai": 7, "impli": 7, "connect": 7, "reason": 7, "why": 7, "suit": 7, "critic": 7, "xi": 7, "lambda_1": 7, "log": 7, "lambda_": 7, "mathrm": 7, "max": 7, "sublead": 7, "typic": 7, "thei": 7, "focuss": 7, "specif": 7, "symmetri": 7, "sector": 7, "target": 7, "associ": 7, "exit": 7, "particular": 7, "plai": 7, "crucial": 7, "role": 7, "techniqu": 7, "scale": 7, "rcc18": [7, 10], "uniqu": 7, "convers": 7, "true": 7, "mai": 7, "give": [7, 11], "rise": 7, "easili": 7, "seen": 7, "transform": 7, "leav": 7, "freedom": [7, 9], "parametr": 7, "canon": 7, "start": [7, 11], "orthonorm": 7, "term": 7, "a_l": 7, "satisfi": 7, "find": 7, "bring": 7, "iter": 7, "procedur": 7, "qr": 7, "docomposit": 7, "initi": 7, "guess": 7, "repeatedli": 7, "perform": 7, "bound": 7, "converg": 7, "i": 7, "construct": 7, "choic": 7, "still": 7, "room": 7, "unitari": 7, "diagon": 7, "a_r": 7, "found": 7, "similar": 7, "final": 7, "mix": 7, "center": 7, "new": 7, "a_c": 7, "obtain": 7, "By": 7, "contrast": 7, "origin": 7, "commonli": 7, "intuit": 7, "lr": 7, "therebi": 7, "relat": 7, "allow": 7, "freeli": 7, "move": 7, "through": 7, "link": 7, "singular": 7, "usv": 7, "absorb": 7, "u": 7, "definit": 7, "residu": 7, "entir": 7, "rm": 7, "tr": 7, "ident": 7, "arriv": 7, "particularli": 7, "straightforwardli": 7, "write": 7, "down": 7, "schmidt": 7, "across": 7, "c_i": 7, "i_l": 7, "i_r": 7, "orthogon": 7, "half": 7, "lattic": 7, "element": 7, "exactli": 7, "bipartit": 7, "sum_i": 7, "2": 7, "enabl": 7, "effici": 7, "sum": 7, "optim": 7, "sens": 7, "maxim": 7, "column": 7, "isometri": 7, "correspondingli": 7, "result": 7, "guarante": 7, "lower": 7, "global": 7, "variat": 7, "cost": 7, "tild": 7, "packag": 7, "mani": [7, 11], "tool": 7, "without": 7, "go": 7, "alreadi": 7, "check": 7, "aspect": 7, "specifi": 7, "virtual": 7, "standard": 7, "tensorkit": [7, 11], "complexspac": 7, "3": 7, "5": 7, "\u2102": 7, "cr": 7, "tensormap": 7, "productspac": 7, "al": 7, "automat": 7, "store": 7, "linearalgebra": 7, "show": 7, "ac": 7, "9999999999999999": 7, "explicitli": 7, "verifi": 7, "network": [7, 10], "diagram": 7, "tensoroper": [7, 11], "macro": 7, "al_id": 7, "conj": 7, "ar_id": 7, "assert": 7, "id": 7, "lh": 7, "rh": 7, "consist": 7, "randn": 7, "expectation_valu": 7, "complexf64": 7, "43028995061224634": 7, "133275703540238im": 7, "encod": 7, "correlation_length": 7, "5366801105432544": 7, "export": 7, "varieti": 7, "algorithm": [7, 11], "next": 7, "statist": [8, 11], "mechan": [8, 11], "mpo": 8, "jordan": 8, "block": 8, "expect": [8, 9], "valu": [8, 9], "polynomi": 9, "ansatz": 9, "entangl": 9, "area": 9, "law": 9, "gaug": 9, "correl": [9, 10], "marek": 10, "ram": 10, "piotr": 10, "czarnik": 10, "lukasz": 10, "cincio": 10, "precis": 10, "extrapol": 10, "function": 10, "asymptot": 10, "state": 10, "applic": [10, 11], "bose": 10, "hubbard": 10, "xxz": 10, "x": 10, "041033": 10, "novemb": 10, "2018": 10, "url": 10, "http": 10, "arxiv": 10, "org": 10, "ab": 10, "1801": 10, "08554": 10, "doi": 10, "10": 10, "1103": 10, "physrevx": 10, "lauren": 10, "vanderstraeten": 10, "jutho": 10, "haegeman": 10, "frank": 10, "verstraet": 10, "matrix": 10, "scipost": 10, "lectur": 10, "page": 10, "007": 10, "januari": 10, "2019": 10, "scipostphyslectnot": 10, "7": 10, "1810": 10, "07006": 10, "21468": 10, "seri": 11, "It": 11, "illustr": 11, "theori": 11, "aim": 11, "hand": 11, "practis": 11, "code": 11, "exampl": 11, "showcas": 11, "softwar": 11, "librari": 11, "been": 11, "develop": 11, "bodi": 11, "get": 11, "infinit": 11}, "objects": {}, "objtypes": {}, "objnames": {}, "titleterms": {"get": 0, "start": 0, "julia": [0, 11], "quantum": 1, "mani": 1, "bodi": 1, "theori": [1, 2], "tensor": [2, 11], "network": [2, 11], "tensorkit": 3, "jl": [3, 4], "tensoroper": 4, "algorithm": 5, "applic": 6, "infinit": 7, "matrix": [7, 8, 9, 11], "product": [7, 8, 9, 11], "state": [7, 9, 11], "content": 7, "thermodynam": 7, "limit": 7, "represent": 7, "normal": 7, "expect": 7, "valu": 7, "correl": 7, "function": 7, "gaug": 7, "revisit": 7, "entangl": 7, "entropi": 7, "truncat": 7, "code": 7, "exampl": 7, "mpskit": 7, "infinitemp": 7, "oper": 8, "refer": 10, "method": 11, "introduct": 11, "other": 11}, "envversion": {"sphinx.domains.c": 2, "sphinx.domains.changeset": 1, "sphinx.domains.citation": 1, "sphinx.domains.cpp": 6, "sphinx.domains.index": 1, "sphinx.domains.javascript": 2, "sphinx.domains.math": 2, "sphinx.domains.python": 3, "sphinx.domains.rst": 2, "sphinx.domains.std": 2, "sphinx.ext.intersphinx": 1, "sphinxcontrib.bibtex": 9, "sphinx": 56}})