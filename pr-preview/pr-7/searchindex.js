Search.setIndex({"docnames": ["1-Introduction/Julia", "1-Introduction/QuantumManyBody", "1-Introduction/TensorNetworks", "2-Tensors/TensorKit", "2-Tensors/TensorOperations", "3-MatrixProductStates/Algorithms", "3-MatrixProductStates/Applications", "3-MatrixProductStates/InfiniteMPS", "3-MatrixProductStates/MatrixProductOperators", "3-MatrixProductStates/MatrixProductStates", "3-MatrixProductStates/TimeEvolution", "References", "intro"], "filenames": ["1-Introduction/Julia.md", "1-Introduction/QuantumManyBody.md", "1-Introduction/TensorNetworks.md", "2-Tensors/TensorKit.ipynb", "2-Tensors/TensorOperations.ipynb", "3-MatrixProductStates/Algorithms.md", "3-MatrixProductStates/Applications.md", "3-MatrixProductStates/InfiniteMPS.md", "3-MatrixProductStates/MatrixProductOperators.md", "3-MatrixProductStates/MatrixProductStates.md", "3-MatrixProductStates/TimeEvolution.md", "References.md", "intro.md"], "titles": ["<span class=\"section-number\">3. </span>Getting started with Julia", "<span class=\"section-number\">1. </span>Quantum Many-Body Theory", "<span class=\"section-number\">2. </span>Tensor Network Theory", "<span class=\"section-number\">4. </span>TensorKit.jl", "<span class=\"section-number\">5. </span>TensorOperations.jl", "<span class=\"section-number\">9. </span>Algorithms", "<span class=\"section-number\">10. </span>Applications", "<span class=\"section-number\">8. </span>Infinite Matrix Product States", "<span class=\"section-number\">7. </span>Matrix Product Operators", "<span class=\"section-number\">6. </span>Matrix Product States", "Time Evolution", "<span class=\"section-number\">11. </span>References", "Tensor Network Methods with Julia"], "terms": {"introduct": 1, "context": [1, 2], "histori": [1, 2], "purpos": [1, 2], "relev": [1, 2], "second": [1, 10], "quantiz": 1, "tensor": [1, 10], "product": [1, 11], "classic": [1, 6], "multi": 2, "linear": [2, 10], "algebra": 2, "graphic": 2, "notat": 2, "comput": 2, "complex": 2, "test": [3, 4], "basic": 5, "simpl": 5, "updat": [5, 10], "trotter": 5, "tebd": [5, 10], "transvers": 6, "field": 6, "Ising": 6, "model": 6, "discuss": 7, "larg": 7, "base": 7, "vhv19": [7, 11], "gaug": [7, 9, 10], "revisit": 7, "expect": [7, 8, 9, 10], "valu": [7, 8, 9, 10], "structur": 7, "factor": [7, 10], "statist": [8, 12], "mechan": [8, 12], "an": [8, 10], "mpo": 8, "jordan": 8, "block": 8, "form": [8, 10], "mp": [9, 10], "polynomi": 9, "ansatz": 9, "entangl": [9, 10], "area": 9, "law": 9, "freedom": 9, "correl": 9, "In": 10, "case": 10, "quantum": [10, 12], "mani": [10, 12], "bodi": [10, 12], "system": 10, "amount": 10, "solv": 10, "depend": 10, "schroding": 10, "equat": 10, "i": 10, "frac": 10, "partial": 10, "t": 10, "ket": 10, "psi": 10, "hat": 10, "h": 10, "given": 10, "hamiltonian": 10, "initi": 10, "condit": 10, "psi_0": 10, "t_0": 10, "For": 10, "independ": 10, "solut": 10, "exp": 10, "By": 10, "approxim": 10, "oper": [10, 12], "network": 10, "languag": 10, "we": 10, "can": 10, "also": 10, "studi": 10, "real": 10, "dynam": 10, "One": 10, "should": 10, "keep": 10, "mind": 10, "gener": 10, "increas": 10, "state": [10, 11], "so": 10, "practic": 10, "onli": 10, "done": 10, "rel": 10, "modest": 10, "quench": 10, "grow": 10, "linearli": 10, "bond": 10, "dimens": 10, "would": 10, "need": 10, "exponenti": 10, "order": 10, "accur": 10, "follow": 10, "evolv": 10, "below": 10, "ar": 10, "some": 10, "method": [10, 11], "thi": [10, 12], "part": 10, "tutori": [10, 12], "focu": 10, "latter": 10, "two": 10, "The": 10, "variat": 10, "principl": 10, "old": 10, "idea": 10, "originali": 10, "develop": [10, 12], "dirac": 10, "frenkel": 10, "1930": 10, "s": 10, "minim": 10, "2": 10, "parametr": 10, "set": 10, "matric": 10, "a_1": 10, "a_2": 10, "dot": 10, "a_n": 10, "where": 10, "n": 10, "size": 10, "finitemp": 10, "unit": 10, "cell": 10, "infinit": [10, 12], "other": 10, "word": 10, "live": 10, "manifold": 10, "determin": 10, "geometr": 10, "problem": 10, "project": 10, "rh": 10, "onto": 10, "d": 10, "dt": 10, "A": 10, "p": 10, "_": 10, "tangent": [10, 11], "space": [10, 11], "As": 10, "consequ": 10, "never": 10, "leav": 10, "term": 10, "make": 10, "sens": 10, "work": 10, "abov": 10, "level": 10, "try": 10, "give": [10, 12], "complic": 10, "non": 10, "one": 10, "favourit": 10, "differ": 10, "scheme": 10, "requir": 10, "invers": 10, "small": 10, "singular": 10, "thu": 10, "numer": 10, "instabl": 10, "instead": 10, "turn": 10, "natur": 10, "free": 10, "wai": 10, "possibl": 10, "mix": 10, "show": 10, "consist": 10, "sum": 10, "act": 10, "ac": 10, "site": 10, "c": 10, "right": 10, "have": [10, 12], "brought": 10, "insight": 10, "each": 10, "seper": 10, "fact": 10, "ani": 10, "_c": 10, "eff": 10, "a_c": 10, "aproxim": 10, "effect": 10, "integr": 10, "exactli": 10, "idt": 10, "text": 10, "start": [10, 12], "from": 10, "first": 10, "accord": 10, "formula": 10, "qr": 10, "result": 10, "new": 10, "get": [10, 12], "a_l": 10, "en": 10, "via": 10, "tild": 10, "absorb": 10, "a_r": 10, "1": 10, "repeat": 10, "At": 10, "end": 10, "chain": 10, "sinc": 10, "do": 10, "left": 10, "sweep": 10, "e": 10, "up": 10, "mathcal": 10, "o": 10, "perform": 10, "revers": 10, "combin": 10, "yield": 10, "becaus": 10, "adjoint": 10, "forward": 10, "like": 10, "until": 10, "criteria": 10, "converg": 10, "obtain": 10, "howev": 10, "costli": 10, "ha": 10, "iter": 10, "anoth": 10, "option": 10, "exploit": 10, "translat": 10, "invari": 10, "demand": 10, "allow": 10, "us": 10, "thing": 10, "around": 10, "find": 10, "newli": 10, "found": 10, "coupl": 10, "nice": 10, "properti": 10, "reminisc": 10, "all": 10, "trivial": 10, "eigenst": 10, "propto": 10, "whole": 10, "pick": 10, "phase": 10, "equal": 10, "dte": 10, "addit": 10, "conserv": 10, "energi": 10, "perhap": 10, "most": 10, "write": 10, "simpli": 10, "contract": 10, "its": 10, "trunctat": 10, "taylor": 10, "seri": [10, 12], "tau": 10, "3": 10, "boil": 10, "down": 10, "implement": 10, "power": 10, "effici": 10, "lowest": 10, "extens": 10, "exampl": [10, 12], "which": 10, "correspond": 10, "sum_i": 10, "_i": 10, "b": 10, "apprixim": 10, "tabl": 10, "matrix": [10, 11], "multiplc": 10, "rememb": 10, "boundari": 10, "upper": 10, "express": 10, "approx": 10, "d_i": 10, "c_i": 10, "b_": 10, "desir": 10, "trick": 10, "involv": 10, "remov": 10, "third": 10, "multipli": 10, "appropri": 10, "\u03c4": 10, "visualis": 10, "extend": 10, "repres": 10, "besid": 10, "simul": 10, "groundstat": 10, "take": 10, "basi": 10, "lim_": 10, "infti": 10, "sqrt": 10, "braket": 10, "orthogon": 10, "ground": 10, "inde": 10, "expand": 10, "eigenbasi": 10, "e_i": 10, "e_0": 10, "e_1": 10, "e_2": 10, "Then": 10, "limit": 10, "slowest": 10, "vanish": 10, "normal": 10, "c_0": 10, "irrelev": 10, "It": [10, 12], "construct": 10, "thermal": 10, "densiti": 10, "rho": 10, "z": 10, "beta": 10, "constant": 10, "here": 10, "constraint": 10, "particula": 10, "ensur": 10, "posit": 10, "semi": 10, "definit": 10, "physic": [10, 11, 12], "note": [10, 11], "d_k": 10, "pure": 10, "desniti": 10, "introduc": 10, "ancilla": 10, "a_k": 10, "immedeatli": 10, "see": 10, "tr": 10, "_a": 10, "bra": 10, "0": 10, "mathbf": 10, "delta": 10, "m": 10, "2m": 10, "Or": 10, "purif": 10, "how": 10, "mpskit": 10, "mpskitmodel": 10, "h\u2080": 10, "transverse_field_is": 10, "j": 10, "g": 10, "creat": 10, "virtual": 10, "50": 10, "optim": 10, "gs": 10, "infinitemp": 10, "env": 10, "find_groundst": 10, "vump": 10, "maxit": 10, "400": 10, "let": 10, "check": 10, "sx_g": 10, "expectation_valu": 10, "\u03c3\u2093": 10, "length": 10, "e_g": 10, "defin": 10, "algorithm": [10, 12], "alg": 10, "ht": 10, "01": 10, "st_t": 10, "timestep": 10, "lauren": 11, "vanderstraeten": 11, "jutho": 11, "haegeman": 11, "frank": 11, "verstraet": 11, "uniform": 11, "scipost": 11, "lectur": 11, "page": 11, "007": 11, "januari": 11, "2019": 11, "url": 11, "http": 11, "org": 11, "scipostphyslectnot": 11, "7": 11, "arxiv": 11, "1810": 11, "07006": 11, "doi": 11, "10": 11, "21468": 11, "applic": 12, "well": 12, "illustr": 12, "theori": 12, "aim": 12, "hand": 12, "practis": 12, "provid": 12, "code": 12, "showcas": 12, "softwar": 12, "librari": 12, "been": 12, "tensorkit": 12, "jl": 12, "tensoroper": 12, "refer": 12}, "objects": {}, "objtypes": {}, "objnames": {}, "titleterms": {"get": 0, "start": 0, "julia": [0, 12], "quantum": 1, "mani": 1, "bodi": 1, "theori": [1, 2], "tensor": [2, 12], "network": [2, 12], "tensorkit": 3, "jl": [3, 4], "tensoroper": 4, "algorithm": 5, "applic": 6, "infinit": 7, "matrix": [7, 8, 9, 12], "product": [7, 8, 9, 12], "state": [7, 9, 12], "oper": 8, "time": 10, "evolut": 10, "tdvp": 10, "mpo": 10, "imaginari": 10, "finit": 10, "temperatur": 10, "out": 10, "box": 10, "code": 10, "refer": 11, "method": 12, "introduct": 12, "other": 12}, "envversion": {"sphinx.domains.c": 2, "sphinx.domains.changeset": 1, "sphinx.domains.citation": 1, "sphinx.domains.cpp": 6, "sphinx.domains.index": 1, "sphinx.domains.javascript": 2, "sphinx.domains.math": 2, "sphinx.domains.python": 3, "sphinx.domains.rst": 2, "sphinx.domains.std": 2, "sphinx.ext.intersphinx": 1, "sphinxcontrib.bibtex": 9, "sphinx": 56}})