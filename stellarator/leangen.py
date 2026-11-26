"""Generate Lean 4 certificates for closed-form field identities.

Given scalar functions of the coordinates written as SymPy expressions (built from + - * /,
natural powers, named rational sub-expressions and one named square root) and an identity
between them and their first partial derivatives, this emits a Lean file proving the identity.
Nothing in the proof is searched for by Lean:

  * each partial derivative is a `HasDerivAt` term assembled rule by rule from the expression;
  * the named sub-expressions are replaced by variables `S`, `F`, ... in Lean, SymPy divides the
    numerator of the identity by the defining relations (`F^2 = R`, `den * S = num`), and Lean
    only has to check a rational-function identity with `field_simp; ring`.

SymPy finds the proof; Lean checks it. A wrong generator gives a file that fails to compile,
not a wrong theorem. What has to be trusted is the generated `def`s and theorem statements:
read those against the paper.

Supported: any number of named rational sub-expressions (not nested in each other) and at most
one square root, whose radicand may use them. Not supported yet: sin/cos/exp, complex values.
"""
import random

import sympy as sp

# ------------------------------------------------------------------ expression tree


class Node:
    def __init__(self, kind, args=(), value=None):
        self.kind, self.args, self.value = kind, tuple(args), value
        self.free = set().union(*[a.free for a in self.args]) if self.args else set()
        if kind == "sym":
            self.free = {value}

    def dep(self, v):
        return v in self.free


def show(n, ren=None):
    """Node -> fully parenthesised Lean term. `ren` maps symbols (and atom names) to names."""
    ren = ren or {}
    k, a = n.kind, n.args
    if k == "sym":
        return ren.get(n.value, str(n.value))
    if k == "atom":  # a named sub-expression: printed expanded unless it has been generalised
        return ren[n.value] if n.value in ren else show(a[0], ren)
    if k == "num":
        return f"({n.value} : ℝ)"
    if k == "neg":
        return f"(-{show(a[0], ren)})"
    if k == "sqrt":
        return ren[n.value] if n.value in ren else f"(Real.sqrt {show(a[0], ren)})"
    if k == "pow":
        return f"({show(a[0], ren)} ^ {n.value})"
    op = {"add": "+", "sub": "-", "mul": "*", "div": "/"}[k]
    return f"({show(a[0], ren)} {op} {show(a[1], ren)})"


class LeanGen:
    def __init__(self, namespace, params, coords, header=""):
        self.ns, self.params, self.coords, self.header = namespace, list(params), list(coords), header
        self.atoms = {}     # Symbol -> (sympy expr in params/coords, Node)
        self.root = None    # (Symbol F, radicand sympy expr in atoms, Node)
        self.defs, self.trees = {}, {}
        self.out = []

    # -- named sub-expressions
    def atom(self, name, expr):
        """A named rational sub-expression of the coordinates and parameters."""
        S = sp.Symbol(name)
        self.atoms[S] = (expr, Node("atom", [self.build(expr)], value=S))
        return S

    def sqrt(self, name, radicand):
        """The (single) named square root; the radicand may use atoms."""
        if self.root is not None:
            raise NotImplementedError("more than one square root")
        F = sp.Symbol(name)
        self.root = (F, radicand, Node("sqrt", [self.build(radicand)], value=F))
        return F

    def build(self, e):
        """SymPy expression -> Node tree. Its shape does not depend on the differentiation variable."""
        b = self.build
        if e.is_Symbol:
            if e in self.atoms:
                return self.atoms[e][1]
            if self.root is not None and e == self.root[0]:
                return self.root[2]
            return Node("sym", value=e)
        if e.is_Integer:
            return Node("neg", [Node("num", value=-int(e))]) if e < 0 else Node("num", value=int(e))
        if e.is_Rational:
            n = Node("div", [Node("num", value=abs(int(e.p))), Node("num", value=int(e.q))])
            return Node("neg", [n]) if e < 0 else n
        if e.is_Add:
            terms = e.as_ordered_terms()
            pos = [t for t in terms if not t.could_extract_minus_sign()]
            neg = [-t for t in terms if t.could_extract_minus_sign()]
            if pos:
                acc = b(pos[0])
                for t in pos[1:]:
                    acc = Node("add", [acc, b(t)])
            else:
                acc, neg = Node("neg", [b(neg[0])]), neg[1:]
            for t in neg:
                acc = Node("sub", [acc, b(t)])
            return acc
        if e.is_Mul:
            coeff, rest = e.as_coeff_Mul()
            if coeff < 0:
                return Node("neg", [b(-e)])
            num, den = [], []
            if coeff.is_Rational and coeff != 1:
                if coeff.p != 1:
                    num.append(Node("num", value=int(coeff.p)))
                if coeff.q != 1:
                    den.append(Node("num", value=int(coeff.q)))
            for f in sp.Mul.make_args(rest):
                base, ex = f.as_base_exp()
                (den if ex.is_negative else num).append(b(base ** (-ex)) if ex.is_negative else b(f))

            def prod(fs):
                acc = fs[0]
                for f in fs[1:]:
                    acc = Node("mul", [acc, f])
                return acc

            n = prod(num) if num else Node("num", value=1)
            return Node("div", [n, prod(den)]) if den else n
        if e.is_Pow:
            base, ex = e.as_base_exp()
            if ex.is_Integer and ex >= 2:
                return Node("pow", [b(base)], value=int(ex))
            if ex.is_negative:
                return Node("div", [Node("num", value=1), b(base ** (-ex))])
        raise NotImplementedError(f"unsupported expression: {e!r} (name square roots with .sqrt)")

    # -- definitions
    def define(self, name, expr, doc=None):
        self.defs[name], self.trees[name] = expr, self.build(expr)
        args = " ".join(str(s) for s in self.params + self.coords)
        if doc:
            self.out.append(f"/-- {doc} -/")
        self.out.append(f"def {name} ({args} : ℝ) : ℝ :=\n  {show(self.trees[name])}\n")

    def call(self, name, ren=None):
        ren = ren or {}
        return f"{name} " + " ".join(ren.get(s, str(s)) for s in self.params + self.coords)

    # -- side conditions
    def _hyp(self, kind, node):
        key = (kind, show(node))
        if key not in self._hyps:
            self._hyps[key] = f"h{len(self._hyps)}"
        return self._hyps[key]

    def _ne(self, n):
        """Proof term that Node n is non-zero at the point; registers hypotheses as needed."""
        if n.kind == "num":
            return "(by norm_num)"
        if n.kind == "mul":
            return f"(mul_ne_zero {self._ne(n.args[0])} {self._ne(n.args[1])})"
        if n.kind == "pow":
            return f"(pow_ne_zero {n.value} {self._ne(n.args[0])})"
        if n.kind == "sqrt":
            return f"(Real.sqrt_pos.mpr {self._hyp('pos', n.args[0])}).ne'"
        if n.kind == "div":
            return f"(div_ne_zero {self._ne(n.args[0])} {self._ne(n.args[1])})"
        return self._hyp("ne", n)

    def _side_conditions(self, n):
        for a in n.args:
            self._side_conditions(a)
        if n.kind == "div":
            self._ne(n.args[1])
        if n.kind == "sqrt":
            self._hyp("pos", n.args[0])

    # -- derivative terms
    def _d(self, n, v):
        """Proof term of `HasDerivAt (fun t => n[v := t]) _ v`, for n depending on v."""
        k, a = n.kind, n.args
        if k == "sym":
            return f"(hasDerivAt_id' {v})"
        if k == "atom":
            return self._d(a[0], v)
        c = show
        if k == "neg":
            return f"({self._d(a[0], v)}).neg"
        if k == "pow":
            return f"({self._d(a[0], v)}).pow {n.value}"
        if k == "sqrt":
            return f"({self._d(a[0], v)}).sqrt ({self._hyp('pos', a[0])}).ne'"
        l, r = a
        ld, rd = l.dep(v), r.dep(v)
        if k == "add":
            return (f"({self._d(l, v)}).add ({self._d(r, v)})" if ld and rd else
                    f"({self._d(l, v)}).add_const {c(r)}" if ld else f"({self._d(r, v)}).const_add {c(l)}")
        if k == "sub":
            return (f"({self._d(l, v)}).sub ({self._d(r, v)})" if ld and rd else
                    f"({self._d(l, v)}).sub_const {c(r)}" if ld else f"({self._d(r, v)}).const_sub {c(l)}")
        if k == "mul":
            return (f"({self._d(l, v)}).mul ({self._d(r, v)})" if ld and rd else
                    f"({self._d(l, v)}).mul_const {c(r)}" if ld else f"({self._d(r, v)}).const_mul {c(l)}")
        if k == "div":
            if rd:
                num = self._d(l, v) if ld else f"(hasDerivAt_const {v} {c(l)})"
                return f"({num}).div ({self._d(r, v)}) {self._ne(r)}"
            return f"({self._d(l, v)}).div_const {c(r)}"
        raise NotImplementedError(k)

    # -- SymPy side: derivatives with the atoms kept as symbols
    def _expand(self, e):
        if self.root is not None:
            e = e.subs(self.root[0], sp.sqrt(self.root[1]))
        return e.subs({S: ex for S, (ex, _) in self.atoms.items()})

    def _diff(self, e, v):
        d = sp.diff(e, v)
        for S, (ex, _) in self.atoms.items():
            d += sp.diff(e, S) * sp.diff(ex, v)
        if self.root is not None:
            F, R, _ = self.root
            dR = sp.diff(R, v) + sum(sp.diff(R, S) * sp.diff(ex, v) for S, (ex, _) in self.atoms.items())
            d += sp.diff(e, F) * dR / (2 * F)
        return d

    # -- identities
    def identity(self, name, build_identity, doc=None, check_points=3):
        """build_identity(val, D) -> SymPy expression that should vanish, where val(name) is a
        component and D(name, coord) its partial derivative. It is evaluated twice: with real
        expressions (to find the proof) and with placeholders (to print the statement)."""
        self._hyps = {}
        E = build_identity(lambda f: self.defs[f], lambda f, v: self._diff(self.defs[f], v))

        # numerical sanity check of the claim itself, in the original variables
        rng = random.Random(0)
        Efull = self._expand(E)
        syms = sorted(Efull.free_symbols, key=str)
        done = tries = 0
        while done < check_points and tries < 200:
            tries += 1
            val = sp.N(Efull.subs({s: sp.Rational(rng.randint(20, 90), 100) for s in syms}), 30)
            if not val.is_real:
                continue  # outside the domain of the square root
            done += 1
            if abs(val) > 1e-20:
                raise ValueError(f"identity {name} is false: value {val} at a test point")
        if done < check_points:
            raise ValueError(f"identity {name}: no real test points found")

        # placeholders for printing the statement
        ph_val = {f: sp.Symbol(f"VAL_{f}") for f in self.defs}
        ph_d = {}

        def D(f, v):
            return ph_d.setdefault((f, v), sp.Symbol(f"D_{f}_{v}"))

        stmt = build_identity(lambda f: ph_val[f], D)
        used = [k for k in ph_d if ph_d[k] in stmt.free_symbols]
        ren = {ph_val[f]: f"({self.call(f)})" for f in self.defs}
        for (f, v), s in ph_d.items():
            ren[s] = f"(deriv (fun t : ℝ => {self.call(f, {v: 't'})}) {v})"
        stmt_lean = show(self.build(stmt), ren)

        # derivative proofs
        lines = []
        for i, (f, v) in enumerate(used):
            tree = self.trees[f]
            term = self._d(tree, v) if tree.dep(v) else f"hasDerivAt_const {v} {show(tree)}"
            lines.append(f"  have d{i} : HasDerivAt (fun t : ℝ => {self.call(f, {v: 't'})}) _ {v} := {term}")
        for f in self.defs:
            self._side_conditions(self.trees[f])
        if used:
            lines.append("  rw [" + ", ".join(f"d{i}.deriv" for i in range(len(used))) + "]")
            lines.append("  clear " + " ".join(f"d{i}" for i in range(len(used))))
        # the inferred derivative values contain unreduced `(fun t => ..) x`, `↑2` and `^ (2 - 1)`
        lines.append("  simp only [Pi.add_apply, Pi.sub_apply, Pi.mul_apply, Pi.pow_apply, Pi.div_apply, "
                     "Pi.neg_apply, Nat.cast_ofNat, Nat.add_one_sub_one, Nat.reduceSub, pow_one, mul_one, one_mul, "
                     + ", ".join(self.defs) + "]")

        # algebra certificate: (multiplier) * numerator = sum of q_i * relation_i
        numer, denom = sp.fraction(sp.together(E))
        numer = sp.expand(numer)
        terms, mult, r = [], sp.Integer(1), numer
        gen_ren = {}
        if self.root is not None:
            F, R, Fnode = self.root
            hpos = self._hyp("pos", Fnode.args[0])
            lines.append(f"  have hFsq := Real.sq_sqrt (le_of_lt {hpos})")
            lines.append(f"  have hFne := (Real.sqrt_pos.mpr {hpos}).ne'")
            q, r = sp.div(r, F**2 - R, F)
            r = sp.expand(r)
            terms.append((sp.expand(q), F**2 - R, "hF"))
        for S, (ex, Snode) in self.atoms.items():
            lines.append(f"  generalize hgen_{S} : {show(Snode)} = {S} at *")
            gen_ren[S] = str(S)
            sn, sd = sp.fraction(sp.together(ex))
            rel = sp.expand(S * sd - sn)
            lines.append(f"  have hrel_{S} : {show(self.build(rel), gen_ren)} = 0 := by\n"
                         f"    rw [← hgen_{S}]; field_simp; ring")
            if r != 0:
                k = max(sp.degree(r, S), 0)  # pdiv multiplies by lc^(deg r - deg rel + 1) = sd^deg r
                q, rem = sp.pdiv(r, rel, S)
                lc = sd**k
                terms = [(sp.expand(tq * lc), trel, th) for tq, trel, th in terms]
                mult *= lc
                terms.append((sp.expand(q), rel, f"hrel_{S}"))
                r = sp.expand(rem)
        if r != 0:
            raise ValueError(f"identity {name}: could not reduce the numerator to zero (remainder {r})")
        if self.root is not None:
            F, R, Fnode = self.root
            lines.append(f"  generalize hgen_{F} : {show(Fnode, gen_ren)} = {F} at *")
            gen_ren[F] = str(F)
            lines.append(f"  have hF : {show(self.build(F**2 - R), gen_ren)} = 0 := by rw [hFsq]; ring")
        live = [(q, rel, h) for q, rel, h in terms if q != 0]
        if live:
            Z = " + ".join(f"{show(self.build(sp.factor(q)), gen_ren)} * {show(self.build(rel), gen_ren)}"
                           for q, rel, _ in live)
            Dl = show(self.build(sp.factor(mult * denom)), gen_ren)
            lines.append(f"  refine (?_ : _ = ({Z}) / {Dl}).trans (by rw [" + ", ".join(h for _, _, h in live) + "]; simp)")
        lines += ["  field_simp", "  ring"]

        hyps = " ".join(f"({h} : {'0 < ' + s if kind == 'pos' else s + ' ≠ 0'})" for (kind, s), h in self._hyps.items())
        args = " ".join(str(s) for s in self.params + self.coords)
        if doc:
            self.out.append(f"/-- {doc} -/")
        self.out.append(f"theorem {name} ({args} : ℝ)\n    {hyps} :\n    {stmt_lean} = 0 := by")
        self.out += lines + [""]

    def render(self):
        body = "\n".join(self.out)
        axioms = "\n".join(f"#print axioms {self.ns}.{l.split()[1]}" for l in self.out if l.startswith("theorem "))
        return (f"import Mathlib\n\n/-!\n{self.header}\n\nGENERATED by stellarator/leangen.py. Do not edit by hand.\n-/\n\n"
                "set_option linter.unusedSimpArgs false\nset_option linter.unusedVariables false\n"
                "set_option linter.unusedTactic false\n\n"
                f"noncomputable section\n\nnamespace {self.ns}\n\n{body}\nend {self.ns}\n\n{axioms}\n")
