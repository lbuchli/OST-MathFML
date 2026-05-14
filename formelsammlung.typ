#import "@preview/fletcher:0.5.8" as fletcher: diagram, node, edge


#set page(flipped: true, columns: 2, margin: 1cm)

#set table(stroke: none)

#set enum(tight: true, spacing: 4pt, numbering: n => box(height: 1em, circle(radius: 0.4em, align(center + horizon, text(8pt)[#n]))))

#table(
    columns: (20%, auto),
    [Matrix addition],
    table(
        columns: (auto, auto),
        $(A + B) + C = A + (B + C)$,
        [Associativity of addition],
        $A + B = B + A$,
        [Commutativity of addition],
    ),
    [Matrix multiplication],
    [
        #block($
            (a_(i j)) = A in RR^(r times n),
            (b_(i j)) = B in RR^(n times c),
            (c_(i j)) = C in RR^(r times c)
        $)
        #stack(dir: ltr,
            diagram(spacing: (4mm, 4mm), {
                let (r, n, c) = ((0, 0), (1, 0), (2, 0))
                node(r, $RR^r$)
                node(n, $RR^n$)
                node(c, $RR^c$)
                edge(r, n, "->", $A$)
                edge(n, c, "->", $B$)
                edge(r, c, "->", $C$, bend: -15deg, shift: -5pt)
            }),
            h(1cm),
            stack(dir: ttb,
                $C = A B <=> c_(i j) = sum^n_(k=1) a_(i k) b_(k j)$,
                v(0.5cm),
                table(
                    columns: (auto, auto),
                    $(A B) C = A (B C)$,
                    [Associativity],
                    $A (B + C) = A B + A C$,
                    [Distributivity]
                )
            )

        )
    ],
    [Matrix transposition],
    [#block($
    dot^T : cases(
    RR^(n times m) -> RR^(m times n),
    a_(i j) |-> a_(j i)
    )
    $)
    #block([$
        (A^T)^T &= A &quad    (lambda A)^T &= lambda A^T \
        (A + B)^T &= A^T + B^T &quad (A B)^T &= B^T A^T \
    $ $M M^T = bold(1) <=> M in RR^(n times n) "is orthogonal"$])
    ],
    [Kroenecker $delta$],
    block($ delta_(i j) = cases(
        0 quad "if" i != j,
        1 quad "if" i = j
    ) $),
    [Image],
    block($ im A = { A x | x in RR^n } quad A in RR^(m times n)$),
    [Kernel],
    block($ ker A = { x | x in RR^n, A x = bold(0) } quad A in RR^(m times n)$),
    [Dimension],
    block([
        $ dim A = n quad A in RR^(m times n) $
        $ underbrace(dim im A, "rank (= number of " \ " independent columns)") + underbrace(dim ker A, "nullity") = n $
    ]),
    [Graph],
    block($
        "graph" : cases(
        (RR^n -> RR^m) -> (RR^n -> RR^(n + m)),
        f |-> (x |-> vec(x, f(x)))
        )
    $),
    [Taylor Series],
    [#block($ sum_(n=0)^infinity frac(f^(n"th")(a), n!) (x - a)^n $)
    where $f^(n"th")(a)$ is the $n$th derivative of $f$
    ],
    [Gradient],
    [
        #block($gradient f(x_1, x_2, ...) = (partial / (partial x_1) f, partial / (partial x_2) f, ...)$)
        #block($
            gradient (f + g) &= gradient f + gradient g &quad gradient (f - g) &= gradient f - gradient g \
            gradient (f dot g) &= g gradient f + f gradient g &quad gradient f/g &= frac(g gradient f - f gradient g, g^2) \

        $)
    ],
    [Gradient Descent],
    table(
        columns: (auto, auto),
        $x_(i+1) = x_i - gamma gradient f(x_i)$, $"("gamma": step size)"$,
        $|x_(i+1) - x_i| < epsilon$, [Stop condition],
    ),
    [Jacobian Matrix],
    stack(dir: ltr,
        $
            f : RR^n -> RR^m quad x_0 in RR^n \
            J_f (x_0) : RR^(m times n) \
            J_f (bold(x_0))_(i j) = frac(partial, partial x_j) f_i (x_0)
        $,
        diagram(spacing: (5mm, 5mm), {
            let (ern, erm, erk, mrn, mrm, mrk) = ((0,0), (0,1), (0,2), (4,0), (4,1), (4,2))
            node(ern, $(RR^n, bold(x)_0)$)
            node(erm, $(RR^m, g(bold(x)_0))$)
            node(erk, $(RR^k, f(g(bold(x)_0)))$)
            node(mrn, $RR^n$)
            node(mrm, $RR^m$)
            node(mrk, $RR^k$)
            edge(ern, erm, "->", $g$, label-side: right)
            edge(erm, erk, "->", $f$, label-side: right)
            edge(ern, erk, "->", $f compose g$, bend: +40deg, shift: 10pt, label-angle: left)
            edge(mrn, mrm, "->", $J_g (bold(x)_0)$, label-side: left)
            edge(mrm, mrk, "->", $J_f (g(bold(x)_0))$, label-side: left)
            edge(mrn, mrk, "->", $J_f (g(bold(x)_0)) J_g (bold(x)_0)$, bend: -40deg, shift: -5pt, label-angle: left)
            edge((1, 1), (1.5, 1), "=>")
        })
    ),
    [Hessian Matrix],
    block($
        f : RR^n -> RR quad bold(x) in RR^n quad H(bold(x)) : RR^(n times n) \
        H(bold(x))_(i j) = partial / (partial bold(x)_i) partial / (partial bold(x)_j) f(bold(x)) = partial / (partial bold(x)_j) partial / (partial bold(x)_i) f(bold(x))
    $),
    [Quadratic Form],
    block($
        M in RR^(n times n) quad q_M (bold(v)) = bold(v)^T M bold(v) \
        H = 1/2 (M + M^T) "upholds" forall bold(v) in RR^n .  q_M (bold(v)) = q_H (bold(v)) \
        H "is symmetric"
    $),
    [Polynomial Quadratic Form],
    block($
        q_M (v) = q_M vec(x_1, x_2, ...) = vec(x_1, x_2, ...)^T M vec(x_1, x_2, ...) -> "multiply out"
    $),
    [Completing Squares],
    table(
        columns: (auto, auto),
    ),
    table.hline(stroke: 0.5pt),
    [Finding stationary points],
    [
        + Solve for $gradient f (bold(x_0)) = bold(0)$
    ],
    [Finding extrema],
    [
        + Find stationary points. Those are the candidates.
        + Calculate the Hessian matrix at those points.
        + Calculate the polynomial quadratic form of the Hessian matrix.
        + Use completing squares one variable at a time.
        + Substitue to get a polynomial of the form $c_1a^2 + c_2b^2 + ...$
        $
        & forall v in RR^n\{0\} . q_H(v) > 0 quad &&->&& "local maximum (positive definite)" \
        & forall v in RR^n\{0\} . q_H(v) < 0 &&->&& "local minimum (negative definite)" \
        & "both, depending on" v &&->&& "saddle point (indefinite)"
        $
    ],

)

/*
== Curves

If all $n$ components of some curve $c : I -> RR^n$ are continuous, the curve is continuous. Limits, derivatives and linearization ($l_t (x) = f'(t) x + f(t) $) can be computed componentwise.

For some curve $c : I -> RR^n$ and some arbitrary bijective $h : I -> I$, $im (c compose h) = im c$.

== Hyper Surfaces

Some function $f : D -> RR$, where
- Differentiable everywhere
- Something with $epsilon$

#table(
    columns: (25%, auto),
    [Partial derivative of $i$th coordinate function],
    block($

    $)
)

*/
