#set page(flipped: true, columns: 2, margin: 1cm)

#set table(stroke: none)

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
    [#block($
    (a_(i j)) = A in RR^(r times c),
    (b_(i j)) = B in RR^(r' times c'),
    (c_(i j)) = C in RR^(r' times c') \
    C = A B quad => quad c_(i j) = sum^(c=r')_(k=1) a_(i k) b_(k j)
    $)
    #table(
        columns: (auto, auto),
        $(A B) C = A (B C)$,
        [Associativity of multiplication],
        $A (B + C) = A B + A C$,
        [Distributivity]
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
        (A + B)^T &= A^T + B^T &quad (A B)^T &= A^T B^T \
    $ $M M^T = bold(1) <=> M in RR^(n times n) "is orthogonal"$])
    ],
    [Kroenecker $delta$],
    block($ delta_(i j) = cases(
        0 quad "if" i != j,
        1 quad "if" i = j
    ) $),
    [Image],
    block($ im A = { A x | x in RR^n }$),
    [Kernel],
    block($ ker A = { x | x in RR^n, A x = bold(0) } $),
    [Dimension],
    block([
        $ dim space (A in RR^(n times n)) = n  $
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
    block($
        f : RR^n -> RR^m quad x_0 in RR^n quad J_f (x_0) : RR^(m times n) \
        J_(i j) = partial_x_i f_j (x_0)
    $),
    [Quadratic Form],
    block($q_M (v) = v^T M v quad M in RR^(n times n)$),
    [Hessian Matrix],
    block($ H_(i j) (x) = partial / (partial x_i) partial / (partial x_j) f(x) $),
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
