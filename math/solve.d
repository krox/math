module math.solve;

/**
 * Numerically solve an equations f(x) = 0.
 * The template type of all functions is assumed to be float/double/real/Floating/...
 */

private import std.math : abs, signbit, isNaN;
private import std.functional : unaryFun;
private import std.algorithm;
private import math.numerics;

/**
 * General purpose method for solving f(x) = 0.
 * The result will be exact to full precision of T (assuming f is smooth).
 */
T solve(T, alias _f)(T a, T b)
{
    alias unaryFun!_f f;

    if(a > b)
        swap(a, b);

    assert(!isNaN(a) && !isNaN(b));

    T fa = f(a);
    T fb = f(b);
    assert(!isNaN(fa) && !isNaN(fb));

    if(fa == 0)
        return a;
    if(fb == 0)
        return b;

    assert(signbit(fa) != signbit(fb));

    return solveSecant!(T,_f)(a, b, fa, fb, 100); // TODO: the limit should depend on T
}

/**
 * Secant method with fallback to bisection if necessary.
 * conditions: a < b and f(a)*f(b) < 0
 */
T solveSecant(T, alias _f)(T a, T b, T fa, T fb, int maxIter)
{
    alias unaryFun!_f f;

    assert(a < b);
    auto s = signbit(fa);
    assert(signbit(fa) != signbit(fb));

    T left = a;
    T right = b;

    // a should be the best guess
    if(abs(fb) < abs(fa))
    {
        swap(a, b);
        swap(fa, fb);
    }

    while(maxIter --> 0)
    {
        // choose new point c
        T c = (b*fa - a*fb) / (fa - fb);   // secant method
        if(!(left < c && c < right)) // outside bracket (or nan) -> fall back to bisection
        {
            c = ieeeMean(left, right);
            assert(left <= c && c <= right);
            if(c == left || c == right) // there is no further floating point number between a and b -> we are done
                return a;
        }

        // evaluate f at new point
        b = a;
        fb = fa;
        a = c;
        fa = f(c);
        assert(!isNaN(fa));
        if(fa == 0)
            return a;

        // update brackets
        if(signbit(fa) == s)
            left = a;
        else
            right = a;
    }

    // TODO: this can be avoided by falling back to bisection not only when
    // when secant method leads outside of bracket, but also if previous step(s)
    // have been bad.
    throw new NumericsException;
}
