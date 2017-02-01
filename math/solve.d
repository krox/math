module math.solve;

/**
 * solving of numerical equations f(x) = 0.
 * So far, only 1-dimensional problems.
 * The template type of all functions is assumed to be float/double/real/Floating/...
 */

private import std.math;
private import std.functional : unaryFun;
private import std.algorithm;
private import math.numerics;

int sign(T)(T x) pure nothrow
{
    if(x < 0)
        return -1;
    if(x > 0)
        return 1;
    return 0;
}

/** bisection method. fully robust but converges only linearly */
T solveBisect(T, alias _f)(T a, T b, T eps = 4 * T.epsilon, ref int maxIter = *new int(int.max))
{
    alias unaryFun!_f f;

    if(maxIter < 2)
        throw new Exception("max iterations reached without convergence");

    T fa = f(a);
    T fb = f(b);
    maxIter -= 2;
    T m, fm;

    if(sign(fa)*sign(fb) > 0)
        throw new NumericsException("invalid bracket");

    while(true)
    {
        if(approxEqual(a, b, eps))
            return b;

        if(maxIter < 1)
            throw new Exception("max iterations reached without convergence");

        m = 0.5*(a+b);
        fm = f(m);
        --maxIter;

        if(sign(fm) == sign(fa))
            { a = m; fa = fm; }
        else if(sign(fm) == sign(fb))
            { b = m; fb = fm; }
        else
            return m;
    }

    return a;
}

/** Dekker method combining secant and bisection. Might converge extremely slow for nasty functions */
T solveDekker(T, alias _f)(T a, T b, T eps = 4 * T.epsilon, ref int maxIter = *new int(int.max))
{
    alias unaryFun!_f f;

    if(maxIter < 2)
        throw new Exception("max iterations reached without convergence");

    /**
    b = current best
    b2 = previous best
    a = contrapoint for a (maybe equal to b2)
    **/

    T fa = f(a);
    T fb = f(b);
    maxIter -= 2;

    if(sign(fa)*sign(fb) > 0)
        throw new NumericsException("invalid bracket");

    if(abs(fa) < abs(fb)) // b should be the best guess
    {
        swap(a, b);
        swap(fa, fb);
    }

    T b2 = a;
    T fb2 = fa;
    T m, fm;

    while(true)
    {
        if(approxEqual(a, b, eps))
            return b;

        if(maxIter < 1)
            throw new Exception("max iterations reached without convergence");

        // evaluate secant method or midpoint
        m = (b*fb2 - b2*fb) / (fb2 - fb);
        if((m >= a && m >= b) || (m <= a && m <= b))
            m = 0.5*(a+b);
        fm = f(m);
        --maxIter;

        // update brackets
        if(sign(fm)*sign(fb) < 0)
            { a = b; fa = fb; }
        b2 = b;
        fb2 = fb;
        b = m;
        fb = fm;
    }
}
