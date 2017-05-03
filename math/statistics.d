module math.statistics;

private import std.math;
private import std.algorithm;
private import std.stdio;
private import std.format;
private import jive.array;
private import jive.internal;

struct Histogram
{
    double low = 0, high = 0;
    int nBins = 0;

    Array!long hist;
    long count = 0;
    long countLow = 0;
    long countHigh = 0;
    double sum = 0;
    double sum2 = 0;
    double min = double.infinity;
    double max = -double.infinity;

    this(double low, double high, int nBins)
    {
        this.low = low;
        this.high = high;
        this.nBins = nBins;
        hist.resize(nBins, 0);
    }

    void add(double x, long n = 1)
    {
        count += n;
        sum += n*x;
        sum2 += n*x*x;
        min = std.algorithm.min(min, x);
        max = std.algorithm.max(max, x);

        if(x < low)
            countLow += n;
        else if(x >= high)
            countHigh += n;
        else
            hist[cast(int)((x-low)/(high-low)*nBins)] += n;
    }

    double avg() const nothrow @property @safe
    {
        return sum / count;
    }

    double var() const nothrow @property @safe
    {
        return sum2/count - sum/count*sum/count;
    }

    /** iterate over center-of-bin / count-of-bin */
    int opApply(int delegate(double, double) dg) const
    {
        int r = 0;
        foreach(i, x; hist)
            if((r = dg(low + (i+0.5)*(high-low)/nBins, x)) != 0)
                break;
        return r;
    }

    void write() const
    {
        if(countLow)
            writefln("<:\t%s", countLow);
        foreach(x, y; this)
            writefln("%s:\t%s", x, y);
        if(countHigh)
            writefln(">:\t%s", countHigh);
        writefln("all:\t%s", count);
        writefln("avg = %s +- %s", avg, sqrt(var));
    }
}

/** put n-dimensional data points in, get average/variance/covariance out */
struct Statistics(size_t n = 1)
{
    double count = 0;
    double[n] sum;
    double[n][n] sum2;

    /** add a new data point */
    void add(Times!(n, double) xs, double w = 1)
    {
        // TODO: figure out a way to do this statically
        if(count == 0)
        {
            foreach(ref x; sum)
                x = 0;
            foreach(ref l; sum2)
                foreach(ref x; l)
                    x = 0;
        }

        count += w;
        foreach(i, x; xs)
            sum[i] += w*x;
        foreach(i, x; xs)
            foreach(j, y; xs)
                sum2[i][j] += w*x*y;
    }

    /** average in dimension i */
    double avg(size_t i = 0)() const nothrow @property @safe
    {
        return sum[i]/count;
    }

    /** variance in dimension i */
    double var(size_t i = 0)() const nothrow @property @safe
    {
        return cov!(i,i);
    }

    /** standard deviation in dimension i */
    double stddev(size_t i = 0)() const nothrow @property @safe
    {
        return sqrt(var!i);
    }

    /** covariance between diemnsions i and j */
    double cov(size_t i = 0, size_t j = 1)() const nothrow @property @safe
    {
        return sum2[i][j]/count - sum[i]/count * sum[j]/count;
    }

    /** correlation between diemnsions i and j */
    double corr(size_t i = 0, size_t j = 1)() const nothrow @property @safe
    {
        return cov!(i,j) / sqrt(var!i*var!j);
    }
}

/**
 * A (mean, variance) tuple. Can be interpreted as a random variable,
 * or a measurement with error. Standard arithmetic is overloaded to propagate
 * the error term, but there are severe limitations:
 *    -  but it assumes there is no correlation between variables,
 * so use carefully. For example "x+x" is not the same as "2*x" (and the latter
 * one is generally the correct one).
 */
struct Var
{
    double mean = double.nan;
    double var = 0;

    /** standard deviation sqrt(variance) */
    double stddev() const pure nothrow @property @safe
    {
        return sqrt(var);
    }

    Var opBinary(string op)(double b) const pure nothrow @safe
    {
        switch(op)
        {
            case "+": return Var(mean + b, var);
            case "-": return Var(mean - b, var);
            case "*": return Var(mean*b, var*b*b);
            case "/": return Var(mean/b, var/(b*b));
            default: assert(false);
        }
    }

    Var opBinaryRight(string op)(double a) const pure nothrow @safe
    {
        switch(op)
        {
            case "+": return Var(a + mean, var);
            case "-": return Var(a - mean, var);
            case "*": return Var(a * mean, a*a*var);
            default: assert(false);
        }
    }

    Var opBinary(string op)(Var b) const pure nothrow @safe
    {
        switch(op)
        {
            case "+": return Var(mean + b.mean, var + b.var);
            case "-": return Var(mean - b.mean, var + b.var);
            case "*": return Var(mean * b.mean, mean*mean*b.var + b.mean*b.mean*var + var*b.var);
            default: assert(false);
        }
    }

    /** returns human readable string "mean +- stddev" */
    string toString() const @property @safe
    {
        return format("%s +- %s", mean, stddev);
    }
}
