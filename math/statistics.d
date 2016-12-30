module math.statistics;

private import std.math;
private import std.algorithm;
private import std.stdio;
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
