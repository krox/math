module math.histogram;

private import std.math;
private import std.algorithm;
private import std.stdio;
private import jive.array;

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
    }
}
