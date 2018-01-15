module math.autodiff;

import std.math;
import std.functional;

/** represents numbers of the form a + b*eps with eps^^2 = 0 */
struct Dual(T)
{
	T a, b;

	this(T a, T b = T(0))
	{
		this.a = a;
		this.b = b;
	}

	Dual opUnary(string op)() const
	{
		     static if(op == "+") return Dual(+a, +b);
		else static if(op == "-") return Dual(-a, -b);
		else static assert(false);
	}

	Dual inverse() const
	{
		return Dual(1/a, -b/(a*a));
	}

	Dual opBinary(string op)(T rhs) const
	{
		     static if(op == "+") return Dual(a + rhs, b);
		else static if(op == "-") return Dual(a - rhs, b);
		else static if(op == "*") return Dual(a * rhs, b * rhs);
		else static if(op == "/") return Dual(a / rhs, b / rhs);
		else static assert(false);
	}

	Dual opBinaryRight(string op)(T lhs) const
	{
		     static if(op == "+") return Dual(lhs + a, b);
		else static if(op == "-") return Dual(lhs - a, -b);
		else static if(op == "*") return Dual(lhs * a, lhs * b);
		else static if(op == "/") return Dual(lhs / a, -lhs*b / (a*a));
		else static assert(false);
	}

	Dual opBinary(string op)(Dual rhs) const
	{
		     static if(op == "+") return Dual(a + rhs.a, b + rhs.b);
		else static if(op == "-") return Dual(a - rhs.a, b - rhs.b);
		else static if(op == "*") return Dual(a * rhs.a, a * rhs.b + b * rhs.a);
		else static if(op == "/") return Dual(a / rhs.a, (b*rhs.a - a*rhs.b)/(rhs.a*rhs.a));
		else static assert(false);
	}

	Dual opOpAssign(string op, S)(S s)
	{
		this = this.opBinary!op(s);
		return this;
	}

	Dual sqr() const
	{
		return Dual(a*a, 2*a*b);
	}

	Dual sqrt() const
	{
		return Dual(a.sqrt, 0.5/a.sqrt*b);
	}

	Dual exp() const
	{
		return Dual(a.exp, a.exp*b);
	}

	Dual log() const
	{
		return Dual(a.log, b/a);
	}

	Dual sin() const
	{
		return Dual(a.sin, a.cos*b);
	}

	Dual cos() const
	{
		return Dual(a.cos, -a.sin*b);
	}
}

/** represents numbers of the form a + b*eps + c*eps*eps with eps^^3 = 0 */
struct Triple(T)
{
	T a, b, c;

	this(T a, T b = T(0), T c = T(0))
	{
		this.a = a;
		this.b = b;
		this.c = c;
	}

	Triple opUnary(string op)() const
	{
		     static if(op == "+") return Triple(+a, +b, +c);
		else static if(op == "-") return Triple(-a, -b, -c);
		else static assert(false);
	}

	Triple inverse() const
	{
		return Triple(1/a, -b/(a*a), (b*b-a*c)/(a*a*a));
	}

	Triple opBinary(string op)(T rhs) const
	{
		     static if(op == "+") return Triple(a + rhs, b, c);
		else static if(op == "-") return Triple(a - rhs, b, c);
		else static if(op == "*") return Triple(a * rhs, b * rhs, c * rhs);
		else static if(op == "/") return Triple(a / rhs, b / rhs, c / rhs);
		else static assert(false);
	}

	Triple opBinaryRight(string op)(T lhs) const
	{
		     static if(op == "+") return Triple(lhs + a, b, c);
		else static if(op == "-") return Triple(lhs - a, -b, -c);
		else static if(op == "*") return Triple(lhs * a, lhs * b, lhs * c);
		else static if(op == "/") return lhs * inverse;
		else static assert(false);
	}

	Triple opBinary(string op)(Triple rhs) const
	{
		     static if(op == "+") return Triple(a + rhs.a, b + rhs.b, c + rhs.c);
		else static if(op == "-") return Triple(a - rhs.a, b - rhs.b, c - rhs.c);
		else static if(op == "*") return Triple(a * rhs.a, a * rhs.b + b * rhs.a, a * rhs.c + b*rhs.b + c*rhs.a);
		else static if(op == "/") return this * rhs.inverse;
		else static assert(false);
	}

	Triple opOpAssign(string op, S)(S s)
	{
		this = this.opBinary!op(s);
		return this;
	}

	Triple sqr() const
	{
		return Triple(a*a, 2*a*b, 2*a*c + b*b);
	}

	Triple sqrt() const
	{
		return Triple(a.sqrt, 0.5/a.sqrt*b, 0.5/a.sqrt*c - a.pow(-1.5)*b*b/8.0);
	}

	Triple exp() const
	{
		return Triple(a.exp, a.exp*b, a.exp*c + 0.5*a.exp*b*b);
	}

	Triple log() const
	{
		return Triple(a.log, b/a, c/a - 0.5*b*b/(a*a));
	}

	Triple sin() const
	{
		return Triple(a.sin, a.cos*b, a.cos*c - 0.5*a.sin*b*b);
	}

	Triple cos() const
	{
		return Triple(a.cos, -a.sin*b, -a.sin*c - 0.5*a.cos*b*b);
	}
}

alias sqrt = std.math.sqrt;
alias exp = std.math.exp;
alias log = std.math.log;
alias sin = std.math.sin;
alias cos = std.math.cos;

Dual!T sqr(T)(Dual!T x) { return x.sqr; }
Dual!T sqrt(T)(Dual!T x) { return x.sqrt; }
Dual!T exp(T)(Dual!T x) { return x.exp; }
Dual!T log(T)(Dual!T x) { return x.log; }
Dual!T sin(T)(Dual!T x) { return x.sin; }
Dual!T cos(T)(Dual!T x) { return x.cos; }

Triple!T sqr(T)(Triple!T x) { return x.sqr; }
Triple!T sqrt(T)(Triple!T x) { return x.sqrt; }
Triple!T exp(T)(Triple!T x) { return x.exp; }
Triple!T log(T)(Triple!T x) { return x.log; }
Triple!T sin(T)(Triple!T x) { return x.sin; }
Triple!T cos(T)(Triple!T x) { return x.cos; }

double[2] diff(alias F)(double x)
{
	auto x_ = Dual!double(x, 1);
	auto y = unaryFun!F(x_);
	return [y.a, y.b];
}

double[3] diff2(alias F)(double x)
{
	auto x_ = Triple!double(x, 1, 1);
	auto y = unaryFun!F(x_);
	return [y.a, y.b, 2*(y.c-y.b)];
}
