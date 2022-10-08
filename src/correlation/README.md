# Correlated variables

Let $X$ and $Z$ be independent Gaussian randomw varaibles distribtued $\sim\mathcal{N}(1,1)$.

Construct the dependent random variable, $Y=\frac{1}{2}(X+Z)$.

Look at the correlation function

$$\rho(X,Y) = \frac{\langle XY \rangle -\langle X \rangle\langle Y \rangle}{\sqrt{(\langle X^2 \rangle -\langle X \rangle^2)(\langle {Y}^2 \rangle -\langle Y \rangle^2)}},$$

The result should be $\rho(X,Y) = 1/\sqrt(2)$. 

This is not what I get out of Signaloid. In fact, I obtain $\rho(X,Y) = 0 $, as if $X$ and $Y$ were independent. I have tried different processors (in particular CO-L+) as I thought that this is the problem that "Autocorrelation Tracking" is made to solve, but this did not change things.

I also find that `Signaloid-Demo-General-C` behaves in a similar way, in that the distribution I get for $c$ is as if we had evaluated $c=(a+b)/(a'-b')$ where $a'$ and $b'$ are distributed identically to $a$ and $b$ but independently of them.
