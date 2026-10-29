mms = raw"""
  # New Farey word formula
  We define Farey words $W_{r/s}$ of the free group $F_2 = \langle f, g \rangle$ via splitting sequence.
  They can also be defined from *Christoffel words* via the *S-sequence* which is defined as follows:
  If $\Omega_{r /s}$ is the Christoffel word defined via the inductive formulae:
  $$
  \begin{align}
  \Omega_{\alpha\oplus\beta} &= \Omega_{\alpha}\Omega_{\beta} &\alpha &< \beta, &\\
  \Omega_{0 /1} &= a; &\Omega_{1 /1} &= ab; & \Omega_{1 /0} &= b
  \end{align}
  $$
  Where $\alpha\oplus\beta$ is 'Farey addition' which is only defined if $\alpha=\frac{r}{s},\beta = \frac{t}{u}\in \mathbb{Q}$ are *Farey Neighbours* :  i.e. if $st-ru = \pm1$ in which case $\alpha\oplus\beta = \frac{{r+t}}{s+u}$
  Each $\Omega$ is of the form:
  $$
  \Omega_{r/s} = a^{k_1} b\,a^{k_2} b \dots \,a^{k_r} b
  $$

  $$
  (\textrm{e.g. } \Omega_{3/8} = a^3 b\,a^3 b\,a^2 b\,)
  $$

  And we define the associated S-sequence by $S[r/s] = \text{S}[\Omega_{r /s}^2] = [k_1, k_2, \dots, k_r,k_{1}, k_{2}, \dots, k_{r}]$. The  Farey word is defined taking an alternating sequence of the generators $f^{\pm}$ and $g ^{\pm}$ and for each of them choosing their sign based on the $k_i$ in the S-sequence

  $$
  W_{r/s} = fg\overset{({k_{1}}\text{ times})}{ \dots f} GFG\overset{({k_{2}}\text{ times})}{\dots F}gf\overset{({k_{2}}\text{ times})}{\dots fg\dots}
  $$

  For instance, for $\alpha=3/8$, $S[3 /8] =\text{S}[({a^3 b\,a^3 b\,a^2 b})^{2}] = [3,3,2,3,3,2]$ so
  $$
  W_{3 /8} = \overset{(3)}{fgf}.\overset{(3)}{GFG}.\overset{(2)}{fg}.\overset{(3)}{FGF}.\overset{(2)}{gfg}.\overset{(2)}{FG}
  $$

  > Suppose $\gamma = [a_{1}, a_{2},a_{3},\dots, a_{n}, l]$ as a continued fraction with $l = 2k$ or $2k+1$.
  >
  > $\gamma$ can be written as: $$\gamma = \alpha\oplus^{2k}\beta$$
  > Where $\alpha = [a_{1},a_{2},a_{3}\dots,a_{n}]$ and :

  > $$
   \beta = \left\{\begin{array}\\
   [a_{1}, a_{2},a_{3},\dots, a_{n}+1] &l = 2k+1 \\
   \\
   [a_{1}, a_{2},a_{3},\dots, a_{n-1}] &l = 2k
   \end{array} \right.
   $$
  >
  > Then if $n$ is odd (so that $\alpha<\beta$ ):
  > $$
   W_{\gamma} = W_{\alpha}^k W_{\beta} {(W_{\alpha}^f)}^{-k}
   $$
  > and if $n$ is even (so $\alpha>\beta$) then
  > $$
   W_{\gamma} = {(W_{\alpha}^F)}^{-k} W_{\beta} W_{\alpha}^k
   $$
  """
mds = md"$mms"
